!-------------------------------------------------------------------------------
!> module GRIDTRANS
!!
!! @par Description
!!          Grid transfer module
!!          Map projection and Terrain-following metrics
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-10-24 (H.Yashiro)  [new] reconstruct from scale_topography
!!
!<
!-------------------------------------------------------------------------------
module scale_gridtrans
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GTRANS_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: GTRANS_MAPF (:,:,:,:) !< map factor
  real(RP), public, allocatable :: GTRANS_ROTC (:,:,:)   !< rotation coefficient

  real(RP), public, allocatable :: GTRANS_GSQRT(:,:,:,:) !< transformation metrics from Z to Xi, {G}^1/2
  real(RP), public, allocatable :: GTRANS_J13G (:,:,:,:) !< (1,3) element of Jacobian matrix * {G}^1/2
  real(RP), public, allocatable :: GTRANS_J23G (:,:,:,:) !< (2,3) element of Jacobian matrix * {G}^1/2
  real(RP), public              :: GTRANS_J33G           !< (3,3) element of Jacobian matrix * {G}^1/2

  real(RP), public, allocatable :: GTRANS_LIMYZ(:,:,:,:) !< flux limiter y-z face
  real(RP), public, allocatable :: GTRANS_LIMXZ(:,:,:,:) !< flux limiter x-z face
  real(RP), public, allocatable :: GTRANS_LIMXY(:,:,:,:) !< flux limiter x-y face

  integer,  public :: I_XYZ = 1 ! at (x,y,z)
  integer,  public :: I_XYW = 2 ! at (x,y,w)
  integer,  public :: I_UYW = 3 ! at (u,y,w)
  integer,  public :: I_XVW = 4 ! at (x,v,w)
  integer,  public :: I_UYZ = 5 ! at (u,y,z)
  integer,  public :: I_XVZ = 6 ! at (x,v,z)
  integer,  public :: I_UVZ = 7 ! at (u,v,z)

  integer,  public :: I_XY  = 1 ! at (x,y)
  integer,  public :: I_UY  = 2 ! at (u,y)
  integer,  public :: I_XV  = 3 ! at (x,v)
  integer,  public :: I_UV  = 4 ! at (u,v)

  integer,  public :: I_FYZ = 1 ! y-z face limiting x-flux
  integer,  public :: I_FXZ = 2 ! x-z face limiting y-flux
  integer,  public :: I_FXY = 3 ! x-y face limiting z-flux

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: GTRANS_mapfactor
  private :: GTRANS_terrainfollowing
  private :: GTRANS_thin_wall
  private :: GTRANS_step_mountain
  private :: GTRANS_write

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private :: GTRANS_OUT_BASENAME  = ''                     !< basename of the output file
  character(len=H_MID),  private :: GTRANS_OUT_TITLE     = 'SCALE-RM GEOMETRICS'  !< title    of the output file
  character(len=H_MID),  private :: GTRANS_OUT_DTYPE     = 'DEFAULT'              !< REAL4 or REAL8

  character(len=H_MID),  private :: GTRANS_TOPO_TYPE     = 'TERRAINFOLLOWING'     !< topographical shceme
  integer,               private :: GTRANS_ThinWall_XDIV = 50                     !< number dividing quarter-cell (x)
  integer,               private :: GTRANS_ThinWall_YDIV = 50                     !< number dividing quarter-cell (y)
  logical,               private :: debug                = .false.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine GTRANS_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_GTRANS / &
       GTRANS_OUT_BASENAME,  &
       GTRANS_OUT_DTYPE,     &
       GTRANS_TOPO_TYPE,     &
       GTRANS_ThinWall_XDIV, &
       GTRANS_ThinWall_YDIV, &
       debug

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[GRIDTRANS] / Categ[ATMOS-RM GRID] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_GTRANS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_GTRANS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_GTRANS)

    allocate( GTRANS_MAPF (IA,JA,2,4) )

    allocate( GTRANS_ROTC (IA,JA,2) )

    allocate( GTRANS_GSQRT(KA,IA,JA,7) )
    allocate( GTRANS_J13G (KA,IA,JA,7) )
    allocate( GTRANS_J23G (KA,IA,JA,7) )

    GTRANS_GSQRT(:,:,:,:) = 1.0_RP
    GTRANS_J13G (:,:,:,:) = 0.0_RP
    GTRANS_J23G (:,:,:,:) = 0.0_RP
    GTRANS_J33G = 1.0_RP

    allocate( GTRANS_LIMYZ(KA,IA,JA,7) )
    allocate( GTRANS_LIMXZ(KA,IA,JA,7) )
    allocate( GTRANS_LIMXY(KA,IA,JA,7) )
    GTRANS_LIMYZ(:,:,:,:) = 1.0_RP
    GTRANS_LIMXZ(:,:,:,:) = 1.0_RP
    GTRANS_LIMXY(:,:,:,:) = 1.0_RP

    ! calc metrics for orthogonal curvelinear coordinate
    call GTRANS_mapfactor

    ! calc coeficient for rotaion of velocity vector
    call GTRANS_rotcoef

    ! calc metrics for terrain-following,step-mountain,thin-wall coordinate
    select case(GTRANS_TOPO_TYPE)
    case ('TERRAINFOLLOWING')
      if( IO_L ) write(IO_FID_LOG,*) '=> Use terrain-following coordinate'
      call GTRANS_terrainfollowing
    case ('STEPMOUNTAIN')
      if( IO_L ) write(IO_FID_LOG,*) '=> Use step mountain method'
      call GTRANS_thin_wall
      call GTRANS_step_mountain
    case ('THINWALL')
      if( IO_L ) write(IO_FID_LOG,*) '=> Use thin-wall approximation'
      call GTRANS_thin_wall
    case default
       write(*,*) 'xxx Not appropriate name for GTRANS_TOPO_TYPE : ', trim(GTRANS_TOPO_TYPE)
       call PRC_MPIstop
    end select

    ! output metrics (for debug)
    call GTRANS_write

    return
  end subroutine GTRANS_setup

  !-----------------------------------------------------------------------------
  !> Calculate map factor
  subroutine GTRANS_mapfactor
    use scale_mapproj, only: &
       MPRJ_mapfactor
    use scale_grid_real, only: &
       REAL_calc_areavol, &
       REAL_LAT,  &
       REAL_LATX, &
       REAL_LATY, &
       REAL_LATXY
    implicit none
    !---------------------------------------------------------------------------

    call MPRJ_mapfactor( REAL_LAT  , GTRANS_MAPF(:,:,1,I_XY), GTRANS_MAPF (:,:,2,I_XY))
    call MPRJ_mapfactor( REAL_LATX , GTRANS_MAPF(:,:,1,I_UY), GTRANS_MAPF (:,:,2,I_UY))
    call MPRJ_mapfactor( REAL_LATY , GTRANS_MAPF(:,:,1,I_XV), GTRANS_MAPF (:,:,2,I_XV))
    call MPRJ_mapfactor( REAL_LATXY, GTRANS_MAPF(:,:,1,I_UV), GTRANS_MAPF (:,:,2,I_UV))

    call REAL_calc_areavol( GTRANS_MAPF(:,:,:,I_XY) )

    return
  end subroutine GTRANS_mapfactor

  !-----------------------------------------------------------------------------
  !> Calculate rotation coeffient
  subroutine GTRANS_rotcoef
    use scale_mapproj, only: &
       MPRJ_rotcoef
    use scale_grid_real, only: &
       REAL_LON,  &
       REAL_LAT
    implicit none
    !---------------------------------------------------------------------------

    call MPRJ_rotcoef( GTRANS_ROTC(:,:,:), & ! [OUT]
                       REAL_LON   (:,:),   & ! [IN]
                       REAL_LAT   (:,:)    ) ! [IN]

    return
  end subroutine GTRANS_rotcoef

  !-----------------------------------------------------------------------------
  !> Calculate G^1/2 & Jacobian
  subroutine GTRANS_terrainfollowing
    use scale_grid, only: &
       GRID_RCDZ, &
       GRID_RCDX, &
       GRID_RCDY, &
       GRID_RFDZ, &
       GRID_RFDX, &
       GRID_RFDY
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_FZ
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: REAL_CZ_U (  KA,IA,JA) !< Z coordinate [m] at (u,y,z)
    real(RP) :: REAL_CZ_V (  KA,IA,JA) !< Z coordinate [m] at (x,v,z)
    real(RP) :: REAL_CZ_UV(  KA,IA,JA) !< Z coordinate [m] at (u,y,z)
    real(RP) :: REAL_FZ_U (0:KA,IA,JA) !< Z coordinate [m] at (u,y,w)
    real(RP) :: REAL_FZ_V (0:KA,IA,JA) !< Z coordinate [m] at (x,v,w)
    real(RP) :: REAL_FZ_UV(0:KA,IA,JA) !< Z coordinate [m] at (u,v,w)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! calc Z-coordinate height at staggered position
    do j = 1, JA
    do i = 1, IA-1
    do k = 1, KA
       REAL_CZ_U(k,i,j) = 0.5_RP * ( REAL_CZ(k,i+1,j) + REAL_CZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA-1
    do k = 0, KA
       REAL_FZ_U(k,i,j) = 0.5_RP * ( REAL_FZ(k,i+1,j) + REAL_FZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
    do k = 1, KA
       REAL_CZ_V(k,i,j) = 0.5_RP * ( REAL_CZ(k,i,j+1) + REAL_CZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
    do k = 0, KA
       REAL_FZ_V(k,i,j) = 0.5_RP * ( REAL_FZ(k,i,j+1) + REAL_FZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA-1
    do k = 1, KA
       REAL_CZ_UV(k,i,j) = 0.25_RP * ( REAL_CZ(k,i+1,j+1) + REAL_CZ(k,i+1,j) &
                                    + REAL_CZ(k,i  ,j+1) + REAL_CZ(k,i  ,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA-1
    do k = 0, KA
       REAL_FZ_UV(k,i,j) = 0.25_RP * ( REAL_FZ(k,i+1,j+1) + REAL_FZ(k,i+1,j) &
                                    + REAL_FZ(k,i  ,j+1) + REAL_FZ(k,i  ,j) )
    enddo
    enddo
    enddo

    ! G^1/2
    do j = JS, JE
    do i = IS, IE
       ! at (x,y,z)
       do k = 1, KA
          GTRANS_GSQRT(k,i,j,I_XYZ) = ( REAL_FZ(k,i,j) - REAL_FZ(k-1,i,j) ) * GRID_RCDZ(k)
       enddo

       ! at (x,y,w)
       do k = 1, KA-1
          GTRANS_GSQRT(k,i,j,I_XYW) = ( REAL_CZ(k+1,i,j) - REAL_CZ(k,i,j) ) * GRID_RFDZ(k)
       enddo
       GTRANS_GSQRT(KA,i,j,I_XYW) = GTRANS_GSQRT(KA-1,i,j,I_XYW)

       ! at (u,y,w)
       do k = 1, KA-1
          GTRANS_GSQRT(k,i,j,I_UYW) = ( REAL_CZ_U(k+1,i,j) - REAL_CZ_U(k,i,j) ) * GRID_RFDZ(k)
       enddo
       GTRANS_GSQRT(KA,i,j,I_UYW) = GTRANS_GSQRT(KA-1,i,j,I_UYW)

       ! at (x,v,w)
       do k = 1, KA-1
          GTRANS_GSQRT(k,i,j,I_XVW) = ( REAL_CZ_V(k+1,i,j) - REAL_CZ_V(k,i,j) ) * GRID_RFDZ(k)
       enddo
       GTRANS_GSQRT(KA,i,j,I_XVW) = GTRANS_GSQRT(KA-1,i,j,I_XVW)

       ! at (u,y,z)
       do k = 1, KA
          GTRANS_GSQRT(k,i,j,I_UYZ) = ( REAL_FZ_U(k,i,j) - REAL_FZ_U(k-1,i,j) ) * GRID_RCDZ(k)
       enddo

       ! at (x,v,z)
       do k = 1, KA
          GTRANS_GSQRT(k,i,j,I_XVZ) = ( REAL_FZ_V(k,i,j) - REAL_FZ_V(k-1,i,j) ) * GRID_RCDZ(k)
       enddo

       ! at (u,v,z)
       do k = 1, KA
          GTRANS_GSQRT(k,i,j,I_UVZ) = ( REAL_FZ_UV(k,i,j) - REAL_FZ_UV(k-1,i,j) ) * GRID_RCDZ(k)
       enddo
    enddo
    enddo

    call COMM_vars8( GTRANS_GSQRT(:,:,:,1), 1 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,2), 2 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,3), 3 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,4), 4 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,5), 5 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,6), 6 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,7), 7 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,1), 1 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,2), 2 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,3), 3 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,4), 4 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,5), 5 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,6), 6 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,7), 7 )

    ! Jacobian * G^1/2
    do j = JS, JE
    do i = IS, IE
    do k = 1,  KA
       GTRANS_J13G(k,i,j,I_XYZ) = -( REAL_CZ_U (k,i  ,j) - REAL_CZ_U (k,i-1,j) ) * GRID_RCDX(i)
       GTRANS_J13G(k,i,j,I_XYW) = -( REAL_FZ_U (k,i  ,j) - REAL_FZ_U (k,i-1,j) ) * GRID_RCDX(i)
       GTRANS_J13G(k,i,j,I_UYW) = -( REAL_FZ   (k,i+1,j) - REAL_FZ   (k,i  ,j) ) * GRID_RFDX(i)
       GTRANS_J13G(k,i,j,I_XVW) = -( REAL_FZ_UV(k,i  ,j) - REAL_FZ_UV(k,i-1,j) ) * GRID_RCDX(i)
       GTRANS_J13G(k,i,j,I_UYZ) = -( REAL_CZ   (k,i+1,j) - REAL_CZ   (k,i  ,j) ) * GRID_RFDX(i)
       GTRANS_J13G(k,i,j,I_XVZ) = -( REAL_CZ_UV(k,i  ,j) - REAL_CZ_UV(k,i-1,j) ) * GRID_RCDX(i)
       GTRANS_J13G(k,i,j,I_UVZ) = -( REAL_CZ_V (k,i+1,j) - REAL_CZ_V (k,i  ,j) ) * GRID_RFDX(i)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 1,  KA
       GTRANS_J23G(k,i,j,I_XYZ) = -( REAL_CZ_V (k,i,j  ) - REAL_CZ_V (k,i,j-1) ) * GRID_RCDY(j)
       GTRANS_J23G(k,i,j,I_XYW) = -( REAL_FZ_V (k,i,j  ) - REAL_FZ_V (k,i,j-1) ) * GRID_RCDY(j)
       GTRANS_J23G(k,i,j,I_XVW) = -( REAL_FZ   (k,i,j+1) - REAL_FZ   (k,i,j  ) ) * GRID_RFDY(j)
       GTRANS_J23G(k,i,j,I_UYW) = -( REAL_FZ_UV(k,i,j  ) - REAL_FZ_UV(k,i,j-1) ) * GRID_RCDY(j)
       GTRANS_J23G(k,i,j,I_XVZ) = -( REAL_CZ   (k,i,j+1) - REAL_CZ   (k,i,j  ) ) * GRID_RFDY(j)
       GTRANS_J23G(k,i,j,I_UYZ) = -( REAL_CZ_UV(k,i,j  ) - REAL_CZ_UV(k,i,j-1) ) * GRID_RCDY(j)
       GTRANS_J23G(k,i,j,I_UVZ) = -( REAL_CZ_U (k,i,j+1) - REAL_CZ_U (k,i,j  ) ) * GRID_RFDY(j)
    enddo
    enddo
    enddo

    GTRANS_J33G = 1.0_RP ! - 1 / G^1/2 * G^1/2

    call COMM_vars8( GTRANS_J13G(:,:,:,I_XYZ),  8 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_XYW),  9 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_UYW), 10 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_XVW), 11 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_UYZ), 12 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_XVZ), 13 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_UVZ), 14 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_XYZ),  8 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_XYW),  9 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_UYW), 10 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_XVW), 11 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_UYZ), 12 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_XVZ), 13 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_UVZ), 14 )

    call COMM_vars8( GTRANS_J23G(:,:,:,I_XYZ), 15 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_XYW), 16 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_UYW), 17 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_XVW), 18 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_UYZ), 19 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_XVZ), 20 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_UVZ), 21 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_XYZ), 15 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_XYW), 16 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_UYW), 17 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_XVW), 18 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_UYZ), 19 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_XVZ), 20 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_UVZ), 21 )

    return
  end subroutine GTRANS_terrainfollowing

  !-----------------------------------------------------------------------------
  subroutine GTRANS_thin_wall
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only: &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY, &
       GRID_FZ, &
       GRID_FX, &
       GRID_FY
    use scale_topography, only : &
       TOPO_Zsfc
    use scale_comm, only : &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: TOPO_ZsfcALL(2*IA,2*JA)        !< doubled resolution of topography
    real(RP) :: TOPO_ZsfcXY (IA,JA)            !< absolute height at (x,y)
    real(RP) :: TOPO_ZsfcUY (IA,JA)            !< absolute height at (u,y)
    real(RP) :: TOPO_ZsfcXV (IA,JA)            !< absolute height at (x,v)
    real(RP) :: TOPO_ZsfcUV (IA,JA)            !< absolute height at (u,v)

    real(RP) :: GTRANS_QLIM (2*KA,2*IA,2*JA,3) !< quarter size flux limiter
    real(RP) :: QDZ(2*KA)                      !< length of control volume of quarter cell
    real(RP) :: QDX(2*IA)                      !< length of control volume of quarter cell
    real(RP) :: QDY(2*JA)                      !< length of control volume of quarter cell
    real(RP) :: AQAF (3)                       !< Area of the part of air on quarter-cell face
    real(RP) :: XSLOPE, YSLOPE
    real(RP) :: Ztop
    real(RP) :: DX_piece, DY_piece
    real(RP) :: DX, DY, DZ

    integer  :: I_QLIMtoLIM(3,7)               !< index when q-Flux lmiter combined

    integer  :: iii, jjj, n
    integer  :: k, i, j, kk, ii, jj
    !---------------------------------------------------------------------------

    ! calc absolute height at staggered position
    ! at (x,y)
    do j = 1, JA
    do i = 1, IA
      TOPO_ZsfcXY(i,j) = TOPO_Zsfc(i,j)
    enddo
    enddo
    ! at (u,y)
    do j = 1, JA
    do i = 1, IA-1
      TOPO_ZsfcUY(i,j) = 0.5_RP * ( TOPO_Zsfc(i,j) + TOPO_Zsfc(i+1,j) )
    enddo
    enddo
    ! at (x,v)
    do j = 1, JA-1
    do i = 1, IA
      TOPO_ZsfcXV(i,j) = 0.5_RP * ( TOPO_Zsfc(i,j) + TOPO_Zsfc(i,j+1) )
    enddo
    enddo
    ! at (u,v)
    do j = 1, JA-1
    do i = 1, IA-1
      TOPO_ZsfcUV(i,j) = 0.25_RP * ( TOPO_Zsfc(i  ,j  ) + TOPO_Zsfc(i  ,j+1) &
                                   + TOPO_Zsfc(i  ,j+1) + TOPO_Zsfc(i+1,j+1) )
    enddo
    enddo

    ! reset topography
    TOPO_Zsfc(:,:) = 0.D0

    call COMM_vars8( TOPO_ZsfcXY(:,:), 1 )
    call COMM_vars8( TOPO_ZsfcUY(:,:), 2 )
    call COMM_vars8( TOPO_ZsfcXV(:,:), 3 )
    call COMM_vars8( TOPO_ZsfcUV(:,:), 4 )
    call COMM_wait ( TOPO_ZsfcXY(:,:), 1 )
    call COMM_wait ( TOPO_ZsfcUY(:,:), 2 )
    call COMM_wait ( TOPO_ZsfcXV(:,:), 3 )
    call COMM_wait ( TOPO_ZsfcUV(:,:), 4 )

    ! all height
    do j = 1, JA
    do i = 1, IA
       ii = (i-1) * 2 + 1
       jj = (j-1) * 2 + 1

       TOPO_ZsfcALL(ii  ,jj  ) = TOPO_ZsfcXY(i,j)
       TOPO_ZsfcALL(ii+1,jj  ) = TOPO_ZsfcUY(i,j)
       TOPO_ZsfcALL(ii  ,jj+1) = TOPO_ZsfcXV(i,j)
       TOPO_ZsfcALL(ii+1,jj+1) = TOPO_ZsfcUV(i,j)
    enddo
    enddo

    ! length of control volume of quarter cell
    do k = 1, KA
       kk = (k-1) * 2 + 1

       QDZ(kk  ) = GRID_CZ(k) - GRID_FZ(k-1)
       QDZ(kk+1) = GRID_FZ(k) - GRID_CZ(k  )
    enddo

    do i = 1, IA-1
       ii = (i-1) * 2 + 1

       QDX(ii  ) = GRID_FX(i  ) - GRID_CX(i)
       QDX(ii+1) = GRID_CX(i+1) - GRID_FX(i)
    enddo

    do j = 1, JA-1
       jj = (j-1) * 2 + 1

       QDY(jj  ) = GRID_FY(j  ) - GRID_CY(j)
       QDY(jj+1) = GRID_CY(j+1) - GRID_FY(j)
    enddo

    ! quarter flux limiter
    do jj = 1, 2*(JA-1)
    do ii = 1, 2*(IA-1)
    do kk = KS, 2*KE
       DX_piece = QDX(ii) / real(GTRANS_ThinWall_XDIV,kind=RP)
       DY_piece = QDY(jj) / real(GTRANS_ThinWall_YDIV,kind=RP)

       AQAF(1:3) = 0.0_RP
       Ztop = sum(QDZ(KS:kk))

       !--- y-z face ---
       YSLOPE = ( TOPO_ZsfcALL(ii,jj+1) - TOPO_ZsfcALL(ii,jj) ) / QDY(jj)

       do jjj = 1, GTRANS_ThinWall_YDIV
          DY = ( real(jjj,kind=RP) - 0.5_RP ) * DY_piece
          DZ = Ztop - TOPO_ZsfcALL(ii,jj) - YSLOPE * DY

          if ( DZ > 0.0_RP ) then
             if ( DZ < QDZ(kk) ) then
                AQAF(I_FYZ) = AQAF(I_FYZ) + DZ      * DY_piece
             else
                AQAF(I_FYZ) = AQAF(I_FYZ) + QDZ(kk) * DY_piece
             endif
          endif
       enddo

       !--- x-z face ---
       XSLOPE = ( TOPO_ZsfcALL(ii+1,jj) - TOPO_ZsfcALL(ii,jj) ) / QDX(ii)

       do iii = 1, GTRANS_ThinWall_XDIV
          DX = ( real(iii,kind=RP) - 0.5_RP ) * DX_piece
          DZ = Ztop - TOPO_ZsfcALL(ii,jj) + XSLOPE * DX

          if ( DZ > 0.0_RP ) then
             if ( DZ < QDZ(kk) ) then
                AQAF(I_FXZ) = AQAF(I_FXZ) + DZ      * DX_piece
             else
                AQAF(I_FXZ) = AQAF(I_FXZ) + QDZ(kk) * DX_piece
             endif
          endif
       enddo

       !--- x-y face ---
       do jjj = 1, GTRANS_ThinWall_YDIV
       do iii = 1, GTRANS_ThinWall_XDIV
          DX = ( real(iii,kind=RP) - 0.5_RP ) * DX_piece
          DY = ( real(jjj,kind=RP) - 0.5_RP ) * DY_piece
          DZ = Ztop - TOPO_ZsfcALL(ii,jj) - XSLOPE * DX - YSLOPE * DY

          if ( DZ > 0.0_RP ) then
             AQAF(I_FXY) = AQAF(I_FXY) + DX_piece * DY_piece
          endif
       enddo
       enddo

       GTRANS_QLIM(kk,ii,jj,I_FYZ) = AQAF(I_FYZ) / ( QDY(jj) * QDZ(kk) )
       GTRANS_QLIM(kk,ii,jj,I_FXZ) = AQAF(I_FXZ) / ( QDX(ii) * QDZ(kk) )
       GTRANS_QLIM(kk,ii,jj,I_FXY) = AQAF(I_FXY) / ( QDY(jj) * QDX(ii) )
    enddo
    enddo
    enddo

    ! index i,j,k
    I_QLIMtoLIM(1:3,I_XYZ) = (/ 1, 1, 1 /)
    I_QLIMtoLIM(1:3,I_XYW) = (/ 1, 1, 0 /)
    I_QLIMtoLIM(1:3,I_UYW) = (/ 0, 1, 0 /)
    I_QLIMtoLIM(1:3,I_XVW) = (/ 1, 0, 0 /)
    I_QLIMtoLIM(1:3,I_UYZ) = (/ 0, 1, 1 /)
    I_QLIMtoLIM(1:3,I_XVZ) = (/ 1, 0, 1 /)
    I_QLIMtoLIM(1:3,I_UVZ) = (/ 0, 0, 1 /)

    do n = 1, 7
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ii = (i-1) * 2 + 1 - I_QLIMtoLIM(1,n)
          jj = (j-1) * 2 + 1 - I_QLIMtoLIM(2,n)
          kk = (k-1) * 2 + 1 - I_QLIMtoLIM(3,n)

         GTRANS_LIMYZ(k,i,j,n) = 0.25_RP * ( GTRANS_QLIM(kk  ,ii,jj,I_FYZ) + GTRANS_QLIM(kk  ,ii  ,jj+1,I_FYZ) &
                                           + GTRANS_QLIM(kk+1,ii,jj,I_FYZ) + GTRANS_QLIM(kk+1,ii  ,jj+1,I_FYZ) )
         GTRANS_LIMXZ(k,i,j,n) = 0.25_RP * ( GTRANS_QLIM(kk  ,ii,jj,I_FXZ) + GTRANS_QLIM(kk  ,ii+1,jj  ,I_FXZ) &
                                           + GTRANS_QLIM(kk+1,ii,jj,I_FXZ) + GTRANS_QLIM(kk+1,ii+1,jj  ,I_FXZ) )
         GTRANS_LIMXY(k,i,j,n) = 0.25_RP * ( GTRANS_QLIM(kk  ,ii,jj,I_FXY) + GTRANS_QLIM(kk  ,ii+1,jj+1,I_FXY) &
                                           + GTRANS_QLIM(kk  ,ii,jj,I_FXY) + GTRANS_QLIM(kk  ,ii+1,jj  ,I_FXY) )

         if ( GTRANS_LIMYZ(k,i,j,n) > 1.D0 ) then
            write(*,*) 'xxx Facter miss! Check!'
            write(*,*) k,i,j,n,GTRANS_LIMYZ(k,i,j,n)
            call PRC_MPIstop
         endif
       enddo
       enddo
       enddo
    enddo

    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_XYZ),  1 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_XYW),  2 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_UYW),  3 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_XVW),  4 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_UYZ),  5 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_XVZ),  6 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_UVZ),  7 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_XYZ),  1 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_XYW),  2 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_UYW),  3 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_XVW),  4 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_UYZ),  5 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_XVZ),  6 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_UVZ),  7 )

    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_XYZ),  8 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_XYW),  9 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_UYW), 10 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_XVW), 11 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_UYZ), 12 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_XVZ), 13 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_UVZ), 14 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_XYZ),  8 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_XYW),  9 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_UYW), 10 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_XVW), 11 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_UYZ), 12 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_XVZ), 13 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_UVZ), 14 )

    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_XYZ), 15 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_XYW), 16 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_UYW), 17 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_XVW), 18 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_UYZ), 19 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_XVZ), 20 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_UVZ), 21 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_XYZ), 15 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_XYW), 16 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_UYW), 17 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_XVW), 18 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_UYZ), 19 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_XVZ), 20 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_UVZ), 21 )

    return
  end subroutine GTRANS_thin_wall

  !-----------------------------------------------------------------------------
  subroutine GTRANS_step_mountain
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i,j,k,n
    !---------------------------------------------------------------------------

    do n = 1, 7
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE

         if ( GTRANS_LIMYZ(k,i,j,n) > 0.0_RP ) then
            GTRANS_LIMYZ(k,i,j,n) = 1.0_RP
         endif

         if ( GTRANS_LIMXZ(k,i,j,n) > 0.0_RP ) then
            GTRANS_LIMXZ(k,i,j,n) = 1.0_RP
         endif

         if ( GTRANS_LIMXY(k,i,j,n) > 0.0_RP ) then
            GTRANS_LIMXY(k,i,j,n) = 1.0_RP
         endif
       enddo
       enddo
       enddo
    enddo

    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_XYZ),  1 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_XYW),  2 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_UYW),  3 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_XVW),  4 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_UYZ),  5 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_XVZ),  6 )
    call COMM_vars8( GTRANS_LIMYZ(:,:,:,I_UVZ),  7 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_XYZ),  1 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_XYW),  2 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_UYW),  3 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_XVW),  4 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_UYZ),  5 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_XVZ),  6 )
    call COMM_wait ( GTRANS_LIMYZ(:,:,:,I_UVZ),  7 )

    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_XYZ),  8 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_XYW),  9 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_UYW), 10 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_XVW), 11 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_UYZ), 12 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_XVZ), 13 )
    call COMM_vars8( GTRANS_LIMXZ(:,:,:,I_UVZ), 14 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_XYZ),  8 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_XYW),  9 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_UYW), 10 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_XVW), 11 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_UYZ), 12 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_XVZ), 13 )
    call COMM_wait ( GTRANS_LIMXZ(:,:,:,I_UVZ), 14 )

    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_XYZ), 15 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_XYW), 16 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_UYW), 17 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_XVW), 18 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_UYZ), 19 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_XVZ), 20 )
    call COMM_vars8( GTRANS_LIMXY(:,:,:,I_UVZ), 21 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_XYZ), 15 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_XYW), 16 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_UYW), 17 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_XVW), 18 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_UYZ), 19 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_XVZ), 20 )
    call COMM_wait ( GTRANS_LIMXY(:,:,:,I_UVZ), 21 )

    return
  end subroutine GTRANS_step_mountain

  !-----------------------------------------------------------------------------
  !> Write metrics
  subroutine GTRANS_write
    use scale_fileio, only: &
       FILEIO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( GTRANS_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output metrics file ***'

       call FILEIO_write( GTRANS_MAPF(:,:,1,I_XY),       GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'MAPF_X_XY', 'Map factor x-dir at XY', 'NIL', 'XY', GTRANS_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( GTRANS_MAPF(:,:,2,I_XY),       GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'MAPF_Y_XY', 'Map factor y-dir at XY', 'NIL', 'XY', GTRANS_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( GTRANS_MAPF(:,:,1,I_UY),       GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'MAPF_X_UY', 'Map factor x-dir at UY', 'NIL', 'UY', GTRANS_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( GTRANS_MAPF(:,:,2,I_UY),       GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'MAPF_Y_UY', 'Map factor y-dir at UY', 'NIL', 'UY', GTRANS_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( GTRANS_MAPF(:,:,1,I_XV),       GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'MAPF_X_XV', 'Map factor x-dir at XV', 'NIL', 'XY', GTRANS_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( GTRANS_MAPF(:,:,2,I_XV),       GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'MAPF_Y_XV', 'Map factor y-dir at XV', 'NIL', 'XY', GTRANS_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( GTRANS_MAPF(:,:,1,I_UV),       GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'MAPF_X_UV', 'Map factor x-dir at UV', 'NIL', 'UY', GTRANS_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( GTRANS_MAPF(:,:,2,I_UV),       GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'MAPF_Y_UV', 'Map factor y-dir at UV', 'NIL', 'UY', GTRANS_OUT_DTYPE  ) ! [IN]

       call FILEIO_write( GTRANS_ROTC(:,:,1),            GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'ROTC_COS',  'Rotation factor (cos)',  'NIL', 'XY', GTRANS_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( GTRANS_ROTC(:,:,2),            GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'ROTC_SIN',  'Rotation factor (sin)',  'NIL', 'XY', GTRANS_OUT_DTYPE  ) ! [IN]

       call FILEIO_write( GTRANS_ROTC(:,:,1),            GTRANS_OUT_BASENAME, GTRANS_OUT_TITLE, & ! [IN]
                          'ROTC_COS',  'Rotation factor (cos)',  'NIL', 'XY', GTRANS_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine GTRANS_write

end module scale_gridtrans
