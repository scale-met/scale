!-------------------------------------------------------------------------------
!> module Atmosphere Grid CartesianC metirc
!!
!! @par Description
!!          Map projection and Terrain-following metrics for the CaresianC grid
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_grid_cartesC_metric
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_GRID_CARTESC_METRIC_setup
  public :: ATMOS_GRID_CARTESC_METRIC_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_METRIC_MAPF (:,:,:,:) !< map factor
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_METRIC_ROTC (:,:,:)   !< rotation coefficient

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_METRIC_GSQRT(:,:,:,:) !< transformation metrics from Z to Xi, {G}^1/2
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_METRIC_J13G (:,:,:,:) !< (1,3) element of Jacobian matrix * {G}^1/2
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_METRIC_J23G (:,:,:,:) !< (2,3) element of Jacobian matrix * {G}^1/2
  real(RP), public              :: ATMOS_GRID_CARTESC_METRIC_J33G           !< (3,3) element of Jacobian matrix * {G}^1/2

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,:) !< flux limiter y-z face
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,:) !< flux limiter x-z face
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,:) !< flux limiter x-y face

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_GRID_CARTESC_METRIC_mapfactor
  private :: ATMOS_GRID_CARTESC_METRIC_terrainfollowing
  private :: ATMOS_GRID_CARTESC_METRIC_thin_wall
  private :: ATMOS_GRID_CARTESC_METRIC_step_mountain
  private :: ATMOS_GRID_CARTESC_METRIC_write

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME  = ''                     !< basename of the output file
  character(len=H_MID),   private :: ATMOS_GRID_CARTESC_METRIC_OUT_TITLE     = 'SCALE-RM GEOMETRICS'  !< title    of the output file
  character(len=H_SHORT), private :: ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE     = 'DEFAULT'              !< REAL4 or REAL8

  character(len=H_SHORT), private :: ATMOS_GRID_CARTESC_METRIC_TOPO_type     = 'TERRAINFOLLOWING'     !< topographical shceme
  integer,                private :: ATMOS_GRID_CARTESC_METRIC_ThinWall_XDIV = 50                     !< number dividing quarter-cell (x)
  integer,                private :: ATMOS_GRID_CARTESC_METRIC_ThinWall_YDIV = 50                     !< number dividing quarter-cell (y)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_GRID_CARTESC_METRIC_setup
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    namelist / PARAM_ATMOS_GRID_CARTESC_METRIC / &
       ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME,  &
       ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE,     &
       ATMOS_GRID_CARTESC_METRIC_TOPO_type,     &
       ATMOS_GRID_CARTESC_METRIC_ThinWall_XDIV, &
       ATMOS_GRID_CARTESC_METRIC_ThinWall_YDIV

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_METRIC_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_GRID_CARTESC_METRIC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_GRID_CARTESC_METRIC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_GRID_CARTESC_METRIC_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_GRID_CARTESC_METRIC. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_GRID_CARTESC_METRIC)

    allocate( ATMOS_GRID_CARTESC_METRIC_MAPF (IA,JA,2,4) )

    allocate( ATMOS_GRID_CARTESC_METRIC_ROTC (IA,JA,2) )

    if ( PRC_TwoD ) then
       allocate( ATMOS_GRID_CARTESC_METRIC_GSQRT(KA,IA,JA,4) )
       allocate( ATMOS_GRID_CARTESC_METRIC_J13G (KA,IA,JA,4) )
       allocate( ATMOS_GRID_CARTESC_METRIC_J23G (KA,IA,JA,4) )
    else
       allocate( ATMOS_GRID_CARTESC_METRIC_GSQRT(KA,IA,JA,7) )
       allocate( ATMOS_GRID_CARTESC_METRIC_J13G (KA,IA,JA,7) )
       allocate( ATMOS_GRID_CARTESC_METRIC_J23G (KA,IA,JA,7) )
    end if

    ATMOS_GRID_CARTESC_METRIC_GSQRT(:,:,:,:) = 1.0_RP
    ATMOS_GRID_CARTESC_METRIC_J13G (:,:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_METRIC_J23G (:,:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_METRIC_J33G = 1.0_RP

    allocate( ATMOS_GRID_CARTESC_METRIC_LIMYZ(KA,IA,JA,7) )
    allocate( ATMOS_GRID_CARTESC_METRIC_LIMXZ(KA,IA,JA,7) )
    allocate( ATMOS_GRID_CARTESC_METRIC_LIMXY(KA,IA,JA,7) )
    ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,:) = 1.0_RP
    ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,:) = 1.0_RP
    ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,:) = 1.0_RP

    ! calc metrics for orthogonal curvelinear coordinate
    call ATMOS_GRID_CARTESC_METRIC_mapfactor

    ! calc coeficient for rotaion of velocity vector
    call ATMOS_GRID_CARTESC_METRIC_rotcoef

    ! calc metrics for terrain-following,step-mountain,thin-wall coordinate
    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_METRIC_setup",*) 'Terrain coordinate type : ', trim(ATMOS_GRID_CARTESC_METRIC_TOPO_type)
    select case(ATMOS_GRID_CARTESC_METRIC_TOPO_type)
    case('TERRAINFOLLOWING')
      LOG_INFO_CONT(*) '=> Terrain-following method'
      call ATMOS_GRID_CARTESC_METRIC_terrainfollowing
    case('STEPMOUNTAIN')
      LOG_INFO_CONT(*) '=> Step-mountain method'
      call ATMOS_GRID_CARTESC_METRIC_thin_wall
      call ATMOS_GRID_CARTESC_METRIC_step_mountain
    case('THINWALL')
      LOG_INFO_CONT(*) '=> Thin-wall approximation method'
      call ATMOS_GRID_CARTESC_METRIC_thin_wall
    case default
       LOG_ERROR("ATMOS_GRID_CARTESC_METRIC_setup",*) 'Unsupported ATMOS_GRID_CARTESC_METRIC_TOPO_type. STOP'
       call PRC_abort
    end select

    ! output metrics (for debug)
    call ATMOS_GRID_CARTESC_METRIC_write

    return
  end subroutine ATMOS_GRID_CARTESC_METRIC_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_GRID_CARTESC_METRIC_finalize
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_METRIC_finalize",*) 'Finalize'

    deallocate( ATMOS_GRID_CARTESC_METRIC_MAPF )

    deallocate( ATMOS_GRID_CARTESC_METRIC_ROTC )

    if ( PRC_TwoD ) then
       deallocate( ATMOS_GRID_CARTESC_METRIC_GSQRT )
       deallocate( ATMOS_GRID_CARTESC_METRIC_J13G  )
       deallocate( ATMOS_GRID_CARTESC_METRIC_J23G  )
    else
       deallocate( ATMOS_GRID_CARTESC_METRIC_GSQRT )
       deallocate( ATMOS_GRID_CARTESC_METRIC_J13G  )
       deallocate( ATMOS_GRID_CARTESC_METRIC_J23G  )
    end if

    deallocate( ATMOS_GRID_CARTESC_METRIC_LIMYZ )
    deallocate( ATMOS_GRID_CARTESC_METRIC_LIMXZ )
    deallocate( ATMOS_GRID_CARTESC_METRIC_LIMXY )

    return
  end subroutine ATMOS_GRID_CARTESC_METRIC_finalize

  !-----------------------------------------------------------------------------
  !> Calculate map factor
  subroutine ATMOS_GRID_CARTESC_METRIC_mapfactor
    use scale_prc_cartesC, only: &
       PRC_TwoD
    use scale_mapprojection, only: &
       MAPPROJECTION_mapfactor
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_RCDX, &
       ATMOS_GRID_CARTESC_RCDY
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_LAT,          &
       ATMOS_GRID_CARTESC_REAL_LATUY,        &
       ATMOS_GRID_CARTESC_REAL_LATXV,        &
       ATMOS_GRID_CARTESC_REAL_LATUV
    use scale_topography, only: &
       TOPOGRAPHY_calc_tan_slope
    implicit none
    !---------------------------------------------------------------------------

    call MAPPROJECTION_mapfactor( IA, 1, IA, JA, 1, JA,   &
         ATMOS_GRID_CARTESC_REAL_LAT  ( :, :), ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,1,I_XY), ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,2,I_XY))
    call MAPPROJECTION_mapfactor( IA, 1, IA, JA, 1, JA, &
         ATMOS_GRID_CARTESC_REAL_LATXV( :,1:), ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,1,I_XV), ATMOS_GRID_CARTESC_METRIC_MAPF (:,:,2,I_XV))
    if ( .not. PRC_TwoD ) then
    call MAPPROJECTION_mapfactor( IA, 1, IA, JA, 1, JA,   &
         ATMOS_GRID_CARTESC_REAL_LATUY(1:, :), ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,1,I_UY), ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,2,I_UY))
    call MAPPROJECTION_mapfactor( IA, 1, IA, JA, 1, JA, &
         ATMOS_GRID_CARTESC_REAL_LATUV(1:,1:), ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,1,I_UV), ATMOS_GRID_CARTESC_METRIC_MAPF (:,:,2,I_UV))
    end if

    call TOPOGRAPHY_calc_tan_slope( IA, IS, IE, JA, JS, JE, &
         ATMOS_GRID_CARTESC_RCDX(:), ATMOS_GRID_CARTESC_RCDY(:), &
         ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,:,I_XY) )

    return
  end subroutine ATMOS_GRID_CARTESC_METRIC_mapfactor

  !-----------------------------------------------------------------------------
  !> Calculate rotation coeffient
  subroutine ATMOS_GRID_CARTESC_METRIC_rotcoef
    use scale_mapprojection, only: &
       MAPPROJECTION_rotcoef
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_LON,  &
       ATMOS_GRID_CARTESC_REAL_LAT
    implicit none
    !---------------------------------------------------------------------------

    call MAPPROJECTION_rotcoef( IA, 1, IA, JA, 1, JA, &
         ATMOS_GRID_CARTESC_REAL_LON   (:,:),   & ! [IN]
         ATMOS_GRID_CARTESC_REAL_LAT   (:,:),   & ! [IN]
         ATMOS_GRID_CARTESC_METRIC_ROTC(:,:,1), & ! [OUT]
         ATMOS_GRID_CARTESC_METRIC_ROTC(:,:,2)  ) ! [OUT]

    return
  end subroutine ATMOS_GRID_CARTESC_METRIC_rotcoef

  !-----------------------------------------------------------------------------
  !> Calculate G^1/2 & Jacobian
  subroutine ATMOS_GRID_CARTESC_METRIC_terrainfollowing
    use scale_prc_cartesC, only: &
       PRC_TwoD
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_RCDZ, &
       ATMOS_GRID_CARTESC_RCDX, &
       ATMOS_GRID_CARTESC_RCDY, &
       ATMOS_GRID_CARTESC_RFDZ, &
       ATMOS_GRID_CARTESC_RFDX, &
       ATMOS_GRID_CARTESC_RFDY
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_CZ,   &
       ATMOS_GRID_CARTESC_REAL_CZUY, &
       ATMOS_GRID_CARTESC_REAL_CZXV, &
       ATMOS_GRID_CARTESC_REAL_CZUV, &
       ATMOS_GRID_CARTESC_REAL_FZ,   &
       ATMOS_GRID_CARTESC_REAL_FZUY, &
       ATMOS_GRID_CARTESC_REAL_FZXV, &
       ATMOS_GRID_CARTESC_REAL_FZUV
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: k, i, j
    integer :: i_start, i_end
    !---------------------------------------------------------------------------

    ! G^1/2
    !$omp parallel do
    do j = 1, JA
    do i = 1, IA
       ! at (x,y,z)
       do k = 1, KA
          ATMOS_GRID_CARTESC_METRIC_GSQRT(k,i,j,I_XYZ) = ( ATMOS_GRID_CARTESC_REAL_FZ(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZ(k-1,i,j) ) * ATMOS_GRID_CARTESC_RCDZ(k)
       enddo

       ! at (x,y,w)
       do k = 1, KA-1
          ATMOS_GRID_CARTESC_METRIC_GSQRT(k,i,j,I_XYW) = ( ATMOS_GRID_CARTESC_REAL_CZ(k+1,i,j) - ATMOS_GRID_CARTESC_REAL_CZ(k,i,j) ) * ATMOS_GRID_CARTESC_RFDZ(k)
       enddo
       ATMOS_GRID_CARTESC_METRIC_GSQRT(KA,i,j,I_XYW) = ATMOS_GRID_CARTESC_METRIC_GSQRT(KA-1,i,j,I_XYW)

       ! at (x,v,w)
       do k = 1, KA-1
          ATMOS_GRID_CARTESC_METRIC_GSQRT(k,i,j,I_XVW) = ( ATMOS_GRID_CARTESC_REAL_CZXV(k+1,i,j) - ATMOS_GRID_CARTESC_REAL_CZXV(k,i,j) ) * ATMOS_GRID_CARTESC_RFDZ(k)
       enddo
       ATMOS_GRID_CARTESC_METRIC_GSQRT(KA,i,j,I_XVW) = ATMOS_GRID_CARTESC_METRIC_GSQRT(KA-1,i,j,I_XVW)

       ! at (x,v,z)
       do k = 1, KA
          ATMOS_GRID_CARTESC_METRIC_GSQRT(k,i,j,I_XVZ) = ( ATMOS_GRID_CARTESC_REAL_FZXV(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZXV(k-1,i,j) ) * ATMOS_GRID_CARTESC_RCDZ(k)
       enddo
    enddo
    enddo
    if ( .not. PRC_TwoD ) then
    !$omp parallel do
    do j = 1, JA
    do i = 1, IA
       ! at (u,y,w)
       do k = 1, KA-1
          ATMOS_GRID_CARTESC_METRIC_GSQRT(k,i,j,I_UYW) = ( ATMOS_GRID_CARTESC_REAL_CZUY(k+1,i,j) - ATMOS_GRID_CARTESC_REAL_CZUY(k,i,j) ) * ATMOS_GRID_CARTESC_RFDZ(k)
       enddo
       ATMOS_GRID_CARTESC_METRIC_GSQRT(KA,i,j,I_UYW) = ATMOS_GRID_CARTESC_METRIC_GSQRT(KA-1,i,j,I_UYW)

       ! at (u,y,z)
       do k = 1, KA
          ATMOS_GRID_CARTESC_METRIC_GSQRT(k,i,j,I_UYZ) = ( ATMOS_GRID_CARTESC_REAL_FZUY(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZUY(k-1,i,j) ) * ATMOS_GRID_CARTESC_RCDZ(k)
       enddo

       ! at (u,v,z)
       do k = 1, KA
          ATMOS_GRID_CARTESC_METRIC_GSQRT(k,i,j,I_UVZ) = ( ATMOS_GRID_CARTESC_REAL_FZUV(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZUV(k-1,i,j) ) * ATMOS_GRID_CARTESC_RCDZ(k)
       enddo
    enddo
    enddo
    end if

    ! Jacobian * G^1/2
    if ( .not. PRC_TwoD ) then
    !$omp parallel do
    do j = 1, JA
    do i = 2, IA
    do k = 1, KA
       ATMOS_GRID_CARTESC_METRIC_J13G(k,i,j,I_XYZ) = -( ATMOS_GRID_CARTESC_REAL_CZUY(k,i,j) - ATMOS_GRID_CARTESC_REAL_CZUY(k,i-1,j) ) * ATMOS_GRID_CARTESC_RCDX(i)
       ATMOS_GRID_CARTESC_METRIC_J13G(k,i,j,I_XYW) = -( ATMOS_GRID_CARTESC_REAL_FZUY(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZUY(k,i-1,j) ) * ATMOS_GRID_CARTESC_RCDX(i)
       ATMOS_GRID_CARTESC_METRIC_J13G(k,i,j,I_XVW) = -( ATMOS_GRID_CARTESC_REAL_FZUV(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZUV(k,i-1,j) ) * ATMOS_GRID_CARTESC_RCDX(i)
       ATMOS_GRID_CARTESC_METRIC_J13G(k,i,j,I_XVZ) = -( ATMOS_GRID_CARTESC_REAL_CZUV(k,i,j) - ATMOS_GRID_CARTESC_REAL_CZUV(k,i-1,j) ) * ATMOS_GRID_CARTESC_RCDX(i)
    enddo
    enddo
    enddo
    !$omp parallel do
    do j = 1, JA
    do i = 1, IA-1
    do k = 1, KA
       ATMOS_GRID_CARTESC_METRIC_J13G(k,i,j,I_UYW) = -( ATMOS_GRID_CARTESC_REAL_FZ  (k,i+1,j) - ATMOS_GRID_CARTESC_REAL_FZ  (k,i,j) ) * ATMOS_GRID_CARTESC_RFDX(i)
       ATMOS_GRID_CARTESC_METRIC_J13G(k,i,j,I_UYZ) = -( ATMOS_GRID_CARTESC_REAL_CZ  (k,i+1,j) - ATMOS_GRID_CARTESC_REAL_CZ  (k,i,j) ) * ATMOS_GRID_CARTESC_RFDX(i)
       ATMOS_GRID_CARTESC_METRIC_J13G(k,i,j,I_UVZ) = -( ATMOS_GRID_CARTESC_REAL_CZXV(k,i+1,j) - ATMOS_GRID_CARTESC_REAL_CZXV(k,i,j) ) * ATMOS_GRID_CARTESC_RFDX(i)
    enddo
    enddo
    enddo
    end if

    !$omp parallel do
    do j = 2, JA-1
    do i = 1, IA
    do k = 1, KA
       ATMOS_GRID_CARTESC_METRIC_J23G(k,i,j,I_XYZ) = -( ATMOS_GRID_CARTESC_REAL_CZXV(k,i,j  ) - ATMOS_GRID_CARTESC_REAL_CZXV(k,i,j-1) ) * ATMOS_GRID_CARTESC_RCDY(j)
       ATMOS_GRID_CARTESC_METRIC_J23G(k,i,j,I_XYW) = -( ATMOS_GRID_CARTESC_REAL_FZXV(k,i,j  ) - ATMOS_GRID_CARTESC_REAL_FZXV(k,i,j-1) ) * ATMOS_GRID_CARTESC_RCDY(j)
       ATMOS_GRID_CARTESC_METRIC_J23G(k,i,j,I_XVW) = -( ATMOS_GRID_CARTESC_REAL_FZ  (k,i,j+1) - ATMOS_GRID_CARTESC_REAL_FZ  (k,i,j  ) ) * ATMOS_GRID_CARTESC_RFDY(j)
       ATMOS_GRID_CARTESC_METRIC_J23G(k,i,j,I_XVZ) = -( ATMOS_GRID_CARTESC_REAL_CZ  (k,i,j+1) - ATMOS_GRID_CARTESC_REAL_CZ  (k,i,j  ) ) * ATMOS_GRID_CARTESC_RFDY(j)
    enddo
    enddo
    enddo
    if ( .not. PRC_TwoD ) then
    !$omp parallel do
    do j = 2, JA-1
    do i = 1, IA
    do k = 1, KA
       ATMOS_GRID_CARTESC_METRIC_J23G(k,i,j,I_UYW) = -( ATMOS_GRID_CARTESC_REAL_FZUV(k,i,j  ) - ATMOS_GRID_CARTESC_REAL_FZUV(k,i,j-1) ) * ATMOS_GRID_CARTESC_RCDY(j)
       ATMOS_GRID_CARTESC_METRIC_J23G(k,i,j,I_UYZ) = -( ATMOS_GRID_CARTESC_REAL_CZUV(k,i,j  ) - ATMOS_GRID_CARTESC_REAL_CZUV(k,i,j-1) ) * ATMOS_GRID_CARTESC_RCDY(j)
       ATMOS_GRID_CARTESC_METRIC_J23G(k,i,j,I_UVZ) = -( ATMOS_GRID_CARTESC_REAL_CZUY(k,i,j+1) - ATMOS_GRID_CARTESC_REAL_CZUY(k,i,j  ) ) * ATMOS_GRID_CARTESC_RFDY(j)
    enddo
    enddo
    enddo
    end if

    ATMOS_GRID_CARTESC_METRIC_J33G = 1.0_RP ! - 1 / G^1/2 * G^1/2

    if ( .not. PRC_TwoD ) then
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XYZ),  8 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XYW),  9 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XVW), 10 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XVZ), 11 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_UYW), 12 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_UYZ), 13 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_UVZ), 14 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XYZ),  8, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XYW),  9, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XVW), 10, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XVZ), 11, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_UYW), 12, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_UYZ), 13, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_UVZ), 14, .false. )
    end if

    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XYZ), 15 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XYW), 16 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XVW), 17 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XVZ), 18 )
    if ( .not. PRC_TwoD ) then
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_UYW), 19 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_UYZ), 20 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_UVZ), 21 )
    end if
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XYZ), 15, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XYW), 16, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XVW), 17, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XVZ), 18, .false. )
    if ( .not. PRC_TwoD ) then
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_UYW), 19, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_UYZ), 20, .false. )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_UVZ), 21, .false. )
    endif

    return
  end subroutine ATMOS_GRID_CARTESC_METRIC_terrainfollowing

  !-----------------------------------------------------------------------------
  subroutine ATMOS_GRID_CARTESC_METRIC_thin_wall
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CZ, &
       ATMOS_GRID_CARTESC_CX, &
       ATMOS_GRID_CARTESC_CY, &
       ATMOS_GRID_CARTESC_FZ, &
       ATMOS_GRID_CARTESC_FX, &
       ATMOS_GRID_CARTESC_FY
    use scale_topography, only : &
       TOPOGRAPHY_Zsfc
    use scale_comm_cartesC, only : &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: TOPO_ZsfcALL(2*IA,2*JA)        !< doubled resolution of topography
    real(RP) :: TOPO_ZsfcXY (IA,JA)            !< absolute height at (x,y)
    real(RP) :: TOPO_ZsfcUY (IA,JA)            !< absolute height at (u,y)
    real(RP) :: TOPO_ZsfcXV (IA,JA)            !< absolute height at (x,v)
    real(RP) :: TOPO_ZsfcUV (IA,JA)            !< absolute height at (u,v)

    real(RP) :: ATMOS_GRID_CARTESC_METRIC_QLIM (2*KA,2*IA,2*JA,3) !< quarter size flux limiter
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
      TOPO_ZsfcXY(i,j) = TOPOGRAPHY_Zsfc(i,j)
    enddo
    enddo
    ! at (u,y)
    do j = 1, JA
    do i = 1, IA-1
      TOPO_ZsfcUY(i,j) = 0.5_RP * ( TOPOGRAPHY_Zsfc(i,j) + TOPOGRAPHY_Zsfc(i+1,j) )
    enddo
    enddo
    ! at (x,v)
    do j = 1, JA-1
    do i = 1, IA
      TOPO_ZsfcXV(i,j) = 0.5_RP * ( TOPOGRAPHY_Zsfc(i,j) + TOPOGRAPHY_Zsfc(i,j+1) )
    enddo
    enddo
    ! at (u,v)
    do j = 1, JA-1
    do i = 1, IA-1
      TOPO_ZsfcUV(i,j) = 0.25_RP * ( TOPOGRAPHY_Zsfc(i  ,j  ) + TOPOGRAPHY_Zsfc(i  ,j+1) &
                                   + TOPOGRAPHY_Zsfc(i  ,j+1) + TOPOGRAPHY_Zsfc(i+1,j+1) )
    enddo
    enddo

    ! reset topography
    TOPOGRAPHY_Zsfc(:,:) = 0.D0

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

       QDZ(kk  ) = ATMOS_GRID_CARTESC_CZ(k) - ATMOS_GRID_CARTESC_FZ(k-1)
       QDZ(kk+1) = ATMOS_GRID_CARTESC_FZ(k) - ATMOS_GRID_CARTESC_CZ(k  )
    enddo

    do i = 1, IA-1
       ii = (i-1) * 2 + 1

       QDX(ii  ) = ATMOS_GRID_CARTESC_FX(i  ) - ATMOS_GRID_CARTESC_CX(i)
       QDX(ii+1) = ATMOS_GRID_CARTESC_CX(i+1) - ATMOS_GRID_CARTESC_FX(i)
    enddo

    do j = 1, JA-1
       jj = (j-1) * 2 + 1

       QDY(jj  ) = ATMOS_GRID_CARTESC_FY(j  ) - ATMOS_GRID_CARTESC_CY(j)
       QDY(jj+1) = ATMOS_GRID_CARTESC_CY(j+1) - ATMOS_GRID_CARTESC_FY(j)
    enddo

    ! quarter flux limiter
    do jj = 1, 2*(JA-1)
    do ii = 1, 2*(IA-1)
    do kk = KS, 2*KE
       DX_piece = QDX(ii) / real(ATMOS_GRID_CARTESC_METRIC_ThinWall_XDIV,kind=RP)
       DY_piece = QDY(jj) / real(ATMOS_GRID_CARTESC_METRIC_ThinWall_YDIV,kind=RP)

       AQAF(1:3) = 0.0_RP
       Ztop = sum(QDZ(KS:kk))

       !--- y-z face ---
       YSLOPE = ( TOPO_ZsfcALL(ii,jj+1) - TOPO_ZsfcALL(ii,jj) ) / QDY(jj)

       do jjj = 1, ATMOS_GRID_CARTESC_METRIC_ThinWall_YDIV
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

       do iii = 1, ATMOS_GRID_CARTESC_METRIC_ThinWall_XDIV
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
       do jjj = 1, ATMOS_GRID_CARTESC_METRIC_ThinWall_YDIV
       do iii = 1, ATMOS_GRID_CARTESC_METRIC_ThinWall_XDIV
          DX = ( real(iii,kind=RP) - 0.5_RP ) * DX_piece
          DY = ( real(jjj,kind=RP) - 0.5_RP ) * DY_piece
          DZ = Ztop - TOPO_ZsfcALL(ii,jj) - XSLOPE * DX - YSLOPE * DY

          if ( DZ > 0.0_RP ) then
             AQAF(I_FXY) = AQAF(I_FXY) + DX_piece * DY_piece
          endif
       enddo
       enddo

       ATMOS_GRID_CARTESC_METRIC_QLIM(kk,ii,jj,I_FYZ) = AQAF(I_FYZ) / ( QDY(jj) * QDZ(kk) )
       ATMOS_GRID_CARTESC_METRIC_QLIM(kk,ii,jj,I_FXZ) = AQAF(I_FXZ) / ( QDX(ii) * QDZ(kk) )
       ATMOS_GRID_CARTESC_METRIC_QLIM(kk,ii,jj,I_FXY) = AQAF(I_FXY) / ( QDY(jj) * QDX(ii) )
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

         ATMOS_GRID_CARTESC_METRIC_LIMYZ(k,i,j,n) = 0.25_RP * ( ATMOS_GRID_CARTESC_METRIC_QLIM(kk  ,ii,jj,I_FYZ) + ATMOS_GRID_CARTESC_METRIC_QLIM(kk  ,ii  ,jj+1,I_FYZ) &
                                           + ATMOS_GRID_CARTESC_METRIC_QLIM(kk+1,ii,jj,I_FYZ) + ATMOS_GRID_CARTESC_METRIC_QLIM(kk+1,ii  ,jj+1,I_FYZ) )
         ATMOS_GRID_CARTESC_METRIC_LIMXZ(k,i,j,n) = 0.25_RP * ( ATMOS_GRID_CARTESC_METRIC_QLIM(kk  ,ii,jj,I_FXZ) + ATMOS_GRID_CARTESC_METRIC_QLIM(kk  ,ii+1,jj  ,I_FXZ) &
                                           + ATMOS_GRID_CARTESC_METRIC_QLIM(kk+1,ii,jj,I_FXZ) + ATMOS_GRID_CARTESC_METRIC_QLIM(kk+1,ii+1,jj  ,I_FXZ) )
         ATMOS_GRID_CARTESC_METRIC_LIMXY(k,i,j,n) = 0.25_RP * ( ATMOS_GRID_CARTESC_METRIC_QLIM(kk  ,ii,jj,I_FXY) + ATMOS_GRID_CARTESC_METRIC_QLIM(kk  ,ii+1,jj+1,I_FXY) &
                                           + ATMOS_GRID_CARTESC_METRIC_QLIM(kk  ,ii,jj,I_FXY) + ATMOS_GRID_CARTESC_METRIC_QLIM(kk  ,ii+1,jj  ,I_FXY) )

         if ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(k,i,j,n) > 1.D0 ) then
            LOG_ERROR("ATMOS_GRID_CARTESC_METRIC_thin_wall",*) 'Facter miss! Check!'
            LOG_ERROR_CONT(*) k,i,j,n,ATMOS_GRID_CARTESC_METRIC_LIMYZ(k,i,j,n)
            call PRC_abort
         endif
       enddo
       enddo
       enddo
    enddo

    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XYZ),  1 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XYW),  2 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UYW),  3 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XVW),  4 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UYZ),  5 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XVZ),  6 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UVZ),  7 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XYZ),  1 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XYW),  2 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UYW),  3 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XVW),  4 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UYZ),  5 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XVZ),  6 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UVZ),  7 )

    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XYZ),  8 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XYW),  9 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UYW), 10 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XVW), 11 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UYZ), 12 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XVZ), 13 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UVZ), 14 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XYZ),  8 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XYW),  9 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UYW), 10 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XVW), 11 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UYZ), 12 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XVZ), 13 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UVZ), 14 )

    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XYZ), 15 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XYW), 16 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UYW), 17 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XVW), 18 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UYZ), 19 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XVZ), 20 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UVZ), 21 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XYZ), 15 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XYW), 16 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UYW), 17 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XVW), 18 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UYZ), 19 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XVZ), 20 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UVZ), 21 )

    return
  end subroutine ATMOS_GRID_CARTESC_METRIC_thin_wall

  !-----------------------------------------------------------------------------
  subroutine ATMOS_GRID_CARTESC_METRIC_step_mountain
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i,j,k,n
    !---------------------------------------------------------------------------

    do n = 1, 7
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE

         if ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(k,i,j,n) > 0.0_RP ) then
            ATMOS_GRID_CARTESC_METRIC_LIMYZ(k,i,j,n) = 1.0_RP
         endif

         if ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(k,i,j,n) > 0.0_RP ) then
            ATMOS_GRID_CARTESC_METRIC_LIMXZ(k,i,j,n) = 1.0_RP
         endif

         if ( ATMOS_GRID_CARTESC_METRIC_LIMXY(k,i,j,n) > 0.0_RP ) then
            ATMOS_GRID_CARTESC_METRIC_LIMXY(k,i,j,n) = 1.0_RP
         endif
       enddo
       enddo
       enddo
    enddo

    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XYZ),  1 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XYW),  2 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UYW),  3 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XVW),  4 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UYZ),  5 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XVZ),  6 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UVZ),  7 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XYZ),  1 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XYW),  2 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UYW),  3 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XVW),  4 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UYZ),  5 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_XVZ),  6 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMYZ(:,:,:,I_UVZ),  7 )

    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XYZ),  8 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XYW),  9 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UYW), 10 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XVW), 11 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UYZ), 12 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XVZ), 13 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UVZ), 14 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XYZ),  8 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XYW),  9 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UYW), 10 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XVW), 11 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UYZ), 12 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_XVZ), 13 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXZ(:,:,:,I_UVZ), 14 )

    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XYZ), 15 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XYW), 16 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UYW), 17 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XVW), 18 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UYZ), 19 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XVZ), 20 )
    call COMM_vars8( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UVZ), 21 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XYZ), 15 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XYW), 16 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UYW), 17 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XVW), 18 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UYZ), 19 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_XVZ), 20 )
    call COMM_wait ( ATMOS_GRID_CARTESC_METRIC_LIMXY(:,:,:,I_UVZ), 21 )

    return
  end subroutine ATMOS_GRID_CARTESC_METRIC_step_mountain

  !-----------------------------------------------------------------------------
  !> Write metrics
  subroutine ATMOS_GRID_CARTESC_METRIC_write
    use scale_prc_cartesC, only: &
       PRC_TwoD
    use scale_const, only: &
       CONST_RADIUS
    use scale_vector, only: &
       VECTR_distance
    use scale_file_cartesC, only: &
       FILE_CARTESC_write
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON, &
       ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT, &
       ATMOS_GRID_CARTESC_REAL_LON,           &
       ATMOS_GRID_CARTESC_REAL_LAT
    use scale_mapprojection, only: &
       MAPPROJECTION_lonlat2xy
    implicit none

    real(RP) :: check_X_XY(IA,JA)
    real(RP) :: check_Y_XY(IA,JA)
    real(RP) :: distance  (IA,JA)

    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("ATMOS_GRID_CARTESC_METRIC_write",*) 'Output metrics file '

       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,1,I_XY),       ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'MAPF_X_XY', 'Map factor x-dir at XY', 'NIL', 'XY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,2,I_XY),       ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'MAPF_Y_XY', 'Map factor y-dir at XY', 'NIL', 'XY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,1,I_UY),       ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'MAPF_X_UY', 'Map factor x-dir at UY', 'NIL', 'UY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,2,I_UY),       ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'MAPF_Y_UY', 'Map factor y-dir at UY', 'NIL', 'UY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,1,I_XV),       ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'MAPF_X_XV', 'Map factor x-dir at XV', 'NIL', 'XV', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,2,I_XV),       ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'MAPF_Y_XV', 'Map factor y-dir at XV', 'NIL', 'XV', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,1,I_UV),       ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'MAPF_X_UV', 'Map factor x-dir at UV', 'NIL', 'UV', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,2,I_UV),       ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'MAPF_Y_UV', 'Map factor y-dir at UV', 'NIL', 'UV', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]

       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_ROTC(:,:,1),            ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'ROTC_COS',  'Rotation factor (cos)',  'NIL', 'XY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_ROTC(:,:,2),            ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'ROTC_SIN',  'Rotation factor (sin)',  'NIL', 'XY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]

       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_ROTC(:,:,1),            ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'ROTC_COS',  'Rotation factor (cos)',  'NIL', 'XY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]

       do j = 1, JA
       do i = 1, IA
          call MAPPROJECTION_lonlat2xy( ATMOS_GRID_CARTESC_REAL_LON(i,j), ATMOS_GRID_CARTESC_REAL_LAT(i,j), check_X_XY(i,j), check_Y_XY(i,j) )
       enddo
       enddo

       call FILE_CARTESC_write( check_X_XY(:,:),     ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'X_XY', 'x at XY for check', 'NIL', 'XY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( check_Y_XY(:,:),     ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'Y_XY', 'y at XY for check', 'NIL', 'XY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]

       do j = 1, JA
       do i = 1, IA
          call VECTR_distance( CONST_RADIUS,       & ! [IN]
                               ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON, & ! [IN]
                               ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT, & ! [IN]
                               ATMOS_GRID_CARTESC_REAL_LON(i,j),      & ! [IN]
                               ATMOS_GRID_CARTESC_REAL_LAT(i,j),      & ! [IN]
                               distance(i,j)       ) ! [OUT]
       enddo
       enddo

       call FILE_CARTESC_write( distance(:,:),               ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                          'distance', 'distance from basepoint', 'm', 'XY', ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]


       ! metrics for terrain-following coordinate
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_GSQRT(:,:,:,I_XYZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'GSQRT_ZXY',  'transformation metrics from Z to Xi, G^1/2 at ZXY', '1', 'ZXY',                                             & ! [IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_GSQRT(:,:,:,I_XYW), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'GSQRT_WXY',  'transformation metrics from Z to Xi, G^1/2 at WXY', '1', 'ZHXY',                                            & ! [IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_GSQRT(:,:,:,I_UYW), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'GSQRT_WUY',  'transformation metrics from Z to Xi, G^1/2 at WUY', '1', 'ZHXHY',                                           & ! [IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_GSQRT(:,:,:,I_XVW), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'GSQRT_WXV',  'transformation metrics from Z to Xi, G^1/2 at WXV', '1', 'ZHXYH',                                           & ! [IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_GSQRT(:,:,:,I_UYZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'GSQRT_ZUY',  'transformation metrics from Z to Xi, G^1/2 at ZUY', '1', 'ZXHY',                                            & ! [IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_GSQRT(:,:,:,I_XVZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'GSQRT_ZXV',  'transformation metrics from Z to Xi, G^1/2 at ZXV', '1', 'ZXYH',                                            & ! [IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_GSQRT(:,:,:,I_UVZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'GSQRT_ZUV',  'transformation metrics from Z to Xi, G^1/2 at ZUV', '1', 'ZXHYH',                                           & ! [IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]

       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XYZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J13G_ZXY',  '(1,3) element of Jacobian matrix * {G}^1/2 at ZXY', '1', 'ZXY',                                             & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XYW), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J13G_WXY',  '(1,3) element of Jacobian matrix * {G}^1/2 at WXY', '1', 'ZHXY',                                            & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_UYW), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J13G_WUY',  '(1,3) element of Jacobian matrix * {G}^1/2 at WUY', '1', 'ZHXHY',                                           & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XVW), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J13G_WXV',  '(1,3) element of Jacobian matrix * {G}^1/2 at WXV', '1', 'ZHXYH',                                           & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_UYZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J13G_ZUY',  '(1,3) element of Jacobian matrix * {G}^1/2 at ZUY', '1', 'ZXHY',                                            & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_XVZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J13G_ZXV',  '(1,3) element of Jacobian matrix * {G}^1/2 at ZXV', '1', 'ZXYH',                                            & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J13G(:,:,:,I_UVZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J13G_ZUV',  '(1,3) element of Jacobian matrix * {G}^1/2 at ZUV', '1', 'ZXHYH',                                           & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]

       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XYZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J23G_ZXY',  '(2,3) element of Jacobian matrix * {G}^1/2 at ZXY', '1', 'ZXY',                                             & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XYW), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J23G_WXY',  '(2,3) element of Jacobian matrix * {G}^1/2 at WXY', '1', 'ZHXY',                                            & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_UYW), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J23G_WUY',  '(2,3) element of Jacobian matrix * {G}^1/2 at WUY', '1', 'ZHXHY',                                           & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XVW), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J23G_WXV',  '(2,3) element of Jacobian matrix * {G}^1/2 at WXV', '1', 'ZHXYH',                                           & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_UYZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J23G_ZUY',  '(2,3) element of Jacobian matrix * {G}^1/2 at ZUY', '1', 'ZXHY',                                            & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_XVZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J23G_ZXV',  '(2,3) element of Jacobian matrix * {G}^1/2 at ZXV', '1', 'ZXYH',                                            & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]
       if ( .not. PRC_TwoD ) &
       call FILE_CARTESC_write( ATMOS_GRID_CARTESC_METRIC_J23G(:,:,:,I_UVZ), ATMOS_GRID_CARTESC_METRIC_OUT_BASENAME, ATMOS_GRID_CARTESC_METRIC_OUT_TITLE, & ! [IN]
                                'J23G_ZUV',  '(2,3) element of Jacobian matrix * {G}^1/2 at ZUV', '1', 'ZXHYH',                                           & ![IN]
                                ATMOS_GRID_CARTESC_METRIC_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_GRID_CARTESC_METRIC_write

end module scale_atmos_grid_cartesC_metric
