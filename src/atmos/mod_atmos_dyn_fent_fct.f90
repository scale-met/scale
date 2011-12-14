!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics FENT + FCT
!!
!! @par Description
!!          Dynamical core for Atmospheric process
!!          Full explicit, no terrain + tracer FCT limiter
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-LES ver.2
!! @li      2011-11-11 (H.Yashiro) [mod] Merge with Y.Miyamoto's
!! @li      2011-12-11 (H.Yashiro) [mod] Use reference state
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_dyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_setup
  public :: ATMOS_DYN
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  integer, parameter :: I_DENS = 1
  integer, parameter :: I_MOMX = 2
  integer, parameter :: I_MOMY = 3
  integer, parameter :: I_MOMZ = 4
  integer, parameter :: I_POTT = 5
  integer, parameter :: I_PRES = 6
  integer, parameter :: I_VELX = 7
  integer, parameter :: I_VELY = 8
  integer, parameter :: I_VELZ = 9

  integer, parameter :: XDIR   = 1
  integer, parameter :: YDIR   = 2
  integer, parameter :: ZDIR   = 3

  ! time settings
  integer, parameter :: RK = 3 ! order of Runge-Kutta scheme

  ! advection settings
  real(8), parameter :: FACT_N =   7.D0 / 6.D0 !  7/6: fourth, 1: second
  real(8), parameter :: FACT_F = - 1.D0 / 6.D0 ! -1/6: fourth, 0: second

  ! numerical filter settings
  integer, parameter :: DF   = 4     ! order of numerical filter

  real(8), save      :: ATMOS_DYN_numerical_diff = 1.D-3 ! nondimensional numerical diffusion
  real(8), save      :: DIFF

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Initialize Dynamical Process
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_grid, only: &
       GRID_DX
    implicit none

    NAMELIST / PARAM_ATMOS_DYN / &
       ATMOS_DYN_numerical_diff

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Dynamics]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_DYN)

    DIFF = ATMOS_DYN_numerical_diff * (-1.D0)**dble( DF/2+1 )


  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN( dens,   momx,   momy,   momz,   pott,   qtrc,  &
                        dens_t, momx_t, momy_t, momz_t, pott_t, qtrc_t )
    use mod_const, only : &
       GRAV   => CONST_GRAV,  &
       Rair   => CONST_Rair,   &
       CPovCV => CONST_CPovCV, &
       Pstd   => CONST_Pstd
    use mod_time, only: &
       TIME_DTSEC_ATMOS_DYN
    use mod_comm, only: &
       COMM_vars
    use mod_grid, only : &
       IA  => GRID_IA, &
       JA  => GRID_JA, &
       KA  => GRID_KA, &
       IS  => GRID_IS, &
       IE  => GRID_IE, &
       JS  => GRID_JS, &
       JE  => GRID_JE, &
       KS  => GRID_KS, &
       KE  => GRID_KE, &
       WS  => GRID_WS, &
       WE  => GRID_WE, &
       RDX => GRID_RDX, &
       RDY => GRID_RDY, &
       RDZ => GRID_RDZ
    use mod_atmos_vars, only: &
       QA => A_QA
    use mod_atmos_refstate, only: &
       REF_dens
    use mod_atmos_boundary, only: &
       DAMP_alphau, &
       DAMP_alphav, &
       DAMP_alphaw, &
       DAMP_alphat, &
       velx_ref,    &
       vely_ref,    &
       velz_ref,    &
       pott_ref
    implicit none

    ! prognostic value
    real(8), intent(in)  :: dens(IA,JA,KA)      ! density [kg/m3]
    real(8), intent(in)  :: momx(IA,JA,KA)      ! momentum (x) [kg/m3 * m/s]
    real(8), intent(in)  :: momy(IA,JA,KA)      ! momentum (y) [kg/m3 * m/s]
    real(8), intent(in)  :: momz(IA,JA,KA)      ! momentum (z) [kg/m3 * m/s]
    real(8), intent(in)  :: pott(IA,JA,KA)      ! potential temperature [K]
    real(8), intent(in)  :: qtrc(IA,JA,KA,QA)   ! tracer mixing ratio   [kg/kg],[1/m3]
    ! prognostic tendency
    real(8), intent(out) :: dens_t(IA,JA,KA)
    real(8), intent(out) :: momx_t(IA,JA,KA)
    real(8), intent(out) :: momy_t(IA,JA,KA)
    real(8), intent(out) :: momz_t(IA,JA,KA)
    real(8), intent(out) :: pott_t(IA,JA,KA)
    real(8), intent(out) :: qtrc_t(IA,JA,KA,QA)

    ! work
    real(8), allocatable :: var(:,:,:,:)       ! work
    real(8), allocatable :: dens_diff(:,:,:)   ! anomary of density

    ! divergence
    real(8), allocatable :: ddiv(:,:,:)        ! density divergence [kg/m3/s]
    real(8)              :: qdiv               ! divergence for any quantity
    ! mass flux
    real(8), allocatable :: mflx_lo  (:,:,:,:) ! rho * vel(x,y,z) @ (u,v,w)-face low  order
    real(8), allocatable :: mflx_hi  (:,:,:,:) ! rho * vel(x,y,z) @ (u,v,w)-face high order
    real(8), allocatable :: mflx_anti(:,:,:,:) ! rho * vel(x,y,z) @ (u,v,w)-face antidiffusive
    real(8), allocatable :: qflx_lo  (:,:,:,:) ! rho * vel(x,y,z) * phi @ (u,v,w)-face low  order
    real(8), allocatable :: qflx_hi  (:,:,:,:) ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
    real(8), allocatable :: qflx_anti(:,:,:,:) ! rho * vel(x,y,z) * phi @ (u,v,w)-face antidiffusive
    ! flux correction
    real(8), allocatable :: rjpls    (:,:,:,:) ! plus  in (x,y,z)-direction
    real(8), allocatable :: rjmns    (:,:,:,:) ! minus in (x,y,z)-direction

    real(8) :: pjpls, pjmns, pjmax, pjmin, var_l
    real(8) :: midvel

    real(8) :: dtrk
    integer :: i, j, k, iq, rko
    !---------------------------------------------------------------------------

    allocate( var(IA,JA,KA,9) )
 
    allocate( dens_diff(IA,JA,KA)           )
    allocate( ddiv     (IA,JA,KA)           )

    allocate( mflx_lo  (IA,JA,KA,XDIR:ZDIR) )
    allocate( mflx_hi  (IA,JA,KA,XDIR:ZDIR) )
    allocate( mflx_anti(IA,JA,KA,XDIR:ZDIR) )
    allocate( qflx_lo  (IA,JA,KA,XDIR:ZDIR) )
    allocate( qflx_hi  (IA,JA,KA,XDIR:ZDIR) )
    allocate( qflx_anti(IA,JA,KA,XDIR:ZDIR) )

    allocate( rjpls    (IA,JA,KA,XDIR:ZDIR) )
    allocate( rjmns    (IA,JA,KA,XDIR:ZDIR) )

    ! copy to work from global vars
    var(:,:,:,I_DENS) = dens(:,:,:)
    var(:,:,:,I_MOMX) = momx(:,:,:)
    var(:,:,:,I_MOMY) = momy(:,:,:)
    var(:,:,:,I_MOMZ) = momz(:,:,:)
    var(:,:,:,I_POTT) = pott(:,:,:)

    do rko = 1, RK
       dtrk = TIME_DTSEC_ATMOS_DYN / (RK - rko + 1)

       !--- calc pressure, velocity & communication

       ! momentum -> velocity
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          var(i,j,k,I_VELX) = 2.D0 * var(i,j,k,I_MOMX) / ( var(i+1,j  ,k,I_DENS)+var(i,j,k,I_DENS) )
          var(i,j,k,I_VELY) = 2.D0 * var(i,j,k,I_MOMY) / ( var(i  ,j+1,k,I_DENS)+var(i,j,k,I_DENS) )
       enddo
       enddo
       enddo

       do k = WS+1, WE-1
       do j = JS,   JE
       do i = IS,   IE
          var(i,j,k,I_VELZ) = 2.D0 * var(i,j,k,I_MOMZ) / ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) )
       enddo
       enddo
       enddo
       var(:,:,WS,I_VELZ) = 0.D0 ! bottom boundary
       var(:,:,WE,I_VELZ) = 0.D0 ! top    boundary

       ! pressure
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          var(i,j,k,I_PRES) = Pstd * ( var(i,j,k,I_DENS) * var(i,j,k,I_POTT) * Rair / Pstd )**CPovCV
       enddo
       enddo
       enddo

       ! fill IHALO & JHALO
       call COMM_vars( var(:,:,:,6:9) )

       !##### continuity equation #####
       do k = KS,   KE
       do j = JS,   JE
       do i = IS-1, IE
          mflx_lo(i,j,k,XDIR) = 0.5D0                                                                &
                              * (     var(i,j,k,I_VELX)  * ( var(i+1,j,k,I_DENS)+var(i,j,k,I_DENS) ) &
                                - abs(var(i,j,k,I_VELX)) * ( var(i+1,j,k,I_DENS)-var(i,j,k,I_DENS) ) )

          mflx_hi(i,j,k,XDIR) = 0.5D0 * var(i,j,k,I_VELX)                              &
                              * ( FACT_N * ( var(i+1,j,k,I_DENS)+var(i  ,j,k,I_DENS) ) &
                                + FACT_F * ( var(i+2,j,k,I_DENS)+var(i-1,j,k,I_DENS) ) )

          mflx_anti(i,j,k,XDIR) = mflx_hi(i,j,k,XDIR) - mflx_lo(i,j,k,XDIR)
       enddo
       enddo
       enddo

       do k = KS,   KE
       do j = JS-1, JE
       do i = IS,   IE
          mflx_lo(i,j,k,YDIR) = 0.5D0                                                                &
                              * (     var(i,j,k,I_VELY)  * ( var(i,j+1,k,I_DENS)+var(i,j,k,I_DENS) ) &
                                - abs(var(i,j,k,I_VELY)) * ( var(i,j+1,k,I_DENS)-var(i,j,k,I_DENS) ) )

          mflx_hi(i,j,k,YDIR) = 0.5D0 * var(i,j,k,I_VELY)                              &
                              * ( FACT_N * ( var(i,j+1,k,I_DENS)+var(i,j  ,k,I_DENS) ) &
                                + FACT_F * ( var(i,j+2,k,I_DENS)+var(i,j-1,k,I_DENS) ) )

          mflx_anti(i,j,k,YDIR) = mflx_hi(i,j,k,YDIR) - mflx_lo(i,j,k,YDIR)
       enddo
       enddo
       enddo

       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          mflx_lo(i,j,k,ZDIR) = 0.5D0                                                                &
                              * (     var(i,j,k,I_VELZ)  * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) ) &
                                - abs(var(i,j,k,I_VELZ)) * ( var(i,j,k+1,I_DENS)-var(i,j,k,I_DENS) ) )

          mflx_hi(i,j,k,ZDIR) = 0.5D0 * var(i,j,k,I_VELZ)                              &
                              * ( FACT_N * ( var(i,j,k+1,I_DENS)+var(i,j,k  ,I_DENS) ) &
                                + FACT_F * ( var(i,j,k+2,I_DENS)+var(i,j,k-1,I_DENS) ) )

          mflx_anti(i,j,k,ZDIR) = mflx_hi(i,j,k,ZDIR) - mflx_lo(i,j,k,ZDIR)
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          mflx_lo(i,j,WS  ,ZDIR) = 0.D0                                          ! bottom boundary
          mflx_lo(i,j,WS+1,ZDIR) = 0.5D0 * var(i,j,WS+1,I_VELZ)                &
                                 * ( var(i,j,WS+2,I_DENS)+var(i,j,WS+1,I_DENS) ) ! just above the bottom boundary
          mflx_lo(i,j,WE-1,ZDIR) = 0.5D0 * var(i,j,WE-1,I_VELZ)                &
                                 * ( var(i,j,WE  ,I_DENS)+var(i,j,WE-1,I_DENS) ) ! just below the top boundary
          mflx_lo(i,j,WE  ,ZDIR) = 0.D0                                          ! top boundary

          mflx_hi(i,j,WS  ,ZDIR) = 0.D0
          mflx_hi(i,j,WS+1,ZDIR) = mflx_lo(i,j,WS+1,ZDIR)
          mflx_hi(i,j,WE-1,ZDIR) = mflx_lo(i,j,WE-1,ZDIR)
          mflx_hi(i,j,WE  ,ZDIR) = 0.D0

          mflx_anti(i,j,WS  ,ZDIR) = 0.D0
          mflx_anti(i,j,WS+1,ZDIR) = 0.D0
          mflx_anti(i,j,WE-1,ZDIR) = 0.D0
          mflx_anti(i,j,WE  ,ZDIR) = 0.D0
       enddo
       enddo

       !--- flux-divergence
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          ddiv(i,j,k) = ( mflx_hi(i,j,k,XDIR)-mflx_hi(i-1,j,  k  ,XDIR) ) * RDX &
                      + ( mflx_hi(i,j,k,YDIR)-mflx_hi(i,  j-1,k  ,YDIR) ) * RDY &
                      + ( mflx_hi(i,j,k,ZDIR)-mflx_hi(i,  j,  k-1,ZDIR) ) * RDZ
       enddo
       enddo
       enddo

       !--- numerical filter
       do k = KS, KE
          dens_diff(:,:,k) = var(:,:,k,I_DENS) - REF_dens(:,:,k)
       enddo

       do k = KS+2, KE-2
       do j = JS,   JE
       do i = IS,   IE
          dens_t(i,j,k) = - ddiv(i,j,k)                                       &
                          - DIFF / dtrk * ( (        dens_diff(i+2,j  ,k  )   &
                                            - 4.D0 * dens_diff(i+1,j  ,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i-1,j  ,k  )   &
                                            +        dens_diff(i-2,j  ,k  ) ) &
                                          + (        dens_diff(i  ,j+2,k  )   &
                                            - 4.D0 * dens_diff(i  ,j+1,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i  ,j-1,k  )   &
                                            +        dens_diff(i  ,j-2,k  ) ) &
                                          + (        dens_diff(i  ,j  ,k+2)   &
                                            - 4.D0 * dens_diff(i  ,j  ,k+1)   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i  ,j  ,k-1)   &
                                            +        dens_diff(i  ,j  ,k-2) ) )
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          k = KS
          dens_t(i,j,k) = - ddiv(i,j,k)                                       &
                          - DIFF / dtrk * ( (        dens_diff(i+2,j  ,k  )   &
                                            - 4.D0 * dens_diff(i+1,j  ,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i-1,j  ,k  )   &
                                            +        dens_diff(i-2,j  ,k  ) ) &
                                          + (        dens_diff(i  ,j+2,k  )   &
                                            - 4.D0 * dens_diff(i  ,j+1,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i  ,j-1,k  )   &
                                            +        dens_diff(i  ,j-2,k  ) ) &
                                          + (        dens_diff(i  ,j  ,k+2)   &
                                            - 3.D0 * dens_diff(i  ,j  ,k+1)   &
                                            + 2.D0 * dens_diff(i  ,j  ,k  ) ) )

          k = KS+1
          dens_t(i,j,k) = - ddiv(i,j,k)                                       &
                          - DIFF / dtrk * ( (        dens_diff(i+2,j  ,k  )   &
                                            - 4.D0 * dens_diff(i+1,j  ,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i-1,j  ,k  )   &
                                            +        dens_diff(i-2,j  ,k  ) ) &
                                          + (        dens_diff(i  ,j+2,k  )   &
                                            - 4.D0 * dens_diff(i  ,j+1,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i  ,j-1,k  )   &
                                            +        dens_diff(i  ,j-2,k  ) ) &
                                          + (        dens_diff(i  ,j  ,k+2)   &
                                            - 4.D0 * dens_diff(i  ,j  ,k+1)   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 3.D0 * dens_diff(i  ,j  ,k-1) ) )

          k = KE-1
          dens_t(i,j,k) = - ddiv(i,j,k)                                       &
                          - DIFF / dtrk * ( (        dens_diff(i+2,j  ,k  )   &
                                            - 4.D0 * dens_diff(i+1,j  ,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i-1,j  ,k  )   &
                                            +        dens_diff(i-2,j  ,k  ) ) &
                                          + (        dens_diff(i  ,j+2,k  )   &
                                            - 4.D0 * dens_diff(i  ,j+1,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i  ,j-1,k  )   &
                                            +        dens_diff(i  ,j-2,k  ) ) &
                                          + (        dens_diff(i  ,j  ,k-2)   &
                                            - 4.D0 * dens_diff(i  ,j  ,k-1)   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 3.D0 * dens_diff(i  ,j  ,k+1) ) )

          k = KE
          dens_t(i,j,k) = - ddiv(i,j,k)                                       &
                          - DIFF / dtrk * ( (        dens_diff(i+2,j  ,k  )   &
                                            - 4.D0 * dens_diff(i+1,j  ,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i-1,j  ,k  )   &
                                            +        dens_diff(i-2,j  ,k  ) ) &
                                          + (        dens_diff(i  ,j+2,k  )   &
                                            - 4.D0 * dens_diff(i  ,j+1,k  )   &
                                            + 6.D0 * dens_diff(i  ,j  ,k  )   &
                                            - 4.D0 * dens_diff(i  ,j-1,k  )   &
                                            +        dens_diff(i  ,j-2,k  ) ) &
                                          + (        dens_diff(i  ,j  ,k-2)   &
                                            - 3.D0 * dens_diff(i  ,j  ,k-1)   &
                                            + 2.D0 * dens_diff(i  ,j  ,k  ) ) )
       enddo
       enddo


       !##### momentum equation #####

       !--- make mass fluxes ( momentum(x) * vel )

       ! at (x, y, layer)
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE+1
          midvel = 0.5D0 * ( var(i,j,k,I_VELX)+var(i-1,j,k,I_VELX) ) ! at center

          qflx_lo(i,j,k,XDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i,j,k,I_MOMX)+var(i-1,j,k,I_MOMX) ) &
                                - abs(midvel) * ( var(i,j,k,I_MOMX)-var(i-1,j,k,I_MOMX) ) )

          qflx_hi(i,j,k,XDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i  ,j,k,I_MOMX)+var(i-1,j,k,I_MOMX) ) &
                                + FACT_F * ( var(i+1,j,k,I_MOMX)+var(i-2,j,k,I_MOMX) ) )

          qflx_anti(i,j,k,XDIR) = qflx_hi(i,j,k,XDIR) - qflx_lo(i,j,k,XDIR)
       enddo
       enddo
       enddo

       ! at (u, v, layer)
       do k = KS,   KE
       do j = JS-1, JE
       do i = IS,   IE
          midvel = 0.5D0 * ( var(i+1,j,k,I_VELY)+var(i,j,k,I_VELY) ) ! at face

          qflx_lo(i,j,k,YDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i,j+1,k,I_MOMX)+var(i,j,k,I_MOMX) ) &
                                - abs(midvel) * ( var(i,j+1,k,I_MOMX)-var(i,j,k,I_MOMX) ) ) 

          qflx_hi(i,j,k,YDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i,j+1,k,I_MOMX)+var(i,j  ,k,I_MOMX) ) &
                                + FACT_F * ( var(i,j+2,k,I_MOMX)+var(i,j-1,k,I_MOMX) ) )

          qflx_anti(i,j,k,YDIR) = qflx_hi(i,j,k,YDIR) - qflx_lo(i,j,k,YDIR)
       enddo
       enddo
       enddo

       ! at (u, y, interface)
       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          midvel = 0.5D0 * ( var(i+1,j,k,I_VELZ)+var(i,j,k,I_VELZ) ) ! at face

          qflx_lo(i,j,k,ZDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i,j,k+1,I_MOMX)+var(i,j,k,I_MOMX) ) &
                                - abs(midvel) * ( var(i,j,k+1,I_MOMX)-var(i,j,k,I_MOMX) ) ) 

          qflx_hi(i,j,k,ZDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i,j,k+1,I_MOMX)+var(i,j,k  ,I_MOMX) ) &
                                + FACT_F * ( var(i,j,k+2,I_MOMX)+var(i,j,k-1,I_MOMX) ) )

          qflx_anti(i,j,k,ZDIR) = qflx_hi(i,j,k,ZDIR) - qflx_lo(i,j,k,ZDIR)
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_lo(i,j,WS  ,ZDIR) = 0.D0                                                     ! bottom boundary
          qflx_lo(i,j,WS+1,ZDIR) = 0.25D0 * ( var(i,j,WS+2,I_MOMX)+var(i,j,WS+1,I_MOMX) ) &
                                 * ( var(i+1,j,WS+2,I_VELZ)+var(i,j,WS+1,I_VELZ) )          ! just above the bottom boundary
          qflx_lo(i,j,WE-1,ZDIR) = 0.25D0 * ( var(i,j,WE  ,I_MOMX)+var(i,j,WE-1,I_MOMX) ) &
                                 * ( var(i+1,j,WE  ,I_VELZ)+var(i,j,WE-1,I_VELZ) )          ! just below the top boundary
          qflx_lo(i,j,WE  ,ZDIR) = 0.D0                                                     ! top boundary

          qflx_hi(i,j,WS  ,ZDIR) = 0.D0
          qflx_hi(i,j,WS+1,ZDIR) = qflx_lo(i,j,WS+1,ZDIR)
          qflx_hi(i,j,WE-1,ZDIR) = qflx_lo(i,j,WE-1,ZDIR)
          qflx_hi(i,j,WE  ,ZDIR) = 0.D0

          qflx_anti(i,j,WS  ,ZDIR) = 0.D0
          qflx_anti(i,j,WS+1,ZDIR) = 0.D0
          qflx_anti(i,j,WE-1,ZDIR) = 0.D0
          qflx_anti(i,j,WE  ,ZDIR) = 0.D0
       enddo
       enddo

       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          !--- update flux-divergence with the monotone scheme
          qdiv = ( qflx_lo(i+1,j,k,XDIR)-qflx_lo(i,j,  k  ,XDIR) ) * RDX &
               + ( qflx_lo(i  ,j,k,YDIR)-qflx_lo(i,j-1,k  ,YDIR) ) * RDY &
               + ( qflx_lo(i  ,j,k,ZDIR)-qflx_lo(i,j,  k-1,ZDIR) ) * RDZ

          !--- first guess of tendency and updated value
          momx_t(i,j,k) = - qdiv                                                  &
                          - (var(i+1,j,k,I_PRES)-var(i,j,k,I_PRES) ) * RDX        &
                          - DAMP_alphau(i,j,k) * (var(i,j,k,I_VELX)-velx_ref(i,j,k) ) &
                          * 0.5D0 * ( var(i+1,j,k,I_DENS)+var(i,j,k,I_DENS) )

          var_l = momx(i,j,k) + momx_t(i,j,k) * dtrk

          !--- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
          pjmax = max( var(i  ,j  ,k  ,I_MOMX), &
                       var(i+1,j  ,k  ,I_MOMX), &
                       var(i-1,j  ,k  ,I_MOMX), &
                       var(i  ,j+1,k  ,I_MOMX), &
                       var(i  ,j-1,k  ,I_MOMX), &
                       var(i  ,j  ,k+1,I_MOMX), &
                       var(i  ,j  ,k-1,I_MOMX)  )

          pjmin = min( var(i  ,j  ,k  ,I_MOMX), &
                       var(i+1,j  ,k  ,I_MOMX), &
                       var(i-1,j  ,k  ,I_MOMX), &
                       var(i  ,j+1,k  ,I_MOMX), &
                       var(i  ,j-1,k  ,I_MOMX), &
                       var(i  ,j  ,k+1,I_MOMX), &
                       var(i  ,j  ,k-1,I_MOMX)  )

          !--- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
          pjpls = max( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i+1,j  ,k  ,XDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j-1,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k-1,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) )
          pjmns = max( 0.D0, qflx_anti(i+1,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j-1,k  ,YDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k-1,ZDIR) )

          !--- incoming fluxes at scalar grid points ---
          if ( pjpls > 0 ) then
             rjpls(i,j,k,XDIR) = (pjmax-var_l) / pjpls * abs(var(i,j,k,I_VELX))

             midvel = 0.25D0 * ( var(i+1,j  ,k  ,I_VELY) &
                               + var(i  ,j  ,k  ,I_VELY) &
                               + var(i+1,j-1,k  ,I_VELY) &
                               + var(i  ,j-1,k  ,I_VELY) )
             rjpls(i,j,k,YDIR) = (pjmax-var_l) / pjpls * abs(midvel)

             midvel = 0.25D0 * ( var(i+1,j  ,k  ,I_VELZ) &
                               + var(i  ,j  ,k  ,I_VELZ) &
                               + var(i+1,j  ,k-1,I_VELZ) &
                               + var(i  ,j  ,k-1,I_VELZ) )
             rjpls(i,j,k,ZDIR) = (pjmax-var_l) / pjpls * abs(midvel)
          else
             rjpls(i,j,k,XDIR:ZDIR) = 0.D0
          endif

          !--- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             rjmns(i,j,k,XDIR) = (var_l-pjmin) / pjmns * abs(var(i,j,k,I_VELX))

             midvel = 0.25D0 * ( var(i+1,j  ,k  ,I_VELY) &
                               + var(i  ,j  ,k  ,I_VELY) &
                               + var(i+1,j-1,k  ,I_VELY) &
                               + var(i  ,j-1,k  ,I_VELY) )
             rjmns(i,j,k,YDIR) = (var_l-pjmin) / pjmns * abs(midvel)

             midvel = 0.25D0 * ( var(i+1,j  ,k  ,I_VELZ) &
                               + var(i  ,j  ,k  ,I_VELZ) &
                               + var(i+1,j  ,k-1,I_VELZ) &
                               + var(i  ,j  ,k-1,I_VELZ) )
             rjmns(i,j,k,ZDIR) = (var_l-pjmin) / pjmns * abs(midvel)
          else
             rjmns(i,j,k,XDIR:ZDIR) = 0.D0
          endif

       enddo
       enddo
       enddo

       !--- [STEP 7S] limit the antidiffusive flux at velocity grid points ---
       do k = KS-1, KE
       do j = JS-1, JE
       do i = IS  , IE+1
          if ( qflx_anti(i,j,k,XDIR) >= 0 ) then
             qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i  ,j,k,XDIR), rjmns(i-1,j,k,XDIR), 1.D0 )
          else
             qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i-1,j,k,XDIR), rjmns(i  ,j,k,XDIR), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,YDIR) >= 0 ) then
             qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j+1,k,YDIR), rjmns(i,j  ,k,YDIR), 1.D0 )
          else
             qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j  ,k,YDIR), rjmns(i,j+1,k,YDIR), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,ZDIR) >= 0 ) then
             qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k+1,ZDIR), rjmns(i,j,k  ,ZDIR), 1.D0 )
          else
             qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k  ,ZDIR), rjmns(i,j,k+1,ZDIR), 1.D0 )
          endif
       enddo
       enddo
       enddo

       !--- modify flux-divergence->tendency with antidiffusive fluxes
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          momx_t(i,j,k) = momx_t(i,j,k) - ( ( qflx_anti(i+1,j,k,XDIR)-qflx_anti(i,j  ,k  ,XDIR) ) * RDX &
                                          + ( qflx_anti(i  ,j,k,YDIR)-qflx_anti(i,j-1,k  ,YDIR) ) * RDY &
                                          + ( qflx_anti(i  ,j,k,ZDIR)-qflx_anti(i,j  ,k-1,ZDIR) ) * RDZ )
       enddo
       enddo
       enddo

       !--- make mass fluxes ( momentum(y) * vel )

       ! at (u, v, layer)
       do k = KS,   KE
       do j = JS,   JE
       do i = IS-1, IE
          midvel = 0.5D0 * ( var(i,j+1,k,I_VELX)+var(i,j,k,I_VELX) ) ! at center

          qflx_lo(i,j,k,XDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i+1,j,k,I_MOMY)+var(i,j,k,I_MOMY) ) &
                                - abs(midvel) * ( var(i+1,j,k,I_MOMY)-var(i,j,k,I_MOMY) ) )

          qflx_hi(i,j,k,XDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i+1,j,k,I_MOMY)+var(i  ,j,k,I_MOMY) ) &
                                + FACT_F * ( var(i+2,j,k,I_MOMY)+var(i-1,j,k,I_MOMY) ) )

          qflx_anti(i,j,k,XDIR) = qflx_hi(i,j,k,XDIR) - qflx_lo(i,j,k,XDIR)
       enddo
       enddo
       enddo
       ! at (x, y, layer)
       do k = KS, KE
       do j = JS, JE+1
       do i = IS, IE
          midvel = 0.5D0 * ( var(i,j,k,I_VELY)+var(i,j-1,k,I_VELY) ) ! at face

          qflx_lo(i,j,k,YDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i,j,k,I_MOMY)+var(i,j-1,k,I_MOMY) ) &
                                - abs(midvel) * ( var(i,j,k,I_MOMY)-var(i,j-1,k,I_MOMY) ) ) 

          qflx_hi(i,j,k,YDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i,j  ,k,I_MOMY)+var(i,j-1,k,I_MOMY) ) &
                                + FACT_F * ( var(i,j+1,k,I_MOMY)+var(i,j-2,k,I_MOMY) ) )

          qflx_anti(i,j,k,YDIR) = qflx_hi(i,j,k,YDIR) - qflx_lo(i,j,k,YDIR)
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          midvel = 0.5D0 * ( var(i,j+1,k,I_VELZ)+var(i,j,k,I_VELZ) ) ! at face

          qflx_lo(i,j,k,ZDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i,j,k+1,I_MOMY)+var(i,j,k,I_MOMY) ) &
                                - abs(midvel) * ( var(i,j,k+1,I_MOMY)-var(i,j,k,I_MOMY) ) ) 

          qflx_hi(i,j,k,ZDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i,j,k+1,I_MOMY)+var(i,j,k  ,I_MOMY) ) &
                                + FACT_F * ( var(i,j,k+2,I_MOMY)+var(i,j,k-1,I_MOMY) ) )

          qflx_anti(i,j,k,ZDIR) = qflx_hi(i,j,k,ZDIR) - qflx_lo(i,j,k,ZDIR)
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_lo(i,j,WS  ,ZDIR) = 0.D0                                                     ! bottom boundary
          qflx_lo(i,j,WS+1,ZDIR) = 0.25D0 * ( var(i,j,WS+2,I_MOMY)+var(i,j,WS+1,I_MOMY) ) &
                                 * ( var(i,j+1,WS+2,I_VELZ)+var(i,j,WS+1,I_VELZ) )          ! just above the bottom boundary
          qflx_lo(i,j,WE-1,ZDIR) = 0.25D0 * ( var(i,j,WE  ,I_MOMY)+var(i,j,WE-1,I_MOMY) ) &
                                 * ( var(i,j+1,WE  ,I_VELZ)+var(i,j,WE-1,I_VELZ) )          ! just below the top boundary
          qflx_lo(i,j,WE  ,ZDIR) = 0.D0                                                     ! top boundary

          qflx_hi(i,j,WS  ,ZDIR) = 0.D0
          qflx_hi(i,j,WS+1,ZDIR) = qflx_lo(i,j,WS+1,ZDIR)
          qflx_hi(i,j,WE-1,ZDIR) = qflx_lo(i,j,WE-1,ZDIR)
          qflx_hi(i,j,WE  ,ZDIR) = 0.D0

          qflx_anti(i,j,WS  ,ZDIR) = 0.D0
          qflx_anti(i,j,WS+1,ZDIR) = 0.D0
          qflx_anti(i,j,WE-1,ZDIR) = 0.D0
          qflx_anti(i,j,WE  ,ZDIR) = 0.D0
       enddo
       enddo

       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          !--- update flux-divergence with the monotone scheme
          qdiv = ( qflx_lo(i,j  ,k,XDIR)-qflx_lo(i-1,j,k  ,XDIR) ) * RDX &
               + ( qflx_lo(i,j+1,k,YDIR)-qflx_lo(i  ,j,k  ,YDIR) ) * RDY &
               + ( qflx_lo(i,j  ,k,ZDIR)-qflx_lo(i  ,j,k-1,ZDIR) ) * RDZ

          !--- first guess of tendency and updated value
          momy_t(i,j,k) = - qdiv                                                  &
                          - (var(i,j+1,k,I_PRES)-var(i,j,k,I_PRES) ) * RDY        &
                          - DAMP_alphav(i,j,k) * (var(i,j,k,I_VELY)-vely_ref(i,j,k) ) &
                          * 0.5D0 * ( var(i,j+1,k,I_DENS)+var(i,j,k,I_DENS) )

          var_l = momy(i,j,k) + dtrk * momy_t(i,j,k)

          !--- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
          pjmax = max( var(i  ,j  ,k  ,I_MOMY), &
                       var(i+1,j  ,k  ,I_MOMY), &
                       var(i-1,j  ,k  ,I_MOMY), &
                       var(i  ,j+1,k  ,I_MOMY), &
                       var(i  ,j-1,k  ,I_MOMY), &
                       var(i  ,j  ,k+1,I_MOMY), &
                       var(i  ,j  ,k-1,I_MOMY)  )

          pjmin = min( var(i  ,j  ,k  ,I_MOMY), &
                       var(i+1,j  ,k  ,I_MOMY), &
                       var(i-1,j  ,k  ,I_MOMY), &
                       var(i  ,j+1,k  ,I_MOMY), &
                       var(i  ,j-1,k  ,I_MOMY), &
                       var(i  ,j  ,k+1,I_MOMY), &
                       var(i  ,j  ,k-1,I_MOMY)  )

          !--- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
          pjpls = max( 0.D0, qflx_anti(i-1,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j+1,k  ,YDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k-1,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) )
          pjmns = max( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i-1,j  ,k  ,XDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j+1,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k-1,ZDIR) )

          !--- incoming fluxes at scalar grid points ---
          if ( pjpls > 0 ) then
             midvel = 0.25D0 * ( var(i  ,j+1,k  ,I_VELX) &
                               + var(i  ,j  ,k  ,I_VELX) &
                               + var(i-1,j+1,k  ,I_VELX) &
                               + var(i-1,j  ,k  ,I_VELX) )
             rjpls(i,j,k,XDIR) = (pjmax-var_l) / pjpls * abs(midvel)

             rjpls(i,j,k,YDIR) = (pjmax-var_l) / pjpls * abs(var(i,j,k,I_VELY))

             midvel = 0.25D0 * ( var(i  ,j+1,k  ,I_VELZ) &
                               + var(i  ,j  ,k  ,I_VELZ) &
                               + var(i  ,j+1,k-1,I_VELZ) &
                               + var(i  ,j  ,k-1,I_VELZ) )
             rjpls(i,j,k,ZDIR) = (pjmax-var_l) / pjpls * abs(midvel)
          else
             rjpls(i,j,k,XDIR:ZDIR) = 0.D0
          endif

          !--- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             midvel = 0.25D0 * ( var(i  ,j+1,k  ,I_VELX) &
                               + var(i  ,j  ,k  ,I_VELX) &
                               + var(i-1,j+1,k  ,I_VELX) &
                               + var(i-1,j  ,k  ,I_VELX) )
             rjmns(i,j,k,XDIR) = (var_l-pjmin) / pjmns * abs(midvel)

             rjmns(i,j,k,YDIR) = (var_l-pjmin) / pjmns * abs(var(i,j,k,I_VELY))

             midvel = 0.25D0 * ( var(i  ,j+1,k  ,I_VELZ) &
                               + var(i  ,j  ,k  ,I_VELZ) &
                               + var(i  ,j+1,k-1,I_VELZ) &
                               + var(i  ,j  ,k-1,I_VELZ) )
             rjmns(i,j,k,ZDIR) = (var_l-pjmin) / pjmns * abs(midvel)
          else
             rjmns(i,j,k,XDIR:ZDIR) = 0.D0
          endif

       enddo
       enddo
       enddo

       !--- [STEP 7S] limit the antidiffusive flux at velocity grid points ---
       do k = KS-1, KE
       do j = JS,   JE+1
       do i = IS-1, IE
          if ( qflx_anti(i,j,k,XDIR) >= 0 ) then
             qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i+1,j,k,XDIR), rjmns(i  ,j,k,XDIR), 1.D0 )
          else
             qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i  ,j,k,XDIR), rjmns(i+1,j,k,XDIR), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,YDIR) >= 0 ) then
             qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j  ,k,YDIR), rjmns(i,j-1,k,YDIR), 1.D0 )
          else
             qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j-1,k,YDIR), rjmns(i,j  ,k,YDIR), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,ZDIR) >= 0 ) then
             qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k+1,ZDIR), rjmns(i,j,k  ,ZDIR), 1.D0 )
          else
             qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k  ,ZDIR), rjmns(i,j,k+1,ZDIR), 1.D0 )
          endif
       enddo
       enddo
       enddo

       !--- modify flux-divergence->tendency with antidiffusive fluxes
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          momy_t(i,j,k) = momy_t(i,j,k) - ( ( qflx_anti(i,j  ,k,XDIR)-qflx_anti(i-1,j,k  ,XDIR) ) * RDX &
                                          + ( qflx_anti(i,j+1,k,YDIR)-qflx_anti(i  ,j,k  ,YDIR) ) * RDY &
                                          + ( qflx_anti(i,j  ,k,ZDIR)-qflx_anti(i  ,j,k-1,ZDIR) ) * RDZ )
       enddo
       enddo
       enddo

       !--- make mass fluxes ( momentum(z) * vel )

       ! at (u, y, interface)
       k = WS ! bottom boundary
       do j = JS,   JE
       do i = IS-1, IE
          midvel = var(i,j,KS,I_VELX)

          qflx_lo(i,j,k,XDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i+1,j,k,I_MOMZ)+var(i,j,k,I_MOMZ) ) &
                                - abs(midvel) * ( var(i+1,j,k,I_MOMZ)-var(i,j,k,I_MOMZ) ) )

          qflx_hi(i,j,k,XDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i+1,j,k,I_MOMZ)+var(i  ,j,k,I_MOMZ) ) &
                                + FACT_F * ( var(i+2,j,k,I_MOMZ)+var(i-1,j,k,I_MOMZ) ) )

          qflx_anti(i,j,k,XDIR) = qflx_hi(i,j,k,XDIR) - qflx_lo(i,j,k,XDIR)
       enddo
       enddo
       k = WE ! top boundary
       do j = JS,   JE
       do i = IS-1, IE
          midvel = var(i,j,KE,I_VELX)

          qflx_lo(i,j,k,XDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i+1,j,k,I_MOMZ)+var(i,j,k,I_MOMZ) ) &
                                - abs(midvel) * ( var(i+1,j,k,I_MOMZ)-var(i,j,k,I_MOMZ) ) )

          qflx_hi(i,j,k,XDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i+1,j,k,I_MOMZ)+var(i  ,j,k,I_MOMZ) ) &
                                + FACT_F * ( var(i+2,j,k,I_MOMZ)+var(i-1,j,k,I_MOMZ) ) )

          qflx_anti(i,j,k,XDIR) = qflx_hi(i,j,k,XDIR) - qflx_lo(i,j,k,XDIR)
       enddo
       enddo

       do k = WS+1, WE-1
       do j = JS,   JE
       do i = IS-1, IE
          midvel = 0.5D0 * ( var(i,j,k+1,I_VELX)+var(i,j,k,I_VELX) ) ! on face

          qflx_lo(i,j,k,XDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i+1,j,k,I_MOMZ)+var(i,j,k,I_MOMZ) ) &
                                - abs(midvel) * ( var(i+1,j,k,I_MOMZ)-var(i,j,k,I_MOMZ) ) )

          qflx_hi(i,j,k,XDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i+1,j,k,I_MOMZ)+var(i  ,j,k,I_MOMZ) ) &
                                + FACT_F * ( var(i+2,j,k,I_MOMZ)+var(i-1,j,k,I_MOMZ) ) )

          qflx_anti(i,j,k,XDIR) = qflx_hi(i,j,k,XDIR) - qflx_lo(i,j,k,XDIR)
       enddo
       enddo
       enddo

       ! at (x, v, interface)
       k = WS ! bottom boundary
       do j = JS-1, JE
       do i = IS,   IE
          midvel = var(i,j,KS,I_VELY)

          qflx_lo(i,j,k,YDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i,j+1,k,I_MOMZ)+var(i,j,k,I_MOMZ) ) &
                                - abs(midvel) * ( var(i,j+1,k,I_MOMZ)-var(i,j,k,I_MOMZ) ) )

          qflx_hi(i,j,k,YDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i,j+1,k,I_MOMZ)+var(i,j  ,k,I_MOMZ) ) &
                                + FACT_F * ( var(i,j+2,k,I_MOMZ)+var(i,j-1,k,I_MOMZ) ) )

          qflx_anti(i,j,k,YDIR) = qflx_hi(i,j,k,YDIR) - qflx_lo(i,j,k,YDIR)
       enddo
       enddo

       k = WE ! top boundary
       do j = JS-1, JE
       do i = IS,   IE
          midvel = var(i,j,KE,I_VELY)

          qflx_lo(i,j,k,YDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i,j+1,k,I_MOMZ)+var(i,j,k,I_MOMZ) ) &
                                - abs(midvel) * ( var(i,j+1,k,I_MOMZ)-var(i,j,k,I_MOMZ) ) )

          qflx_hi(i,j,k,YDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i,j+1,k,I_MOMZ)+var(i,j  ,k,I_MOMZ) ) &
                                + FACT_F * ( var(i,j+2,k,I_MOMZ)+var(i,j-1,k,I_MOMZ) ) )

          qflx_anti(i,j,k,YDIR) = qflx_hi(i,j,k,YDIR) - qflx_lo(i,j,k,YDIR)
       enddo
       enddo

       do k = WS+1, WE-1
       do j = JS-1, JE
       do i = IS,   IE
          midvel = 0.5D0 * ( var(i,j,k+1,I_VELY)+var(i,j,k,I_VELY) ) ! on face

          qflx_lo(i,j,k,YDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i,j+1,k,I_MOMZ)+var(i,j,k,I_MOMZ) ) &
                                - abs(midvel) * ( var(i,j+1,k,I_MOMZ)-var(i,j,k,I_MOMZ) ) )

          qflx_hi(i,j,k,YDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i,j+1,k,I_MOMZ)+var(i,j  ,k,I_MOMZ) ) &
                                + FACT_F * ( var(i,j+2,k,I_MOMZ)+var(i,j-1,k,I_MOMZ) ) )

          qflx_anti(i,j,k,YDIR) = qflx_hi(i,j,k,YDIR) - qflx_lo(i,j,k,YDIR)
       enddo
       enddo
       enddo

       ! at (x, y, layer)
       do k = KS+1, KE-1
       do j = JS,   JE
       do i = IS,   IE
          midvel = 0.5D0 * ( var(i,j,k,I_VELZ)+var(i,j,k-1,I_VELZ) ) ! at center

          qflx_lo(i,j,k,ZDIR) = 0.5D0                                                     &
                              * (     midvel  * ( var(i,j,k,I_MOMZ)+var(i,j,k-1,I_MOMZ) ) &
                                - abs(midvel) * ( var(i,j,k,I_MOMZ)-var(i,j,k-1,I_MOMZ) ) )

          qflx_hi(i,j,k,ZDIR) = 0.5D0 * midvel                                         &
                              * ( FACT_N * ( var(i,j,k  ,I_MOMZ)+var(i,j,k-1,I_MOMZ) ) &
                                + FACT_F * ( var(i,j,k+1,I_MOMZ)+var(i,j,k-2,I_MOMZ) ) )

          qflx_anti(i,j,k,ZDIR) = qflx_hi(i,j,k,ZDIR) - qflx_lo(i,j,k,ZDIR)
       enddo
       enddo
       enddo
       qflx_lo  (:,:,KS,ZDIR) = 0.D0 ! bottom cell center
       qflx_lo  (:,:,KE,ZDIR) = 0.D0 ! top    cell center
       qflx_hi  (:,:,KS,ZDIR) = 0.D0 ! bottom cell center
       qflx_hi  (:,:,KE,ZDIR) = 0.D0 ! top    cell center
       qflx_anti(:,:,KS,ZDIR) = 0.D0 ! bottom cell center
       qflx_anti(:,:,KE,ZDIR) = 0.D0 ! top    cell center

       do k = WS+1, WE-1
       do j = JS,   JE
       do i = IS,   IE
          !--- update flux-divergence with the monotone scheme
          qdiv = ( qflx_lo(i,j,k  ,XDIR)-qflx_lo(i-1,j  ,k,XDIR) ) * RDX &
               + ( qflx_lo(i,j,k  ,YDIR)-qflx_lo(i  ,j-1,k,YDIR) ) * RDY &
               + ( qflx_lo(i,j,k+1,ZDIR)-qflx_lo(i  ,j  ,k,ZDIR) ) * RDZ

          !--- first guess of tendency and updated value
          momz_t(i,j,k) = - qdiv                                                      &
                          - (var(i,j,k+1,I_PRES)-var(i,j,k,I_PRES) ) * RDZ            &
                          - GRAV * 0.5D0 * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) )  &
                          - DAMP_alphaw(i,j,k) * (var(i,j,k,I_VELZ)-velz_ref(i,j,k) ) &
                          * 0.5D0 * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) )

          var_l = momz(i,j,k) + dtrk * momz_t(i,j,k)

          ! --- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
          pjmax = max( var(i  ,j  ,k  ,I_MOMZ), &
                       var(i+1,j  ,k  ,I_MOMZ), &
                       var(i-1,j  ,k  ,I_MOMZ), &
                       var(i  ,j+1,k  ,I_MOMZ), &
                       var(i  ,j-1,k  ,I_MOMZ), &
                       var(i  ,j  ,k+1,I_MOMZ), &
                       var(i  ,j  ,k-1,I_MOMZ)  )

          pjmin = min( var(i  ,j  ,k  ,I_MOMZ), &
                       var(i+1,j  ,k  ,I_MOMZ), &
                       var(i-1,j  ,k  ,I_MOMZ), &
                       var(i  ,j+1,k  ,I_MOMZ), &
                       var(i  ,j-1,k  ,I_MOMZ), &
                       var(i  ,j  ,k+1,I_MOMZ), &
                       var(i  ,j  ,k-1,I_MOMZ)  )

          ! --- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
          pjpls = max( 0.D0, qflx_anti(i-1,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j-1,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k+1,ZDIR) )
          pjmns = max( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i-1,j  ,k  ,XDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j-1,k  ,YDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k+1,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) )

          ! --- incoming fluxes at scalar grid points ---
          if ( pjpls > 0 ) then
             midvel = 0.25D0 * ( var(i  ,j  ,k+1,I_VELX) &
                               + var(i  ,j  ,k  ,I_VELX) &
                               + var(i-1,j  ,k+1,I_VELX) &
                               + var(i-1,j  ,k  ,I_VELX) )
             rjpls(i,j,k,XDIR) = (pjmax-var_l) / pjpls * abs(midvel)

             midvel = 0.25D0 * ( var(i  ,j  ,k+1,I_VELY) &
                               + var(i  ,j  ,k  ,I_VELY) &
                               + var(i  ,j-1,k+1,I_VELY) &
                               + var(i  ,j-1,k  ,I_VELY) )
             rjpls(i,j,k,YDIR) = (pjmax-var_l) / pjpls * abs(midvel)

             rjpls(i,j,k,ZDIR) = (pjmax-var_l) / pjpls * abs(var(i,j,k,I_VELZ))
          else
             rjpls(i,j,k,XDIR:ZDIR) = 0.D0
          endif

          ! --- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             midvel = 0.25D0 * ( var(i  ,j  ,k+1,I_VELX) &
                               + var(i  ,j  ,k  ,I_VELX) &
                               + var(i-1,j  ,k+1,I_VELX) &
                               + var(i-1,j  ,k  ,I_VELX) )
             rjmns(i,j,k,XDIR) = (var_l-pjmin) / pjmns * abs(midvel)

             midvel = 0.25D0 * ( var(i  ,j  ,k+1,I_VELY) &
                               + var(i  ,j  ,k  ,I_VELY) &
                               + var(i  ,j-1,k+1,I_VELY) &
                               + var(i  ,j-1,k  ,I_VELY) )
             rjmns(i,j,k,YDIR) = (var_l-pjmin) / pjmns * abs(midvel)

             rjmns(i,j,k,ZDIR) = (var_l-pjmin) / pjmns * abs(var(i,j,k,I_VELZ))
           else
             rjmns(i,j,k,XDIR:ZDIR) = 0.D0
         endif

       enddo
       enddo
       enddo

       ! --- [STEP 7S] limit the antidiffusive flux at velocity grid points ---
       do k = WS+1, WE
       do j = JS-1, JE
       do i = IS-1, IE
          if ( qflx_anti(i,j,k,XDIR) >= 0 ) then
             qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i+1,j,k,XDIR), rjmns(i  ,j,k,XDIR), 1.D0 )
          else
             qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i  ,j,k,XDIR), rjmns(i+1,j,k,XDIR), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,YDIR) >= 0 ) then
             qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j+1,k,YDIR), rjmns(i,j  ,k,YDIR), 1.D0 )
          else
             qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j  ,k,YDIR), rjmns(i,j+1,k,YDIR), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,ZDIR) >= 0 ) then
             qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k  ,ZDIR), rjmns(i,j,k-1,ZDIR), 1.D0 )
          else
             qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k-1,ZDIR), rjmns(i,j,k  ,ZDIR), 1.D0 )
          endif
       enddo
       enddo
       enddo
        
       !--- modify flux-divergence->tendency with antidiffusive fluxes
       do k = WS+1, WE-1
       do j = JS,   JE
       do i = IS,   IE
          momz_t(i,j,k) = momz_t(i,j,k) - ( ( qflx_anti(i,j,k  ,XDIR)-qflx_anti(i-1,j  ,k,XDIR) ) * RDX &
                                          + ( qflx_anti(i,j,k  ,YDIR)-qflx_anti(i  ,j-1,k,YDIR) ) * RDY &
                                          + ( qflx_anti(i,j,k+1,ZDIR)-qflx_anti(i  ,j  ,k,ZDIR) ) * RDZ )
       enddo
       enddo
       enddo

       !##### Thermodynamic Equation #####

       ! -- make mass fluxes ( scalar value * mass flux )

       do k = KS,   KE
       do j = JS,   JE
       do i = IS-1, IE
          qflx_lo(i,j,k,XDIR) = 0.5D0                                                                  &
                              * (     mflx_lo(i,j,k,XDIR)  * ( var(i+1,j,k,I_POTT)+var(i,j,k,I_POTT) ) &
                                - abs(mflx_lo(i,j,k,XDIR)) * ( var(i+1,j,k,I_POTT)-var(i,j,k,I_POTT) ) )

          qflx_hi(i,j,k,XDIR) = 0.5D0 * mflx_hi(i,j,k,XDIR)                            &
                              * ( FACT_N * ( var(i+1,j,k,I_POTT)+var(i  ,j,k,I_POTT) ) &
                                + FACT_F * ( var(i+2,j,k,I_POTT)+var(i-1,j,k,I_POTT) ) )

          qflx_anti(i,j,k,XDIR) = qflx_hi(i,j,k,XDIR) - qflx_lo(i,j,k,XDIR)
       enddo
       enddo
       enddo

       do k = KS,   KE
       do j = JS-1, JE
       do i = IS,   IE
          qflx_lo(i,j,k,YDIR) = 0.5D0                                                                  &
                              * (     mflx_lo(i,j,k,YDIR)  * ( var(i,j+1,k,I_POTT)+var(i,j,k,I_POTT) ) &
                                - abs(mflx_lo(i,j,k,YDIR)) * ( var(i,j+1,k,I_POTT)-var(i,j,k,I_POTT) ) )

          qflx_hi(i,j,k,YDIR) = 0.5D0 * mflx_hi(i,j,k,YDIR)                            &
                              * ( FACT_N * ( var(i,j+1,k,I_POTT)+var(i,j  ,k,I_POTT) ) &
                                + FACT_F * ( var(i,j+2,k,I_POTT)+var(i,j-1,k,I_POTT) ) )

          qflx_anti(i,j,k,YDIR) = qflx_hi(i,j,k,YDIR) - qflx_lo(i,j,k,YDIR)
       enddo
       enddo
       enddo

       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          qflx_lo(i,j,k,ZDIR) = 0.5D0                                                                  &
                              * (     mflx_lo(i,j,k,ZDIR)  * ( var(i,j,k+1,I_POTT)+var(i,j,k,I_POTT) ) &
                                - abs(mflx_lo(i,j,k,ZDIR)) * ( var(i,j,k+1,I_POTT)-var(i,j,k,I_POTT) ) )

          qflx_hi(i,j,k,ZDIR) = 0.5D0 * mflx_hi(i,j,k,ZDIR)                              &
                              * ( FACT_N * ( var(i,j,k+1,I_POTT)+var(i,j,k  ,I_POTT) ) &
                                + FACT_F * ( var(i,j,k+2,I_POTT)+var(i,j,k-1,I_POTT) ) )

          qflx_anti(i,j,k,ZDIR) = qflx_hi(i,j,k,ZDIR) - qflx_lo(i,j,k,ZDIR)
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_lo(i,j,WS  ,ZDIR) = 0.D0                                          ! bottom boundary
          qflx_lo(i,j,WS+1,ZDIR) = 0.5D0 * mflx_lo(i,j,WS+1,ZDIR)              &
                                 * ( var(i,j,WS+2,I_POTT)+var(i,j,WS+1,I_POTT) ) ! just above the bottom boundary
          qflx_lo(i,j,WE-1,ZDIR) = 0.5D0 * mflx_lo(i,j,WE-1,ZDIR)              &
                                 * ( var(i,j,WE  ,I_POTT)+var(i,j,WE-1,I_POTT) ) ! just below the top boundary
          qflx_lo(i,j,WE  ,ZDIR) = 0.D0                                          ! top boundary

          qflx_hi(i,j,WS  ,ZDIR) = 0.D0 
          qflx_hi(i,j,WS+1,ZDIR) = qflx_lo(i,j,WS+1,ZDIR)
          qflx_hi(i,j,WE-1,ZDIR) = qflx_lo(i,j,WE-1,ZDIR)
          qflx_hi(i,j,WE  ,ZDIR) = 0.D0 

          qflx_anti(i,j,WS  ,ZDIR) = 0.D0
          qflx_anti(i,j,WS+1,ZDIR) = 0.D0
          qflx_anti(i,j,WE-1,ZDIR) = 0.D0
          qflx_anti(i,j,WE  ,ZDIR) = 0.D0
       enddo
       enddo

       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          ! -- update flux-divergence with the monotone scheme
          qdiv = ( qflx_lo(i,j,k,XDIR)-qflx_lo(i-1,j,  k  ,XDIR) ) * RDX &
               + ( qflx_lo(i,j,k,YDIR)-qflx_lo(i,  j-1,k  ,YDIR) ) * RDY &
               + ( qflx_lo(i,j,k,ZDIR)-qflx_lo(i,  j,  k-1,ZDIR) ) * RDZ

          !--- first guess of tendency and updated value
          pott_t(i,j,k) = ( - qdiv + var(i,j,k,I_POTT)*ddiv(i,j,k) ) / var(i,j,k,I_DENS) &
                          - DAMP_alphat(i,j,k) * (var(i,j,k,I_POTT)-pott_ref(i,j,k) )

          var_l = pott(i,j,k) + dtrk * pott_t(i,j,k)

          ! --- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
          pjmax = max( var(i  ,j  ,k  ,I_POTT), &
                       var(i+1,j  ,k  ,I_POTT), &
                       var(i-1,j  ,k  ,I_POTT), &
                       var(i  ,j+1,k  ,I_POTT), &
                       var(i  ,j-1,k  ,I_POTT), &
                       var(i  ,j  ,k+1,I_POTT), &
                       var(i  ,j  ,k-1,I_POTT)  )

          pjmin = min( var(i  ,j  ,k  ,I_POTT), &
                       var(i+1,j  ,k  ,I_POTT), &
                       var(i-1,j  ,k  ,I_POTT), &
                       var(i  ,j+1,k  ,I_POTT), &
                       var(i  ,j-1,k  ,I_POTT), &
                       var(i  ,j  ,k+1,I_POTT), &
                       var(i  ,j  ,k-1,I_POTT)  )

          ! --- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
          pjpls = max( 0.D0, qflx_anti(i-1,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j-1,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k-1,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) )
          pjmns = max( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i-1,j  ,k  ,XDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j-1,k  ,YDIR) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k-1,ZDIR) )

          ! --- incoming fluxes ---
          if ( pjpls > 0 ) then
             rjpls(i,j,k,XDIR) = (pjmax-var_l) / pjpls * abs((mflx_lo(i,j,k,XDIR)+mflx_lo(i-1,j  ,k  ,XDIR)) * 0.5D0)
             rjpls(i,j,k,YDIR) = (pjmax-var_l) / pjpls * abs((mflx_lo(i,j,k,YDIR)+mflx_lo(i  ,j-1,k  ,YDIR)) * 0.5D0)
             rjpls(i,j,k,ZDIR) = (pjmax-var_l) / pjpls * abs((mflx_lo(i,j,k,ZDIR)+mflx_lo(i  ,j  ,k-1,ZDIR)) * 0.5D0)
          else
             rjpls(i,j,k,XDIR:ZDIR) = 0.D0
          endif

          ! --- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             rjmns(i,j,k,XDIR) = (var_l-pjmin) / pjmns * abs((mflx_lo(i,j,k,XDIR)+mflx_lo(i-1,j  ,k  ,XDIR)) * 0.5D0)
             rjmns(i,j,k,YDIR) = (var_l-pjmin) / pjmns * abs((mflx_lo(i,j,k,YDIR)+mflx_lo(i  ,j-1,k  ,YDIR)) * 0.5D0)
             rjmns(i,j,k,ZDIR) = (var_l-pjmin) / pjmns * abs((mflx_lo(i,j,k,ZDIR)+mflx_lo(i  ,j  ,k-1,ZDIR)) * 0.5D0)
          else
             rjmns(i,j,k,XDIR:ZDIR) = 0.D0
          endif

       enddo
       enddo
       enddo

       ! --- [STEP 7S] limit the antidiffusive flux ---
       do k = KS-1, KE 
       do j = JS-1, JE 
       do i = IS-1, IE
          if ( qflx_anti(i,j,k,XDIR) >= 0 ) then
             qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i+1,j,k,XDIR), rjmns(i  ,j,k,XDIR), 1.D0 )
          else
             qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i  ,j,k,XDIR), rjmns(i+1,j,k,XDIR), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,YDIR) >= 0 ) then
             qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j+1,k,YDIR), rjmns(i,j  ,k,YDIR), 1.D0 )
          else
             qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j  ,k,YDIR), rjmns(i,j+1,k,YDIR), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,ZDIR) >= 0 ) then
             qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k+1,ZDIR), rjmns(i,j,k  ,ZDIR), 1.D0 )
          else
             qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k  ,ZDIR), rjmns(i,j,k+1,ZDIR), 1.D0 )
          endif
       enddo
       enddo
       enddo

       !--- modify advection->tendency with antidiffusive fluxes
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          pott_t(i,j,k) = pott_t(i,j,k)                                                 &
                        - ( ( qflx_anti(i,j,k,XDIR)-qflx_anti(i-1,j  ,k  ,XDIR) ) * RDX &
                          + ( qflx_anti(i,j,k,YDIR)-qflx_anti(i  ,j-1,k  ,YDIR) ) * RDY &
                          + ( qflx_anti(i,j,k,ZDIR)-qflx_anti(i  ,j  ,k-1,ZDIR) ) * RDZ &
                          ) / var(i,j,k,I_DENS)
       enddo
       enddo
       enddo

       if ( rko == RK ) then ! do only at last step

          !##### advection of scalar quantity #####

          do iq = 1, QA

             ! -- make mass fluxes ( scalar value * mass flux )

             do k = KS,   KE
             do j = JS,   JE
             do i = IS-1, IE
                qflx_lo(i,j,k,XDIR) = 0.5D0                                                            &
                                    * (     mflx_lo(i,j,k,XDIR)  * ( qtrc(i+1,j,k,iq)+qtrc(i,j,k,iq) ) &
                                      - abs(mflx_lo(i,j,k,XDIR)) * ( qtrc(i+1,j,k,iq)-qtrc(i,j,k,iq) ) )

                qflx_hi(i,j,k,XDIR) = 0.5D0 * mflx_hi(i,j,k,XDIR)                      &
                                    * ( FACT_N * ( qtrc(i+1,j,k,iq)+qtrc(i  ,j,k,iq) ) &
                                      + FACT_F * ( qtrc(i+2,j,k,iq)+qtrc(i-1,j,k,iq) ) )

                qflx_anti(i,j,k,XDIR) = qflx_hi(i,j,k,XDIR) - qflx_lo(i,j,k,XDIR)
             enddo
             enddo
             enddo

             do k = KS,   KE
             do j = JS-1, JE
             do i = IS,   IE
                qflx_lo(i,j,k,YDIR) = 0.5D0                                                            &
                                    * (     mflx_lo(i,j,k,YDIR)  * ( qtrc(i,j+1,k,iq)+qtrc(i,j,k,iq) ) &
                                      - abs(mflx_lo(i,j,k,YDIR)) * ( qtrc(i,j+1,k,iq)-qtrc(i,j,k,iq) ) )

                qflx_hi(i,j,k,YDIR) = 0.5D0 * mflx_hi(i,j,k,YDIR)                      &
                                    * ( FACT_N * ( qtrc(i,j+1,k,iq)+qtrc(i,j  ,k,iq) ) &
                                      + FACT_F * ( qtrc(i,j+2,k,iq)+qtrc(i,j-1,k,iq) ) )

                qflx_anti(i,j,k,YDIR) = qflx_hi(i,j,k,YDIR) - qflx_lo(i,j,k,YDIR)
             enddo
             enddo
             enddo

             do k = WS+2, WE-2
             do j = JS,   JE
             do i = IS,   IE
                qflx_lo(i,j,k,ZDIR) = 0.5D0                                                            &
                                    * (     mflx_lo(i,j,k,ZDIR)  * ( qtrc(i,j,k+1,iq)+qtrc(i,j,k,iq) ) &
                                      - abs(mflx_lo(i,j,k,ZDIR)) * ( qtrc(i,j,k+1,iq)-qtrc(i,j,k,iq) ) )

                qflx_hi(i,j,k,ZDIR) = 0.5D0 * mflx_hi(i,j,k,ZDIR)                      &
                                    * ( FACT_N * ( qtrc(i,j,k+1,iq)+qtrc(i,j,k  ,iq) ) &
                                      + FACT_F * ( qtrc(i,j,k+2,iq)+qtrc(i,j,k-1,iq) ) )

                qflx_anti(i,j,k,ZDIR) = qflx_hi(i,j,k,ZDIR) - qflx_lo(i,j,k,ZDIR)
             enddo
             enddo
             enddo

             do j = JS, JE
             do i = IS, IE
                qflx_lo(i,j,WS  ,ZDIR) = 0.D0                                    ! bottom boundary
                qflx_lo(i,j,WS+1,ZDIR) = 0.5D0 * mflx_lo(i,j,WS+1,ZDIR)        &
                                       * ( qtrc(i,j,WS+2,iq)+qtrc(i,j,WS+1,iq) ) ! just above the bottom boundary
                qflx_lo(i,j,WE-1,ZDIR) = 0.5D0 * mflx_lo(i,j,WE-1,ZDIR)        &
                                       * ( qtrc(i,j,WE  ,iq)+qtrc(i,j,WE-1,iq) ) ! just below the top boundary
                qflx_lo(i,j,WE  ,ZDIR) = 0.D0                                    ! top boundary

                qflx_hi(i,j,WS  ,ZDIR) = 0.D0 
                qflx_hi(i,j,WS+1,ZDIR) = qflx_lo(i,j,WS+1,ZDIR)
                qflx_hi(i,j,WE-1,ZDIR) = qflx_lo(i,j,WE-1,ZDIR)
                qflx_hi(i,j,WE  ,ZDIR) = 0.D0 

                qflx_anti(i,j,WS  ,ZDIR) = 0.D0
                qflx_anti(i,j,WS+1,ZDIR) = 0.D0
                qflx_anti(i,j,WE-1,ZDIR) = 0.D0
                qflx_anti(i,j,WE  ,ZDIR) = 0.D0
             enddo
             enddo

             do k = KS, KE
             do j = JS, JE
             do i = IS, IE
                ! -- update flux-divergence with the monotone scheme
                qdiv = ( qflx_lo(i,j,k,XDIR)-qflx_lo(i-1,j,  k  ,XDIR) ) * RDX &
                     + ( qflx_lo(i,j,k,YDIR)-qflx_lo(i,  j-1,k  ,YDIR) ) * RDY &
                     + ( qflx_lo(i,j,k,ZDIR)-qflx_lo(i,  j,  k-1,ZDIR) ) * RDZ

                !--- first guess of tendency and updated value
                qtrc_t(i,j,k,iq) = ( - qdiv + qtrc(i,j,k,iq)*ddiv(i,j,k) ) / var(i,j,k,I_DENS)

                var_l = qtrc(i,j,k,iq) + dtrk * qtrc_t(i,j,k,iq)

                ! --- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
                pjmax = max( qtrc(i  ,j  ,k  ,iq), &
                             qtrc(i+1,j  ,k  ,iq), &
                             qtrc(i-1,j  ,k  ,iq), &
                             qtrc(i  ,j+1,k  ,iq), &
                             qtrc(i  ,j-1,k  ,iq), &
                             qtrc(i  ,j  ,k+1,iq), &
                             qtrc(i  ,j  ,k-1,iq)  )

                pjmin = min( qtrc(i  ,j  ,k  ,iq), &
                             qtrc(i+1,j  ,k  ,iq), &
                             qtrc(i-1,j  ,k  ,iq), &
                             qtrc(i  ,j+1,k  ,iq), &
                             qtrc(i  ,j-1,k  ,iq), &
                             qtrc(i  ,j  ,k+1,iq), &
                             qtrc(i  ,j  ,k-1,iq)  )

                ! --- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
                pjpls = max( 0.D0, qflx_anti(i-1,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) &
                      + max( 0.D0, qflx_anti(i  ,j-1,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) &
                      + max( 0.D0, qflx_anti(i  ,j  ,k-1,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) )
                pjmns = max( 0.D0, qflx_anti(i  ,j  ,k  ,XDIR) ) - min( 0.D0, qflx_anti(i-1,j  ,k  ,XDIR) ) &
                      + max( 0.D0, qflx_anti(i  ,j  ,k  ,YDIR) ) - min( 0.D0, qflx_anti(i  ,j-1,k  ,YDIR) ) &
                      + max( 0.D0, qflx_anti(i  ,j  ,k  ,ZDIR) ) - min( 0.D0, qflx_anti(i  ,j  ,k-1,ZDIR) )

                ! --- incoming fluxes ---
                if ( pjpls > 0 ) then
                   rjpls(i,j,k,XDIR) = (pjmax-var_l) / pjpls * abs((mflx_lo(i,j,k,XDIR)+mflx_lo(i-1,j  ,k  ,XDIR)) * 0.5D0)
                   rjpls(i,j,k,YDIR) = (pjmax-var_l) / pjpls * abs((mflx_lo(i,j,k,YDIR)+mflx_lo(i  ,j-1,k  ,YDIR)) * 0.5D0)
                   rjpls(i,j,k,ZDIR) = (pjmax-var_l) / pjpls * abs((mflx_lo(i,j,k,ZDIR)+mflx_lo(i  ,j  ,k-1,ZDIR)) * 0.5D0)
                else
                   rjpls(i,j,k,XDIR:ZDIR) = 0.D0
                endif

                ! --- outgoing fluxes at scalar grid points ---
                if ( pjmns > 0 ) then
                   rjmns(i,j,k,XDIR) = (var_l-pjmin) / pjmns * abs((mflx_lo(i,j,k,XDIR)+mflx_lo(i-1,j  ,k  ,XDIR)) * 0.5D0)
                   rjmns(i,j,k,YDIR) = (var_l-pjmin) / pjmns * abs((mflx_lo(i,j,k,YDIR)+mflx_lo(i  ,j-1,k  ,YDIR)) * 0.5D0)
                   rjmns(i,j,k,ZDIR) = (var_l-pjmin) / pjmns * abs((mflx_lo(i,j,k,ZDIR)+mflx_lo(i  ,j  ,k-1,ZDIR)) * 0.5D0)
                else
                   rjmns(i,j,k,XDIR:ZDIR) = 0.D0
                endif

             enddo
             enddo
             enddo

             ! --- [STEP 7S] limit the antidiffusive flux ---
             do k = KS-1, KE 
             do j = JS-1, JE 
             do i = IS-1, IE
                if ( qflx_anti(i,j,k,XDIR) >= 0 ) then
                   qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i+1,j,k,XDIR), rjmns(i  ,j,k,XDIR), 1.D0 )
                else
                   qflx_anti(i,j,k,XDIR) = qflx_anti(i,j,k,XDIR) * min( rjpls(i  ,j,k,XDIR), rjmns(i+1,j,k,XDIR), 1.D0 )
                endif

                if ( qflx_anti(i,j,k,YDIR) >= 0 ) then
                   qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j+1,k,YDIR), rjmns(i,j  ,k,YDIR), 1.D0 )
                else
                   qflx_anti(i,j,k,YDIR) = qflx_anti(i,j,k,YDIR) * min( rjpls(i,j  ,k,YDIR), rjmns(i,j+1,k,YDIR), 1.D0 )
                endif

                if ( qflx_anti(i,j,k,ZDIR) >= 0 ) then
                   qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k+1,ZDIR), rjmns(i,j,k  ,ZDIR), 1.D0 )
                else
                   qflx_anti(i,j,k,ZDIR) = qflx_anti(i,j,k,ZDIR) * min( rjpls(i,j,k  ,ZDIR), rjmns(i,j,k+1,ZDIR), 1.D0 )
                endif
             enddo
             enddo
             enddo

             !--- modify advection->tendency with antidiffusive fluxes
             do k = KS, KE
             do j = JS, JE
             do i = IS, IE
                qtrc_t(i,j,k,iq) = qtrc_t(i,j,k,iq)                                              &
                                 - ( ( qflx_anti(i,j,k,XDIR)-qflx_anti(i-1,j  ,k  ,XDIR) ) * RDX &
                                   + ( qflx_anti(i,j,k,YDIR)-qflx_anti(i  ,j-1,k  ,YDIR) ) * RDY &
                                   + ( qflx_anti(i,j,k,ZDIR)-qflx_anti(i  ,j  ,k-1,ZDIR) ) * RDZ &
                                   ) / var(i,j,k,I_DENS)
             enddo
             enddo
             enddo

          enddo ! scalar quantities loop

       else

          !##### RK time integration #####
          do k = KS, KE
          do j = JS, JE
          do i = IS, IE
             var(i,j,k,I_DENS) = dens(i,j,k) + dtrk * dens_t(i,j,k)
             var(i,j,k,I_POTT) = pott(i,j,k) + dtrk * pott_t(i,j,k)
          enddo
          enddo
          enddo

          do k = WS-1, WE+1
          do j = JS, JE
          do i = IS, IE
             var(i,j,k,I_MOMX) = momx(i,j,k) + dtrk * momx_t(i,j,k)
             var(i,j,k,I_MOMY) = momy(i,j,k) + dtrk * momy_t(i,j,k)
             var(i,j,k,I_MOMZ) = momz(i,j,k) + dtrk * momz_t(i,j,k)
          enddo
          enddo
          enddo

          ! fill IHALO & JHALO
          call COMM_vars( var(:,:,:,1:5) )

       endif

    enddo ! RK loop

    return
  end subroutine ATMOS_DYN

end module mod_atmos_dyn
