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
  integer, private, parameter :: I_DENS =  1
  integer, private, parameter :: I_MOMX =  2
  integer, private, parameter :: I_MOMY =  3
  integer, private, parameter :: I_MOMZ =  4
  integer, private, parameter :: I_RHOT =  5
  integer, private, parameter :: I_PRES =  6
  integer, private, parameter :: I_VELX =  7
  integer, private, parameter :: I_VELY =  8
  integer, private, parameter :: I_VELZ =  9
  integer, private, parameter :: I_POTT = 10

  integer, private, parameter :: XDIR   = 1
  integer, private, parameter :: YDIR   = 2
  integer, private, parameter :: ZDIR   = 3

  ! time settings
  integer, private, parameter :: RK = 3                 ! order of Runge-Kutta scheme

  ! advection settings
  real(8), private, parameter :: FACT_N =   7.D0 / 6.D0 !  7/6: fourth, 1: second
  real(8), private, parameter :: FACT_F = - 1.D0 / 6.D0 ! -1/6: fourth, 0: second

  ! numerical filter settings
  integer, private, parameter :: DF     = 4             ! order of numerical filter

  real(8), private, save      :: ATMOS_DYN_numerical_diff = 1.D-3 ! nondimensional numerical diffusion
  real(8), private, save      :: DIFF

  ! work
  real(8), private, allocatable, save :: var      (:,:,:,:) ! work
  real(8), private, allocatable, save :: rhot_t   (:,:,:)   ! tendency  of rho * theta
  real(8), private, allocatable, save :: dens_diff(:,:,:)   ! deviation of rho
  real(8), private, allocatable, save :: rhot_diff(:,:,:)   ! deviation of rho * theta
  ! mass flux
  real(8), private, allocatable, save :: mflx_hi  (:,:,:,:) ! rho * vel(x,y,z)       @ (u,v,w)-face high order
  real(8), private, allocatable, save :: mflx_lo  (:,:,:,:) ! rho * vel(x,y,z)       @ (u,v,w)-face low  order
  real(8), private, allocatable, save :: qflx_hi  (:,:,:,:) ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
  real(8), private, allocatable, save :: qflx_lo  (:,:,:,:) ! rho * vel(x,y,z) * phi @ (u,v,w)-face low  order
  real(8), private, allocatable, save :: qflx_anti(:,:,:,:) ! rho * vel(x,y,z) * phi @ (u,v,w)-face antidiffusive
  ! flux correction term
  real(8), private, allocatable, save :: rjpls    (:,:,:,:) ! plus  in (x,y,z)-direction
  real(8), private, allocatable, save :: rjmns    (:,:,:,:) ! minus in (x,y,z)-direction

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
       IA  => GRID_IA, &
       JA  => GRID_JA, &
       KA  => GRID_KA
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

    allocate( var(IA,JA,KA,I_DENS:I_POTT) )

    allocate( rhot_t(IA,JA,KA)    )

    allocate( dens_diff(IA,JA,KA) )
    allocate( rhot_diff(IA,JA,KA) )

    allocate( mflx_lo  (IA,JA,KA,XDIR:ZDIR) )
    allocate( mflx_hi  (IA,JA,KA,XDIR:ZDIR) )
    allocate( qflx_lo  (IA,JA,KA,XDIR:ZDIR) )
    allocate( qflx_hi  (IA,JA,KA,XDIR:ZDIR) )
    allocate( qflx_anti(IA,JA,KA,XDIR:ZDIR) )

    allocate( rjpls    (IA,JA,KA,XDIR:ZDIR) )
    allocate( rjmns    (IA,JA,KA,XDIR:ZDIR) )

    return
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
       REF_dens, &
       REF_pott
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
    real(8), intent(in)    :: dens(IA,JA,KA)      ! density [kg/m3]
    real(8), intent(in)    :: momx(IA,JA,KA)      ! momentum (x) [kg/m3 * m/s]
    real(8), intent(in)    :: momy(IA,JA,KA)      ! momentum (y) [kg/m3 * m/s]
    real(8), intent(in)    :: momz(IA,JA,KA)      ! momentum (z) [kg/m3 * m/s]
    real(8), intent(in)    :: pott(IA,JA,KA)      ! potential temperature [K]
    real(8), intent(in)    :: qtrc(IA,JA,KA,QA)   ! tracer mixing ratio   [kg/kg],[1/m3]
    ! prognostic tendency
    real(8), intent(inout) :: dens_t(IA,JA,KA)
    real(8), intent(inout) :: momx_t(IA,JA,KA)
    real(8), intent(inout) :: momy_t(IA,JA,KA)
    real(8), intent(inout) :: momz_t(IA,JA,KA)
    real(8), intent(inout) :: pott_t(IA,JA,KA)
    real(8), intent(inout) :: qtrc_t(IA,JA,KA,QA)

    real(8) :: dtrk ! dt for each RK step
    real(8) :: ddiv ! density divergence [kg/m3/s]
    real(8) :: qdiv ! flux divergence for any quantity

    ! FCT work
    real(8) :: pjpls, pjmns, pjmax, pjmin, var_l

    integer :: i, j, k, iq, rko
    !---------------------------------------------------------------------------

    rhot_t(:,:,:) = 0.D0
    dtrk = 0.D0

    do rko = 1, RK

       !##### RK time integration for previous RK step #####
       ! at first time, this section is used for copying to work from global vars
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          var(i,j,k,I_DENS) = dens(i,j,k) + dtrk * dens_t(i,j,k)

          var(i,j,k,I_RHOT) = pott(i,j,k) * dens(i,j,k) + dtrk * rhot_t(i,j,k)

          var(i,j,k,I_MOMX) = momx(i,j,k) + dtrk * momx_t(i,j,k)
          var(i,j,k,I_MOMY) = momy(i,j,k) + dtrk * momy_t(i,j,k)
       enddo
       enddo
       enddo

       do k = WS, WE
       do j = JS, JE
       do i = IS, IE
          var(i,j,k,I_MOMZ) = momz(i,j,k) + dtrk * momz_t(i,j,k)
       enddo
       enddo
       enddo

       ! fill IHALO & JHALO
       call COMM_vars( var(:,:,:,I_DENS:I_RHOT) )


       !##### calc pressure, velocity & communication
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

       ! pressure, potential temperature
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          var(i,j,k,I_PRES) = Pstd * ( var(i,j,k,I_RHOT) * Rair / Pstd )**CPovCV
          var(i,j,k,I_POTT) = var(i,j,k,I_RHOT) / var(i,j,k,I_DENS)
       enddo
       enddo
       enddo

       ! fill IHALO & JHALO
       call COMM_vars( var(:,:,:,I_PRES:I_POTT) )



       !##### Update RK dt #####
       dtrk = TIME_DTSEC_ATMOS_DYN / (RK - rko + 1)


       !##### continuity equation #####
       ! at (u, y, layer)
       do k = KS,   KE
       do j = JS,   JE
       do i = IS-1, IE
          mflx_hi(i,j,k,XDIR) = 0.5D0 * var(i,j,k,I_VELX)                              &
                              * ( FACT_N * ( var(i+1,j,k,I_DENS)+var(i  ,j,k,I_DENS) ) &
                                + FACT_F * ( var(i+2,j,k,I_DENS)+var(i-1,j,k,I_DENS) ) )
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do k = KS,   KE
       do j = JS-1, JE
       do i = IS,   IE
          mflx_hi(i,j,k,YDIR) = 0.5D0 * var(i,j,k,I_VELY)                              &
                              * ( FACT_N * ( var(i,j+1,k,I_DENS)+var(i,j  ,k,I_DENS) ) &
                                + FACT_F * ( var(i,j+2,k,I_DENS)+var(i,j-1,k,I_DENS) ) )
       enddo
       enddo
       enddo
       ! at (x, y, interface)
       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          mflx_hi(i,j,k,ZDIR) = 0.5D0 * var(i,j,k,I_VELZ)                              &
                              * ( FACT_N * ( var(i,j,k+1,I_DENS)+var(i,j,k  ,I_DENS) ) &
                                + FACT_F * ( var(i,j,k+2,I_DENS)+var(i,j,k-1,I_DENS) ) )
       enddo
       enddo
       enddo
       do j = JS, JE
       do i = IS, IE
          mflx_hi(i,j,WS  ,ZDIR) = 0.D0
          mflx_hi(i,j,WS+1,ZDIR) = 0.5D0 * var(i,j,WS+1,I_VELZ)                &
                                 * ( var(i,j,WS+2,I_DENS)+var(i,j,WS+1,I_DENS) ) ! just above the bottom boundary
          mflx_hi(i,j,WE-1,ZDIR) = 0.5D0 * var(i,j,WE-1,I_VELZ)                &
                                 * ( var(i,j,WE  ,I_DENS)+var(i,j,WE-1,I_DENS) ) ! just below the top boundary
          mflx_hi(i,j,WE  ,ZDIR) = 0.D0
       enddo
       enddo

       !--- use difference for numerical filter
       do k = KS, KE
          dens_diff(:,:,k) = var(:,:,k,I_DENS) - REF_dens(:,:,k)
       enddo

       !--- calc tendency with numerical filter
       do k = KS+2, KE-2
       do j = JS,   JE
       do i = IS,   IE
          dens_t(i,j,k) = - ( ( mflx_hi(i,j,k,XDIR)-mflx_hi(i-1,j,  k  ,XDIR) ) * RDX   &
                            + ( mflx_hi(i,j,k,YDIR)-mflx_hi(i,  j-1,k  ,YDIR) ) * RDY   &
                            + ( mflx_hi(i,j,k,ZDIR)-mflx_hi(i,  j,  k-1,ZDIR) ) * RDZ ) & ! density divergence
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
                                            +        dens_diff(i  ,j  ,k-2) ) ) ! numerical filter
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          k = KS
          dens_t(i,j,k) = - ( ( mflx_hi(i,j,k,XDIR)-mflx_hi(i-1,j,  k  ,XDIR) ) * RDX   &
                            + ( mflx_hi(i,j,k,YDIR)-mflx_hi(i,  j-1,k  ,YDIR) ) * RDY   &
                            + ( mflx_hi(i,j,k,ZDIR)-mflx_hi(i,  j,  k-1,ZDIR) ) * RDZ ) &
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
          dens_t(i,j,k) = - ( ( mflx_hi(i,j,k,XDIR)-mflx_hi(i-1,j,  k  ,XDIR) ) * RDX   &
                            + ( mflx_hi(i,j,k,YDIR)-mflx_hi(i,  j-1,k  ,YDIR) ) * RDY   &
                            + ( mflx_hi(i,j,k,ZDIR)-mflx_hi(i,  j,  k-1,ZDIR) ) * RDZ ) &
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
          dens_t(i,j,k) = - ( ( mflx_hi(i,j,k,XDIR)-mflx_hi(i-1,j,  k  ,XDIR) ) * RDX   &
                            + ( mflx_hi(i,j,k,YDIR)-mflx_hi(i,  j-1,k  ,YDIR) ) * RDY   &
                            + ( mflx_hi(i,j,k,ZDIR)-mflx_hi(i,  j,  k-1,ZDIR) ) * RDZ ) &
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
          dens_t(i,j,k) = - ( ( mflx_hi(i,j,k,XDIR)-mflx_hi(i-1,j,  k  ,XDIR) ) * RDX   &
                            + ( mflx_hi(i,j,k,YDIR)-mflx_hi(i,  j-1,k  ,YDIR) ) * RDY   &
                            + ( mflx_hi(i,j,k,ZDIR)-mflx_hi(i,  j,  k-1,ZDIR) ) * RDZ ) &
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

       !--- x-direction momentum equation

       ! at (x, y, layer)
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE+1
          qflx_hi(i,j,k,XDIR) = 0.25D0 * ( var(i,j,k,I_VELX)+var(i-1,j,k,I_VELX) )     &
                              * ( FACT_N * ( var(i  ,j,k,I_MOMX)+var(i-1,j,k,I_MOMX) ) &
                                + FACT_F * ( var(i+1,j,k,I_MOMX)+var(i-2,j,k,I_MOMX) ) )
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       do k = KS,   KE
       do j = JS-1, JE
       do i = IS,   IE
          qflx_hi(i,j,k,YDIR) = 0.25D0 * ( var(i+1,j,k,I_VELY)+var(i,j,k,I_VELY) )     &
                              * ( FACT_N * ( var(i,j+1,k,I_MOMX)+var(i,j  ,k,I_MOMX) ) &
                                + FACT_F * ( var(i,j+2,k,I_MOMX)+var(i,j-1,k,I_MOMX) ) )
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          qflx_hi(i,j,k,ZDIR) = 0.25D0 * ( var(i+1,j,k,I_VELZ)+var(i,j,k,I_VELZ) )     &
                              * ( FACT_N * ( var(i,j,k+1,I_MOMX)+var(i,j,k  ,I_MOMX) ) &
                                + FACT_F * ( var(i,j,k+2,I_MOMX)+var(i,j,k-1,I_MOMX) ) )
       enddo
       enddo
       enddo
       do j = JS, JE
       do i = IS, IE
          qflx_hi(i,j,WS  ,ZDIR) = 0.D0                                                     ! bottom boundary
          qflx_hi(i,j,WS+1,ZDIR) = 0.25D0 * ( var(i+1,j,WS+2,I_VELZ)+var(i,j,WS+1,I_VELZ) ) &
                                 * ( var(i,j,WS+2,I_MOMX)+var(i,j,WS+1,I_MOMX) )            ! just above the bottom boundary
          qflx_hi(i,j,WE-1,ZDIR) = 0.25D0 * ( var(i+1,j,WE  ,I_VELZ)+var(i,j,WE-1,I_VELZ) ) &
                                 * ( var(i,j,WE  ,I_MOMX)+var(i,j,WE-1,I_MOMX) )            ! just below the top boundary
          qflx_hi(i,j,WE  ,ZDIR) = 0.D0                                                     ! top boundary
       enddo
       enddo

       !--- calc tendency with numerical filter
       do k = KS+2, KE-2
       do j = JS,   JE
       do i = IS,   IE
          momx_t(i,j,k) = - ( ( qflx_hi(i+1,j,k,XDIR)-qflx_hi(i,j,  k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i  ,j,k,YDIR)-qflx_hi(i,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i  ,j,k,ZDIR)-qflx_hi(i,j,  k-1,ZDIR) ) * RDZ ) &    ! flux
                          - ( var(i+1,j,k,I_PRES)-var(i,j,k,I_PRES) ) * RDX                & ! Pressure gradient
                          - ( DAMP_alphau(i,j,k) * ( var(i,j,k,I_VELX)-velx_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i+1,j,k,I_DENS)+var(i,j,k,I_DENS) )          ) & ! Damping factor
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMX)   &
                                            +        var(i-2,j  ,k  ,I_MOMX) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMX)   &
                                            +        var(i  ,j-2,k  ,I_MOMX) ) &
                                          + (        var(i  ,j  ,k+2,I_MOMX)   &
                                            - 4.D0 * var(i  ,j  ,k+1,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j  ,k-1,I_MOMX)   &
                                            +        var(i  ,j  ,k-2,I_MOMX) ) ) ! numerical filter
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          k = KS
          momx_t(i,j,k) = - ( ( qflx_hi(i+1,j,k,XDIR)-qflx_hi(i,j,  k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i  ,j,k,YDIR)-qflx_hi(i,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i  ,j,k,ZDIR)-qflx_hi(i,j,  k-1,ZDIR) ) * RDZ ) &
                          - ( var(i+1,j,k,I_PRES)-var(i,j,k,I_PRES) ) * RDX                &
                          - ( DAMP_alphau(i,j,k) * ( var(i,j,k,I_VELX)-velx_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i+1,j,k,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMX)   &
                                            +        var(i-2,j  ,k  ,I_MOMX) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMX)   &
                                            +        var(i  ,j-2,k  ,I_MOMX) ) &
                                          + (        var(i  ,j  ,k+2,I_MOMX)   &
                                            - 3.D0 * var(i  ,j  ,k+1,I_MOMX)   &
                                            + 2.D0 * var(i  ,j  ,k  ,I_MOMX) ) )

          k = KS+1
          momx_t(i,j,k) = - ( ( qflx_hi(i+1,j,k,XDIR)-qflx_hi(i,j,  k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i  ,j,k,YDIR)-qflx_hi(i,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i  ,j,k,ZDIR)-qflx_hi(i,j,  k-1,ZDIR) ) * RDZ ) &
                          - ( var(i+1,j,k,I_PRES)-var(i,j,k,I_PRES) ) * RDX                &
                          - ( DAMP_alphau(i,j,k) * ( var(i,j,k,I_VELX)-velx_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i+1,j,k,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMX)   &
                                            +        var(i-2,j  ,k  ,I_MOMX) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMX)   &
                                            +        var(i  ,j-2,k  ,I_MOMX) ) &
                                          + (        var(i  ,j  ,k+2,I_MOMX)   &
                                            - 4.D0 * var(i  ,j  ,k+1,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 3.D0 * var(i  ,j  ,k-1,I_MOMX) ) )

          k = KE-1
          momx_t(i,j,k) = - ( ( qflx_hi(i+1,j,k,XDIR)-qflx_hi(i,j,  k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i  ,j,k,YDIR)-qflx_hi(i,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i  ,j,k,ZDIR)-qflx_hi(i,j,  k-1,ZDIR) ) * RDZ ) &
                          - ( var(i+1,j,k,I_PRES)-var(i,j,k,I_PRES) ) * RDX                &
                          - ( DAMP_alphau(i,j,k) * ( var(i,j,k,I_VELX)-velx_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i+1,j,k,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMX)   &
                                            +        var(i-2,j  ,k  ,I_MOMX) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMX)   &
                                            +        var(i  ,j-2,k  ,I_MOMX) ) &
                                          + (        var(i  ,j  ,k-2,I_MOMX)   &
                                            - 4.D0 * var(i  ,j  ,k-1,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 3.D0 * var(i  ,j  ,k+1,I_MOMX) ) )

          k = KE
          momx_t(i,j,k) = - ( ( qflx_hi(i+1,j,k,XDIR)-qflx_hi(i,j,  k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i  ,j,k,YDIR)-qflx_hi(i,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i  ,j,k,ZDIR)-qflx_hi(i,j,  k-1,ZDIR) ) * RDZ ) &
                          - ( var(i+1,j,k,I_PRES)-var(i,j,k,I_PRES) ) * RDX                &
                          - ( DAMP_alphau(i,j,k) * ( var(i,j,k,I_VELX)-velx_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i+1,j,k,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMX)   &
                                            +        var(i-2,j  ,k  ,I_MOMX) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMX)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMX)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMX)   &
                                            +        var(i  ,j-2,k  ,I_MOMX) ) &
                                          + (        var(i  ,j  ,k-2,I_MOMX)   &
                                            - 3.D0 * var(i  ,j  ,k-1,I_MOMX)   &
                                            + 2.D0 * var(i  ,j  ,k  ,I_MOMX) ) )
       enddo
       enddo

       !--- y-direction momentum equation

       ! at (u, v, layer)
       do k = KS,   KE
       do j = JS,   JE
       do i = IS-1, IE
          qflx_hi(i,j,k,XDIR) = 0.25D0 * ( var(i,j+1,k,I_VELX)+var(i,j,k,I_VELX) )     &
                              * ( FACT_N * ( var(i+1,j,k,I_MOMY)+var(i  ,j,k,I_MOMY) ) &
                                + FACT_F * ( var(i+2,j,k,I_MOMY)+var(i-1,j,k,I_MOMY) ) )
       enddo
       enddo
       enddo
       ! at (x, y, layer)
       do k = KS, KE
       do j = JS, JE+1
       do i = IS, IE
          qflx_hi(i,j,k,YDIR) = 0.25D0 * ( var(i,j,k,I_VELY)+var(i,j-1,k,I_VELY) )     &
                              * ( FACT_N * ( var(i,j  ,k,I_MOMY)+var(i,j-1,k,I_MOMY) ) &
                                + FACT_F * ( var(i,j+1,k,I_MOMY)+var(i,j-2,k,I_MOMY) ) )
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          qflx_hi(i,j,k,ZDIR) = 0.25D0 * ( var(i,j+1,k,I_VELZ)+var(i,j,k,I_VELZ) )     &
                              * ( FACT_N * ( var(i,j,k+1,I_MOMY)+var(i,j,k  ,I_MOMY) ) &
                                + FACT_F * ( var(i,j,k+2,I_MOMY)+var(i,j,k-1,I_MOMY) ) )
       enddo
       enddo
       enddo
       do j = JS, JE
       do i = IS, IE
          qflx_hi(i,j,WS  ,ZDIR) = 0.D0                                                     ! bottom boundary
          qflx_hi(i,j,WS+1,ZDIR) = 0.25D0 * ( var(i,j+1,WS+2,I_VELZ)+var(i,j,WS+1,I_VELZ) ) &
                                 * ( var(i,j,WS+2,I_MOMY)+var(i,j,WS+1,I_MOMY) )            ! just above the bottom boundary
          qflx_hi(i,j,WE-1,ZDIR) = 0.25D0 * ( var(i,j+1,WE  ,I_VELZ)+var(i,j,WE-1,I_VELZ) ) &
                                  * ( var(i,j,WE  ,I_MOMY)+var(i,j,WE-1,I_MOMY) )           ! just below the top boundary
          qflx_hi(i,j,WE  ,ZDIR) = 0.D0                                                     ! top boundary
       enddo
       enddo

       !--- calc tendency with numerical filter
       do k = KS+2, KE-2
       do j = JS,   JE
       do i = IS,   IE
          momy_t(i,j,k) = - ( ( qflx_hi(i,j  ,k,XDIR)-qflx_hi(i-1,j,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j+1,k,YDIR)-qflx_hi(i  ,j,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j  ,k,ZDIR)-qflx_hi(i  ,j,k-1,ZDIR) ) * RDZ ) &    ! flux
                          - ( var(i,j+1,k,I_PRES)-var(i,j,k,I_PRES) ) * RDY                & ! Pressure gradient
                          - ( DAMP_alphav(i,j,k) * ( var(i,j,k,I_VELY)-vely_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i,j+1,k,I_DENS)+var(i,j,k,I_DENS) )          ) & ! Damping factor
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMY)   &
                                            +        var(i-2,j  ,k  ,I_MOMY) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMY)   &
                                            +        var(i  ,j-2,k  ,I_MOMY) ) &
                                          + (        var(i  ,j  ,k+2,I_MOMY)   &
                                            - 4.D0 * var(i  ,j  ,k+1,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j  ,k-1,I_MOMY)   &
                                            +        var(i  ,j  ,k-2,I_MOMY) ) ) ! numerical filter
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          k = KS
          momy_t(i,j,k) = - ( ( qflx_hi(i,j  ,k,XDIR)-qflx_hi(i-1,j,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j+1,k,YDIR)-qflx_hi(i  ,j,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j  ,k,ZDIR)-qflx_hi(i  ,j,k-1,ZDIR) ) * RDZ ) &
                          - ( var(i,j+1,k,I_PRES)-var(i,j,k,I_PRES) ) * RDY                &
                          - ( DAMP_alphav(i,j,k) * ( var(i,j,k,I_VELY)-vely_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i,j+1,k,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMY)   &
                                            +        var(i-2,j  ,k  ,I_MOMY) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMY)   &
                                            +        var(i  ,j-2,k  ,I_MOMY) ) &
                                          + (        var(i  ,j  ,k+2,I_MOMY)   &
                                            - 3.D0 * var(i  ,j  ,k+1,I_MOMY)   &
                                            + 2.D0 * var(i  ,j  ,k  ,I_MOMY) ) )

          k = KS+1
          momy_t(i,j,k) = - ( ( qflx_hi(i,j  ,k,XDIR)-qflx_hi(i-1,j,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j+1,k,YDIR)-qflx_hi(i  ,j,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j  ,k,ZDIR)-qflx_hi(i  ,j,k-1,ZDIR) ) * RDZ ) &
                          - ( var(i,j+1,k,I_PRES)-var(i,j,k,I_PRES) ) * RDY                &
                          - ( DAMP_alphav(i,j,k) * ( var(i,j,k,I_VELY)-vely_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i,j+1,k,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMY)   &
                                            +        var(i-2,j  ,k  ,I_MOMY) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMY)   &
                                            +        var(i  ,j-2,k  ,I_MOMY) ) &
                                          + (        var(i  ,j  ,k+2,I_MOMY)   &
                                            - 4.D0 * var(i  ,j  ,k+1,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 3.D0 * var(i  ,j  ,k-1,I_MOMY) ) )

          k = KE-1
          momy_t(i,j,k) = - ( ( qflx_hi(i,j  ,k,XDIR)-qflx_hi(i-1,j,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j+1,k,YDIR)-qflx_hi(i  ,j,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j  ,k,ZDIR)-qflx_hi(i  ,j,k-1,ZDIR) ) * RDZ ) &
                          - ( var(i,j+1,k,I_PRES)-var(i,j,k,I_PRES) ) * RDY                &
                          - ( DAMP_alphav(i,j,k) * ( var(i,j,k,I_VELY)-vely_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i,j+1,k,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMY)   &
                                            +        var(i-2,j  ,k  ,I_MOMY) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMY)   &
                                            +        var(i  ,j-2,k  ,I_MOMY) ) &
                                          + (        var(i  ,j  ,k-2,I_MOMY)   &
                                            - 4.D0 * var(i  ,j  ,k-1,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 3.D0 * var(i  ,j  ,k+1,I_MOMY) ) )

          k = KE
          momy_t(i,j,k) = - ( ( qflx_hi(i,j  ,k,XDIR)-qflx_hi(i-1,j,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j+1,k,YDIR)-qflx_hi(i  ,j,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j  ,k,ZDIR)-qflx_hi(i  ,j,k-1,ZDIR) ) * RDZ ) &
                          - ( var(i,j+1,k,I_PRES)-var(i,j,k,I_PRES) ) * RDY                &
                          - ( DAMP_alphav(i,j,k) * ( var(i,j,k,I_VELY)-vely_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i,j+1,k,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMY)   &
                                            +        var(i-2,j  ,k  ,I_MOMY) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMY)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMY)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMY)   &
                                            +        var(i  ,j-2,k  ,I_MOMY) ) &
                                          + (        var(i  ,j  ,k-2,I_MOMY)   &
                                            - 3.D0 * var(i  ,j  ,k-1,I_MOMY)   &
                                            + 2.D0 * var(i  ,j  ,k  ,I_MOMY) ) )
       enddo
       enddo

       !--- z-direction momentum equation

       ! at (u, y, interface)
       k = WS ! bottom boundary
       do j = JS,   JE
       do i = IS-1, IE
          qflx_hi(i,j,k,XDIR) = 0.5D0 * var(i,j,KS,I_VELX)                             &
                              * ( FACT_N * ( var(i+1,j,k,I_MOMZ)+var(i  ,j,k,I_MOMZ) ) &
                                + FACT_F * ( var(i+2,j,k,I_MOMZ)+var(i-1,j,k,I_MOMZ) ) )
       enddo
       enddo
       k = WE ! top boundary
       do j = JS,   JE
       do i = IS-1, IE
          qflx_hi(i,j,k,XDIR) = 0.5D0 * var(i,j,KE,I_VELX)                             &
                              * ( FACT_N * ( var(i+1,j,k,I_MOMZ)+var(i  ,j,k,I_MOMZ) ) &
                                + FACT_F * ( var(i+2,j,k,I_MOMZ)+var(i-1,j,k,I_MOMZ) ) )
       enddo
       enddo
       do k = WS+1, WE-1
       do j = JS,   JE
       do i = IS-1, IE
          qflx_hi(i,j,k,XDIR) = 0.25D0 * ( var(i,j,k+1,I_VELX)+var(i,j,k,I_VELX) )     &
                              * ( FACT_N * ( var(i+1,j,k,I_MOMZ)+var(i  ,j,k,I_MOMZ) ) &
                                + FACT_F * ( var(i+2,j,k,I_MOMZ)+var(i-1,j,k,I_MOMZ) ) )
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       k = WS ! bottom boundary
       do j = JS-1, JE
       do i = IS,   IE
          qflx_hi(i,j,k,YDIR) = 0.5D0 * var(i,j,KS,I_VELY)                             &
                              * ( FACT_N * ( var(i,j+1,k,I_MOMZ)+var(i,j  ,k,I_MOMZ) ) &
                                + FACT_F * ( var(i,j+2,k,I_MOMZ)+var(i,j-1,k,I_MOMZ) ) )
       enddo
       enddo
       k = WE ! top boundary
       do j = JS-1, JE
       do i = IS,   IE
          qflx_hi(i,j,k,YDIR) = 0.5D0 * var(i,j,KE,I_VELY)                             &
                              * ( FACT_N * ( var(i,j+1,k,I_MOMZ)+var(i,j  ,k,I_MOMZ) ) &
                                + FACT_F * ( var(i,j+2,k,I_MOMZ)+var(i,j-1,k,I_MOMZ) ) )
       enddo
       enddo
       do k = WS+1, WE-1
       do j = JS-1, JE
       do i = IS,   IE
          qflx_hi(i,j,k,YDIR) = 0.25D0 * ( var(i,j,k+1,I_VELY)+var(i,j,k,I_VELY) )     &
                              * ( FACT_N * ( var(i,j+1,k,I_MOMZ)+var(i,j  ,k,I_MOMZ) ) &
                                + FACT_F * ( var(i,j+2,k,I_MOMZ)+var(i,j-1,k,I_MOMZ) ) )
       enddo
       enddo
       enddo
       ! at (x, y, layer)
       do k = KS+1, KE-1
       do j = JS,   JE
       do i = IS,   IE
          qflx_hi(i,j,k,ZDIR) = 0.25D0 * ( var(i,j,k,I_VELZ)+var(i,j,k-1,I_VELZ) )      &
                              * ( FACT_N * ( var(i,j,k  ,I_MOMZ)+var(i,j,k-1,I_MOMZ) ) &
                                + FACT_F * ( var(i,j,k+1,I_MOMZ)+var(i,j,k-2,I_MOMZ) ) )
       enddo
       enddo
       enddo
       qflx_hi  (:,:,KS,ZDIR) = 0.D0 ! bottom cell center
       qflx_hi  (:,:,KE,ZDIR) = 0.D0 ! top    cell center

       !--- calc tendency with numerical filter
       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          momz_t(i,j,k) = - ( ( qflx_hi(i,j,k  ,XDIR)-qflx_hi(i-1,j  ,k,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j,k  ,YDIR)-qflx_hi(i  ,j-1,k,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j,k+1,ZDIR)-qflx_hi(i  ,j  ,k,ZDIR) ) * RDZ ) &    ! flux
                          - ( var(i,j,k+1,I_PRES)-var(i,j,k,I_PRES) ) * RDZ                & ! Pressure gradient
                          - ( DAMP_alphaw(i,j,k) * ( var(i,j,k,I_VELZ)-velz_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) )          ) & ! Damping factor
                          - GRAV * 0.5D0 * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) )       & ! Gravity
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMZ)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMZ)   &
                                            +        var(i-2,j  ,k  ,I_MOMZ) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMZ)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMZ)   &
                                            +        var(i  ,j-2,k  ,I_MOMZ) ) &
                                          + (        var(i  ,j  ,k+2,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j  ,k+1,I_MOMZ)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j  ,k-1,I_MOMZ)   &
                                            +        var(i  ,j  ,k-2,I_MOMZ) ) ) ! numerical filter
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          k = WS+1
          momz_t(i,j,k) = - ( ( qflx_hi(i,j,k  ,XDIR)-qflx_hi(i-1,j  ,k,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j,k  ,YDIR)-qflx_hi(i  ,j-1,k,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j,k+1,ZDIR)-qflx_hi(i  ,j  ,k,ZDIR) ) * RDZ ) &
                          - ( var(i,j,k+1,I_PRES)-var(i,j,k,I_PRES) ) * RDZ                &
                          - ( DAMP_alphaw(i,j,k) * ( var(i,j,k,I_VELZ)-velz_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - GRAV * 0.5D0 * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) )       &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMZ)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMZ)   &
                                            +        var(i-2,j  ,k  ,I_MOMZ) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMZ)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMZ)   &
                                            +        var(i  ,j-2,k  ,I_MOMZ) ) &
                                          + (        var(i  ,j  ,k+2,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j  ,k+1,I_MOMZ)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMZ)   &
                                            - 3.D0 * var(i  ,j  ,k-1,I_MOMZ) ) )

          k = WE-1
          momz_t(i,j,k) = - ( ( qflx_hi(i,j,k  ,XDIR)-qflx_hi(i-1,j  ,k,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j,k  ,YDIR)-qflx_hi(i  ,j-1,k,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j,k+1,ZDIR)-qflx_hi(i  ,j  ,k,ZDIR) ) * RDZ ) &
                          - ( var(i,j,k+1,I_PRES)-var(i,j,k,I_PRES) ) * RDZ                &
                          - ( DAMP_alphaw(i,j,k) * ( var(i,j,k,I_VELZ)-velz_ref(i,j,k) )   &
                            * 0.5D0 * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) )          ) &
                          - GRAV * 0.5D0 * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) )       &
                          - DIFF / dtrk * ( (        var(i+2,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i+1,j  ,k  ,I_MOMZ)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i-1,j  ,k  ,I_MOMZ)   &
                                            +        var(i-2,j  ,k  ,I_MOMZ) ) &
                                          + (        var(i  ,j+2,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j+1,k  ,I_MOMZ)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j-1,k  ,I_MOMZ)   &
                                            +        var(i  ,j-2,k  ,I_MOMZ) ) &
                                          + (        var(i  ,j  ,k-2,I_MOMZ)   &
                                            - 4.D0 * var(i  ,j  ,k-1,I_MOMZ)   &
                                            + 6.D0 * var(i  ,j  ,k  ,I_MOMZ)   &
                                            - 3.D0 * var(i  ,j  ,k+1,I_MOMZ) ) )
       enddo
       enddo

       !##### Thermodynamic Equation #####
       ! at (u, y, layer)
       do k = KS,   KE
       do j = JS,   JE
       do i = IS-1, IE
          qflx_hi(i,j,k,XDIR) = 0.5D0 * mflx_hi(i,j,k,XDIR)                            &
                              * ( FACT_N * ( var(i+1,j,k,I_POTT)+var(i  ,j,k,I_POTT) ) &
                                + FACT_F * ( var(i+2,j,k,I_POTT)+var(i-1,j,k,I_POTT) ) )
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do k = KS,   KE
       do j = JS-1, JE
       do i = IS,   IE
          qflx_hi(i,j,k,YDIR) = 0.5D0 * mflx_hi(i,j,k,YDIR)                            &
                              * ( FACT_N * ( var(i,j+1,k,I_POTT)+var(i,j  ,k,I_POTT) ) &
                                + FACT_F * ( var(i,j+2,k,I_POTT)+var(i,j-1,k,I_POTT) ) )
       enddo
       enddo
       enddo
       ! at (x, y, interface)
       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          qflx_hi(i,j,k,ZDIR) = 0.5D0 * mflx_hi(i,j,k,ZDIR)                            &
                              * ( FACT_N * ( var(i,j,k+1,I_POTT)+var(i,j,k  ,I_POTT) ) &
                                + FACT_F * ( var(i,j,k+2,I_POTT)+var(i,j,k-1,I_POTT) ) )
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_hi(i,j,WS  ,ZDIR) = 0.D0                                          ! bottom boundary
          qflx_hi(i,j,WS+1,ZDIR) = 0.5D0 * mflx_hi(i,j,WS+1,ZDIR)              &
                                 * ( var(i,j,WS+2,I_POTT)+var(i,j,WS+1,I_POTT) ) ! just above the bottom boundary
          qflx_hi(i,j,WE-1,ZDIR) = 0.5D0 * mflx_hi(i,j,WE-1,ZDIR)              &
                                 * ( var(i,j,WE  ,I_POTT)+var(i,j,WE-1,I_POTT) ) ! just below the top boundary
          qflx_hi(i,j,WE  ,ZDIR) = 0.D0                                          ! top boundary
       enddo
       enddo


       !--- use difference for numerical filter
       do k = KS, KE
          rhot_diff(:,:,k) = var(:,:,k,I_RHOT) - REF_pott(:,:,k)*var(:,:,k,I_DENS)
       enddo

       !--- calc tendency with numerical filter
       do k = KS+2, KE-2
       do j = JS,   JE
       do i = IS,   IE
          rhot_t(i,j,k) = - ( ( qflx_hi(i,j,k,XDIR)-qflx_hi(i-1,j  ,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j,k,YDIR)-qflx_hi(i  ,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j,k,ZDIR)-qflx_hi(i  ,j  ,k-1,ZDIR) ) * RDZ ) &   ! flux
                          - DAMP_alphat(i,j,k) * ( var(i,j,k,I_RHOT)-pott_ref(i,j,k)*var(i,j,k,I_DENS) ) & ! Damping factor
                          - DIFF / dtrk * ( (        rhot_diff(i+2,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i+1,j  ,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i-1,j  ,k  )   &
                                            +        rhot_diff(i-2,j  ,k  ) ) &
                                          + (        rhot_diff(i  ,j+2,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j+1,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j-1,k  )   &
                                            +        rhot_diff(i  ,j-2,k  ) ) &
                                          + (        rhot_diff(i  ,j  ,k+2)   &
                                            - 4.D0 * rhot_diff(i  ,j  ,k+1)   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j  ,k-1)   &
                                            +        rhot_diff(i  ,j  ,k-2) ) ) ! numerical filter
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          k = KS
          rhot_t(i,j,k) = - ( ( qflx_hi(i,j,k,XDIR)-qflx_hi(i-1,j  ,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j,k,YDIR)-qflx_hi(i  ,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j,k,ZDIR)-qflx_hi(i  ,j  ,k-1,ZDIR) ) * RDZ ) &
                          - DAMP_alphat(i,j,k) * ( var(i,j,k,I_RHOT)-pott_ref(i,j,k)*var(i,j,k,I_DENS) ) &
                          - DIFF / dtrk * ( (        rhot_diff(i+2,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i+1,j  ,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i-1,j  ,k  )   &
                                            +        rhot_diff(i-2,j  ,k  ) ) &
                                          + (        rhot_diff(i  ,j+2,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j+1,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j-1,k  )   &
                                            +        rhot_diff(i  ,j-2,k  ) ) &
                                          + (        rhot_diff(i  ,j  ,k+2)   &
                                            - 3.D0 * rhot_diff(i  ,j  ,k+1)   &
                                            + 2.D0 * rhot_diff(i  ,j  ,k  ) ) )

          k = KS+1
          rhot_t(i,j,k) = - ( ( qflx_hi(i,j,k,XDIR)-qflx_hi(i-1,j  ,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j,k,YDIR)-qflx_hi(i  ,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j,k,ZDIR)-qflx_hi(i  ,j  ,k-1,ZDIR) ) * RDZ ) &
                          - DAMP_alphat(i,j,k) * ( var(i,j,k,I_RHOT)-pott_ref(i,j,k)*var(i,j,k,I_DENS) ) &
                          - DIFF / dtrk * ( (        rhot_diff(i+2,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i+1,j  ,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i-1,j  ,k  )   &
                                            +        rhot_diff(i-2,j  ,k  ) ) &
                                          + (        rhot_diff(i  ,j+2,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j+1,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j-1,k  )   &
                                            +        rhot_diff(i  ,j-2,k  ) ) &
                                          + (        rhot_diff(i  ,j  ,k+2)   &
                                            - 4.D0 * rhot_diff(i  ,j  ,k+1)   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 3.D0 * rhot_diff(i  ,j  ,k-1) ) )

          k = KE-1
          rhot_t(i,j,k) = - ( ( qflx_hi(i,j,k,XDIR)-qflx_hi(i-1,j  ,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j,k,YDIR)-qflx_hi(i  ,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j,k,ZDIR)-qflx_hi(i  ,j  ,k-1,ZDIR) ) * RDZ ) &
                          - DAMP_alphat(i,j,k) * ( var(i,j,k,I_RHOT)-pott_ref(i,j,k)*var(i,j,k,I_DENS) ) &
                          - DIFF / dtrk * ( (        rhot_diff(i+2,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i+1,j  ,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i-1,j  ,k  )   &
                                            +        rhot_diff(i-2,j  ,k  ) ) &
                                          + (        rhot_diff(i  ,j+2,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j+1,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j-1,k  )   &
                                            +        rhot_diff(i  ,j-2,k  ) ) &
                                          + (        rhot_diff(i  ,j  ,k-2)   &
                                            - 4.D0 * rhot_diff(i  ,j  ,k-1)   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 3.D0 * rhot_diff(i  ,j  ,k+1) ) )

          k = KE
          rhot_t(i,j,k) = - ( ( qflx_hi(i,j,k,XDIR)-qflx_hi(i-1,j  ,k  ,XDIR) ) * RDX   &
                            + ( qflx_hi(i,j,k,YDIR)-qflx_hi(i  ,j-1,k  ,YDIR) ) * RDY   &
                            + ( qflx_hi(i,j,k,ZDIR)-qflx_hi(i  ,j  ,k-1,ZDIR) ) * RDZ ) &
                          - DAMP_alphat(i,j,k) * ( var(i,j,k,I_RHOT)-pott_ref(i,j,k)*var(i,j,k,I_DENS) ) &
                          - DIFF / dtrk * ( (        rhot_diff(i+2,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i+1,j  ,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i-1,j  ,k  )   &
                                            +        rhot_diff(i-2,j  ,k  ) ) &
                                          + (        rhot_diff(i  ,j+2,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j+1,k  )   &
                                            + 6.D0 * rhot_diff(i  ,j  ,k  )   &
                                            - 4.D0 * rhot_diff(i  ,j-1,k  )   &
                                            +        rhot_diff(i  ,j-2,k  ) ) &
                                          + (        rhot_diff(i  ,j  ,k-2)   &
                                            - 3.D0 * rhot_diff(i  ,j  ,k-1)   &
                                            + 2.D0 * rhot_diff(i  ,j  ,k  ) ) )
       enddo
       enddo

    enddo ! RK loop

    ! tendency rho * theta -> theta
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       pott_t(i,j,k) = rhot_t(i,j,k) / var(i,j,k,I_DENS)
    enddo
    enddo
    enddo

    !##### advection of scalar quantity #####

    ! calc low-order mass flux and high-low difference (at last step of RK)
    do k = KS,   KE
    do j = JS,   JE
    do i = IS-1, IE
       mflx_lo(i,j,k,XDIR) = 0.5D0 * (     var(i,j,k,I_VELX)  * ( var(i+1,j,k,I_DENS)+var(i,j,k,I_DENS) ) &
                                     - abs(var(i,j,k,I_VELX)) * ( var(i+1,j,k,I_DENS)-var(i,j,k,I_DENS) ) )
    enddo
    enddo
    enddo

    do k = KS,   KE
    do j = JS-1, JE
    do i = IS,   IE
       mflx_lo(i,j,k,YDIR) = 0.5D0 * (     var(i,j,k,I_VELY)  * ( var(i,j+1,k,I_DENS)+var(i,j,k,I_DENS) ) &
                                     - abs(var(i,j,k,I_VELY)) * ( var(i,j+1,k,I_DENS)-var(i,j,k,I_DENS) ) )
    enddo
    enddo
    enddo

    do k = WS+2, WE-2
    do j = JS,   JE
    do i = IS,   IE
       mflx_lo(i,j,k,ZDIR) = 0.5D0 * (     var(i,j,k,I_VELZ)  * ( var(i,j,k+1,I_DENS)+var(i,j,k,I_DENS) ) &
                                     - abs(var(i,j,k,I_VELZ)) * ( var(i,j,k+1,I_DENS)-var(i,j,k,I_DENS) ) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       mflx_lo(i,j,WS  ,ZDIR) = 0.D0                                                                         ! bottom boundary
       mflx_lo(i,j,WS+1,ZDIR) = 0.5D0 * var(i,j,WS+1,I_VELZ) * ( var(i,j,WS+2,I_DENS)+var(i,j,WS+1,I_DENS) ) ! just above the bottom boundary
       mflx_lo(i,j,WE-1,ZDIR) = 0.5D0 * var(i,j,WE-1,I_VELZ) * ( var(i,j,WE  ,I_DENS)+var(i,j,WE-1,I_DENS) ) ! just below the top boundary
       mflx_lo(i,j,WE  ,ZDIR) = 0.D0                                                                         ! top boundary
    enddo
    enddo

    do iq = 1, QA

       !--- calc mass flux * phi
       do k = KS,   KE
       do j = JS,   JE
       do i = IS-1, IE
          qflx_lo(i,j,k,XDIR) = 0.5D0 * (     mflx_lo(i,j,k,XDIR)  * ( qtrc(i+1,j,k,iq)+qtrc(i,j,k,iq) ) &
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
          qflx_lo(i,j,k,YDIR) = 0.5D0 * (     mflx_lo(i,j,k,YDIR)  * ( qtrc(i,j+1,k,iq)+qtrc(i,j,k,iq) ) &
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
          qtrc_t(i,j,k,iq) = - qdiv

          var_l = qtrc(i,j,k,iq) + dtrk * qtrc_t(i,j,k,iq) / var(i,j,k,I_DENS)

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
          qtrc_t(i,j,k,iq) = qtrc_t(i,j,k,iq)                                             &
                           - ( ( qflx_anti(i,j,k,XDIR)-qflx_anti(i-1,j  ,k  ,XDIR) ) * RDX &
                             + ( qflx_anti(i,j,k,YDIR)-qflx_anti(i  ,j-1,k  ,YDIR) ) * RDY &
                             + ( qflx_anti(i,j,k,ZDIR)-qflx_anti(i  ,j  ,k-1,ZDIR) ) * RDZ &
                             ) / var(i,j,k,I_DENS)
       enddo
       enddo
       enddo

    enddo ! scalar quantities loop

    return
  end subroutine ATMOS_DYN

end module mod_atmos_dyn
