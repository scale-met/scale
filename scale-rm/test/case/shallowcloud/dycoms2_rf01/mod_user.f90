!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-20 (S.Nishizawa)   [new] split from dynamical core
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  use scale_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

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
  real(RP), private, allocatable :: MOMZ_LS(:,:)
  real(RP), private, allocatable :: MOMZ_LS_DZ(:,:)
  real(RP), private, allocatable :: QV_LS(:,:)
  logical,  private, save        :: MOMZ_LS_FLG(6)
  real(RP), private, allocatable :: U_GEOS(:)
  real(RP), private, allocatable :: V_GEOS(:)

  real(RP), private, save :: USER_LSsink_D = 0.0_RP         ! large scale sinking parameter [1/s]
  real(RP), private, save :: USER_LSsink_bottom = 0.0_RP    ! large scale sinking parameter [m]
  ! flag
  integer,  private, save :: USER_LS_FLG = 0 !-- 0->no force, 1->dycoms, 2->rico

  real(RP), private, allocatable :: Q_rate(:,:,:,:)
  real(RP), private, save :: corioli

  real(RP), private, save :: F0    = 70.00_RP  ! Upward [J/m2/s]
  real(RP), private, save :: F1    = 22.00_RP  ! [K/m**-1/3]
  real(RP), private, save :: Dval  = 3.75E-6_RP ! divergence of large scale horizontal winds [1/s]
  real(RP), private, parameter :: kappa = 85.00_RP  ! scaling factor for LWP [m2/kg]
  real(RP), private, parameter :: a     =  1.00_RP  ! [K/m**-1/3]
  logical, private :: first = .true.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       PI => CONST_PI
    use scale_grid, only: &
       CZ => GRID_CZ, &
       FZ => GRID_FZ
    implicit none

    namelist / PARAM_USER / &
       USER_LSsink_D,          &
       USER_LSsink_bottom,     &
       USER_LS_FLG

    real(RP) :: zovzb
    real(RP) :: cr, xi, ar, yc, xz

    integer :: k
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    allocate( MOMZ_LS(KA,2) )
    allocate( MOMZ_LS_DZ(KA,2) )
    allocate( QV_LS(KA,2) )
    allocate( V_GEOS(KA) )
    allocate( U_GEOS(KA) )
    allocate( Q_rate( KA,IA,JA,QA ) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_USER)

    if( USER_LS_FLG == 0 ) then  ! no large scale sinking

       MOMZ_LS(:,:) = 0.0_RP
       MOMZ_LS_DZ(:,:) = 0.0_RP
       MOMZ_LS_FLG( : ) = .false.
       QV_LS(:,:) = 0.0_RP
       Q_rate( :,:,:,: ) = 0.0_RP
       V_GEOS(:) = 0.0_RP
       U_GEOS(:) = 0.0_RP
       corioli = 0.0_RP

    elseif( USER_LS_FLG == 1 ) then ! DYCOMS

       do k = KS, KE
          if ( CZ(k) < 0.0_RP ) then
             MOMZ_LS(k,1) = 0.0_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          else if ( CZ(k) < USER_LSsink_bottom ) then
             zovzb = CZ(k) / USER_LSsink_bottom
             MOMZ_LS(k,1) = - 0.5_RP * USER_LSsink_D * USER_LSsink_bottom &
                  * ( zovzb - sin(PI*zovzb)/PI )
             MOMZ_LS_DZ(k,1) = - 0.5_RP * USER_LSsink_d &
                  * ( 1.0_RP - cos(PI*zovzb) )
          else
             MOMZ_LS(k,1) = - USER_LSsink_D * ( CZ(k) - USER_LSsink_bottom * 0.5_RP )
             MOMZ_LS_DZ(k,1) = - USER_LSsink_D
          end if
       enddo
       do k = KS-1, KE
          if ( FZ(k) < 0.0_RP ) then
             MOMZ_LS(k,2) = 0.0_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          else if ( FZ(k) < USER_LSsink_bottom ) then
             zovzb = FZ(k) / USER_LSsink_bottom
             MOMZ_LS(k,2) = - 0.5_RP * USER_LSsink_D * USER_LSsink_bottom &
                  * ( zovzb - sin(PI*zovzb)/PI )
             MOMZ_LS_DZ(k,2) = - 0.5_RP * USER_LSsink_D &
                  * ( 1.0_RP - cos(PI*zovzb) )
          else
             MOMZ_LS(k,2) = - USER_LSsink_D * ( FZ(k) - USER_LSsink_bottom * 0.5_RP )
             MOMZ_LS_DZ(k,2) = - USER_LSsink_D
          end if
       enddo
       MOMZ_LS_FLG( : ) = .true.
       QV_LS(:,:) = 0.0_RP
       Q_rate( :,:,:,: ) = 0.0_RP

       do k = KS-1, KE
          V_GEOS(k) = -5.5_RP
          U_GEOS(k) =  7.0_RP
       enddo
       corioli = 7.6E-5_RP

    elseif( USER_LS_FLG == 2 ) then ! RICO

       do k = KS, KE
          if ( CZ(k) < 0.0_RP ) then
             MOMZ_LS(k,1) = 0.0_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          else if( CZ(k) < 2000.0_RP ) then
             MOMZ_LS(k,1) = - 5.0E-3_RP / 2260.0_RP * CZ(k)
             MOMZ_LS_DZ(k,1) = - 5.0E-3_RP / 2260.0_RP
             !--- fitted by circle
          else if( CZ(k) < 2520.0_RP ) then
             xi   = 2260.0_RP
             ar   = -5.0E-3_RP/xi
             cr   = ar* ( 2520.0_RP-xi )/( 1.0_RP - sqrt( 1.0_RP + ar*ar ) )
             xz   = CZ(k)
             yc   = ar * xi + cr
             MOMZ_LS(k,1) = ar * xi + cr  &
                  - sqrt( cr*cr - ( xz - xi - cr / ar * ( 1.0_RP - sqrt( 1.0_RP + ar*ar ) ) )**2.0_RP )
             MOMZ_LS_DZ(k,1) = - ( xz - 2520._RP )/( MOMZ_LS(k,1) - yc )
          else
             MOMZ_LS(k,1) = - 5.0E-3_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          end if
          if( CZ(k) < 2980.0_RP ) then
             QV_LS(k,1) = ( -1.0_RP + 1.3456_RP / 2980.0_RP  &
                  * CZ(k) ) * 1.E-3_RP / 86400.0_RP  !--- [kg/kg/s]
          else
             QV_LS(k,1) = 4.0_RP * 1.E-6_RP * 1.E-3_RP !--- [kg/kg/s]
          endif
       enddo
       do k = KS-1, KE
          if ( FZ(k) < 0.0_RP ) then
             MOMZ_LS(k,2) = 0.0_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          else if( FZ(k) < 2000.0_RP ) then
             MOMZ_LS(k,2) = - 5.0E-3_RP / 2260.0_RP * FZ(k)
             MOMZ_LS_DZ(k,2) = - 5.0E-3_RP / 2260.0_RP
          else if( FZ(k) < 2520.0_RP ) then
             xi   = 2260.0_RP
             ar   = -5.0E-3_RP/xi
             cr   = ar* ( 2520.0_RP-xi )/( 1.0_RP - sqrt( 1.0_RP + ar*ar ) )
             xz   = FZ(k)
             yc   = ar * xi + cr
             MOMZ_LS(k,2) = ar * xi + cr  &
                  - sqrt( cr*cr - ( xz - xi - cr / ar * ( 1.0_RP - sqrt( 1.0_RP + ar*ar ) ) )**2.0_RP )
             MOMZ_LS_DZ(k,2) = - ( xz - 2520._RP )/( MOMZ_LS(k,2) - yc )
          else
             MOMZ_LS(k,2) = - 5.0E-3_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          end if
       enddo
       MOMZ_LS_FLG( : ) = .false.
       MOMZ_LS_FLG( I_RHOT ) = .true.
       MOMZ_LS_FLG( I_QTRC ) = .true.
       Q_rate( :,:,:,: ) = 0.0_RP

       do k = KS-1, KE
          V_GEOS(k) = -3.8_RP
          U_GEOS(k) = -9.9_RP + 2.0E-3_RP * CZ(k)
       enddo
       corioli = 4.5E-5_RP

    endif

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_step
    use mod_atmos_vars, only: &
         DENS,    &
         MOMZ,    &
         MOMX,    &
         MOMY,    &
         RHOT,    &
         QTRC,    &
         MOMZ_tp, &
         MOMX_tp, &
         MOMY_tp, &
         RHOT_tp, &
         RHOQ_tp
    use scale_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ, &
         CZ   => GRID_CZ, &
         FZ   => GRID_FZ
    use scale_const, only: &
         Rdry  => CONST_Rdry,  &
         CPdry => CONST_CPdry, &
         CVdry => CONST_CVdry, &
         P00   => CONST_PRE00
    use scale_time, only: &
       dtrd =>  TIME_DTSEC_ATMOS_PHY_RD
    use scale_atmos_hydrometeor, only: &
       I_QV, &
       QLS, &
       QLE
    use mod_admin_time, only: &
       do_phy_rd => TIME_DOATMOS_PHY_RD
    implicit none

    real(RP) :: ratesum

    real(RP) :: WORK(KA,IA,JA)
    real(RP) :: VELX(KA,IA,JA), VELY(KA,IA,JA)
    integer :: k, i, j, iq
    integer :: IIS, IIE, JJS, JJE

    real(RP) :: CDZ(KA)
    real(RP) :: Zi(IA,JA) ! Cloud top height [m]
    real(RP) :: dZ, dZ_CBRT
    integer  :: k_cldtop

    real(RP) :: QTOT ! Qv + Qc + Qr [kg/kg]
    real(RP) :: Qbelow, Qabove ! scaled LWP (above/below layer)
    real(RP) :: dQ, QWSUM
    real(RP) :: TEMP_t(KA,IA,JA)
    real(RP) :: flux_rad(KA,IA,JA)
    real(RP) :: RovCP, RovCV

    integer :: k2
    !---------------------------------------------------------------------------

    RovCP = Rdry / CPdry
    RovCV = Rdry / CVdry

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       if ( MOMZ_LS_FLG(I_MOMZ) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             WORK(k,i,j) = MOMZ(k,i,j) * 2.0_RP / ( DENS(k+1,i,j) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-2
             MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) &
                  - MOMZ_LS(k,2) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RCDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             MOMZ_tp(KE-1,i,j) = MOMZ_tp(KE-1,i,j) &
                  - MOMZ_LS(KE-1,2) * (           - WORK(KE-1,i,j) ) * RCDZ(KE-1)
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_MOMX) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             WORK(k,i,j) = MOMX(k,i,j) * 2.0_RP / ( DENS(k,i+1,j) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE+1
          do k = KS, KE
             VELY(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             MOMX_tp(k,i,j) = MOMX_tp(k,i,j) &
                  + 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                  * ( - CORIOLI * V_GEOS(k) &
                      + CORIOLI * 0.25_RP &
                      * ( VELY(k,i,j)+VELY(k,i+1,j)+VELY(k,i,j-1)+VELY(k,i+1,j-1) ) &
                    ) &
                  - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             MOMX_tp(KE,i,j) = MOMX_tp(KE,i,j) &
                  + 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                  *  ( - CORIOLI * V_GEOS(KE) &
                       + CORIOLI * 0.25_RP &
                       * ( VELY(KE,i,j)+VELY(KE,i+1,j)+VELY(KE,i,j-1)+VELY(KE,i+1,j-1) ) &
                     ) &
                  - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1)
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_MOMY) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             WORK(k,i,j) = MOMY(k,i,j) * 2.0_RP / ( DENS(k,i,j+1) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE+1
          do i = IIS-1, IIE
          do k = KS, KE
             VELX(k,i,j) = MOMX(k,i,j) * 2.0_RP / ( DENS(k,i+1,j)+DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             MOMY_tp(k,i,j) = MOMY_tp(k,i,j) &
                  + 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )  &
                  * ( + CORIOLI * U_GEOS(k) &
                      - CORIOLI * 0.25_RP &
                      * ( VELX(k,i,j)+VELX(k,i,j+1)+VELX(k,i-1,j)+VELX(k,i-1,j+1) ) &
                    ) &
                  - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE+1
          do i = IIS-1, IIE
             MOMY_tp(KE,i,j) = MOMY_tp(KE,i,j) &
                  + 0.5_RP * ( DENS(KE,i,j+1)+DENS(KE,i,j) ) &
                  * ( + CORIOLI * U_GEOS(KE) &
                      - CORIOLI * 0.25_RP  &
                      * ( VELX(KE,i,j)+VELX(KE,i,j+1)+VELX(KE,i-1,j)+VELX(KE,i-1,j+1) ) &
                    ) &
                  - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1)
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_RHOT) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             WORK(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             RHOT_tp(k,i,j) = RHOT_tp(k,i,j) &
                  - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             RHOT_tp(KE,i,j) = RHOT_tp(KE,i,j) &
                  - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1)
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_QTRC) ) then

          !##### advection of scalar quantity #####
          if( USER_LS_FLG == 2 ) then ! RICO
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE
                Q_rate( k,i,j,: ) = 0.0_RP
                ratesum = 0.0_RP
                do iq = I_QV, QLE
                   ratesum = ratesum + QTRC(k,i,j,iq)
                enddo
                do iq = I_QV, QLE
                   Q_rate( k,i,j,iq ) = QTRC(k,i,j,iq) / ratesum
                enddo
             enddo
             enddo
             enddo
          endif


          do iq = 1, QA
             !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE-1
                RHOQ_tp(k,i,j,iq) = RHOQ_tp(k,i,j,iq) &
                     + QV_LS(k,1) * Q_rate(k,i,j,iq) &
                     - MOMZ_LS(k,1) * ( QTRC(k+1,i,j,iq) - QTRC(k,i,j,iq) ) * RFDZ(k)
             enddo
             enddo
             enddo
             !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
                RHOQ_tp(KE,i,j,iq) = RHOQ_tp(KE,i,j,iq) &
                     + QV_LS(KE,1) * Q_rate(KE,i,j,iq) &
                     - MOMZ_LS(KE,1) * ( QTRC(KE,i,j,iq) - QTRC(KE-1,i,j,iq) ) * RFDZ(KE-1)
             enddo
             enddo
          end do
       end if

    enddo
    enddo


    if( do_phy_rd ) then

      if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Parametarized Radiation (DYCOMS-II)'

      flux_rad    (:,:,:) = 0.0_RP

      Zi(:,:) = 0.0_RP ! for history

      do j = JS, JE
      do i = IS, IE
       do k = KS, KE
          CDZ(k) = FZ(k) - FZ(k-1)
       enddo

       ! diagnose cloud top
       k_cldtop = -1
       do k = KS, KE
          QTOT = 0.0_RP
          do iq = I_QV, QLE
             QTOT = QTOT + QTRC(k,i,j,iq)
          enddo
          if( QTOT < 8.E-3_RP ) exit ! above cloud

          k_cldtop = k
       enddo

       if( k_cldtop == -1 ) cycle ! no cloud

       Zi(i,j) = CZ(k_cldtop)

       do k = KS-1, KE
          Qbelow = 0.0_RP
          Qabove = 0.0_RP

          do k2 = KS, KE
             QWSUM = 0.0_RP
             do iq = QLS, QLE
                QWSUM = QWSUM + QTRC(k2,i,j,iq)
             enddo
             dQ = kappa * CDZ(k2) * DENS(k2,i,j) * QWSUM
             if ( k2 <= k ) then ! below layer
                Qbelow = Qbelow + dQ
             else                ! above layer
                Qabove = Qabove + dQ
             endif
          enddo

          flux_rad(k,i,j) = F0 * exp( -Qabove ) + F1 * exp( -Qbelow )
        enddo

        do k = k_cldtop, KE
           dZ      = FZ(k)-CZ(k_cldtop)
           dZ_CBRT = dZ**(1.0_RP/3.0_RP)

           QTOT = 0.0_RP
           do iq = I_QV, QLE
              QTOT = QTOT + QTRC(k,i,j,iq)
           enddo

           flux_rad(k,i,j) = flux_rad(k,i,j) &
                           + a * DENS(k_cldtop,i,j)*( 1.0_RP-QTOT ) * CPdry * Dval &
                           * ( 0.250_RP * dZ * dZ_CBRT + CZ(k_cldtop) * dZ_CBRT )
        enddo

        do k = KS, KE
          TEMP_t(k,i,j) = - ( flux_rad(k,i,j) - flux_rad(k-1,i,j) ) / CPdry * RCDZ(k)
          RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + ( 1.0_RP - RovCP ) &
                         * ( P00/(RHOT(k,i,j)*Rdry) )**RovCV * TEMP_t(k,i,j)
        enddo


       enddo
       enddo

    endif

    return
  end subroutine USER_step

end module mod_user
