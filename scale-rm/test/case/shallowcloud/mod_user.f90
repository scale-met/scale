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
  real(RP),               private :: USER_LSsink_D      = 0.0_RP ! large scale sinking parameter [1/s]
  real(RP),               private :: USER_LSsink_bottom = 0.0_RP ! large scale sinking parameter [m]
  character(len=H_SHORT), private :: USER_LS_TYPE       = 'NONE'
                                                        ! 'DYCOMS2_RF01'
                                                        ! 'DYCOMS2_RF02'
                                                        ! 'RICO'
                                                        ! 'BOMEX'

  real(RP), private, allocatable :: MOMZ_LS    (:,:)
  real(RP), private, allocatable :: MOMZ_LS_DZ (:,:)
  real(RP), private, allocatable :: QV_LS      (:,:)
  real(RP), private, allocatable :: Q_rate     (:,:,:,:)
  real(RP), private, allocatable :: U_GEOS     (:)
  real(RP), private, allocatable :: V_GEOS     (:)
  logical,  private              :: MOMZ_LS_FLG(6)
  real(RP), private              :: corioli

  !--- for radiation flux (DYCOMS-II)
  real(RP), private              :: F0    =   70.00_RP ! Upward [J/m2/s]
  real(RP), private              :: F1    =   22.00_RP ! [K/m**-1/3]
  real(RP), private              :: Dval  = 3.75E-6_RP ! divergence of large scale horizontal winds [1/s]
  real(RP), private, parameter   :: kappa =   85.00_RP ! scaling factor for LWP [m2/kg]
  real(RP), private, parameter   :: a     =    1.00_RP ! [K/m**-1/3]

  !--- for surface flux (RICO)
  real(RP), private              :: Cm_min      =    0.0E-5_RP ! minimum bulk coef. of u,v,w
  real(RP), private, parameter   :: Cm_max      =    2.5E-3_RP ! maximum bulk coef. of u,v,w
  real(RP), private              :: U_minM      =       0.0_RP ! minimum U_abs for u,v,w
  real(RP), private, parameter   :: U_maxM      =     100.0_RP ! maximum U_abs for u,v,w
  real(RP), private              :: U_minE      =       0.0_RP ! minimum U_abs for u,v,w
  real(RP), private, parameter   :: U_maxE      =     100.0_RP ! maximum U_abs for u,v,w
  real(RP), private              :: U_minH      =       0.0_RP ! minimum U_abs for u,v,w
  real(RP), private, parameter   :: U_maxH      =     100.0_RP ! maximum U_abs for u,v,w
  real(RP), private              :: Cm_const    =  0.001229_RP ! constant bulk coef. of u,v,w
  real(RP), private              :: Ce_const    =  0.001133_RP ! constant bulk coef. of u,v,w
  real(RP), private              :: Ch_const    =  0.001094_RP ! constant bulk coef. of u,v,w
  real(RP), private              :: FIXED_PTSST =     298.5_RP ! fixed PT of surface [K]
  real(RP), private              :: FIXED_SST   =     299.8_RP ! fixed Temperature of surface [K]
  real(RP), private, parameter   :: pres_sfc    =  1015.4E2_RP ! fixed surface pressure

  !--- for radiation flux (BOMEX)
  real(RP), private, allocatable :: dQrad(:)
  real(RP), private, parameter   :: USER_Ustar  = 0.28_RP      ! constant frection velocity [m/s]
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config
    implicit none
    !---------------------------------------------------------------------------

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
       USER_LSsink_D,      &
       USER_LSsink_bottom, &
       USER_LS_TYPE

    real(RP) :: zovzb
    real(RP) :: cr, xi, ar, yc, xz

    integer  :: k
    integer  :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    allocate( MOMZ_LS   (KA,2)        )
    allocate( MOMZ_LS_DZ(KA,2)        )
    allocate( QV_LS     (KA,2)        )
    allocate( Q_rate    (KA,IA,JA,QA) )
    allocate( V_GEOS    (KA)          )
    allocate( U_GEOS    (KA)          )
    allocate( dQrad     (KA)          )

    if ( USER_LS_TYPE == 'NONE' ) then  ! no large scale sinking

       MOMZ_LS    (:,:)     = 0.0_RP
       MOMZ_LS_DZ (:,:)     = 0.0_RP
       QV_LS      (:,:)     = 0.0_RP
       Q_rate     (:,:,:,:) = 0.0_RP
       V_GEOS     (:)       = 0.0_RP
       U_GEOS     (:)       = 0.0_RP
       MOMZ_LS_FLG(:)       = .false.
       corioli              = 0.0_RP

    elseif( USER_LS_TYPE == 'DYCOMS2_RF01' ) then ! DYCOMS RF-01

       do k = KS, KE
          if    ( CZ(k) < 0.0_RP ) then
             MOMZ_LS   (k,1) = 0.0_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          elseif( CZ(k) < USER_LSsink_bottom ) then
             zovzb = CZ(k) / USER_LSsink_bottom

             MOMZ_LS   (k,1) = - 0.5_RP * USER_LSsink_D * USER_LSsink_bottom * ( zovzb - sin(PI*zovzb)/PI )
             MOMZ_LS_DZ(k,1) = - 0.5_RP * USER_LSsink_D * ( 1.0_RP - cos(PI*zovzb) )
          else
             MOMZ_LS   (k,1) = - USER_LSsink_D * ( CZ(k) - USER_LSsink_bottom * 0.5_RP )
             MOMZ_LS_DZ(k,1) = - USER_LSsink_D
          endif
       end do

       do k = KS-1, KE
          if    ( FZ(k) < 0.0_RP ) then
             MOMZ_LS   (k,2) = 0.0_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          elseif( FZ(k) < USER_LSsink_bottom ) then
             zovzb = FZ(k) / USER_LSsink_bottom

             MOMZ_LS   (k,2) = - 0.5_RP * USER_LSsink_D * USER_LSsink_bottom * ( zovzb - sin(PI*zovzb)/PI )
             MOMZ_LS_DZ(k,2) = - 0.5_RP * USER_LSsink_D * ( 1.0_RP - cos(PI*zovzb) )
          else
             MOMZ_LS   (k,2) = - USER_LSsink_D * ( FZ(k) - USER_LSsink_bottom * 0.5_RP )
             MOMZ_LS_DZ(k,2) = - USER_LSsink_D
          endif
       end do
       QV_LS (:,:)     = 0.0_RP
       Q_rate(:,:,:,:) = 0.0_RP

       do k = KS-1, KE
          V_GEOS(k) = -5.5_RP
          U_GEOS(k) =  7.0_RP
       end do

       MOMZ_LS_FLG(:) = .true.
       corioli        = 7.6E-5_RP

    elseif( USER_LS_TYPE == 'DYCOMS2_RF02' ) then ! DYCOMS RF-02

       do k = KS, KE
          if    ( CZ(k) < 0.0_RP ) then
             MOMZ_LS   (k,1) = 0.0_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          elseif( CZ(k) < USER_LSsink_bottom ) then
             zovzb = CZ(k) / USER_LSsink_bottom

             MOMZ_LS   (k,1) = - 0.5_RP * USER_LSsink_D * USER_LSsink_bottom * ( zovzb - sin(PI*zovzb)/PI )
             MOMZ_LS_DZ(k,1) = - 0.5_RP * USER_LSsink_D * ( 1.0_RP - cos(PI*zovzb) )
          else
             MOMZ_LS   (k,1) = - USER_LSsink_D * ( CZ(k) - USER_LSsink_bottom * 0.5_RP )
             MOMZ_LS_DZ(k,1) = - USER_LSsink_D
          endif
       end do

       do k = KS-1, KE
          if    ( FZ(k) < 0.0_RP ) then
             MOMZ_LS   (k,2) = 0.0_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          elseif( FZ(k) < USER_LSsink_bottom ) then
             zovzb = FZ(k) / USER_LSsink_bottom

             MOMZ_LS   (k,2) = - 0.5_RP * USER_LSsink_D * USER_LSsink_bottom * ( zovzb - sin(PI*zovzb)/PI )
             MOMZ_LS_DZ(k,2) = - 0.5_RP * USER_LSsink_D * ( 1.0_RP - cos(PI*zovzb) )
          else
             MOMZ_LS   (k,2) = - USER_LSsink_D * ( FZ(k) - USER_LSsink_bottom * 0.5_RP )
             MOMZ_LS_DZ(k,2) = - USER_LSsink_D
          endif
       end do

       QV_LS (:,:)     = 0.0_RP
       Q_rate(:,:,:,:) = 0.0_RP

       do k = KS-1, KE
          U_GEOS(k) =  3.0_RP + 4.3 * CZ(k)*1.E-3_RP
          V_GEOS(k) = -9.0_RP + 5.6 * CZ(k)*1.E-3_RP
       end do

       MOMZ_LS_FLG(:) = .true.
       corioli        = 7.6E-5_RP

    elseif( USER_LS_TYPE == 'RICO' ) then ! RICO

       xi = 2260.0_RP
       ar = -5.E-3_RP / xi
       cr = ar * ( 2520.0_RP-xi ) / ( 1.0_RP - sqrt( 1.0_RP+ar*ar ) )
       yc = ar * xi + cr

       do k = KS, KE
          if   ( CZ(k) < 0.0_RP ) then
             MOMZ_LS   (k,1) =    0.0_RP
             MOMZ_LS_DZ(k,1) =    0.0_RP
          elseif( CZ(k) < 2000.0_RP ) then
             MOMZ_LS   (k,1) = -5.E-3_RP / 2260.0_RP * CZ(k)
             MOMZ_LS_DZ(k,1) = -5.E-3_RP / 2260.0_RP
          elseif( CZ(k) < 2520.0_RP ) then ! fitted by circle
             xz = CZ(k)

             MOMZ_LS   (k,1) = yc - sqrt( cr*cr - ( xz - xi - cr/ar*( 1.0_RP - sqrt( 1.0_RP+ar*ar ) ) )**2 )
             MOMZ_LS_DZ(k,1) = -( xz - 2520._RP ) / ( MOMZ_LS(k,1) - yc )
          else
             MOMZ_LS   (k,1) = -5.E-3_RP
             MOMZ_LS_DZ(k,1) =    0.0_RP
          endif

          if( CZ(k) < 2980.0_RP ) then
             QV_LS(k,1) = ( -1.0_RP + 1.3456_RP / 2980.0_RP * CZ(k) ) * 1.E-3_RP / 86400.0_RP  ! [kg/kg/s]
          else
             QV_LS(k,1) = 4.0_RP * 1.E-6_RP * 1.E-3_RP ! [kg/kg/s]
          endif
       end do

       do k = KS-1, KE
          if    ( FZ(k) < 0.0_RP ) then
             MOMZ_LS   (k,2) =    0.0_RP
             MOMZ_LS_DZ(k,2) =    0.0_RP
          elseif( FZ(k) < 2000.0_RP ) then
             MOMZ_LS   (k,2) = -5.E-3_RP / 2260.0_RP * FZ(k)
             MOMZ_LS_DZ(k,2) = -5.E-3_RP / 2260.0_RP
          elseif( FZ(k) < 2520.0_RP ) then
             xz = FZ(k)

             MOMZ_LS   (k,2) = yc - sqrt( cr*cr - ( xz - xi - cr/ar*( 1.0_RP - sqrt( 1.0_RP+ar*ar ) ) )**2 )
             MOMZ_LS_DZ(k,2) = -( xz - 2520._RP ) / ( MOMZ_LS(k,2) - yc )
          else
             MOMZ_LS   (k,2) = -5.E-3_RP
             MOMZ_LS_DZ(k,2) =    0.0_RP
          endif
       end do

       Q_rate(:,:,:,:) = 0.0_RP

       do k = KS-1, KE
          V_GEOS(k) = -3.8_RP
          U_GEOS(k) = -9.9_RP + 2.E-3_RP * CZ(k)
       end do

       MOMZ_LS_FLG(:)      = .false.
       MOMZ_LS_FLG(I_RHOT) = .true.
       MOMZ_LS_FLG(I_QTRC) = .true.
       corioli             = 4.5E-5_RP

    elseif( USER_LS_TYPE == 'BOMEX' ) then ! BOMEX

       do k = KS, KE
          if   ( CZ(k) < 0.0_RP ) then
             MOMZ_LS   (k,1) =    0.0_RP
             MOMZ_LS_DZ(k,1) =    0.0_RP
          elseif( CZ(k) < 1500.0_RP ) then
             MOMZ_LS   (k,1) = -6.5E-3_RP / 1500.0_RP * CZ(k)
             MOMZ_LS_DZ(k,1) = -6.5E-3_RP / 1500.0_RP
          elseif( CZ(k) < 2100.0_RP ) then
             MOMZ_LS   (k,1) = - 6.5E-3_RP + 6.5E-3_RP / ( 2100.0_RP - 1500.0_RP ) * ( CZ(k) - 1500.0_RP )
             MOMZ_LS_DZ(k,1) =   6.5E-3_RP / ( 2100.0_RP - 1500.0_RP )
          else
             MOMZ_LS   (k,1) = 0.0_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          endif

          if( CZ(k) < 300.0_RP ) then
             QV_LS(k,1) = -1.2E-8_RP
          elseif( CZ(k) < 500.0_RP ) then
             QV_LS(k,1) = - 1.2E-8_RP + 1.2E-8_RP * ( CZ(k) - 300.0_RP ) / 200.0_RP  !--- [kg/kg/s]
          else
             QV_LS(k,1) = 0.0_RP !--- [kg/kg/s]
          endif
       end do

       do k = KS-1, KE
          if    ( FZ(k) < 0.0_RP ) then
             MOMZ_LS   (k,2) =    0.0_RP
             MOMZ_LS_DZ(k,2) =    0.0_RP
          elseif( FZ(k) < 1500.0_RP ) then
             MOMZ_LS   (k,2) = -6.5E-3_RP / 1500.0_RP * FZ(k)
             MOMZ_LS_DZ(k,2) = -6.5E-3_RP / 1500.0_RP
          elseif( FZ(k) < 2100.0_RP ) then
             MOMZ_LS(k,2) = - 6.5E-3_RP + 6.5E-3_RP / ( 2100.0_RP - 1500.0_RP ) * ( FZ(k) - 1500.0_RP )
             MOMZ_LS_DZ(k,1) = 6.5E-3_RP / ( 2100.0_RP - 1500.0_RP )
          else
             MOMZ_LS(k,2) = 0.0_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          end if
       end do

       Q_rate(:,:,:,:) = 0.0_RP

       do k = KS-1, KE
          V_GEOS(k) = 0.0_RP
          U_GEOS(k) = -10.0_RP + 1.8E-3_RP * CZ(k)
       end do

       MOMZ_LS_FLG(:)      = .false.
       MOMZ_LS_FLG(I_MOMX) = .true.
       MOMZ_LS_FLG(I_MOMY) = .true.
       MOMZ_LS_FLG(I_RHOT) = .true.
       MOMZ_LS_FLG(I_QTRC) = .true.
       corioli             = 3.76E-5_RP

       do k = KS, KE
         if( CZ(k) < 1500.0_RP ) then
           dQrad(k) = - 2.315E-5_RP
         elseif( CZ(k) < 2500.0_RP ) then
           dQrad(k) = - 2.315E-5_RP + 2.315E-5_RP / 1000.0_RP * ( CZ(k) - 1500.0_RP )
         else
           dQrad(k) = 0.0_RP
         endif
       enddo

    else
       write(*,*) 'xxx Not appropriate type for USER_LS_TYPE. STOP.', trim(USER_LS_TYPE)
       call PRC_MPIstop
    endif

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    call USER_step

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
    use scale_const, only: &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       Rdry  => CONST_Rdry,  &
       LHV   => CONST_LHV0,  &
       P00   => CONST_PRE00
    use scale_atmos_hydrometeor, only: &
       I_QV, &
       QHS,  &
       QHE
    use scale_grid, only: &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ, &
       CDZ  => GRID_CDZ,  &
       CZ   => GRID_CZ,   &
       FZ   => GRID_FZ
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd => ATMOS_THERMODYN_qd, &
       THERMODYN_cp => ATMOS_THERMODYN_cp, &
       THERMODYN_r  => ATMOS_THERMODYN_r
    use scale_atmos_saturation, only : &
       moist_pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq
    use scale_history, only: &
       HIST_in
    use mod_admin_time, only: &
       do_phy_rd => TIME_DOATMOS_PHY_RD, &
       do_phy_sf => TIME_DOATMOS_PHY_SF
    use mod_atmos_vars, only: &
       DENS,    &
       MOMZ,    &
       MOMX,    &
       MOMY,    &
       RHOT,    &
       QTRC,    &
       DENS_tp, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       RHOQ_tp
    implicit none

    real(RP) :: WORK(KA,IA,JA)
    real(RP) :: VELX(KA,IA,JA)
    real(RP) :: VELY(KA,IA,JA)
    real(RP) :: ratesum

    real(RP) :: TEMP_t  (KA,IA,JA) ! tendency rho*theta     [K*kg/m3/s]
    real(RP) :: flux_rad(KA,IA,JA)
    real(RP) :: Zi      (IA,JA)    ! Cloud top height [m]
    real(RP) :: dZ, dZ_CBRT
    integer  :: k_cldtop
    real(RP) :: Qbelow, Qabove ! scaled LWP (above/below layer)
    real(RP) :: dQ, QWSUM

    real(RP) :: SHFLX    (IA,JA)
    real(RP) :: LHFLX    (IA,JA)
    real(RP) :: SFLX_MOMZ(IA,JA)
    real(RP) :: SFLX_MOMX(IA,JA)
    real(RP) :: SFLX_MOMY(IA,JA)
    real(RP) :: SFLX_POTT(IA,JA)
    real(RP) :: SFLX_QV  (IA,JA)
    real(RP) :: drhot, LWPT, LWPT_a
    real(RP) :: QHYD, PRES
    real(RP) :: qv_evap
    real(RP) :: Uabs  ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Cm, Ch, Ce

    real(RP) :: q(QA)
    real(RP) :: qdry
    real(RP) :: qtot
    real(RP) :: Rtot
    real(RP) :: CPtot
    real(RP) :: RovCP
    real(RP) :: RovCV
    real(RP) :: CPovCV

    integer  :: k, i, j, iq, k2
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Large scale sinking'

    RovCP = Rdry / CPdry
    RovCV = Rdry / CVdry

    if ( MOMZ_LS_FLG(I_MOMZ) ) then

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE-1
          WORK(k,i,j) = MOMZ(k,i,j) * 2.0_RP / ( DENS(k+1,i,j)+DENS(k,i,j) )
       end do
       end do
       end do

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE-2
          MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) - MOMZ_LS(k,2) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RCDZ(k)
       end do
       end do
       end do

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
          MOMZ_tp(KE-1,i,j) = MOMZ_tp(KE-1,i,j) - MOMZ_LS(KE-1,2) * ( - WORK(KE-1,i,j) ) * RCDZ(KE-1)
       end do
       end do

    endif

    if ( MOMZ_LS_FLG(I_MOMX) .OR. MOMZ_LS_FLG(I_MOMY) ) then

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS,   JE+1
       do i = IS-1, IE
       do k = KS, KE
          VELX(k,i,j) = MOMX(k,i,j) * 2.0_RP / ( DENS(k,i+1,j)+DENS(k,i,j) )
       end do
       end do
       end do

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS-1, JE
       do i = IS,   IE+1
       do k = KS, KE
          VELY(k,i,j) = MOMY(k,i,j) * 2.0_RP / ( DENS(k,i,j+1)+DENS(k,i,j) )
       end do
       end do
       end do

    endif

    if ( MOMZ_LS_FLG(I_MOMX) ) then

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE-1
          MOMX_tp(k,i,j) = MOMX_tp(k,i,j) - MOMZ_LS(k,1) * ( VELX(k+1,i,j) - VELX(k,i,j) ) * RFDZ(k) &
                         + 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                         * ( - CORIOLI * V_GEOS(k)                                                             &
                             + CORIOLI * 0.25_RP * ( VELY(k,i,j)+VELY(k,i+1,j)+VELY(k,i,j-1)+VELY(k,i+1,j-1) ) )
       end do
       end do
       end do

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
          MOMX_tp(KE,i,j) = MOMX_tp(KE,i,j) - MOMZ_LS(KE,1) * ( VELX(KE,i,j) - VELX(KE-1,i,j) ) * RFDZ(KE-1) &
                          + 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                          * ( - CORIOLI * V_GEOS(KE)                                                                &
                              + CORIOLI * 0.25_RP * ( VELY(KE,i,j)+VELY(KE,i+1,j)+VELY(KE,i,j-1)+VELY(KE,i+1,j-1) ) )
       end do
       end do

    endif

    if ( MOMZ_LS_FLG(I_MOMY) ) then

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE-1
          MOMY_tp(k,i,j) = MOMY_tp(k,i,j) - MOMZ_LS(k,1) * ( VELY(k+1,i,j) - VELY(k,i,j) ) * RFDZ(k) &
                         + 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) ) &
                         * ( + CORIOLI * U_GEOS(k)                                                             &
                             - CORIOLI * 0.25_RP * ( VELX(k,i,j)+VELX(k,i,j+1)+VELX(k,i-1,j)+VELX(k,i-1,j+1) ) )
       end do
       end do
       end do

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
          MOMY_tp(KE,i,j) = MOMY_tp(KE,i,j) - MOMZ_LS(KE,1) * ( VELY(KE,i,j) - VELY(KE-1,i,j) ) * RFDZ(KE-1) &
                          + 0.5_RP * ( DENS(KE,i,j+1)+DENS(KE,i,j) ) &
                          * ( + CORIOLI * U_GEOS(KE)                                                                &
                              - CORIOLI * 0.25_RP * ( VELX(KE,i,j)+VELX(KE,i,j+1)+VELX(KE,i-1,j)+VELX(KE,i-1,j+1) ) )
       end do
       end do

    endif

    if ( MOMZ_LS_FLG(I_RHOT) ) then

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          WORK(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       end do
       end do
       end do

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE-1
          RHOT_tp(k,i,j) = RHOT_tp(k,i,j) - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k)
       end do
       end do
       end do

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JS, JE
       do i = IS, IE
          RHOT_tp(KE,i,j) = RHOT_tp(KE,i,j) - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1)
       end do
       end do

    endif

    if ( MOMZ_LS_FLG(I_QTRC) ) then

       !##### advection of scalar quantity #####
       if ( USER_LS_TYPE == 'RICO' .OR. USER_LS_TYPE == 'BOMEX' ) then ! RICO or BOMEX
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             Q_rate(k,i,j,:) = 0.0_RP

             ratesum = QTRC(k,i,j,I_QV)
             do iq = QHS, QHE
                ratesum = ratesum + QTRC(k,i,j,iq)
             end do

             do iq = I_QV, QHE
                Q_rate(k,i,j,iq) = QTRC(k,i,j,iq) / ratesum
             end do
          end do
          end do
          end do
       endif

       do iq = 1, QA
          if ( .not. TRACER_ADVC(iq) ) cycle

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE-1
             RHOQ_tp(k,i,j,iq) = RHOQ_tp(k,i,j,iq) - MOMZ_LS(k,1) * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) * RFDZ(k) &
                               + QV_LS(k,1) * Q_rate(k,i,j,iq)
          end do
          end do
          end do

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JS, JE
          do i = IS, IE
             RHOQ_tp(KE,i,j,iq) = RHOQ_tp(KE,i,j,iq) - MOMZ_LS(KE,1) * ( QTRC(KE,i,j,iq)-QTRC(KE-1,i,j,iq) ) * RFDZ(KE-1) &
                                + QV_LS(KE,1) * Q_rate(KE,i,j,iq)
          end do
          end do
       end do

    endif

    if ( do_phy_rd ) then
       if ( USER_LS_TYPE == 'DYCOMS2_RF01' .OR. USER_LS_TYPE == 'DYCOMS2_RF02' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Parametarized Radiation (DYCOMS-II)'

          flux_rad(:,:,:) = 0.0_RP

          do j = JS, JE
          do i = IS, IE
             ! diagnose cloud top
             Zi(i,j) = 0.0_RP
             k_cldtop = -1
             do k = KS, KE
                qtot = QTRC(k,i,j,I_QV)
                do iq = QHS, QHE
                   qtot = qtot + QTRC(k,i,j,iq)
                end do
                if( qtot < 8.E-3_RP ) exit ! above cloud
                k_cldtop = k
             end do

             if( k_cldtop == -1 ) cycle ! no cloud

             Zi(i,j) = CZ(k_cldtop)

             do k = KS-1, KE
                Qbelow = 0.0_RP
                Qabove = 0.0_RP

                do k2 = KS, KE
                   QWSUM = 0.0_RP
                   do iq = QHS, QHE
                      QWSUM = QWSUM + QTRC(k2,i,j,iq)
                   end do
                   dQ = kappa * CDZ(k2) * DENS(k2,i,j) * QWSUM

                   if ( k2 <= k ) then ! below layer
                      Qbelow = Qbelow + dQ
                   else                ! above layer
                      Qabove = Qabove + dQ
                   endif
                end do

                flux_rad(k,i,j) = F0 * exp( -Qabove ) &
                                + F1 * exp( -Qbelow )
             end do

             do k = k_cldtop, KE
                dZ      = FZ(k) - CZ(k_cldtop)
                dZ_CBRT = dZ**(1.0_RP/3.0_RP)

                qtot = QTRC(k,i,j,I_QV)
                do iq = QHS, QHE
                   qtot = qtot + QTRC(k,i,j,iq)
                end do

                flux_rad(k,i,j) = flux_rad(k,i,j) + a * DENS(k_cldtop,i,j) * ( 1.0_RP-qtot ) * CPdry * Dval &
                                                  * ( 0.250_RP * dZ * dZ_CBRT + CZ(k_cldtop) * dZ_CBRT )
             end do
          end do
          end do

          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             TEMP_t (k,i,j) = -( flux_rad(k,i,j) - flux_rad(k-1,i,j) ) / CPdry * RCDZ(k)
             drhot  = ( 1.0_RP - RovCP ) * ( P00/(RHOT(k,i,j)*Rdry) )**RovCV * TEMP_t(k,i,j)

             RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + drhot
          end do
          end do
          end do

          call HIST_in( flux_rad(:,:,:), 'RADFLUX_DYCOMS', 'net radiation flux', 'W/m2', nohalo=.true. )
          call HIST_in( Zi      (:,:),   'Zi',             'cloud top height'  , 'm'   , nohalo=.true. )

       elseif( USER_LS_TYPE == 'RICO' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Parametarized Radiation (RICO)'

          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             do iq = 1, QA
                q(iq) = QTRC(KS,i,j,iq)
             enddo
             call THERMODYN_qd( qdry,  q, TRACER_MASS )
             call THERMODYN_cp( CPtot, q, TRACER_CP, qdry )
             call THERMODYN_r ( Rtot,  q, TRACER_R,  qdry )
             CPovCV = CPTOT / ( CPTOT - RTOT )
             RovCP  = RTOT / CPTOT

             QHYD = 0.0_RP
             do iq = QHS, QHE
                QHYD = QHYD + QTRC(k,i,j,iq)
             enddo

             PRES   = P00 * ( RHOT(k,i,j) * RTOT / P00 )**CPovCV
             LWPT   = RHOT(k,i,j) / DENS(k,i,j) - ( LHV / CPdry * QHYD ) * ( P00 / PRES )**RovCP
             LWPT_a = LWPT - 2.5_RP / 86400.0_RP
             drhot  = ( LWPT_a + ( LHV / CPdry * QHYD ) * ( P00 / PRES )**RovCP ) * DENS(k,i,j) - RHOT(k,i,j)

             RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + drhot
          enddo
          enddo
          enddo

       elseif( USER_LS_TYPE == 'BOMEX' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Parametarized Radiation (BOMEX)'

          do j = JS, JE
          do i = IS, IE

             do k = KS, KE
                  RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + dQrad(k) * DENS(k,i,j)
             enddo

          enddo
          enddo

       endif
    endif

    if ( do_phy_sf ) then
       if ( USER_LS_TYPE == 'RICO' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Surface flux (RICO)'

          do j = JS-1, JE
          do i = IS-1, IE

             ! at cell center

             !--- absolute velocity
             Uabs = sqrt( ( MOMZ(KS,i,j)                  )**2 &
                        + ( MOMX(KS,i-1,j) + MOMX(KS,i,j) )**2 &
                        + ( MOMY(KS,i,j-1) + MOMY(KS,i,j) )**2 ) / DENS(KS,i,j) * 0.5_RP

             !--- Bulk coef. at w, theta, and qv points
             Cm = Cm_const
             Ch = Ch_const
             Ce = Ce_const

             !--- saturation at surface
             qtot = 0.0_RP
             do iq = 1, QA
                q(iq) = QTRC(KS,i,j,iq)
                qtot  = QTRC(KS,i,j,iq) * TRACER_MASS(iq)
             enddo

             call THERMODYN_qd( qdry,  q, TRACER_MASS )
             call THERMODYN_cp( CPtot, q, TRACER_CP, qdry )
             call THERMODYN_r ( Rtot,  q, TRACER_R,  qdry )

             CPovCV = CPtot / ( CPtot - Rtot )
             pres   = P00 * ( RHOT(KS,i,j) * Rtot / P00 )**CPovCV

             call moist_pres2qsat_liq( qv_evap, FIXED_SST, pres_sfc )

             ! flux
             SFLX_MOMZ(i,j) = 0.0_RP
             SFLX_POTT(i,j) = Ch * min(max(Uabs,U_minH),U_maxH) * ( FIXED_PTSST * DENS(KS,i,j) - RHOT(KS,i,j) )
             SFLX_QV  (i,j) = Ce * min(max(Uabs,U_minE),U_maxE) * DENS(KS,i,j) * ( qv_evap - qtot )

             ! at (u, y, layer)
             Uabs = sqrt( &
                  + ( 2.0_RP *   MOMX(KS,i,j) )**2 &
                  + ( 0.5_RP * ( MOMY(KS,i,j-1) + MOMY(KS,i,j) + MOMY(KS,i+1,j-1) + MOMY(KS,i+1,j) ) )**2 &
                  ) / ( DENS(KS,i,j) + DENS(KS,i+1,j) )
             SFLX_MOMX(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMX(KS,i,j)

             ! at (x, v, layer)
             Uabs = sqrt( &
                  + ( 0.5_RP * ( MOMX(KS,i-1,j) + MOMX(KS,i,j) + MOMX(KS,i-1,j+1) + MOMX(KS,i,j+1) ) )**2 &
                  + ( 2.0_RP *   MOMY(KS,i,j) )**2 &
                ) / ( DENS(KS,i,j) + DENS(KS,i,j+1) )
             SFLX_MOMY(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMY(KS,i,j)

             SHFLX(i,j) = SFLX_POTT(i,j) * CPdry
             LHFLX(i,j) = SFLX_QV  (i,j) * LHV
          enddo
          enddo

          call HIST_in( SHFLX(:,:), 'SHFLX', 'sensible heat flux', 'W/m2' )
          call HIST_in( LHFLX(:,:), 'LHFLX', 'latent heat flux',   'W/m2' )

          !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IE
             RHOT_tp(KS,i,j)      = RHOT_tp(KS,i,j) + ( SFLX_POTT(i,j) &
                                                      + SFLX_QV(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
                                                      ) * RCDZ(KS)
             DENS_tp(KS,i,j)      = DENS_tp(KS,i,j)      + SFLX_QV  (i,j) * RCDZ(KS)
             MOMZ_tp(KS,i,j)      = MOMZ_tp(KS,i,j)      + SFLX_MOMZ(i,j) * RFDZ(KS)
             MOMX_tp(KS,i,j)      = MOMX_tp(KS,i,j)      + SFLX_MOMX(i,j) * RCDZ(KS)
             MOMY_tp(KS,i,j)      = MOMY_tp(KS,i,j)      + SFLX_MOMY(i,j) * RCDZ(KS)
             RHOQ_tp(KS,i,j,I_QV) = RHOQ_tp(KS,i,j,I_QV) + SFLX_QV  (i,j) * RCDZ(KS)
          enddo
          enddo

       elseif ( USER_LS_TYPE == 'BOMEX' ) then

          if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Surface flux (BOMEX)'

          do j = JS-1, JE
          do i = IS-1, IE

             ! at cell center

             !--- saturation at surface
             qtot = 0.0_RP
             do iq = 1, QA
                q(iq) = QTRC(KS,i,j,iq)
                qtot  = QTRC(KS,i,j,iq) * TRACER_MASS(iq)
             enddo

             call THERMODYN_qd( qdry,  q, TRACER_MASS )
             call THERMODYN_cp( CPtot, q, TRACER_CP, qdry )
             call THERMODYN_r ( Rtot,  q, TRACER_R,  qdry )

             CPovCV = CPtot / ( CPtot - Rtot )
             pres   = P00 * ( RHOT(KS,i,j) * Rtot / P00 )**CPovCV

             call moist_pres2qsat_liq( qv_evap, FIXED_SST, pres_sfc )

             ! flux
             SFLX_MOMZ(i,j) = 0.0_RP
             SFLX_POTT(i,j) = 8.0E-3_RP
             SFLX_QV  (i,j) = 5.2E-5_RP

             ! at (u, y, layer)
             Uabs = sqrt( &
                  + ( 2.0_RP *   MOMX(KS,i,j) )**2 &
                  + ( 0.5_RP * ( MOMY(KS,i,j-1) + MOMY(KS,i,j) + MOMY(KS,i+1,j-1) + MOMY(KS,i+1,j) ) )**2 &
                  ) / ( DENS(KS,i,j) + DENS(KS,i+1,j) )
             Cm = min( max(USER_Ustar**2 / Uabs**2, Cm_min), Cm_max )
             SFLX_MOMX(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMX(KS,i,j)

             ! at (x, v, layer)
             Uabs = sqrt( &
                  + ( 0.5_RP * ( MOMX(KS,i-1,j) + MOMX(KS,i,j) + MOMX(KS,i-1,j+1) + MOMX(KS,i,j+1) ) )**2 &
                  + ( 2.0_RP *   MOMY(KS,i,j) )**2 &
                ) / ( DENS(KS,i,j) + DENS(KS,i,j+1) )
             Cm = min( max(USER_Ustar**2 / Uabs**2, Cm_min), Cm_max )
             SFLX_MOMY(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMY(KS,i,j)

             SHFLX(i,j) = SFLX_POTT(i,j) * CPdry
             LHFLX(i,j) = SFLX_QV  (i,j) * LHV
          enddo
          enddo

          call HIST_in( SHFLX(:,:), 'SHFLX', 'sensible heat flux', 'W/m2' )
          call HIST_in( LHFLX(:,:), 'LHFLX', 'latent heat flux',   'W/m2' )

          !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IE
             RHOT_tp(KS,i,j)      = RHOT_tp(KS,i,j) + ( SFLX_POTT(i,j) &
                                                      + SFLX_QV(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
                                                      ) * RCDZ(KS)
             DENS_tp(KS,i,j)      = DENS_tp(KS,i,j)      + SFLX_QV  (i,j) * RCDZ(KS)
             MOMZ_tp(KS,i,j)      = MOMZ_tp(KS,i,j)      + SFLX_MOMZ(i,j) * RFDZ(KS)
             MOMX_tp(KS,i,j)      = MOMX_tp(KS,i,j)      + SFLX_MOMX(i,j) * RCDZ(KS)
             MOMY_tp(KS,i,j)      = MOMY_tp(KS,i,j)      + SFLX_MOMY(i,j) * RCDZ(KS)
             RHOQ_tp(KS,i,j,I_QV) = RHOQ_tp(KS,i,j,I_QV) + SFLX_QV  (i,j) * RCDZ(KS)
          enddo
          enddo

       endif

    endif

    return
  end subroutine USER_step

end module mod_user
