!-------------------------------------------------------------------------------
!> Module ideal initial condition
!!
!! @par Description
!!          This module calculates initial condition of ideal test cases
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_ideal_init
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio

  use scale_const, only: &
     pi    => CONST_PI,     &
     g     => CONST_GRAV,   &
     a     => CONST_RADIUS, &
     omega => CONST_OHM,    &
     Rd    => CONST_Rdry,   &
     Rv    => CONST_Rvap,   &
     Cp    => CONST_CPdry,  &
     PRE00 => CONST_PRE00
  use mod_adm, only: &
     ADM_KNONE, &
     ADM_lall,    &
     ADM_lall_pl, &
     ADM_gall,    &
     ADM_gall_pl, &
     ADM_kall,    &
     ADM_kmin,    &
     ADM_kmax
  use mod_runconf, only: &
     TRC_vmax
  use dcmip_initial_conditions_test_1_2_3, only: &
     test2_steady_state_mountain, &
     test2_schaer_mountain,       &
     test3_gravity_wave
  use baroclinic_wave, only: &
     baroclinic_wave_test
  use Terminator, only: &
     initial_value_Terminator
  use supercell, only: &
     supercell_init,   &
     supercell_test
  use tropical_cyclone, only: &
     tropical_cyclone_test
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: dycore_input
  public :: tracer_input

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public :: DCTEST_type = ''
  character(len=H_SHORT), public :: DCTEST_case = ''

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: hs_init
  private :: jbw_init
  private :: jbw_moist_init
  private :: sc_init
  private :: tc_init
  private :: tracer_init
  private :: mountwave_init
  private :: gravwave_init
  private :: tomita_init

  private :: tomita_2004
  private :: eta_vert_coord_NW
  private :: steady_state
  private :: geo2prs
  private :: ps_estimation
  private :: perturbation
  private :: conv_vxvyvz
  private :: diag_pressure
  private :: simpson

  private :: Sp_Unit_East
  private :: Sp_Unit_North

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  ! physical parameters configurations
  real(RP), private :: Kap                    ! temporal value
  real(RP), private :: d2r                    ! Degree to Radian
  real(RP), private :: r2d                    ! Radian to Degree
  real(RP), private :: zero = 0.0_RP          ! zero

  ! for Jablonowski
  real(RP), private :: clat = 40.0_RP              ! perturbation center: latitude [deg]
  real(RP), private :: clon = 20.0_RP              ! perturbation center: longitude [deg]
  real(RP), private :: etaT = 0.2_RP              ! threashold of vertical profile
  real(RP), private :: eta0 = 0.252_RP            ! threashold of vertical profile
  real(RP), private :: t0 = 288.0_RP               ! [K]
  real(RP), private :: delT = 4.8E+5_RP             ! [K]
  real(RP), private :: ganma = 0.005_RP           ! [K m^-1]
  real(RP), private :: u0 = 35.0_RP                ! [m s^-1]
  real(RP), private :: uP = 1.0_RP                 ! [m s^-1]
  real(RP), private :: p0 = 1.E+5_RP                ! [Pa]
  logical, private, parameter :: message = .false.
  integer, private, parameter :: itrmax = 100       ! # of iteration maximum

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine dycore_input( &
       DIAG_var )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: DIAG_var(ADM_gall,ADM_kall,ADM_lall,6+TRC_VMAX)

    character(len=H_SHORT) :: init_type   = ''
    character(len=H_SHORT) :: test_case   = ''
    real(RP)               :: eps_geo2prs = 1.E-2_RP
    logical                :: nicamcore   = .true.
    logical                :: chemtracer  = .false.
    logical                :: prs_rebuild = .false.

    namelist / DYCORETESTPARAM / &
       init_type,   &
       test_case,   &
       eps_geo2prs, &
       nicamcore,   &
       chemtracer,  &
       prs_rebuild

    integer :: ierr
    !---------------------------------------------------------------------------

    Kap = Rd / Cp
    d2r = pi/180.0_RP
    r2d = 180.0_RP/pi

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[dycoretest]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=DYCORETESTPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** DYCORETESTPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist DYCORETESTPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=DYCORETESTPARAM)

    DCTEST_type = init_type
    DCTEST_case = test_case

    if( IO_L ) write(IO_FID_LOG,*) '*** type: ', trim(init_type)
    select case(init_type)
!    case ('DCMIP2012-11','DCMIP2012-12','DCMIP2012-13' &
!          'DCMIP2012-200','DCMIP2012-21','DCMIP2012-22')

!       if( IO_L ) write(IO_FID_LOG,*) '*** test case: ', trim(test_case)
!       call IDEAL_init_DCMIP2012( ADM_gall, ADM_kall, ADM_lall, init_type, DIAG_var(:,:,:,:) )

    case ('Heldsuarez')

       call hs_init ( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )

    case ('Jablonowski')

       if( IO_L ) write(IO_FID_LOG,*) '*** test case   : ', trim(test_case)
       if( IO_L ) write(IO_FID_LOG,*) '*** eps_geo2prs = ', eps_geo2prs
       if( IO_L ) write(IO_FID_LOG,*) '*** nicamcore   = ', nicamcore
       call jbw_init( ADM_gall, ADM_kall, ADM_lall, test_case, eps_geo2prs, nicamcore, DIAG_var(:,:,:,:) )

    case ('Jablonowski-Moist')

       if( IO_L ) write(IO_FID_LOG,*) '*** test case   : ', trim(test_case)
       if( IO_L ) write(IO_FID_LOG,*) '*** nicamcore   = ', nicamcore
       if( IO_L ) write(IO_FID_LOG,*) '*** chemtracer  = ', chemtracer
       call jbw_moist_init( ADM_gall, ADM_kall, ADM_lall, test_case, chemtracer, DIAG_var(:,:,:,:) )

    case ('Supercell')

       if( IO_L ) write(IO_FID_LOG,*) '*** test case   : ', trim(test_case)
       if( IO_L ) write(IO_FID_LOG,*) '*** nicamcore   = ', nicamcore
       call sc_init( ADM_gall, ADM_kall, ADM_lall, test_case, prs_rebuild, DIAG_var(:,:,:,:) )

    case ('Tropical-Cyclone')

       if( IO_L ) write(IO_FID_LOG,*) '*** nicamcore   = ', nicamcore
       call tc_init( ADM_gall, ADM_kall, ADM_lall, prs_rebuild, DIAG_var(:,:,:,:) )

    case ('Traceradvection')

       if( IO_L ) write(IO_FID_LOG,*) '*** test case: ', trim(test_case)
       call tracer_init( ADM_gall, ADM_kall, ADM_lall, test_case, DIAG_var(:,:,:,:) )

    case ('Mountainwave')

       if( IO_L ) write(IO_FID_LOG,*) '*** test case: ', trim(test_case)
       call mountwave_init( ADM_gall, ADM_kall, ADM_lall, test_case, DIAG_var(:,:,:,:) )

    case ('Gravitywave')

       call gravwave_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )

    case ('Tomita2004')

       call tomita_init( ADM_gall, ADM_kall, ADM_lall, DIAG_var(:,:,:,:) )

    case default

       write(*,*) 'xxx [dycore_input] Invalid init_type. STOP.'
       call PRC_MPIstop

    end select

    return
  end subroutine dycore_input

  !-----------------------------------------------------------------------------
  subroutine tracer_input( &
       TRC_var )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_random, only: &
       RANDOM_get
    use mod_grd, only: &
       GRD_LAT, &
       GRD_LON, &
       GRD_s
    implicit none

    real(RP), intent(out) :: TRC_var(ADM_gall,ADM_kall,ADM_lall,TRC_VMAX)

    real(RP) :: random(ADM_gall,ADM_kall,ADM_lall)
    integer  :: deg

    integer  :: g, k, l, nq, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    do nq = 1, TRC_VMAX
       call RANDOM_get( random(:,:,:) )

       if ( nq == 1 ) then ! vapor (dummy for thermodynamics)
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             TRC_var(g,k,l,nq) = 0.0_RP
          enddo
          enddo
          enddo
       elseif( nq == 2 ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             deg = nint( GRD_s(g,k0,l,GRD_LON) / D2R )
             if ( mod(deg,10) == 0 ) then
                TRC_var(g,k,l,nq) = real(ADM_kall-k+1,kind=RP)
             else
                TRC_var(g,k,l,nq) = 0.0_RP
             endif
          enddo
          enddo
          enddo
       elseif( nq == 3 ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             deg = nint( GRD_s(g,k0,l,GRD_LAT) / D2R )
             if ( mod(deg,10) == 0 ) then
                TRC_var(g,k,l,nq) = real(ADM_kall-k+1,kind=RP)
             else
                TRC_var(g,k,l,nq) = 0.0_RP
             endif
          enddo
          enddo
          enddo
       else
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             TRC_var(g,k,l,nq) = random(g,k,l) * exp(real(nq,kind=RP)/10.0_RP)
          enddo
          enddo
          enddo
       endif
    enddo

    return
  end subroutine tracer_input

  !-----------------------------------------------------------------------------
  subroutine hs_init( &
       ijdim,   &
       kdim,    &
       lall,    &
       DIAG_var )
    use mod_grd, only: &
       GRD_LAT, &
       GRD_LON, &
       GRD_s,   &
       GRD_Z,   &
       GRD_vz
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: lall
    real(RP), intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(RP) :: pre(kdim), tem(kdim), dz(kdim)
    real(RP) :: pre_sfc, tem_sfc
    real(RP) :: pre_save

    real(RP), parameter :: deltaT  = 60.0_RP
    real(RP), parameter :: deltaTh = 10.0_RP
!   real(RP), parameter :: eps_hs  = 1.E-7_RP
    real(RP), parameter :: eps_hs  = max( 1.E-10_RP, (10.0_RP)**(-RP_PREC+1) )

    real(RP) :: f, df
    real(RP) :: lat, lon

    integer  :: n, k, l, k0, itr
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.0_RP

    k0 = ADM_KNONE

    do l = 1, lall
    do n = 1, ijdim

       dz(ADM_kmin) = GRD_vz(n,ADM_kmin,l,GRD_Z)
       do k = ADM_kmin+1, ADM_kmax+1
          dz(k) = GRD_vz(n,k,l,GRD_Z) - GRD_vz(n,k-1,l,GRD_Z)
       enddo

       lat = GRD_s(n,k0,l,GRD_LAT)
       lon = GRD_s(n,k0,l,GRD_LON)

       pre_sfc = PRE00
!       tem_sfc = 300.0_RP
       tem_sfc = 315.0_RP - deltaT*sin(lat)**2

       !---< from ground surface to lowermost atmosphere >---
       k = ADM_kmin
       ! first guess
       pre(k) = pre_sfc
       tem(k) = tem_sfc

       ! Newton-Lapson
       do itr = 1, itrmax
          pre_save = pre(k) ! save

          f  = log(pre(k)/pre_sfc) / dz(k) + g / ( Rd * 0.5_RP * (tem(k)+tem_sfc) )
          df = 1.0_RP / (pre(k)*dz(k))

          pre(k) = pre(k) - f / df
!          tem(k) = 300.0_RP * ( pre(k)/PRE00 )**(Rd/Cp)
          tem(k) = ( 315.0_RP - deltaT*sin(lat)**2 - deltaTh*log(pre(k)/PRE00)*cos(lat)**2 ) &
                 * ( pre(k)/PRE00 )**(Rd/Cp)
          tem(k) = max( 200.0_RP, tem(k) )

!          if( abs(pre_save-pre(k)) <= eps_hs ) exit
          if( abs(pre_save/pre(k)-1.0_RP) <= eps_hs ) exit
       enddo

       if ( itr > itrmax ) then
          write(*,*) 'xxx [dycore_input/hs_init] iteration not converged!', &
                     k, pre_save-pre(k), pre(k), pre_sfc, tem(k), tem_sfc
          stop
       endif

       !---< from lowermost to uppermost atmosphere >---
       do k = ADM_kmin+1, ADM_kmax+1

          ! first guess
          pre(k) = pre(k-1)
          tem(k) = 300.0_RP * ( pre(k)/PRE00 )**(Rd/Cp)
          tem(k) = max( 200.0_RP, tem(k) )

          ! Newton-Lapson
          do itr = 1, itrmax
             pre_save = pre(k) ! save

             f  = log(pre(k)/pre(k-1)) / dz(k) + g / ( Rd * 0.5_RP * (tem(k)+tem(k-1)) )
             df = 1.0_RP / (pre(k)*dz(k))

             pre(k) = pre(k) - f / df
!             tem(k) = 300.0_RP * ( pre(k)/PRE00 )**(Rd/Cp)
             tem(k) = ( 315.0_RP - deltaT*sin(lat)**2 - deltaTh*log(pre(k)/PRE00)*cos(lat)**2 ) &
                    * ( pre(k)/PRE00 )**(Rd/Cp)
             tem(k) = max( 200.0_RP, tem(k) )

!             if( abs(pre_save-pre(k)) <= eps_hs ) exit
             if( abs(pre_save/pre(k)-1.0_RP) <= eps_hs ) exit
          enddo

          if ( itr > itrmax ) then
             write(*,*) 'xxx [dycore_input/hs_init] iteration not converged!', &
                        k, pre_save-pre(k), pre(k), pre(k-1), tem(k), tem(k-1)
             stop
          endif
       enddo

       DIAG_var(n,ADM_kmin-1,l,1) = pre_sfc ! tentative
       DIAG_var(n,ADM_kmin-1,l,2) = tem_sfc ! tentative
       do k = ADM_kmin, ADM_kmax+1
          DIAG_var(n,k,l,1) = pre(k)
          DIAG_var(n,k,l,2) = tem(k)
       enddo

    enddo
    enddo

    return
  end subroutine hs_init

  !-----------------------------------------------------------------------------
  subroutine jbw_init( &
       ijdim,        &
       kdim,         &
       lall,         &
       test_case,    &
       eps_geo2prs,  &
       nicamcore,    &
       DIAG_var      )
    use scale_process, only: &
       PRC_MPIstop
    use mod_grd, only: &
       GRD_LAT, &
       GRD_LON, &
       GRD_s,   &
       GRD_Z,   &
       GRD_ZH,  &
       GRD_vz
    implicit none

    integer,          intent(in)  :: ijdim
    integer,          intent(in)  :: kdim
    integer,          intent(in)  :: lall
    character(len=*), intent(in)  :: test_case
    real(RP),         intent(in)  :: eps_geo2prs
    logical,          intent(in)  :: nicamcore
    real(RP),         intent(out) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(RP) :: lat, lon               ! latitude, longitude on Icosahedral grid
    real(RP) :: eta(kdim,2), geo(kdim) ! eta & geopotential in ICO-grid field
    real(RP) :: prs(kdim),   tmp(kdim) ! pressure & temperature in ICO-grid field
    real(RP) :: wix(kdim),   wiy(kdim) ! zonal/meridional wind components in ICO-grid field

    real(RP) :: z_local (kdim)
    real(RP) :: vx_local(kdim)
    real(RP) :: vy_local(kdim)
    real(RP) :: vz_local(kdim)
    real(RP) :: ps

    logical  :: signal    ! if true, continue iteration
    logical  :: pertb     ! if true, with perturbation
    logical  :: psgm      ! if true, PS Gradient Method
    logical  :: eta_limit ! if true, value of eta is limited upto 1.0
    logical  :: logout    ! log output switch for Pressure Convert

    integer  :: n, k, l, k0, itr
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.0_RP

    k0 = ADM_KNONE

    eta_limit = .true.
    psgm = .false.
    logout = .true.

    select case( trim(test_case) )
    case ('1', '4-1')  ! with perturbation
       if( IO_L ) write(IO_FID_LOG,*) "Jablonowski Initialize - case 1: with perturbation (no rebalance)"
       pertb = .true.
    case ('2', '4-2')  ! without perturbation
       if( IO_L ) write(IO_FID_LOG,*) "Jablonowski Initialize - case 2: without perturbation (no rebalance)"
       pertb = .false.
    case ('3')  ! with perturbation (PS Distribution Method)
       if( IO_L ) write(IO_FID_LOG,*) "Jablonowski Initialize - PS Distribution Method: with perturbation"
       if( IO_L ) write(IO_FID_LOG,*) "### DO NOT INPUT ANY TOPOGRAPHY ###"
       pertb = .true.
       psgm = .true.
       eta_limit = .false.
    case ('4')  ! without perturbation (PS Distribution Method)
       if( IO_L ) write(IO_FID_LOG,*) "Jablonowski Initialize - PS Distribution Method: without perturbation"
       if( IO_L ) write(IO_FID_LOG,*) "### DO NOT INPUT ANY TOPOGRAPHY ###"
       pertb = .false.
       psgm = .true.
       eta_limit = .false.
    case default
       if( IO_L ) write(IO_FID_LOG,*) "Unknown test_case: '"//trim(test_case)//"' specified."
       if( IO_L ) write(IO_FID_LOG,*) "Force changed to case 1 (with perturbation)"
       pertb = .true.
    end select
    if( IO_L ) write(IO_FID_LOG,*) " | eps for geo2prs: ", eps_geo2prs
    if( IO_L ) write(IO_FID_LOG,*) " | nicamcore switch for geo2prs: ", nicamcore

    do l = 1, lall
    do n = 1, ijdim
       z_local(ADM_kmin-1) = GRD_vz(n,2,l,GRD_ZH)
       do k = ADM_kmin, ADM_kmax+1
          z_local(k) = GRD_vz(n,k,l,GRD_Z)
       enddo

       lat = GRD_s(n,k0,l,GRD_LAT)
       lon = GRD_s(n,k0,l,GRD_LON)

       signal = .true.

       ! iteration
       do itr = 1, itrmax

          if ( itr == 1 ) then
             eta(:,:) = 1.E-7_RP ! This initial value is recommended by Jablonowsky.
          else
             call eta_vert_coord_NW( kdim, itr, z_local, tmp, geo, eta_limit, eta, signal )
          endif

          call steady_state( kdim, lat, eta, wix, wiy, tmp, geo )

          if( .NOT. signal ) exit
       enddo

       if ( itr > itrmax ) then
          write(*         ,*) 'ETA ITERATION ERROR: NOT CONVERGED', n, l
          if( IO_L ) write(IO_FID_LOG,*) 'ETA ITERATION ERROR: NOT CONVERGED', n, l
          call PRC_MPIstop
       endif

       if (psgm) then
          call ps_estimation ( kdim, lat, eta(:,1), tmp, geo, wix, ps, nicamcore )
          call geo2prs ( kdim, ps, lat, tmp, geo, wix, prs, eps_geo2prs, nicamcore, logout )
       else
          call geo2prs ( kdim, p0, lat, tmp, geo, wix, prs, eps_geo2prs, nicamcore, logout )
       endif
       logout = .false.

       call conv_vxvyvz ( kdim, lat, lon, wix, wiy, vx_local, vy_local, vz_local )

       do k = 1, kdim
          DIAG_var(n,k,l,1) = prs(k)
          DIAG_var(n,k,l,2) = tmp(k)
          DIAG_var(n,k,l,3) = vx_local(k)
          DIAG_var(n,k,l,4) = vy_local(k)
          DIAG_var(n,k,l,5) = vz_local(k)
       enddo

    enddo
    enddo

    if (pertb) call perturbation( ijdim, kdim, lall, 5, DIAG_var(:,:,:,1:5) )

    if( IO_L ) write(IO_FID_LOG,*) " |            Vertical Coordinate used in JBW initialization              |"
    if( IO_L ) write(IO_FID_LOG,*) " |------------------------------------------------------------------------|"
    do k = 1, kdim
       if( IO_L ) write(IO_FID_LOG,'(3X,"(k=",I3,") HGT:",F8.2," [m]",2X,"PRS: ",F9.2," [Pa]",2X,"GH: ",F8.2," [m]",2X,"ETA: ",F9.5)') &
       k, z_local(k), prs(k), geo(k)/g, eta(k,1)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) " |------------------------------------------------------------------------|"

    return
  end subroutine jbw_init

  !-----------------------------------------------------------------------------
  subroutine jbw_moist_init( &
       ijdim,       &
       kdim,        &
       lall,        &
       test_case,   &
       chemtracer,  &
       DIAG_var     )
    use scale_process, only: &
       PRC_MPIstop
    use mod_grd, only: &
       GRD_LAT,  &
       GRD_LON,  &
       GRD_s,    &
       GRD_afac, &
       GRD_bfac, &
       GRD_Z,    &
       GRD_ZH,   &
       GRD_vz
    use mod_vmtr, only : &
       VMTR_getIJ_PHI
    use mod_runconf, only: &
       I_QV,      &
       NQW_MAX,   &
       NCHEM_MAX, &
       NCHEM_STR, &
       NCHEM_END
    use mod_bndcnd, only : &
      BNDCND_thermo
    implicit none

    integer,          intent(in)  :: ijdim
    integer,          intent(in)  :: kdim
    integer,          intent(in)  :: lall
    character(len=*), intent(in)  :: test_case
    logical,          intent(in)  :: chemtracer
    real(RP),         intent(out) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(DP) :: DP_lon          ! latitude, longitude on Icosahedral grid
    real(DP) :: DP_lat          ! latitude, longitude on Icosahedral grid
    real(DP) :: DP_p            ! dummy variable
    real(DP) :: DP_z
    real(DP) :: DP_wix          ! zonal/meridional wind components in ICO-grid field
    real(DP) :: DP_wiy          ! zonal/meridional wind components in ICO-grid field
    real(DP) :: DP_tmp          ! presssure and temperature in ICO-grid field
    real(DP) :: DP_thetav       ! Virtual potential temperature in ICO-grid field
    real(DP) :: DP_phis
    real(DP) :: DP_ps
    real(DP) :: DP_rho          ! density and water vapor mixing ratio in ICO-grid field
    real(DP) :: DP_q            ! density and water vapor mixing ratio in ICO-grid field
    real(DP) :: cl, cl2

    real(RP) :: lat, lon
    real(RP) :: vx_local(kdim)
    real(RP) :: vy_local(kdim)
    real(RP) :: vz_local(kdim)
    real(RP) :: z       (kdim)
    real(RP) :: prs     (kdim)
    real(RP) :: tmp     (kdim)
    real(RP) :: wix     (kdim)
    real(RP) :: wiy     (kdim)
    real(RP) :: rho     (kdim)
    real(RP) :: q       (kdim)
    real(RP) :: thetav  (kdim)
    real(RP) :: ps

    real(RP) :: dpdz
    real(RP) :: rho0(ijdim,kdim)

    real(RP) :: VMTR_PHI   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_PHI_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer,  parameter :: deep    = 0      ! deep atmosphere (1 = yes or 0 = no)
    integer             :: moist            ! include moisture (1 = yes or 0 = no)
    integer             :: pertt            ! type of perturbation (0 = exponential, 1 = stream function)
    real(DP), parameter :: Xfact   = 1.0_DP ! Earth scaling parameter
    integer,  parameter :: zcoords = 1      ! 1 if z is specified, 0 if p is specified

    integer :: n, k, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    call VMTR_getIJ_PHI( VMTR_PHI, VMTR_PHI_pl )

    DIAG_var(:,:,:,:) = 0.0_RP

    DP_p = 0.0_DP

    select case( trim(test_case) )
    case ('1')  ! perturbation: exponential / with moisture
       if( IO_L ) write(IO_FID_LOG,*) "Moist Baroclinic Wave Initialize - case 1: perturbation: exponential / with moisture"
       moist = 1
       pertt = 0
    case ('2')  ! perturbation: stream function / with moisture
       if( IO_L ) write(IO_FID_LOG,*) "Moist Baroclinic Wave Initialize - case 2: perturbation: stream function / with moisture"
       moist = 1
       pertt = 1
    case ('3')  ! perturbation: exponential / without moisture
       if( IO_L ) write(IO_FID_LOG,*) "Moist Baroclinic Wave Initialize - case 3: perturbation: exponential / without moisture"
       moist = 0
       pertt = 0
    case ('4')  ! perturbation: stream function / without moisture
       if( IO_L ) write(IO_FID_LOG,*) "Moist Baroclinic Wave Initialize - case 4: perturbation: stream function / without moisture"
       moist = 0
       pertt = 1
    case ('5')  ! no perturbation / without moisture
       if( IO_L ) write(IO_FID_LOG,*) "Moist Baroclinic Wave Initialize - case 5: no perturbation / with moisture"
       moist = 1
       pertt = -99
    case ('6')  ! no perturbation / without moisture
       if( IO_L ) write(IO_FID_LOG,*) "Moist Baroclinic Wave Initialize - case 6: no perturbation / without moisture"
       moist = 0
       pertt = -99
    case default
       if( IO_L ) write(IO_FID_LOG,*) "xxx Invalid test_case: '"//trim(test_case)//"' specified."
       if( IO_L ) write(IO_FID_LOG,*) 'STOP.'
       call PRC_MPIstop
    end select
    if( IO_L ) write(IO_FID_LOG,*) "Chemical Tracer: ", chemtracer
    if( IO_L ) write(IO_FID_LOG,*) "### DO NOT INPUT ANY TOPOGRAPHY ###"

    if ( moist == 1 ) then
       if ( NQW_MAX < 3 ) then
          write(*         ,*) 'NQW_MAX is not enough! requires more than 3.', NQW_MAX
          if( IO_L ) write(IO_FID_LOG,*) 'NQW_MAX is not enough! requires more than 3.', NQW_MAX
          call PRC_MPIstop
       endif
    endif

    do l = 1, lall
    do n = 1, ijdim
       lat    = GRD_s(n,k0,l,GRD_LAT)
       lon    = GRD_s(n,k0,l,GRD_LON)
       DP_lat = real(GRD_s(n,k0,l,GRD_LAT),kind=DP)
       DP_lon = real(GRD_s(n,k0,l,GRD_LON),kind=DP)

       do k = 1, kdim
          DP_z = real(GRD_vz(n,k,l,GRD_Z),kind=DP)
          call baroclinic_wave_test( deep,      & ! [IN]
                                     moist,     & ! [IN]
                                     pertt,     & ! [IN]
                                     Xfact,     & ! [IN]
                                     DP_lon,    & ! [IN]
                                     DP_lat,    & ! [IN]
                                     DP_p,      & ! [INOUT]
                                     DP_z,      & ! [INOUT]
                                     zcoords,   & ! [IN]
                                     DP_wix,    & ! [OUT] Zonal wind (m s^-1)
                                     DP_wiy,    & ! [OUT] Meridional wind (m s^-1)
                                     DP_tmp,    & ! [OUT] Temperature (K)
                                     DP_thetav, & ! [OUT] Virtual potential Temperature (K)
                                     DP_phis,   & ! [OUT] Surface Geopotential (m^2 s^-2)
                                     DP_ps,     & ! [OUT] Surface Pressure (Pa)
                                     DP_rho,    & ! [OUT] density (kg m^-3)
                                     DP_q       ) ! [OUT] water vapor mixing ratio (kg/kg)

          z     (k) = GRD_vz(n,k,l,GRD_Z)
          wix   (k) = real( DP_wix   , kind=RP )
          wiy   (k) = real( DP_wiy   , kind=RP )
          tmp   (k) = real( DP_tmp   , kind=RP )
          thetav(k) = real( DP_thetav, kind=RP )
          rho   (k) = real( DP_rho   , kind=RP )
          q     (k) = real( max( DP_q / (1.0_DP+DP_q), 0.0_DP ), kind=RP ) ! fix negative value
       enddo
       ps = real( DP_ps, kind=RP )

       prs(1) = ps
       k      = ADM_kmin
       prs(k) = ps - rho(k) * g * ( GRD_vz(n,k,l,GRD_Z) - GRD_vz(n,k,l,GRD_ZH) )
       do k = ADM_kmin+1, ADM_kmax
          dpdz   = -g * ( rho(k)   * GRD_afac(k) &
                        + rho(k-1) * GRD_bfac(k) )
          prs(k) = prs(k-1) + dpdz * ( z(k)- z(k-1) )
       enddo

       do k = 1, kdim
          tmp(k) = prs(k) / ( rho(k) * ( rd * (1.0_RP-q(k)) + rv * q(k) ) )
       enddo
       rho0(n,:) = rho(:)

       call conv_vxvyvz  ( kdim, lat, lon, wix, wiy, vx_local, vy_local, vz_local )

       do k = 1, kdim
          DIAG_var(n,k,l,1     ) = prs     (k)
          DIAG_var(n,k,l,2     ) = tmp     (k)
          DIAG_var(n,k,l,3     ) = vx_local(k)
          DIAG_var(n,k,l,4     ) = vy_local(k)
          DIAG_var(n,k,l,5     ) = vz_local(k)
          DIAG_var(n,k,l,6+I_QV) = q       (k)
       enddo

       if ( chemtracer ) then
          if ( NCHEM_MAX /= 2 ) then
             write(*         ,*) 'NCHEM_MAX is not enough! requires 2.', NCHEM_MAX
             if( IO_L ) write(IO_FID_LOG,*) 'NCHEM_MAX is not enough! requires 2.', NCHEM_MAX
             call PRC_MPIstop
          endif

          call initial_value_Terminator( DP_lat*r2d, DP_lon*r2d, cl, cl2 )

          ! Todo : the mixing ratios are dry
          ! i.e. the ratio between the density of the species and the density of dry air.
          ! calc density at reference level
          DIAG_var(n,:,l,6+NCHEM_STR) = real(cl ,kind=RP) * ( 1.0_RP - q(:) )
          DIAG_var(n,:,l,6+NCHEM_END) = real(cl2,kind=RP) * ( 1.0_RP - q(:) )
       endif

    enddo
    call BNDCND_thermo(     &
         ijdim,             &
         rho0,              &
         DIAG_var(:,:,l,1), &
         DIAG_var(:,:,l,2), &
         VMTR_PHI(:,:,l)    )
    enddo

    return
  end subroutine jbw_moist_init

  !-----------------------------------------------------------------------------
  subroutine sc_init( &
       ijdim,       &
       kdim,        &
       lall,        &
       test_case,   &
       prs_rebuild, &
       DIAG_var     )
    use scale_process, only: &
       PRC_MPIstop
    use mod_grd, only: &
       GRD_LAT, &
       GRD_LON, &
       GRD_s,   &
       GRD_Z,   &
       GRD_vz
    use mod_runconf, only: &
       I_QV,    &
       NQW_MAX
    implicit none

    integer,          intent(in)  :: ijdim
    integer,          intent(in)  :: kdim
    integer,          intent(in)  :: lall
    character(len=*), intent(in)  :: test_case
    logical,          intent(in)  :: prs_rebuild
    real(RP),         intent(out) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(DP) :: DP_lon          ! latitude, longitude on Icosahedral grid
    real(DP) :: DP_lat          ! latitude, longitude on Icosahedral grid
    real(DP) :: DP_p            ! dummy variable
    real(DP) :: DP_z
    real(DP) :: DP_wix          ! zonal/meridional wind components in ICO-grid field
    real(DP) :: DP_wiy          ! zonal/meridional wind components in ICO-grid field
    real(DP) :: DP_tmp          ! presssure and temperature in ICO-grid field
    real(DP) :: DP_thetav       ! Virtual potential temperature in ICO-grid field
    real(DP) :: DP_ps
    real(DP) :: DP_rho          ! density and water vapor mixing ratio in ICO-grid field
    real(DP) :: DP_q            ! density and water vapor mixing ratio in ICO-grid field

    real(RP) :: lat, lon
    real(RP) :: vx_local(kdim)
    real(RP) :: vy_local(kdim)
    real(RP) :: vz_local(kdim)
    real(RP) :: z       (kdim)
    real(RP) :: prs     (kdim)
    real(RP) :: tmp     (kdim)
    real(RP) :: wix     (kdim)
    real(RP) :: wiy     (kdim)
    real(RP) :: rho     (kdim)
    real(RP) :: q       (kdim)
    real(RP) :: thetav  (kdim)
    real(RP) :: ps

    real(RP), parameter :: Mvap = 0.61_RP  ! Ratio of molar mass of dry air/water (imported from supercell_test.f90)
    real(RP)            :: RdovRv
    real(RP)            :: Mvap2           ! Ratio of molar mass of dry air/water based on NICAM CONSTANTs

    integer,  parameter :: zcoords = 1     ! 1 if z is specified, 0 if p is specified
    integer             :: pert            ! type of perturbation (0 = no perturbation, 1 = perturbation)
    logical,  parameter :: prs_dry = .false.

    integer :: n, k, l, k0
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.0_RP

    k0 = ADM_KNONE

    DP_p   = 0.0_DP
    RdovRv = Rd / Rv
    Mvap2  = ( 1.0_RP - RdovRv ) / RdovRv

    select case( trim(test_case) )
    case ('1')  ! with perturbation
       if( IO_L ) write(IO_FID_LOG,*) "Super-Cell Initialize - case 1: with perturbation"
       pert = 1
    case ('2')  ! without perturbation
       if( IO_L ) write(IO_FID_LOG,*) "Super-Cell Initialize - case 2: no perturbation"
       pert = 0
    case default
       if( IO_L ) write(IO_FID_LOG,*) "xxx Invalid test_case: '"//trim(test_case)//"' specified."
       if( IO_L ) write(IO_FID_LOG,*) 'STOP.'
       call PRC_MPIstop
    end select
    if( IO_L ) write(IO_FID_LOG,*) "### DO NOT INPUT ANY TOPOGRAPHY ###"

    if ( NQW_MAX < 3 ) then
       write(*         ,*) 'NQW_MAX is not enough! requires more than 3.', NQW_MAX
       if( IO_L ) write(IO_FID_LOG,*) 'NQW_MAX is not enough! requires more than 3.', NQW_MAX
       call PRC_MPIstop
    endif

    call supercell_init

    do l = 1, lall
    do n = 1, ijdim
       lat    = GRD_s(n,k0,l,GRD_LAT)
       lon    = GRD_s(n,k0,l,GRD_LON)
       DP_lat = real(GRD_s(n,k0,l,GRD_LAT),kind=DP)
       DP_lon = real(GRD_s(n,k0,l,GRD_LON),kind=DP)

       do k = 1, kdim
          DP_z = real(GRD_vz(n,k,l,GRD_Z),kind=DP)
          call supercell_test( DP_lon,    & ! [IN]
                               DP_lat,    & ! [IN]
                               DP_p,      & ! [INOUT]
                               DP_z,      & ! [INOUT]
                               zcoords,   & ! [IN]
                               DP_wix,    & ! [OUT] Zonal wind (m s^-1)
                               DP_wiy,    & ! [OUT] Meridional wind (m s^-1)
                               DP_tmp,    & ! [OUT] Temperature (K)
                               DP_thetav, & ! [OUT] Virtual potential Temperature (K)
                               DP_ps,     & ! [OUT] Surface Pressure (Pa)
                               DP_rho,    & ! [OUT] density (kg m^-3)
                               DP_q,      & ! [OUT] water vapor mixing ratio (kg/kg)
                               pert       ) ! [IN]  perturbation switch


          z     (k) = GRD_vz(n,k,l,GRD_Z)
          wix   (k) = real( DP_wix             , kind=RP )
          wiy   (k) = real( DP_wiy             , kind=RP )
          tmp   (k) = real( DP_tmp             , kind=RP )
          thetav(k) = real( DP_thetav          , kind=RP )
          rho   (k) = real( DP_rho             , kind=RP )
          q     (k) = real( max( DP_q, 0.0_DP ), kind=RP ) ! fix negative value

          ! force zero for q upper tropopause (12km)
          if ( z(k) > 12000.0_RP ) q(k) = 0.0_RP
       enddo
       ps = real( DP_ps, kind=RP )

       ! Re-Evaluation of temperature from virtual temperature
       tmp(:) = tmp(:) * ( (1.0+RP + Mvap * q(:) ) / ( 1.0_RP + Mvap2 * q(:) ) )

       call diag_pressure( kdim, z, rho, tmp, q, ps, prs, prs_rebuild, prs_dry )
       call conv_vxvyvz ( kdim, lat, lon, wix, wiy, vx_local, vy_local, vz_local )

       do k = 1, kdim
          DIAG_var(n,k,l,1     ) = prs     (k)
          DIAG_var(n,k,l,2     ) = tmp     (k)
          DIAG_var(n,k,l,3     ) = vx_local(k)
          DIAG_var(n,k,l,4     ) = vy_local(k)
          DIAG_var(n,k,l,5     ) = vz_local(k)
          DIAG_var(n,k,l,6+I_QV) = q       (k)
       enddo

    enddo
    enddo

    return
  end subroutine sc_init

  !-----------------------------------------------------------------------------
  subroutine tc_init( &
       ijdim,       &
       kdim,        &
       lall,        &
       prs_rebuild, &
       DIAG_var     )
    use scale_process, only: &
       PRC_MPIstop
    use mod_grd, only: &
       GRD_LAT, &
       GRD_LON, &
       GRD_s,   &
       GRD_Z,   &
       GRD_vz
    use mod_runconf, only: &
       I_QV,    &
       NQW_MAX
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: lall
    logical,  intent(in)  :: prs_rebuild
    real(RP), intent(out) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(DP) :: DP_lon          ! latitude, longitude on Icosahedral grid
    real(DP) :: DP_lat          ! latitude, longitude on Icosahedral grid
    real(DP) :: DP_p            ! dummy variable
    real(DP) :: DP_z
    real(DP) :: DP_wix          ! zonal/meridional wind components in ICO-grid field
    real(DP) :: DP_wiy          ! zonal/meridional wind components in ICO-grid field
    real(DP) :: DP_tmp          ! presssure and temperature in ICO-grid field
    real(DP) :: DP_thetav       ! Virtual potential temperature in ICO-grid field
    real(DP) :: DP_phis
    real(DP) :: DP_ps
    real(DP) :: DP_rho          ! density and water vapor mixing ratio in ICO-grid field
    real(DP) :: DP_q            ! density and water vapor mixing ratio in ICO-grid field

    real(RP) :: lat, lon
    real(RP) :: vx_local(kdim)
    real(RP) :: vy_local(kdim)
    real(RP) :: vz_local(kdim)
    real(RP) :: z       (kdim)
    real(RP) :: prs     (kdim)
    real(RP) :: tmp     (kdim)
    real(RP) :: wix     (kdim)
    real(RP) :: wiy     (kdim)
    real(RP) :: rho     (kdim)
    real(RP) :: q       (kdim)
    real(RP) :: thetav  (kdim)
    real(RP) :: ps

    real(RP), parameter :: Mvap = 0.608_RP ! Ratio of molar mass of dry air/water (imported from tropical_cyclone_test.f90)
    real(RP)            :: RdovRv
    real(RP)            :: Mvap2           ! Ratio of molar mass of dry air/water based on NICAM CONSTANTs

    integer,  parameter :: zcoords = 1     ! 1 if z is specified, 0 if p is specified
    logical,  parameter :: prs_dry = .false.

    integer :: n, k, l, k0
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.0_RP

    k0 = ADM_KNONE

    DP_p   = 0.0_DP
    RdovRv = Rd / Rv
    Mvap2  = ( 1.0_RP - RdovRv ) / RdovRv

    if( IO_L ) write(IO_FID_LOG,*) "### DO NOT INPUT ANY TOPOGRAPHY ###"
    if ( NQW_MAX < 3 ) then
       write(*         ,*) 'NQW_MAX is not enough! requires more than 3.', NQW_MAX
       if( IO_L ) write(IO_FID_LOG,*) 'NQW_MAX is not enough! requires more than 3.', NQW_MAX
       call PRC_MPIstop
    endif

    do l = 1, lall
    do n = 1, ijdim
       lat    = GRD_s(n,k0,l,GRD_LAT)
       lon    = GRD_s(n,k0,l,GRD_LON)
       DP_lat = real(GRD_s(n,k0,l,GRD_LAT),kind=DP)
       DP_lon = real(GRD_s(n,k0,l,GRD_LON),kind=DP)

       do k = 1, kdim
          DP_z = real(GRD_vz(n,k,l,GRD_Z),kind=DP)
          call tropical_cyclone_test( DP_lon,    & ! [IN]
                                      DP_lat,    & ! [IN]
                                      DP_p,      & ! [INOUT]
                                      DP_z,      & ! [INOUT]
                                      zcoords,   & ! [IN]
                                      DP_wix,    & ! [OUT] Zonal wind (m s^-1)
                                      DP_wiy,    & ! [OUT] Meridional wind (m s^-1)
                                      DP_tmp,    & ! [OUT] Temperature (K)
                                      DP_thetav, & ! [OUT] Virtual potential Temperature (K)
                                      DP_phis,   & ! [OUT] Surface Geopotential (m^2 s^-2)
                                      DP_ps,     & ! [OUT] Surface Pressure (Pa)
                                      DP_rho,    & ! [OUT] density (kg m^-3)
                                      DP_q       ) ! [OUT] water vapor mixing ratio (kg/kg)

          z     (k) = GRD_vz(n,k,l,GRD_Z)
          wix   (k) = real( DP_wix             , kind=RP )
          wiy   (k) = real( DP_wiy             , kind=RP )
          tmp   (k) = real( DP_tmp             , kind=RP )
          thetav(k) = real( DP_thetav          , kind=RP )
          rho   (k) = real( DP_rho             , kind=RP )
          q     (k) = real( max( DP_q, 0.0_DP ), kind=RP ) ! fix negative value
       enddo
       ps = real( DP_ps, kind=RP )

       ! Re-Evaluation of temperature from virtual temperature
       tmp(:) = tmp(:) * ( (1.0+RP + Mvap * q(:) ) / ( 1.0_RP + Mvap2 * q(:) ) )

       call diag_pressure( kdim, z, rho, tmp, q, ps, prs, prs_rebuild, prs_dry )
       call conv_vxvyvz ( kdim, lat, lon, wix, wiy, vx_local, vy_local, vz_local )

       do k = 1, kdim
          DIAG_var(n,k,l,1     ) = prs     (k)
          DIAG_var(n,k,l,2     ) = tmp     (k)
          DIAG_var(n,k,l,3     ) = vx_local(k)
          DIAG_var(n,k,l,4     ) = vy_local(k)
          DIAG_var(n,k,l,5     ) = vz_local(k)
          DIAG_var(n,k,l,6+I_QV) = q       (k)
       enddo

    enddo
    enddo

    return
  end subroutine tc_init

  !-----------------------------------------------------------------------------
  subroutine tracer_init( &
       ijdim,      &
       kdim,       &
       lall,       &
       test_case,  &
       DIAG_var    )
    use scale_process, only: &
       PRC_MPIstop
    use mod_grd, only: &
       GRD_LAT, &
       GRD_LON, &
       GRD_s,   &
       GRD_gz,  &
       GRD_Z,   &
       GRD_ZH,  &
       GRD_vz
    use mod_runconf, only: &
       NCHEM_STR
    use mod_chemvar, only: &
       chemvar_getid
    use dcmip_initial_conditions_test_1_2_3, only: &
       test1_advection_deformation, &
       test1_advection_hadley,      &
       test1_advection_orography
    implicit none

    integer,          intent(in)    :: ijdim
    integer,          intent(in)    :: kdim
    integer,          intent(in)    :: lall
    character(len=*), intent(in)    :: test_case
    real(RP),         intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(RP) :: lon      ! longitude            [rad]
    real(RP) :: lat      ! latitude             [rad]
    real(RP) :: z(kdim)  ! Height               [m]
    real(RP) :: p(kdim)  ! pressure             [Pa]
    real(RP) :: u(kdim)  ! zonal      wind      [m/s]
    real(RP) :: v(kdim)  ! meridional wind      [m/s]
    real(RP) :: w(kdim)  ! vertical   wind      [m/s]
    real(RP) :: t(kdim)  ! temperature          [K]
    real(RP) :: phis     ! surface geopotential [m2/s2], not in use
    real(RP) :: ps       ! surface pressure     [Pa]   , not in use
    real(RP) :: rho      ! density              [kg/m3], not in use
    real(RP) :: q        ! specific humidity    [kg/kg], not in use
    real(RP) :: q1(kdim) ! passive tracer       [kg/kg]
    real(RP) :: q2(kdim) ! passive tracer       [kg/kg]
    real(RP) :: q3(kdim) ! passive tracer       [kg/kg]
    real(RP) :: q4(kdim) ! passive tracer       [kg/kg]
    real(RP) :: vx(kdim)
    real(RP) :: vy(kdim)
    real(RP) :: vz(kdim)

    real(DP) :: DP_lon  ! longitude            [rad]
    real(DP) :: DP_lat  ! latitude             [rad]
    real(DP) :: DP_z    ! Height               [m]
    real(DP) :: DP_p    ! pressure             [Pa]
    real(DP) :: DP_u    ! zonal      wind      [m/s]
    real(DP) :: DP_v    ! meridional wind      [m/s]
    real(DP) :: DP_w    ! vertical   wind      [m/s]
    real(DP) :: DP_t    ! temperature          [K]
    real(DP) :: DP_phis ! surface geopotential [m2/s2], not in use
    real(DP) :: DP_ps   ! surface pressure     [Pa]   , not in use
    real(DP) :: DP_rho  ! density              [kg/m3], not in use
    real(DP) :: DP_q    ! specific humidity    [kg/kg], not in use
    real(DP) :: DP_q1   ! passive tracer       [kg/kg]
    real(DP) :: DP_q2   ! passive tracer       [kg/kg]
    real(DP) :: DP_q3   ! passive tracer       [kg/kg]
    real(DP) :: DP_q4   ! passive tracer       [kg/kg]

    logical, parameter :: hybrid_eta = .false. ! dont use hybrid sigma-p (eta) coordinate
    integer, parameter :: zcoords    = 1       ! if zcoords = 1, then we use z and output p
    integer, parameter :: cfv        = 2       ! if cfv = 2 then our velocities follow
                                               ! Gal-Chen coordinates and we need to specify w
    real(DP)           :: hyam       = 0.0_DP  ! dont use hybrid sigma-p (eta) coordinate
    real(DP)           :: hybm       = 0.0_DP  ! dont use hybrid sigma-p (eta) coordinate
    real(DP)           :: DP_gc                ! bar{z} for Gal-Chen coordinate

    integer :: I_pasv1, I_pasv2
    integer :: I_pasv3, I_pasv4
    integer :: n, k, l, k0
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.0_RP

    k0 = ADM_KNONE

    I_pasv1 = 6 + NCHEM_STR + chemvar_getid( "passive001" ) - 1
    I_pasv2 = 6 + NCHEM_STR + chemvar_getid( "passive002" ) - 1
    I_pasv3 = 6 + NCHEM_STR + chemvar_getid( "passive003" ) - 1
    I_pasv4 = 6 + NCHEM_STR + chemvar_getid( "passive004" ) - 1

    select case(test_case)
    case ('1', '1-1') ! DCMIP 2012 Test 1-1: 3D Deformational Flow

       do l = 1, lall
       do n = 1, ijdim
          z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
          do k = ADM_kmin, ADM_kmax+1
             z(k) = GRD_vz(n,k,l,GRD_Z)
          enddo
          p(:) = 0.0_RP

          lat = GRD_s(n,k0,l,GRD_LAT)
          lon = GRD_s(n,k0,l,GRD_LON)

          do k = 1, kdim
             DP_lon = real(lon ,kind=DP)
             DP_lat = real(lat ,kind=DP)
             DP_p   = real(p(k),kind=DP)
             DP_z   = real(z(k),kind=DP)

             call test1_advection_deformation( DP_lon , & ! [IN]
                                               DP_lat , & ! [IN]
                                               DP_p   , & ! [INOUT]
                                               DP_z   , & ! [IN]
                                               zcoords, & ! [IN]
                                               DP_u   , & ! [OUT]
                                               DP_v   , & ! [OUT]
                                               DP_w   , & ! [OUT]
                                               DP_t   , & ! [OUT]
                                               DP_phis, & ! [OUT]
                                               DP_ps  , & ! [OUT]
                                               DP_rho , & ! [OUT]
                                               DP_q   , & ! [OUT]
                                               DP_q1  , & ! [OUT]
                                               DP_q2  , & ! [OUT]
                                               DP_q3  , & ! [OUT]
                                               DP_q4    ) ! [OUT]

             p(k)  = real(DP_p   ,kind=RP)
             u(k)  = real(DP_u   ,kind=RP)
             v(k)  = real(DP_v   ,kind=RP)
             w(k)  = real(DP_w   ,kind=RP)
             t(k)  = real(DP_t   ,kind=RP)
             phis  = real(DP_phis,kind=RP)
             ps    = real(DP_ps  ,kind=RP)
             rho   = real(DP_rho ,kind=RP)
             q     = real(DP_q   ,kind=RP)
             q1(k) = real(DP_q1  ,kind=RP)
             q2(k) = real(DP_q2  ,kind=RP)
             q3(k) = real(DP_q3  ,kind=RP)
             q4(k) = real(DP_q4  ,kind=RP)
          enddo

          call conv_vxvyvz( kdim, lat, lon, u(:), v(:), vx(:), vy(:), vz(:) )

          do k = 1, kdim
             DIAG_var(n,k,l,1) = p (k)
             DIAG_var(n,k,l,2) = t (k)
             DIAG_var(n,k,l,3) = vx(k)
             DIAG_var(n,k,l,4) = vy(k)
             DIAG_var(n,k,l,5) = vz(k)
             DIAG_var(n,k,l,6) = w (k)

             DIAG_var(n,k,l,I_pasv1) = q1(k)
             DIAG_var(n,k,l,I_pasv2) = q2(k)
             DIAG_var(n,k,l,I_pasv3) = q3(k)
             DIAG_var(n,k,l,I_pasv4) = q4(k)
          enddo
       enddo
       enddo

    case ('2', '1-2') ! DCMIP 2012 Test 1-2: Hadley-like Meridional Circulation

       do l = 1, lall
       do n = 1, ijdim
          z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
          do k = ADM_kmin, ADM_kmax+1
             z(k) = GRD_vz(n,k,l,GRD_Z)
          enddo
          p(:) = 0.0_RP

          lat = GRD_s(n,k0,l,GRD_LAT)
          lon = GRD_s(n,k0,l,GRD_LON)

          do k = 1, kdim
             DP_lon = real(lon ,kind=DP)
             DP_lat = real(lat ,kind=DP)
             DP_p   = real(p(k),kind=DP)
             DP_z   = real(z(k),kind=DP)

             call test1_advection_hadley( DP_lon , & ! [IN]
                                          DP_lat , & ! [IN]
                                          DP_p   , & ! [INOUT]
                                          DP_z   , & ! [IN]
                                          zcoords, & ! [IN]
                                          DP_u   , & ! [OUT]
                                          DP_v   , & ! [OUT]
                                          DP_w   , & ! [OUT]
                                          DP_t   , & ! [OUT]
                                          DP_phis, & ! [OUT]
                                          DP_ps  , & ! [OUT]
                                          DP_rho , & ! [OUT]
                                          DP_q   , & ! [OUT]
                                          DP_q1    ) ! [OUT]

             p(k)  = real(DP_p   ,kind=RP)
             u(k)  = real(DP_u   ,kind=RP)
             v(k)  = real(DP_v   ,kind=RP)
             w(k)  = real(DP_w   ,kind=RP)
             t(k)  = real(DP_t   ,kind=RP)
             phis  = real(DP_phis,kind=RP)
             ps    = real(DP_ps  ,kind=RP)
             rho   = real(DP_rho ,kind=RP)
             q     = real(DP_q   ,kind=RP)
             q1(k) = real(DP_q1  ,kind=RP)
             q2(k) = 0.0_RP
             q3(k) = 0.0_RP
             q4(k) = 0.0_RP
          enddo

          call conv_vxvyvz( kdim, lat, lon, u(:), v(:), vx(:), vy(:), vz(:) )

          do k = 1, kdim
             DIAG_var(n,k,l,1) = p (k)
             DIAG_var(n,k,l,2) = t (k)
             DIAG_var(n,k,l,3) = vx(k)
             DIAG_var(n,k,l,4) = vy(k)
             DIAG_var(n,k,l,5) = vz(k)
             DIAG_var(n,k,l,6) = w (k)

             DIAG_var(n,k,l,I_pasv1) = q1(k)
             DIAG_var(n,k,l,I_pasv2) = q2(k)
             DIAG_var(n,k,l,I_pasv3) = q3(k)
             DIAG_var(n,k,l,I_pasv4) = q4(k)
          enddo
       enddo
       enddo

    case ('3', '1-3') ! DCMIP 2012 Test 1-3: Horizontal advection of thin cloud-like tracers in the presence of orography

       do l = 1, lall
       do n = 1, ijdim
          z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
          do k = ADM_kmin, ADM_kmax+1
             z(k) = GRD_vz(n,k,l,GRD_Z)
          enddo
          p(:) = 0.0_RP

          lat = GRD_s(n,k0,l,GRD_LAT)
          lon = GRD_s(n,k0,l,GRD_LON)

          do k = 1, kdim
             DP_gc = real(GRD_gz(k),kind=DP)

             DP_lon = real(lon ,kind=DP)
             DP_lat = real(lat ,kind=DP)
             DP_p   = real(p(k),kind=DP)
             DP_z   = real(z(k),kind=DP)

             call test1_advection_orography( DP_lon    , & ! [IN]
                                             DP_lat    , & ! [IN]
                                             DP_p      , & ! [INOUT]
                                             DP_z      , & ! [IN]
                                             zcoords   , & ! [IN]
                                             cfv,        & ! [IN]
                                             hybrid_eta, & ! [IN]
                                             hyam,       & ! [IN]
                                             hybm,       & ! [IN]
                                             DP_gc     , & ! [IN]
                                             DP_u      , & ! [OUT]
                                             DP_v      , & ! [OUT]
                                             DP_w      , & ! [OUT]
                                             DP_t      , & ! [OUT]
                                             DP_phis   , & ! [OUT]
                                             DP_ps     , & ! [OUT]
                                             DP_rho    , & ! [OUT]
                                             DP_q      , & ! [OUT]
                                             DP_q1     , & ! [OUT]
                                             DP_q2     , & ! [OUT]
                                             DP_q3     , & ! [OUT]
                                             DP_q4       ) ! [OUT]

             p(k)  = real(DP_p   ,kind=RP)
             u(k)  = real(DP_u   ,kind=RP)
             v(k)  = real(DP_v   ,kind=RP)
             w(k)  = real(DP_w   ,kind=RP)
             t(k)  = real(DP_t   ,kind=RP)
             phis  = real(DP_phis,kind=RP)
             ps    = real(DP_ps  ,kind=RP)
             rho   = real(DP_rho ,kind=RP)
             q     = real(DP_q   ,kind=RP)
             q1(k) = real(DP_q1  ,kind=RP)
             q2(k) = real(DP_q2  ,kind=RP)
             q3(k) = real(DP_q3  ,kind=RP)
             q4(k) = real(DP_q4  ,kind=RP)
          enddo

          call conv_vxvyvz( kdim, lat, lon, u(:), v(:), vx(:), vy(:), vz(:) )

          do k = 1, kdim
             DIAG_var(n,k,l,1) = p (k)
             DIAG_var(n,k,l,2) = t (k)
             DIAG_var(n,k,l,3) = vx(k)
             DIAG_var(n,k,l,4) = vy(k)
             DIAG_var(n,k,l,5) = vz(k)
             DIAG_var(n,k,l,6) = w (k)

             DIAG_var(n,k,l,I_pasv1) = q1(k)
             DIAG_var(n,k,l,I_pasv2) = q2(k)
             DIAG_var(n,k,l,I_pasv3) = q3(k)
             DIAG_var(n,k,l,I_pasv4) = q4(k)
          enddo
       enddo
       enddo

    case default
       write(*         ,*) "Unknown test_case: ", trim(test_case)," specified. STOP"
       if( IO_L ) write(IO_FID_LOG,*) "Unknown test_case: ", trim(test_case)," specified. STOP"
       call PRC_MPIstop
    end select

    return
  end subroutine tracer_init

  !-----------------------------------------------------------------------------
  subroutine mountwave_init( &
       ijdim,     &
       kdim,      &
       lall,      &
       test_case, &
       DIAG_var   )
    use scale_process, only: &
       PRC_MPIstop
    use mod_grd, only: &
       GRD_LAT, &
       GRD_LON, &
       GRD_s,   &
       GRD_vz,  &
       GRD_Z,   &
       GRD_ZH
    use mod_runconf, only: &
       NCHEM_STR
    use mod_chemvar, only: &
       chemvar_getid
    implicit none

    integer,          intent(in)    :: ijdim
    integer,          intent(in)    :: kdim
    integer,          intent(in)    :: lall
    character(len=*), intent(in)    :: test_case
    real(RP),         intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    real(RP) :: lon      ! longitude            [rad]
    real(RP) :: lat      ! latitude             [rad]
    real(RP) :: z(kdim)  ! Height               [m]
    real(RP) :: p(kdim)  ! pressure             [Pa]
    real(RP) :: u(kdim)  ! zonal      wind      [m/s]
    real(RP) :: v(kdim)  ! meridional wind      [m/s]
    real(RP) :: w(kdim)  ! vertical   wind      [m/s]
    real(RP) :: t(kdim)  ! temperature          [K]
    real(RP) :: phis     ! surface geopotential [m2/s2], not in use
    real(RP) :: ps       ! surface pressure     [Pa]   , not in use
    real(RP) :: rho      ! density              [kg/m3], not in use
    real(RP) :: q(kdim)  ! specific humidity    [kg/kg], not in use
    real(RP) :: vx(kdim)
    real(RP) :: vy(kdim)
    real(RP) :: vz(kdim)

    real(DP) :: DP_lon  ! longitude            [rad]
    real(DP) :: DP_lat  ! latitude             [rad]
    real(DP) :: DP_z    ! Height               [m]
    real(DP) :: DP_p    ! pressure             [Pa]
    real(DP) :: DP_u    ! zonal      wind      [m/s]
    real(DP) :: DP_v    ! meridional wind      [m/s]
    real(DP) :: DP_w    ! vertical   wind      [m/s]
    real(DP) :: DP_t    ! temperature          [K]
    real(DP) :: DP_phis ! surface geopotential [m2/s2], not in use
    real(DP) :: DP_ps   ! surface pressure     [Pa]   , not in use
    real(DP) :: DP_rho  ! density              [kg/m3], not in use
    real(DP) :: DP_q    ! specific humidity    [kg/kg], not in use

    integer :: I_pasv1

    real(DP) :: hyam       = 0.0_DP
    real(DP) :: hybm       = 0.0_DP
    logical  :: hybrid_eta = .false.

    integer, parameter :: zcoords = 1
    integer  :: shear

    integer :: n, l, k, k0
    !---------------------------------------------------------------------------

    DIAG_var(:,:,:,:) = 0.0_RP

    k0 = ADM_KNONE

    I_pasv1 = 6 + chemvar_getid( "passive001" ) + NCHEM_STR - 1

    select case( test_case )
    case ('0', '2-0') ! DCMIP: TEST CASE 2-0 - Steady-State Atmosphere at Rest in the Presence of Orography

       do l = 1, lall
       do n = 1, ijdim
          z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
          do k = ADM_kmin, ADM_kmax+1
             z(k) = GRD_vz(n,k,l,GRD_Z)
          enddo
          p(:) = 0.0_RP

          lat = GRD_s(n,k0,l,GRD_LAT)
          lon = GRD_s(n,k0,l,GRD_LON)

          do k = 1, kdim
             DP_lon = real(lon ,kind=DP)
             DP_lat = real(lat ,kind=DP)
             DP_p   = real(p(k),kind=DP)
             DP_z   = real(z(k),kind=DP)

             call test2_steady_state_mountain( DP_lon    , & ! [IN]
                                               DP_lat    , & ! [IN]
                                               DP_p      , & ! [INOUT]
                                               DP_z      , & ! [IN]
                                               zcoords   , & ! [IN]
                                               hybrid_eta, & ! [IN]
                                               hyam,       & ! [IN]
                                               hybm,       & ! [IN]
                                               DP_u      , & ! [OUT]
                                               DP_v      , & ! [OUT]
                                               DP_w      , & ! [OUT]
                                               DP_t      , & ! [OUT]
                                               DP_phis   , & ! [OUT]
                                               DP_ps     , & ! [OUT]
                                               DP_rho    , & ! [OUT]
                                               DP_q        ) ! [OUT]

             p(k)  = real(DP_p   ,kind=RP)
             u(k)  = real(DP_u   ,kind=RP)
             v(k)  = real(DP_v   ,kind=RP)
             w(k)  = real(DP_w   ,kind=RP)
             t(k)  = real(DP_t   ,kind=RP)
             phis  = real(DP_phis,kind=RP)
             ps    = real(DP_ps  ,kind=RP)
             rho   = real(DP_rho ,kind=RP)
             q(k)  = real(DP_q   ,kind=RP)
          enddo

          call conv_vxvyvz( kdim, lat, lon, u(:), v(:), vx(:), vy(:), vz(:) )

          do k = 1, kdim
             DIAG_var(n,k,l,1) = p (k)
             DIAG_var(n,k,l,2) = t (k)
             DIAG_var(n,k,l,3) = vx(k)
             DIAG_var(n,k,l,4) = vy(k)
             DIAG_var(n,k,l,5) = vz(k)
             DIAG_var(n,k,l,6) = w (k)

             DIAG_var(n,k,l,I_pasv1) = q(k)
          enddo
       enddo
       enddo

    case ('1', '2-1') ! DCMIP: TEST CASE 2-1 - Non-hydrostatic Mountain Waves over a Schaer-type Mountain
       shear = 0 ! constant u

       do l = 1, lall
       do n = 1, ijdim
          z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
          do k = ADM_kmin, ADM_kmax+1
             z(k) = GRD_vz(n,k,l,GRD_Z)
          enddo
          p(:) = 0.0_RP

          lat = GRD_s(n,k0,l,GRD_LAT)
          lon = GRD_s(n,k0,l,GRD_LON)

          do k = 1, kdim
             DP_lon = real(lon ,kind=DP)
             DP_lat = real(lat ,kind=DP)
             DP_p   = real(p(k),kind=DP)
             DP_z   = real(z(k),kind=DP)

             call test2_steady_state_mountain( DP_lon    , & ! [IN]
                                               DP_lat    , & ! [IN]
                                               DP_p      , & ! [INOUT]
                                               DP_z      , & ! [IN]
                                               zcoords   , & ! [IN]
                                               hybrid_eta, & ! [IN]
                                               hyam,       & ! [IN]
                                               hybm,       & ! [IN]
                                               DP_u      , & ! [OUT]
                                               DP_v      , & ! [OUT]
                                               DP_w      , & ! [OUT]
                                               DP_t      , & ! [OUT]
                                               DP_phis   , & ! [OUT]
                                               DP_ps     , & ! [OUT]
                                               DP_rho    , & ! [OUT]
                                               DP_q        ) ! [OUT]

             p(k)  = real(DP_p   ,kind=RP)
             u(k)  = real(DP_u   ,kind=RP)
             v(k)  = real(DP_v   ,kind=RP)
             w(k)  = real(DP_w   ,kind=RP)
             t(k)  = real(DP_t   ,kind=RP)
             phis  = real(DP_phis,kind=RP)
             ps    = real(DP_ps  ,kind=RP)
             rho   = real(DP_rho ,kind=RP)
             q(k)  = real(DP_q   ,kind=RP)
          enddo

          call conv_vxvyvz( kdim, lat, lon, u(:), v(:), vx(:), vy(:), vz(:) )

          do k = 1, kdim
             DIAG_var(n,k,l,1) = p (k)
             DIAG_var(n,k,l,2) = t (k)
             DIAG_var(n,k,l,3) = vx(k)
             DIAG_var(n,k,l,4) = vy(k)
             DIAG_var(n,k,l,5) = vz(k)
             DIAG_var(n,k,l,6) = w (k)

             DIAG_var(n,k,l,I_pasv1) = q(k)
          enddo
       enddo
       enddo

    case ('2', '2-2') ! DCMIP: TEST CASE 2-2 - Non-hydrostatic Mountain Waves over a Schaer-type Mountain
       shear = 1 ! sheared u

       do l = 1, lall
       do n = 1, ijdim
          z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
          do k = ADM_kmin, ADM_kmax+1
             z(k) = GRD_vz(n,k,l,GRD_Z)
          enddo
          p(:) = 0.0_RP

          lat = GRD_s(n,k0,l,GRD_LAT)
          lon = GRD_s(n,k0,l,GRD_LON)

          do k = 1, kdim
             DP_lon = real(lon ,kind=DP)
             DP_lat = real(lat ,kind=DP)
             DP_p   = real(p(k),kind=DP)
             DP_z   = real(z(k),kind=DP)

             call test2_steady_state_mountain( DP_lon    , & ! [IN]
                                               DP_lat    , & ! [IN]
                                               DP_p      , & ! [INOUT]
                                               DP_z      , & ! [IN]
                                               zcoords   , & ! [IN]
                                               hybrid_eta, & ! [IN]
                                               hyam,       & ! [IN]
                                               hybm,       & ! [IN]
                                               DP_u      , & ! [OUT]
                                               DP_v      , & ! [OUT]
                                               DP_w      , & ! [OUT]
                                               DP_t      , & ! [OUT]
                                               DP_phis   , & ! [OUT]
                                               DP_ps     , & ! [OUT]
                                               DP_rho    , & ! [OUT]
                                               DP_q        ) ! [OUT]

             p(k)  = real(DP_p   ,kind=RP)
             u(k)  = real(DP_u   ,kind=RP)
             v(k)  = real(DP_v   ,kind=RP)
             w(k)  = real(DP_w   ,kind=RP)
             t(k)  = real(DP_t   ,kind=RP)
             phis  = real(DP_phis,kind=RP)
             ps    = real(DP_ps  ,kind=RP)
             rho   = real(DP_rho ,kind=RP)
             q(k)  = real(DP_q   ,kind=RP)
          enddo

          call conv_vxvyvz( kdim, lat, lon, u(:), v(:), vx(:), vy(:), vz(:) )

          do k = 1, kdim
             DIAG_var(n,k,l,1) = p (k)
             DIAG_var(n,k,l,2) = t (k)
             DIAG_var(n,k,l,3) = vx(k)
             DIAG_var(n,k,l,4) = vy(k)
             DIAG_var(n,k,l,5) = vz(k)
             DIAG_var(n,k,l,6) = w (k)

             DIAG_var(n,k,l,I_pasv1) = q(k)
          enddo
       enddo
       enddo

    case default
       write(*         ,*) "Unknown test_case: ", trim(test_case)," specified. STOP"
       if( IO_L ) write(IO_FID_LOG,*) "Unknown test_case: ", trim(test_case)," specified. STOP"
       call PRC_MPIstop
    end select

    return
  end subroutine mountwave_init

  !-----------------------------------------------------------------------------
  subroutine gravwave_init( &
     ijdim,   &
     kdim,    &
     lall,    &
     DIAG_var )
    use mod_grd, only: &
       GRD_LAT, &
       GRD_LON, &
       GRD_s,   &
       GRD_Z,   &
       GRD_ZH,  &
       GRD_vz
    implicit none

    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim
    integer, intent(in)    :: lall
    real(RP), intent(inout) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    integer, parameter :: zcoords = 1

    real(RP) :: lon     ! longitude            [rad]
    real(RP) :: lat     ! latitude             [rad]
    real(RP) :: z(kdim) ! Height               [m]
    real(RP) :: p(kdim) ! pressure             [Pa]
    real(RP) :: u(kdim) ! zonal      wind      [m/s]
    real(RP) :: v(kdim) ! meridional wind      [m/s]
    real(RP) :: w(kdim) ! vertical   wind      [m/s]
    real(RP) :: t(kdim) ! temperature          [K]
    real(RP) :: phis    ! surface geopotential [m2/s2], not in use
    real(RP) :: ps      ! surface pressure     [Pa]   , not in use
    real(RP) :: rho     ! density              [kg/m3], not in use
    real(RP) :: q       ! specific humidity    [kg/kg], not in use
    real(RP) :: vx(kdim)
    real(RP) :: vy(kdim)
    real(RP) :: vz(kdim)

    real(DP) :: DP_lon  ! longitude            [rad]
    real(DP) :: DP_lat  ! latitude             [rad]
    real(DP) :: DP_z    ! Height               [m]
    real(DP) :: DP_p    ! pressure             [Pa]
    real(DP) :: DP_u    ! zonal      wind      [m/s]
    real(DP) :: DP_v    ! meridional wind      [m/s]
    real(DP) :: DP_w    ! vertical   wind      [m/s]
    real(DP) :: DP_t    ! temperature          [K]
    real(DP) :: DP_phis ! surface geopotential [m2/s2], not in use
    real(DP) :: DP_ps   ! surface pressure     [Pa]   , not in use
    real(DP) :: DP_rho  ! density              [kg/m3], not in use
    real(DP) :: DP_q    ! specific humidity    [kg/kg], not in use

    integer :: n, k, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    DIAG_var(:,:,:,:) = 0.0_RP

    do l = 1, lall
    do n = 1, ijdim
       z(ADM_kmin-1) = GRD_vz(n,ADM_kmin,l,GRD_ZH)
       do k = ADM_kmin, ADM_kmax+1
          z(k) = GRD_vz(n,k,l,GRD_Z)
       enddo
       p(:) = 0.0_RP

       lat = GRD_s(n,k0,l,GRD_LAT)
       lon = GRD_s(n,k0,l,GRD_LON)

       do k = 1, kdim
          DP_lon = real(lon ,kind=DP)
          DP_lat = real(lat ,kind=DP)
          DP_p   = real(p(k),kind=DP)
          DP_z   = real(z(k),kind=DP)

          call test3_gravity_wave( DP_lon , & ! [IN]
                                   DP_lat , & ! [IN]
                                   DP_p   , & ! [INOUT]
                                   DP_z   , & ! [IN]
                                   zcoords, & ! [IN]
                                   DP_u   , & ! [OUT]
                                   DP_v   , & ! [OUT]
                                   DP_w   , & ! [OUT]
                                   DP_t   , & ! [OUT]
                                   DP_phis, & ! [OUT]
                                   DP_ps  , & ! [OUT]
                                   DP_rho , & ! [OUT]
                                   DP_q     ) ! [OUT]

          p(k)  = real(DP_p   ,kind=RP)
          u(k)  = real(DP_u   ,kind=RP)
          v(k)  = real(DP_v   ,kind=RP)
          w(k)  = real(DP_w   ,kind=RP)
          t(k)  = real(DP_t   ,kind=RP)
          phis  = real(DP_phis,kind=RP)
          ps    = real(DP_ps  ,kind=RP)
          rho   = real(DP_rho ,kind=RP)
          q     = real(DP_q   ,kind=RP)
       enddo

       call conv_vxvyvz( kdim, lat, lon, u, v, vx, vy, vz )

       do k = 1, kdim
          DIAG_var(n,k,l,1) = p (k)
          DIAG_var(n,k,l,2) = t (k)
          DIAG_var(n,k,l,3) = vx(k)
          DIAG_var(n,k,l,4) = vy(k)
          DIAG_var(n,k,l,5) = vz(k)
          DIAG_var(n,k,l,6) = w (k)
       enddo
    enddo
    enddo

    return
  end subroutine gravwave_init

  !-----------------------------------------------------------------------------
  subroutine tomita_init( &
       ijdim,      &
       kdim,       &
       lall,       &
       DIAG_var    )
    use scale_vector, only: &
       VECTR_xyz2latlon
    use mod_grd, only: &
       GRD_x,          &
       GRD_XDIR,       &
       GRD_YDIR,       &
       GRD_ZDIR,       &
       GRD_vz,         &
       GRD_Z,          &
       GRD_ZH
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: lall
    real(RP), intent(out) :: DIAG_var(ijdim,kdim,lall,6+TRC_VMAX)

    ! work paramters
    real(RP) :: lat, lon                 ! latitude, longitude on Icosahedral grid
    real(RP) :: prs(kdim),   tmp(kdim)   ! pressure & temperature in ICO-grid field
    real(RP) :: wix(kdim),   wiy(kdim)   ! zonal/meridional wind components in ICO-grid field

    real(RP) :: z_local (kdim)
    real(RP) :: vx_local(kdim)
    real(RP) :: vy_local(kdim)
    real(RP) :: vz_local(kdim)

    integer :: n, l, k, K0
    logical :: logout
    !---------------------------------------------------------------------------

    K0 = ADM_KNONE
    logout = .true.

    DIAG_var(:,:,:,:) = 0.0_RP

    if( IO_L ) write(IO_FID_LOG,*) "Qian98 Like Mountain Wave Exp. (Tomita and Satoh 2004)"

    do l = 1, lall
    do n = 1, ijdim
       z_local(1) = GRD_vz(n,2,l,GRD_ZH)
       do k = 2, kdim
          z_local(k) = GRD_vz(n,k,l,GRD_Z)
       enddo

       call VECTR_xyz2latlon( GRD_x(n,K0,l,GRD_XDIR), & ! [IN]
                              GRD_x(n,K0,l,GRD_YDIR), & ! [IN]
                              GRD_x(n,K0,l,GRD_ZDIR), & ! [IN]
                              lat,                    & ! [OUT]
                              lon                     ) ! [OUT]

       call tomita_2004( kdim, lat, z_local, wix, wiy, tmp, prs, logout )
       logout = .false.
       call conv_vxvyvz ( kdim, lat, lon, wix, wiy, vx_local, vy_local, vz_local )

       do k = 1, kdim
          DIAG_var(n,k,l,1) = prs(k)
          DIAG_var(n,k,l,2) = tmp(k)
          DIAG_var(n,k,l,3) = vx_local(k)
          DIAG_var(n,k,l,4) = vy_local(k)
          DIAG_var(n,k,l,5) = vz_local(k)
       enddo

    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) " |            Vertical Coordinate used in JBW initialization              |"
    if( IO_L ) write(IO_FID_LOG,*) " |------------------------------------------------------------------------|"
    do k = 1, kdim
       if( IO_L ) write(IO_FID_LOG,'(3X,"(k=",I3,") HGT:",F8.2," [m]",2X,"PRS: ",F9.2," [Pa]")') &
       k, z_local(k), prs(k)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) " |------------------------------------------------------------------------|"

    return
  end subroutine tomita_init

  !-----------------------------------------------------------------------------
  ! estimation ps distribution by using topography
  subroutine tomita_2004( &
      kdim,     &  !--- IN : # of z dimension
      lat,      &  !--- IN : latitude
      z_local,  &  !--- IN : z vertical coordinate
      wix,      &  !--- INOUT : zonal wind field
      wiy,      &  !--- INOUT : meridional wind field
      tmp,      &  !--- INOUT : temperature
      prs,      &  !--- INOUT : pressure
      logout    )  !--- IN : log output switch
    implicit none

    integer, intent(in) :: kdim
    real(RP), intent(in) :: lat
    real(RP), intent(in) :: z_local(kdim)
    real(RP), intent(inout) :: wix(kdim)
    real(RP), intent(inout) :: wiy(kdim)
    real(RP), intent(inout) :: tmp(kdim)
    real(RP), intent(inout) :: prs(kdim)
    logical, intent(in) :: logout

    integer :: k
    real(RP) :: g1, g2, Gphi, Gzero, Pphi
    real(RP), parameter :: N = 0.0187_RP        ! Brunt-Vaisala Freq.
    real(RP), parameter :: prs0 = 1.E+5_RP      ! pressure at the equator [Pa]
    real(RP), parameter :: ux0 = 40.0_RP        ! zonal wind at the equator [ms-1]
    real(RP) :: N2                              ! Square of Brunt-Vaisala Freq.
    real(RP) :: work
    !-----

    if (logout) then
       if( IO_L ) write(IO_FID_LOG, '("| == Tomita 2004 Mountain Wave Exp.  ( Z levels:", I4, ")")') kdim
       if( IO_L ) write(IO_FID_LOG, '("| -- Brunt-Vaisala Freq.:", F20.10)') N
       if( IO_L ) write(IO_FID_LOG, '("| -- Earth Angular Velocity:", F20.10)') omega
       if( IO_L ) write(IO_FID_LOG, '("| -- Earth Radius:", F20.10)') a
       if( IO_L ) write(IO_FID_LOG, '("| -- Earth Gravity Accel.:", F20.10)') g
    endif

    N2 = N**2.0_RP
    work = (N2*a) / (4.0_RP*g*Kap)

    g1 =   2.0_RP * ( 3.0_RP + 4.0_RP*cos(2.0_RP*zero) + cos(4.0_RP*zero) ) * (ux0**4.0_RP)   &
        +  8.0_RP * ( 3.0_RP + 4.0_RP*cos(2.0_RP*zero) + cos(4.0_RP*zero) ) * (ux0**3.0_RP)*a*omega   &
        +  8.0_RP * ( 3.0_RP + 4.0_RP*cos(2.0_RP*zero) + cos(4.0_RP*zero) ) * (ux0**2.0_RP)*(a**2.0_RP)*(omega**2.0_RP)   &
        - 16.0_RP * ( 1.0_RP + cos(2.0_RP*zero) ) * (ux0**2.0_RP)*a*g   &
        - 32.0_RP * ( 1.0_RP + cos(2.0_RP*zero) ) * ux0*(a**2.0_RP)*g*omega   &
        + 16.0_RP * (a**2.0_RP) * (g**2.0_RP)
    g2 = (ux0**4.0_RP) + 4.0_RP*a*omega*(ux0**3.0_RP) + 4.0_RP*(a**2.0_RP)*(omega**2.0_RP)*(ux0**2.0_RP)
    Gzero = ( g1 / g2 )**work

    g1 =   2.0_RP * ( 3.0_RP + 4.0_RP*cos(2.0_RP*lat) + cos(4.0_RP*lat) ) * (ux0**4.0_RP)   &
        +  8.0_RP * ( 3.0_RP + 4.0_RP*cos(2.0_RP*lat) + cos(4.0_RP*lat) ) * (ux0**3.0_RP)*a*omega   &
        +  8.0_RP * ( 3.0_RP + 4.0_RP*cos(2.0_RP*lat) + cos(4.0_RP*lat) ) * (ux0**2.0_RP)*(a**2.0_RP)*(omega**2.0_RP)   &
        - 16.0_RP * ( 1.0_RP + cos(2.0_RP*lat) ) * (ux0**2.0_RP)*a*g   &
        - 32.0_RP * ( 1.0_RP + cos(2.0_RP*lat) ) * ux0*(a**2.0_RP)*g*omega   &
        + 16.0_RP * (a**2.0_RP) * (g**2.0_RP)
    g2 = (ux0**4.0_RP) + 4.0_RP*a*omega*(ux0**3.0_RP) + 4.0_RP*(a**2.0_RP)*(omega**2.0_RP)*(ux0**2.0_RP)
    Gphi = ( g1 / g2 )**work

    Pphi = prs0 * ( Gzero / Gphi )

    do k = 1, kdim
       wix(k) = ux0 * cos(lat)
       prs(k) = Pphi * exp( (-1.0_RP*N2*z_local(k)) / (g*Kap) )
       tmp(k) = (g * Kap * ( g - (wix(k)**2.0_RP)/a - 2.0_RP*omega*wix(k)*cos(lat) )) / (N2*Rd)
    enddo

    return
  end subroutine tomita_2004

  !-----------------------------------------------------------------------------
  ! eta vertical coordinate by Newton Method
  subroutine eta_vert_coord_NW( &
      kdim,      &
      itr,       &
      z,         &
      tmp,       &
      geo,       &
      eta_limit, &
      eta,       &
      signal     )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    integer,  intent(in)    :: kdim        ! # of z dimension
    integer,  intent(in)    :: itr         ! iteration number
    real(RP), intent(in)    :: z  (kdim)   ! z-height vertical coordinate
    real(RP), intent(in)    :: tmp(kdim)   ! guessed temperature
    real(RP), intent(in)    :: geo(kdim)   ! guessed geopotential
    logical,  intent(in)    :: eta_limit   ! eta limitation flag
    real(RP), intent(inout) :: eta(kdim,2) ! eta level vertical coordinate
    logical,  intent(inout) :: signal      ! iteration signal

    real(RP) :: diff(kdim)
    real(RP) :: F(kdim), Feta(kdim)
    real(RP) :: criteria
    integer  :: k
    !---------------------------------------------------------------------------

    criteria = max( EPS * 10.0_RP, 1.E-14_RP )

    do k = 1, kdim
       F   (k) = -g*z(k) + geo(k)
       Feta(k) = -1.0_RP * ( Rd/eta(k,1) ) * tmp(k)

       eta(k,2) = eta(k,1) - ( F(k)/Feta(k) )

       if (eta_limit) then ! [add] for PSDM (2013/12/20 R.Yoshida)
          eta(k,2) = min(eta(k,2),1.0_RP) ! not allow over 1.0 for eta
       endif

       eta(k,2) = max(eta(k,2),EPS) ! not allow over 1.0 for eta

       diff(k) = abs( eta(k,2) - eta(k,1) )
    enddo

    eta(:,1) = eta(:,2)
    if (message) then
       if( IO_L ) write(IO_FID_LOG,'(A,I4,A,ES20.10,A,ES20.10)') &
                  " | Eta  ",itr,": -- MAX: ",maxval(diff(:))," MIN: ",minval(diff(:))
       if( IO_L ) write(IO_FID_LOG,'(A,I4,A,ES20.10,A,ES20.10)') &
                  " | Diff ",itr,": -- MAX: ",maxval(diff(:))," MIN: ",minval(diff(:))
    endif

    if ( maxval(diff(:)) < criteria ) then
       signal = .false.
    else
       if (message) then
          if( IO_L ) write(IO_FID_LOG,*) "| Iterating : ", itr, "criteria = ", criteria
       endif
    endif

    return
  end subroutine eta_vert_coord_NW

  !-----------------------------------------------------------------------------
  ! calculation of steady state
  subroutine steady_state( &
      kdim, &  !--- IN : # of z dimension
      lat,  &  !--- IN : latitude information
      eta,  &  !--- IN : eta level vertical coordinate
      wix,  &  !--- INOUT : zonal wind component
      wiy,  &  !--- INOUT : meridional wind component
      tmp,  &  !--- INOUT : mean temperature
      geo   )  !--- INOUT : mean geopotential height
    !
    implicit none
    integer :: k
    integer, intent(in) :: kdim
    real(RP), intent(in) :: lat
    real(RP), intent(in) :: eta(kdim,2)
    real(RP), intent(inout) :: wix(kdim)
    real(RP), intent(inout) :: wiy(kdim)
    real(RP), intent(inout) :: tmp(kdim)
    real(RP), intent(inout) :: geo(kdim)
    real(RP) :: eta_v
    real(RP) :: work1, work2
    !
    ! ---------- horizontal mean
    work1 = pi/2.0_RP
    work2 = Rd*ganma/g
    do k = 1, kdim
       eta_v = (eta(k,1) - eta0)*(work1)
       wix(k) = u0 * (cos(eta_v))**1.5_RP * (sin(2.0_RP*lat))**2.0_RP
       !
       !if( etaS >= eta(k,1) .and. eta(k,1) >= etaT ) then  ! not allow over 1.0 for eta
       if( eta(k,1) >= etaT ) then
          tmp(k) = t0 * eta(k,1)**work2
          geo(k) = t0*g/ganma * ( 1.0_RP - eta(k,1)**work2 )
       elseif( eta(k,1) < etaT ) then
          tmp(k) = t0 * eta(k,1)**work2 + delT*(etaT - eta(k,1))**5.0_RP

          geo(k) = t0*g/ganma * ( 1.0_RP - eta(k,1)**work2 ) - Rd * delT *                                          &
                 ( ( log(eta(k,1)/etaT) + 137.0_RP/60.0_RP )*etaT**5.0_RP - 5.0_RP*(etaT**4.0_RP)*eta(k,1)          &
                   +  5.0_RP*(etaT**3.0_RP)*(eta(k,1)**2.0_RP) - (10.0_RP/3.0_RP)*(etaT**2.0_RP)*(eta(k,1)**3.0_RP) &
                   + (5.0_RP/4.0_RP)*etaT*(eta(k,1)**4.0_RP) - (1.0_RP/5.0_RP)*(eta(k,1)**5.0_RP)                   )
       else
          if( IO_L ) write(IO_FID_LOG,'(A)') "|-- ETA BOUNDARY ERROR: [steady state calc.]"
          if( IO_L ) write(IO_FID_LOG,'("|-- (",I3,")  eta: ",F10.4)') k, eta(k,1)
          stop
       endif
       !else
       !   if( IO_L ) write(IO_FID_LOG,'(A)') "|-- OVER 1.0 for eta: [steady state calc.]"
       !   stop
       !endif
       !
    enddo

    ! ---------- meridional distribution for temeperature and geopotential
    work1 = pi/2.0_RP
    work2 = 3.0_RP/4.0_RP * ( pi*u0 / Rd )
    do k = 1, kdim
       eta_v = (eta(k,1) - eta0)*(work1)
       tmp(k) = tmp(k)                                             &
              + work2*eta(k,1) * sin(eta_v) * (cos(eta_v))**0.5_RP &
              * ( ( -2.0_RP * (sin(lat))**6.0_RP * (cos(lat)**2.0_RP + 1.0_RP/3.0_RP) + 10.0_RP/63.0_RP )   &
                * 2.0_RP*u0*(cos(eta_v))**1.5_RP                                                            &
                + ( 8.0_RP/5.0_RP * (cos(lat))**3.0_RP * ((sin(lat))**2.0_RP + 2.0_RP/3.0_RP) - pi/4.0_RP ) &
                * a*omega                                                                                   )
       geo(k) = geo(k)                   &
              + u0*(cos(eta_v))**1.5_RP  &
              * ( ( -2.0_RP * (sin(lat))**6.0_RP * (cos(lat)**2.0_RP + 1.0_RP/3.0_RP) + 10.0_RP/63.0_RP )   &
                * u0*(cos(eta_v))**1.5_RP                                                                   &
                + ( 8.0_RP/5.0_RP * (cos(lat))**3.0_RP * ((sin(lat))**2.0_RP + 2.0_RP/3.0_RP) - pi/4.0_RP ) &
                * a*omega                                                                                   )
    enddo

    wiy(:) = 0.0_RP

    return
  end subroutine steady_state

  !-----------------------------------------------------------------------------
  ! convert geopotential height to pressure
  subroutine geo2prs( &
      kdim,         &  !--- IN : # of z dimension
      ps,           &  !--- IN : surface pressure
      lat,          &  !--- IN : latitude
      tmp,          &  !--- IN : temperature
      geo,          &  !--- IN : geopotential height at full height
      wix,          &  !--- IN : zonal wind
      prs,          &  !--- IN : pressure
      eps_geo2prs,  &  !--- IN : eps
      nicamcore,    &  !--- IN : nicamcore switch
      logout        )  !--- IN : switch of log output
    !
    implicit none
    integer, intent(in) :: kdim
    real(RP), intent(in) :: ps
    real(RP), intent(in) :: lat
    real(RP), intent(in) :: tmp(kdim)
    real(RP), intent(in) :: geo(kdim)
    real(RP), intent(in) :: wix(kdim)
    real(RP), intent(inout) :: prs(kdim)
    real(RP), intent(in) :: eps_geo2prs
    logical, intent(in) :: nicamcore
    logical, intent(in) :: logout

    integer :: i, k
    integer, parameter :: limit = 400
    real(RP) :: dz, uave, diff
    real(RP) :: f_cf(3)
    real(RP) :: pp(kdim)
    logical :: iteration = .false.
    logical :: do_iter = .true.
    !-----

    pp(1) = ps

    ! first guess (upward: trapezoid)
    do k=2, kdim
       dz = (geo(k) - geo(k-1))/g
       if (nicamcore) then
          uave = (wix(k) + wix(k-1)) * 0.5_RP
          f_cf(1) = 2.0_RP*omega*uave*cos(lat) + (uave**2.0_RP)/a
       else
          f_cf(1) = 0.0_RP
       endif

       pp(k) = pp(k-1) * ( 1.0_RP + dz*(f_cf(1) - g ) / ( 2.0_RP*Rd*tmp(k-1) ) ) &
                       / ( 1.0_RP - dz*(f_cf(1) - g ) / ( 2.0_RP*Rd*tmp(k  ) ) )
    enddo
    prs(:) = pp(:)

    ! iteration (simpson)
    if (iteration) then
    do i=1, limit
       prs(1) = ps
       do k=3, kdim        ! upward
          pp(k) = simpson( prs(k), prs(k-1), prs(k-2), tmp(k), tmp(k-1), tmp(k-2), &
                  wix(k), wix(k-1), wix(k-2), geo(k), geo(k-2), lat, .false., nicamcore )
       enddo
       prs(:) = pp(:)
       do k=kdim-2, 1, -1  ! downward
          pp(k) = simpson( prs(k+2), prs(k+1), prs(k), tmp(k+2), tmp(k+1), tmp(k), &
                  wix(k+2), wix(k+1), wix(k), geo(k+2), geo(k), lat, .true., nicamcore )
       enddo
       prs(:) = pp(:)
       diff = pp(1) - ps

       if ( abs(diff) < eps_geo2prs ) then
          do_iter = .false.
          exit
       endif
    enddo
    else
       !if( IO_L ) write(IO_FID_LOG,*) 'ETA ITERATION SKIPPED'
       do_iter = .false.
    endif

    if (do_iter) then
       if( IO_L ) write(IO_FID_LOG,*) 'ETA ITERATION ERROR: NOT CONVERGED at GEO2PRS', diff
       stop
    endif

    ! fininalize
    prs(1) = ps
    pp(1) = ps
    do k=3, kdim        ! upward
       pp(k) = simpson( prs(k), prs(k-1), prs(k-2), tmp(k), tmp(k-1), tmp(k-2), &
               wix(k), wix(k-1), wix(k-2), geo(k), geo(k-2), lat, .false., nicamcore )
    enddo
    prs(:) = pp(:)

    if (logout) then
       if (iteration) then
          if( IO_L ) write(IO_FID_LOG, *) " | diff (guess - ps) : ", diff, "[Pa]  --  itr times: ", (i-1)
       else
          if( IO_L ) write(IO_FID_LOG, *) " | no iteration in geo2prs"
       endif
    endif
    if (message) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A)') " | ----- Pressure (Final Guess) -----"
       do k = 1, kdim
          if( IO_L ) write(IO_FID_LOG, '(" | K(",I3,") -- ",F20.13)') k, prs(k)
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    endif

    return
  end subroutine geo2prs

  !-----------------------------------------------------------------------------
  ! estimation ps distribution by using topography
  subroutine ps_estimation( &
      kdim,   &  !--- IN : # of z dimension
      lat,  &  !--- IN : latitude information
      eta,  &  !--- IN : eta coordinate
      tmp,  &  !--- IN : temperature
      geo,  &  !--- IN : geopotential height at full height
      wix,  &  !--- IN : zonal wind speed
      ps,   &  !--- OUT : surface pressure
      nicamcore )  !--- IN : nicamcore switch
    implicit none

    integer, intent(in) :: kdim
    real(RP), intent(in) :: lat
    real(RP), intent(in) :: eta(kdim)
    real(RP), intent(in) :: tmp(kdim)
    real(RP), intent(in) :: geo(kdim)
    real(RP), intent(in) :: wix(kdim)
    real(RP), intent(out) :: ps
    logical, intent(in) :: nicamcore

    real(RP), parameter :: lat0 = 0.691590985442682_RP
    real(RP) :: cs32ev, f1, f2
    real(RP) :: eta_v, tmp0, tmp1
    real(RP) :: ux1, ux2, hgt0, hgt1
    real(RP) :: dz, uave
    real(RP) :: f_cf(3)
    real(RP), parameter :: eta1 = 1.0_RP
    !-----

    eta_v = (eta1 - eta0)*(pi*0.5_RP)

    ! temperature at bottom of eta-grid
    tmp0 = t0                                                   &
           + (3.0_RP/4.0_RP * (pi*u0/Rd))*eta1 * sin(eta_v) * (cos(eta_v))**0.5_RP                &
           * ( ( -2.0_RP * (sin(lat0))**6.0_RP * (cos(lat0)**2.0_RP + 1.0_RP/3.0_RP) + 10.0_RP/63.0_RP )   &
                * 2.0_RP*u0*(cos(eta_v))**1.5_RP                   &
                + ( 8.0_RP/5.0_RP * (cos(lat0))**3.0_RP * ((sin(lat0))**2.0_RP + 2.0_RP/3.0_RP) - pi/4.0_RP ) &
                * a*omega  )
    tmp1 = tmp(1)

    ! wind speed at bottom of eta-grid
    ux1 = (u0 * cos(eta_v)**1.5_RP) * (sin(2.0_RP*lat0))**2.0_RP
    ux2 = wix(1)

    ! topography calculation (imported from mod_grd.f90)
    cs32ev = ( cos( (1.0_RP-0.252_RP) * pi * 0.5_RP ) )**1.5_RP
    f1 = 10.0_RP/63.0_RP - 2.0_RP * sin(lat)**6 * ( cos(lat)**2 + 1.0_RP/3.0_RP )
    f2 = 1.6_RP * cos(lat)**3 * ( sin(lat)**2 + 2.0_RP/3.0_RP ) - 0.25_RP * pi
    hgt1 = -1.0_RP * u0 * cs32ev * ( f1*u0*cs32ev + f2*a*omega ) / g
    hgt0 = 0.0_RP

    ! ps estimation
    dz = hgt1 - hgt0
    if (nicamcore) then
       uave = (ux1 + ux2) * 0.5_RP
       f_cf(1) = 2.0_RP*omega*uave*cos(lat) + (uave**2.0_RP)/a
    else
       f_cf(1) = 0.0_RP
    endif
    ps = p0 * ( 1.0_RP + dz*(f_cf(1) - g)/(2.0_RP*Rd*tmp0) ) &
            / ( 1.0_RP - dz*(f_cf(1) - g)/(2.0_RP*Rd*tmp1) )

    return
  end subroutine ps_estimation

  !-----------------------------------------------------------------------------
  ! setting perturbation
  subroutine perturbation( &
       ijdim,   &
       kdim,    &
       lall,    &
       vmax,    &
       DIAG_var )
    use mod_adm, only: &
       K0 => ADM_KNONE
    use scale_vector, only: &
       VECTR_xyz2latlon
    use mod_grd, only: &
       GRD_x,    &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR
    implicit none

    integer, intent(in) :: ijdim, kdim, lall, vmax
    real(RP), intent(inout) :: DIAG_var(ijdim,kdim,lall,vmax)

    integer, parameter :: ID_vx  = 3
    integer, parameter :: ID_vy  = 4
    integer, parameter :: ID_vz  = 5
    integer :: n, k, l
    real(RP) :: lat, lon
    real(RP) :: r, rr, rbyrr, cla, clo
    real(RP) :: ptb_wix(kdim), ptb_wiy(kdim)
    real(RP) :: ptb_vx(kdim), ptb_vy(kdim), ptb_vz(kdim)

    cla = clat * d2r
    clo = clon * d2r

    do l = 1, lall
    do n = 1, ijdim
       call VECTR_xyz2latlon( GRD_x(n,K0,l,GRD_XDIR), & ! [IN]
                              GRD_x(n,K0,l,GRD_YDIR), & ! [IN]
                              GRD_x(n,K0,l,GRD_ZDIR), & ! [IN]
                              lat,                    & ! [OUT]
                              lon                     ) ! [OUT]
       r = a * acos( sin(cla)*sin(lat) + cos(cla)*cos(lat)*cos(lon-clo) )
       rr = a / 10.0_RP
       rbyrr = r/rr
       do k = 1, kdim
          ptb_wix(k) = uP * exp( -1.0_RP*rbyrr**2.0_RP )
          ptb_wiy(k) = 0.0_RP
       enddo

       call conv_vxvyvz( kdim, lat, lon, ptb_wix, ptb_wiy, ptb_vx, ptb_vy, ptb_vz )
       do k = 1, kdim
          DIAG_var(n,k,l,ID_vx) = DIAG_var(n,k,l,ID_vx) + ptb_vx(k)
          DIAG_var(n,k,l,ID_vy) = DIAG_var(n,k,l,ID_vy) + ptb_vy(k)
          DIAG_var(n,k,l,ID_vz) = DIAG_var(n,k,l,ID_vz) + ptb_vz(k)
       enddo
    enddo
    enddo
    !
    return
  end subroutine perturbation

  !-----------------------------------------------------------------------------
  subroutine conv_vxvyvz( &
      kdim, &  !--- IN : # of z dimension
      lat,  &  !--- IN : latitude information
      lon,  &  !--- IN : longitude information
      wix,  &  !--- IN : zonal wind component on latlon
      wiy,  &  !--- IN : meridional wind component on latlon
      vx1d, &  !--- INOUT : horizontal-x component on absolute system for horizontal wind
      vy1d, &  !--- INOUT : horizontal-y component on absolute system for horizontal wind
      vz1d  )  !--- INOUT : vertical component on absolute system for horizontal wind
    !
    implicit none
    integer, intent(in) :: kdim
    real(RP), intent(in)    :: lat
    real(RP), intent(in)    :: lon
    real(RP), intent(in)    :: wix(kdim)
    real(RP), intent(in)    :: wiy(kdim)
    real(RP), intent(inout) :: vx1d(kdim)
    real(RP), intent(inout) :: vy1d(kdim)
    real(RP), intent(inout) :: vz1d(kdim)
    !
    integer :: k
    real(RP) :: unit_east(3), unit_north(3)
    !
    ! imported from NICAM/nhm/mkinit/prg_mkinit_ncep.f90 (original written by H.Miura)
    ! *** compute vx, vy, vz as 1-dimensional variables
    do k = 1, kdim
       unit_east  = Sp_Unit_East( lon )
       unit_north = Sp_Unit_North( lon, lat )
       !
       vx1d(k) = unit_east(1) * wix(k) + unit_north(1) * wiy(k)
       vy1d(k) = unit_east(2) * wix(k) + unit_north(2) * wiy(k)
       vz1d(k) = unit_east(3) * wix(k) + unit_north(3) * wiy(k)
    enddo
    !
    return
  end subroutine conv_vxvyvz

  !-----------------------------------------------------------------------------
  subroutine diag_pressure( &
      kdim,        &  !--- IN : # of z dimension
      z,           &  !--- IN : height information
      rho,         &  !--- IN : density
      t,           &  !--- IN : temperature
      q,           &  !--- IN : water vapor content
      ps,          &  !--- IN : surface pressure
      prs,         &  !--- OUT: pressure
      prs_rebuild, &  !--- IN : rebuild switch
      prs_dry      )  !--- IN : dry condition switch
    !
    implicit none
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: z(:)
    real(RP), intent(in)  :: rho(:)
    real(RP), intent(in)  :: t(:)
    real(RP), intent(in)  :: q(:)
    real(RP), intent(in)  :: ps
    real(RP), intent(out) :: prs(:)
    logical,  intent(in)  :: prs_rebuild
    logical,  intent(in)  :: prs_dry
    !
    real(RP) :: dz, R0, R1
    integer  :: k
    !-----
    if ( prs_dry ) then
       do k=1, kdim
          prs(k) = rho(k) * t(k) * Rd
       enddo

       if ( prs_rebuild ) then ! rebuild
          do k=2, kdim
             dz = z(k) - z(k-1)
             prs(k) = prs(k-1) * ( 1.0_RP - dz*g/(2.0_RP*Rd*t(k-1)) ) &
                               / ( 1.0_RP + dz*g/(2.0_RP*Rd*t(k  )) )
          enddo
       endif
    else
       do k=1, kdim
          prs(k) = rho(k) * t(k) * ( (1.0_RP-q(k))*Rd + q(k)*Rv )
       enddo

       if ( prs_rebuild ) then ! rebuild
          do k=2, kdim
             dz = z(k) - z(k-1)
             R0 = ( 1.0_RP - q(k-1) )*Rd + q(k-1)*Rv
             R1 = ( 1.0_RP - q(k  ) )*Rd + q(k  )*Rv
             prs(k) = prs(k-1) * ( 1.0_RP - dz*g/(2.0_RP*R0*t(k-1)) ) &
                               / ( 1.0_RP + dz*g/(2.0_RP*R1*t(k  )) )
          enddo
       endif
    endif

    return
  end subroutine diag_pressure

  !-----------------------------------------------------------------------------
  function simpson( &
      pin1,        &  !--- IN : pressure (top)
      pin2,        &  !--- IN : pressure (middle)
      pin3,        &  !--- IN : pressure (bottom)
      t1,          &  !--- IN : temperature (top)
      t2,          &  !--- IN : temperature (middle)
      t3,          &  !--- IN : temperature (bottom)
      u1,          &  !--- IN : zonal wind (top)
      u2,          &  !--- IN : zonal wind (middle)
      u3,          &  !--- IN : zonal wind (bottom)
      geo1,        &  !--- IN : geopotential (top)
      geo3,        &  !--- IN : geopotential (bottom)
      lat,         &  !--- IN : latitude
      downward,    &  !--- IN : downward switch
      nicamcore )  &  !--- IN : nicamcore switch
      result (pout)
    !
    implicit none
    real(RP), intent(in) :: pin1, pin2, pin3
    real(RP), intent(in) :: t1, t2, t3
    real(RP), intent(in) :: u1, u2, u3
    real(RP), intent(in) :: geo1, geo3, lat
    logical, intent(in) :: downward
    logical, intent(in) :: nicamcore
    !
    real(RP) :: dz, pout
    real(RP) :: f_cf(3), rho(3)
    !---------------------------------------------------------------------------

    dz = (geo1-geo3) / g * 0.5_RP

    if (nicamcore) then
       f_cf(1) = 2.0_RP*omega*u1*cos(lat) + (u1**2.0_RP)/a
       f_cf(2) = 2.0_RP*omega*u2*cos(lat) + (u2**2.0_RP)/a
       f_cf(3) = 2.0_RP*omega*u3*cos(lat) + (u3**2.0_RP)/a
    else
       f_cf(:) = 0.0_RP
    endif
    rho(1) = pin1 / ( Rd*t1 )
    rho(2) = pin2 / ( Rd*t2 )
    rho(3) = pin3 / ( Rd*t3 )

    if (downward) then
       pout = pin1 - ( (1.0_RP/3.0_RP) * rho(1) * ( f_cf(1) - g ) &
                     + (4.0_RP/3.0_RP) * rho(2) * ( f_cf(2) - g ) &
                     + (1.0_RP/3.0_RP) * rho(3) * ( f_cf(3) - g ) ) * dz
    else
       pout = pin3 + ( (1.0_RP/3.0_RP) * rho(1) * ( f_cf(1) - g ) &
                     + (4.0_RP/3.0_RP) * rho(2) * ( f_cf(2) - g ) &
                     + (1.0_RP/3.0_RP) * rho(3) * ( f_cf(3) - g ) ) * dz
    endif

    return
  end function simpson

  !-----------------------------------------------------------------------------
  function Sp_Unit_East( lon ) result( unit_east )
    implicit none

    real(RP), intent(in) :: lon ! [rad]
    real(RP)             :: unit_east(3)
    !---------------------------------------------------------------------------

    unit_east(1) = -sin(lon) ! x-direction
    unit_east(2) =  cos(lon) ! y-direction
    unit_east(3) = 0.0_RP      ! z-direction

    return
  end function Sp_Unit_East

  !-----------------------------------------------------------------------------
  function Sp_Unit_North( lon, lat ) result( unit_north )
    implicit none

    real(RP), intent(in) :: lon, lat ! [rad]
    real(RP)             :: unit_north(3)
    !---------------------------------------------------------------------------

    unit_north(1) = -sin(lat) * cos(lon) ! x-direction
    unit_north(2) = -sin(lat) * sin(lon) ! y-direction
    unit_north(3) =  cos(lat)            ! z-direction

    return
  end function Sp_Unit_North

end module mod_ideal_init
