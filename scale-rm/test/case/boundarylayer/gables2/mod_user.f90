!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-08-28 (S.Nishizawa)   [new]
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
  use scale_index
  use scale_land_grid_index
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
  integer,  parameter :: K_ini = 10
  real(RP), parameter :: z_ini(K_ini) = &
       (/   0.0_RP,  200.0_RP,  850.0_RP,  900.0_RP,  1000.0_RP, &
         2000.0_RP, 3500.0_RP, 4000.0_RP, 5000.0_RP, 40000.0_RP /)
  real(RP), parameter :: pt_ini(K_ini) = &
       (/ 288.0_RP, 286.0_RP, 286.0_RP, 288.0_RP, 292.0_RP, &
          300.0_RP, 310.0_RP, 312.0_RP, 316.0_RP, 800.0_RP /)
  real(RP), parameter :: qv_ini(K_ini) = &
       (/ 0.0025_RP, 0.0025_RP, 0.0025_RP, 0.0025_RP, 0.005_RP, &
          0.0030_RP, 0.0020_RP, 0.0015_RP, 0.0000_RP, 0.000_RP /)

  real(RP), parameter :: BETA = 0.025_RP

  real(RP) :: Ps = 0.972E5_RP
  real(RP) :: Ug =  3.0_RP
  real(RP) :: Vg = -9.0_RP
  real(RP) :: lat = 37.6_RP


  integer, parameter :: QA = 3
  character(len=H_SHORT) :: QNAME(QA)
  character(len=H_MID)   :: QDESC(QA)
  character(len=H_SHORT) :: QUNIT(QA)

  data QNAME / &
                 'QV', &
                 'QC', &
                 'QI'  /

  data QDESC / &
                 'Ratio of Water Vapor mass to total mass (Specific humidity)',   &
                 'Ratio of Cloud Water mass to total mass', &
                 'Ratio of Cloud Ice mass to total mass'   /

  data QUNIT / &
                 'kg/kg', &
                 'kg/kg', &
                 'kg/kg'  /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config
    use scale_atmos_hydrometeor, only: &
         ATMOS_HYDROMETEOR_regist, &
         I_QV, &
         I_QC, &
         I_QI
    implicit none
    integer :: QS

    call ATMOS_HYDROMETEOR_regist( QS, & ! (out)
                                  1, 1, 1,  & ! (in)
                                  QNAME, QDESC, QUNIT ) ! (in)

    I_QV = QS
    I_QC = QS + 1
    I_QI = QS + 2

    return
  end subroutine USER_config
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_USER / &
       Ps, &
       Ug, &
       Vg, &
       lat

    integer :: ierr
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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    use scale_atmos_hydrostatic, only: &
       buildrho => ATMOS_HYDROSTATIC_buildrho
    use scale_atmos_hydrometeor, only: &
       I_QV
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: TEMP(KA)
    real(RP) :: PRES(KA)
    real(RP) :: PT(KA)
    real(RP) :: QV(KA)
    real(RP) :: QC(KA)
    real(RP) :: Ts

    integer :: i, j
    !---------------------------------------------------------------------------

    call interporate(PT, pt_ini)
    call interporate(QV, qv_ini)

    QC = 0.0_RP

    call buildrho( &
         RHO, & ! (out)
         TEMP, & ! (out)
         PRES, & ! (out)
         PT,   & ! (in)
         QV,   & ! (in)
         QC,   & ! (in)
         Ts, & ! (out)
         Ps, & ! (in)
         pt_ini(1), & ! (in)
         qv_ini(1), & ! (in)
         0.0_RP ) ! (in)

    MOMZ = 0.0_RP
    QTRC = 0.0_RP

    do j = 1, JA
    do i = 1, IA
       DENS(:,i,j) = RHO
       MOMX(:,i,j) = RHO * Ug
       MOMY(:,i,j) = RHO * Vg
       RHOT(:,i,j) = RHO * PT
       QTRC(:,i,j,I_QV) = QV
    end do
    end do

    call set_tg( 0.0_RP )

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
         MOMX,    &
         MOMY,    &
         MOMX_tp, &
         MOMY_tp
    use scale_const, only: &
         OHM   => CONST_OHM, &
         PI    => CONST_PI
    use scale_time, only: &
         NOWSEC => TIME_NOWSEC, &
         dt     => TIME_DTSEC
    use mod_land_vars, only: &
         I_WaterCritical, &
         LAND_WATER, &
         LAND_PROPERTY
    implicit none
    real(RP) :: f, fdt
    real(RP) :: dudt, dvdt
    real(RP) :: rden
    integer :: k, i, j
    !---------------------------------------------------------------------------

    f = 2.0_RP * OHM * sin(lat*PI/180.0_RP)
    fdt = f * dt

    rden = 1.0_RP / (1.0_RP + fdt**2)

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dudt =   f * ( MOMY(k,i,j)/DENS(k,i,j) - Vg )
       dvdt = - f * ( MOMX(k,i,j)/DENS(k,i,j) - Ug )

       MOMX_tp(k,i,j) = MOMX_tp(k,i,j) &
                      +(dudt + fdt*dvdt) * rden * DENS(k,i,j)
       MOMY_tp(k,i,j) = MOMY_tp(k,i,j) &
                      +(dvdt + fdt*dudt) * rden * DENS(k,i,j)
    end do
    end do
    end do

    call set_tg( NOWSEC )

    LAND_WATER(LKS,:,:) = BETA * LAND_PROPERTY(:,:,I_WaterCritical)

    return
  end subroutine USER_step

  subroutine interporate( d_out, d_in )
    use scale_grid, only: &
         CZ   => GRID_CZ
    implicit none
    real(RP), intent(out) :: d_out(KA)
    real(RP), intent(in ) :: d_in(K_ini)
    integer :: k, kk

    do k = KS, KE
       if ( CZ(k) >= z_ini(K_ini) ) then
          d_out(k) = d_in(K_ini)
       else if ( CZ(k) < z_ini(1) ) then
          d_out(k) = d_in(1)
       else
          do kk = 1, K_ini-1
             if ( z_ini(kk) <= CZ(k) .and. CZ(k) < z_ini(kk+1) ) then
                d_out(k) = d_in(kk) &
                     + ( d_in(kk+1)-d_in(kk) ) * ( CZ(k)-z_ini(kk) ) / ( z_ini(kk+1)-z_ini(kk) )
                exit
             end if
          end do
       end if
    end do

    return
  end subroutine interporate

  subroutine set_tg( tsec )
!    use mod_atmos_phy_sf_vars, only: &
!       T_sf => ATMOS_PHY_SF_SFC_TEMP
    use mod_land_vars, only: &
       T_sf => LAND_SFC_TEMP
    implicit none
    real(RP), intent(in)  :: tsec
    real(RP) :: time
    real(RP) :: Tg

    time = tsec / 3600.0_RP + 16.0_RP

    if (time <= 17.4_RP) then
       Tg = -10.0_RP - 25.0_RP * cos(time*0.22_RP + 0.2_RP)
    else if (time <= 30.0_RP) then
       Tg = -0.54_RP * time + 15.2_RP
    else if (time <= 41.9_RP) then
       Tg = -7.0_RP - 25.0_RP * cos(time*0.21_RP + 1.8_RP)
    else if (time <= 53.3_RP) then
       Tg = -0.37_RP * time + 18.0_RP
    else if (time <= 65.6_RP) then
       Tg = -4.0_RP - 25.0_RP * cos(time * 0.22_RP + 2.5_RP)
    else
       Tg = 4.4_RP
    end if

    T_sf = Tg + 273.0_RP

    return
  end subroutine set_tg

end module mod_user
