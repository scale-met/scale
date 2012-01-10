!-------------------------------------------------------------------------------
!> Program Warm BUbble Test for SCALE-LES ver.3
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-LES ver.2
!!
!<
!-------------------------------------------------------------------------------
program warmbubble
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_setup
  use mod_process, only: &
     PRC_setup,    &
     PRC_MPIstart, &
     PRC_MPIstop
  use mod_const, only: &
     CONST_setup
  use mod_time, only: &
     TIME_setup,    &
     TIME_rapstart, &
     TIME_rapend,   &
     TIME_rapreport
  use mod_grid, only: &
     GRID_setup
  use mod_comm, only: &
     COMM_setup
  use mod_fileio, only: &
     FIO_setup, &
     FIO_finalize
  use mod_atmos_vars, only: &
     ATMOS_vars_setup, &
     ATMOS_vars_restart_write
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  !########## Initial setup ##########

  ! setup standard I/O
  call IO_setup

  ! start MPI
  call PRC_MPIstart

  ! setup process
  call PRC_setup

  ! setup constants
  call CONST_setup

  ! setup time
  call TIME_setup

  ! setup file I/O
  call FIO_setup

  ! setup horisontal/veritical grid system
  call GRID_setup

  ! setup mpi communication
  call COMM_setup

  ! setup atmosphere
  call ATMOS_vars_setup


  !########## main ##########

  call TIME_rapstart('Main')

  ! make initial state (restart)
  call MKEXP_warmbubble

  ! output restart
  call ATMOS_vars_restart_write

  call TIME_rapend('Main')


  !########## Finalize ##########
  call TIME_rapreport

  call FIO_finalize
  ! stop MPI
  call PRC_MPIstop

  stop
  !=============================================================================
contains

  !-----------------------------------------------------------------------------
  !> Make initial state for cold bubble experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_warmbubble
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       CPdry  => CONST_CPdry,  &
       CPovR  => CONST_CPovR,  &
       RovCP  => CONST_RovCP,  &
       CVovCP => CONST_CVovCP, &
       EPSvap => CONST_EPSvap, &
       Pstd   => CONST_Pstd
    use mod_grid, only : &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KA => GRID_KA, &
       IS => GRID_IS, &
       IE => GRID_IE, &
       JS => GRID_JS, &
       JE => GRID_JE, &
       KS => GRID_KS, &
       KE => GRID_KE, &
       GRID_CX, &
       GRID_CY, &
       GRID_CZ
    use mod_atmos_vars, only: &
       QA => A_QA,     &
       I_QV,           &
       ATMOS_vars_get, &
       ATMOS_vars_put
    implicit none

    real(8) :: ENV_THETA = 300.D0 ! Potential Temperature of environment [K]
    real(8) :: ENV_RH    = 80.D0  ! Relative Humidity of environment [%]
    real(8) :: XC_BBL = 18.D3     ! center location [m]: x
    real(8) :: YC_BBL = 18.D3     ! center location [m]: y
    real(8) :: ZC_BBL = 3.D3      ! center location [m]: z
    real(8) :: XR_BBL = 4.D3      ! bubble radius   [m]: x
    real(8) :: YR_BBL = 4.D3      ! bubble radius   [m]: y
    real(8) :: ZR_BBL = 2.D3      ! bubble radius   [m]: z

    NAMELIST / PARAM_MKEXP_warmbubble / &
       ENV_THETA, &
       ENV_RH,    &
       XC_BBL,    &
       YC_BBL,    &
       ZC_BBL,    &
       XR_BBL,    &
       YR_BBL,    &
       ZR_BBL

    real(8) :: dens(IA,JA,KA)      ! density     [kg/m3]
    real(8) :: momx(IA,JA,KA)      ! momentum(x) [kg/m3 * m/s]
    real(8) :: momy(IA,JA,KA)      ! momentum(y) [kg/m3 * m/s]
    real(8) :: momz(IA,JA,KA)      ! momentum(z) [kg/m3 * m/s]
    real(8) :: rhot(IA,JA,KA)      ! rho * theta [kg/m3 * K]
    real(8) :: qtrc(IA,JA,KA,QA)   ! tracer mixing ratio [kg/kg],[1/m3]

    real(8) :: pres(IA,JA,KA)    ! pressure [Pa]
    real(8) :: temp(IA,JA,KA)    ! temperature [K]
    real(8) :: pott(IA,JA,KA)    ! potential temperature [K]

    real(8) :: dist
    real(8) :: rh(IA,JA,KA)
    real(8) :: psat, qsat

    integer :: i, j, k
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[WARMBUBBLE]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_warmbubble,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_warmbubble. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_warmbubble)

    call ATMOS_vars_get( dens, momx, momy, momz, rhot, qtrc )

    momx(:,:,:)   = 0.D0
    momy(:,:,:)   = 0.D0
    momz(:,:,:)   = 0.D0
    qtrc(:,:,:,:) = 0.D0

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       temp(i,j,k) = ENV_THETA - GRAV / CPdry * GRID_CZ(k)

       dist = ( (GRID_CX(i)-XC_BBL)/XR_BBL )**2.D0 &
            + ( (GRID_CY(j)-YC_BBL)/YR_BBL )**2.D0 &
            + ( (GRID_CZ(k)-ZC_BBL)/ZR_BBL )**2.D0

       if ( dist > 1.D0 ) then ! out of warm bubble
          rh(i,j,k) = ENV_RH / ( 1.D0 + GRID_CZ(k) )
       else
          rh(i,j,k) = 100.D0 ! 100%, saturated
       endif

       call moist_psat_water0( temp(i,j,k), psat )

       pres(i,j,k) = Pstd * ( temp(i,j,k)/ENV_THETA )**CPovR - rh(i,j,k)*1.D-2 * psat

       qsat = EPSvap * psat / ( pres(i,j,k) - ( 1.D0-EPSvap )*psat )


       qtrc(i,j,k,I_QV) = rh(i,j,k)*1.D-2 * qsat

       pott(i,j,k) = temp(i,j,k) * ( Pstd/pres(i,j,k) )**RovCP
       dens(i,j,k) = Pstd / Rdry / pott(i,j,k) * ( pres(i,j,k)/Pstd )**CVovCP
       rhot(i,j,k) = dens(i,j,k) * pott(i,j,k)
    enddo
    enddo
    enddo

    call ATMOS_vars_put( dens, momx, momy, momz, rhot, qtrc  )

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine MKEXP_warmbubble

  subroutine moist_psat_water0( t, psat )
    ! psat : Clasius-Clapeyron: based on CPV, CPL constant
    use mod_const, only : &
       Rvap  => CONST_Rvap,  &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL,    &
       LH0   => CONST_LH00,  &
       PSAT0 => CONST_PSAT0, &
       T00   => CONST_TEM00
    implicit none

    real(8), intent(in)  :: t
    real(8), intent(out) :: psat

    real(8)              :: Tmin = 10.D0
    !---------------------------------------------------------------------------

    psat = PSAT0 * ( max(t,Tmin)/T00 ) ** ( ( CPvap-CL )/Rvap ) &
         * exp ( LH0/Rvap * ( 1.0D0/T00 - 1.0D0/max(t,Tmin) ) )

    return
  end subroutine moist_psat_water0

end program warmbubble
