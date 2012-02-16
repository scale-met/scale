!-------------------------------------------------------------------------------
!> Program Diffusion Test for SCALE-LES ver.3
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-02-06 (Y.Miyamoto) [new]
!!
!!
!<
!-------------------------------------------------------------------------------
program diffusion
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
  call MKEXP_diffusion

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
  !> Make initial state for diffusion experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_diffusion
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       PI     => CONST_PI,     &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       CPdry  => CONST_CPdry,  &
       CVdry  => CONST_CVdry,  &
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

    real(8) :: ENV_THETA  = 300.D0 ! Potential Temperature of environment [K]
    real(8) :: ENV_RH     =  0.D0  ! Relative Humidity of environment [%]
    real(8) :: ENV_DENS   =  1.D0  ! environment density [kg m-3]
    real(8) :: AMP_DIST   =  0.1D0 ! amplitude of the potential temperature disturbance [K]

    NAMELIST / PARAM_MKEXP_DIFFUSION / &
       ENV_THETA,  &
       ENV_RH,     &
       ENV_DENS,   &
       AMP_DIST 

    real(8) :: dens(KA,IA,JA)      ! density     [kg/m3]
    real(8) :: momx(KA,IA,JA)      ! momentum(x) [kg/m3 * m/s]
    real(8) :: momy(KA,IA,JA)      ! momentum(y) [kg/m3 * m/s]
    real(8) :: momz(KA,IA,JA)      ! momentum(z) [kg/m3 * m/s]
    real(8) :: rhot(KA,IA,JA)      ! rho * theta [kg/m3 * K]
    real(8) :: qtrc(KA,IA,JA,QA)   ! tracer mixing ratio [kg/kg],[1/m3]

    real(8) :: pres(KA,IA,JA)    ! pressure [Pa]
    real(8) :: temp(KA,IA,JA)    ! temperature [K]
    real(8) :: pott(KA,IA,JA)    ! potential temperature [K]

    real(8) :: rh(KA,IA,JA)
    real(8) :: rndm(KA,IA,JA)
    real(8) :: RovP, dist, dz

    integer :: i, j, k
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[DIFFUSION]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_diffusion,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_diffusion. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_diffusion)

    call ATMOS_vars_get( dens, momx, momy, momz, rhot, qtrc )

    pott(:,:,:)   = ENV_THETA
    rh  (:,:,:)   = ENV_RH
    dens(:,:,:)   = ENV_DENS
    momx(:,:,:)   = 0.D0
    momy(:,:,:)   = 0.D0
    momz(:,:,:)   = 0.D0
    qtrc(:,:,:,:) = 0.D0
    RovP = Rdry / (Pstd)**CPovR
    call random_number(rndm)

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       pres(k,i,j) = ( dens(k,i,j) * Rdry * pott(k,i,j) )**( CPdry/CVdry ) * ( Pstd )**( -Rdry/CVdry )
       temp(k,i,j) = pres(k,i,j) / dens(k,i,j) * Rdry

       pott(k,i,j) = pott(k,i,j) + AMP_DIST * rndm(k,i,j)
       rhot(k,i,j) = dens(k,i,j) * pott(k,i,j)

    enddo
    enddo
    enddo

    call ATMOS_vars_put( dens, momx, momy, momz, rhot, qtrc  )

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*) 

    return
  end subroutine MKEXP_diffusion

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

end program diffusion
