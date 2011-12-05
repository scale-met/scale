!-------------------------------------------------------------------------------
!> Program Cold BUbble Test for SCALE-LES ver.3
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
program coldbubble
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

  ! setup horisontal/veritical grid system
  call GRID_setup

  ! setup file I/O
  call FIO_setup

  ! setup atmosphere
  call ATMOS_vars_setup


  !########## main ##########

  call TIME_rapstart('Main')

  ! make initial state (restart)
  call MKEXP_coldbubble

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
  subroutine MKEXP_coldbubble
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       PI     => CONST_PI,    &
       GRAV   => CONST_GRAV,  &
       Rair   => CONST_Rair,  &
       CPair  => CONST_CPair, &
       CPovR  => CONST_CPovR, &
       RovCP  => CONST_RovCP, &
       CVovCP => CONST_CVovCP, &
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
       QA => A_QA,                   &
       ATMOS_vars_get, &
       ATMOS_vars_putDMP
    implicit none

    real(8) :: ENV_THETA = 300.D0 ! Potential Temperature of environment
    real(8) :: XC_BBL = 18.D3 ! center location [m]: x
    real(8) :: YC_BBL = 18.D3 ! center location [m]: y
    real(8) :: ZC_BBL = 3.D3  ! center location [m]: z
    real(8) :: XR_BBL = 4.D3  ! bubble radius   [m]: x
    real(8) :: YR_BBL = 4.D3  ! bubble radius   [m]: y
    real(8) :: ZR_BBL = 2.D3  ! bubble radius   [m]: z

    NAMELIST / PARAM_MKEXP_coldbubble / &
       ENV_THETA, &
       XC_BBL,    &
       YC_BBL,    &
       ZC_BBL,    &
       XR_BBL,    &
       YR_BBL,    &
       ZR_BBL

    real(8) :: dens(IA,JA,KA)    ! density [kg/m**3]
    real(8) :: momx(IA,JA,KA)    ! momentum (x) [kg/m**3 * m/s]
    real(8) :: momy(IA,JA,KA)    ! momentum (y) [kg/m**3 * m/s]
    real(8) :: momz(IA,JA,KA)    ! momentum (z) [kg/m**3 * m/s]
    real(8) :: lwpt(IA,JA,KA)    ! liquid water potential temperature [K]

    real(8) :: qtrc(IA,JA,KA,QA) ! tracer mixing ratio   [kg/kg],[1/m3]

    real(8) :: pres(IA,JA,KA)    ! pressure [Pa]
    real(8) :: velx(IA,JA,KA)    ! velocity (x) [m/s]
    real(8) :: vely(IA,JA,KA)    ! velocity (y) [m/s]
    real(8) :: velz(IA,JA,KA)    ! velocity (z) [m/s]
    real(8) :: temp(IA,JA,KA)    ! temperature [K]

    real(8) :: dist

    integer :: i, j, k
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[COLDBUBBLE]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_coldbubble,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_coldbubble. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_coldbubble)


    call ATMOS_vars_get( dens, momx, momy, momz, lwpt, qtrc, &
                         pres, velx, vely, velz, temp        )

    momx(:,:,:)   = 0.D0
    momy(:,:,:)   = 0.D0
    momz(:,:,:)   = 0.D0
    qtrc(:,:,:,:) = 0.D0

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       temp(i,j,k) = ENV_THETA - GRAV / CPair * GRID_CZ(k)

       pres(i,j,k) = Pstd * ( temp(i,j,k)/ENV_THETA )**CPovR

       dist = ( (GRID_CX(i)-XC_BBL)/XR_BBL )**2.D0 &
            + ( (GRID_CY(j)-YC_BBL)/YR_BBL )**2.D0 &
            + ( (GRID_CZ(k)-ZC_BBL)/ZR_BBL )**2.D0

       if ( dist > 1.D0 ) then ! out of cold bubble
          lwpt(i,j,k) = ENV_THETA
       else
          lwpt(i,j,k) = ENV_THETA &
                      - 15.D0 * dcos( 0.5D0*PI*sqrt(dist) )**2.D0 &
                      * ( Pstd/pres(i,j,k) )**RovCP
       endif

       dens(i,j,k) = Pstd / Rair / lwpt(i,j,k) &
                   * ( pres(i,j,k)/Pstd )**CVovCP
    enddo
    enddo
    enddo

    call ATMOS_vars_putDMP( dens, momx, momy, momz, lwpt, qtrc  )

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

  end subroutine MKEXP_coldbubble

end program coldbubble
