!-------------------------------------------------------------------------------
!> Program Lamb wave Test for SCALE-LES ver.3
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-01-31 (Y.Miyamoto) [new] for Lamb wave test
!!
!!
!<
!-------------------------------------------------------------------------------
program lambwave
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
  call MKEXP_lambwave

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
  subroutine MKEXP_lambwave
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

    real(8) :: ENV_THETA = 300.D0  ! environmental potential temperature [K]
    real(8) :: ZC_BBL    =   3.D3  ! center location [m]: z
    real(8) :: XC_BBL    =  18.D3  ! center location [m]: x
    real(8) :: YC_BBL    =  18.D3  ! center location [m]: y
    real(8) :: ZR_BBL    =   2.D3  ! bubble radius   [m]: z
    real(8) :: XR_BBL    =   4.D3  ! bubble radius   [m]: x
    real(8) :: YR_BBL    =   4.D3  ! bubble radius   [m]: y
    real(8) :: DELTA_P   =   4.D3  ! peak of pressure perturbation [Pa]
    real(8) :: VER_MODE  =   4.D3  ! vertical mode [-]

    NAMELIST / PARAM_MKEXP_LAMBWAVE / &
       ENV_THETA, &
       ZC_BBL,    &
       XC_BBL,    &
       YC_BBL,    &
       ZR_BBL,    &
       XR_BBL,    &
       YR_BBL,    &
       DELTA_P,   &
       VER_MODE

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
    real(8) :: dist, radi
    real(8) :: DENS_Z0, dhyd, dgrd, dz, tt, pp, dd, d1, d2

    integer :: k, i, j, n
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[LAMBWAVE]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_lambwave,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_lambwave. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_lambwave)

    call ATMOS_vars_get( dens, momx, momy, momz, rhot, qtrc )

    temp(:,:,:)   = ENV_THETA
    rh  (:,:,:)   = 0.D0
    momx(:,:,:)   = 0.D0
    momy(:,:,:)   = 0.D0
    momz(:,:,:)   = 0.D0
    qtrc(:,:,:,:) = 0.D0
    radi = dsqrt( XR_BBL**2.D0 + YR_BBL**2.D0 )

    tt = ENV_THETA - GRAV / CPdry * GRID_CZ(KS)
    pp = Pstd * ( tt/ENV_THETA )**CPovR
    DENS_Z0 = Pstd / Rdry / ENV_THETA * ( pp/Pstd )**CVovCP

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       if ( k == KS ) then
          dens(k,i,j) = DENS_Z0
       else
          dz = GRID_CZ(k) - GRID_CZ(k-1)
          dhyd = 0.D0
          d1 = 0.D0
          d2 = dens(k-1,i,j)
          n = 0
          do while ( dabs(d2-d1) > 1.D-10 )
             n = n + 1
             d1 = d2
             dhyd = - ( Pstd**( -RovCP )*Rdry*pott(k  ,i,j)*d1            )**( CPdry/CVdry ) / dz - 0.5D0*GRAV*d1 &   
                    + ( Pstd**( -RovCP )*Rdry*pott(k-1,i,j)*dens(k-1,i,j) )**( CPdry/CVdry ) / dz - 0.5D0*GRAV*dens(k-1,i,j)
             dgrd = - ( Pstd**( -RovCP )*Rdry*pott(k,i,j) )**( CPdry/CVdry ) *CPdry/CVdry/dz * d1**( Rdry/CVdry ) - 0.5D0*GRAV 
             d2 = d1 - dhyd / dgrd
          end do
          dens(k,i,j) = d2
          if ( n < 100 ) write(IO_FID_LOG,*) 'iteration converged',n,dhyd,d2,d1
       end if
                       
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       pres(k,i,j) = dens(k,i,j) * Rdry * temp(k,i,j)

       dist = dsqrt( ( GRID_CX(i)-XC_BBL )**2.D0 &
                   + ( GRID_CY(j)-YC_BBL )**2.D0 )

       if ( dist <= radi ) then
          pres(k,i,j) = pres(k,i,j) &
                      + DELTA_P * 0.5D0 * ( 1.D0 + dcos( PI*dist/radi ) ) &
                      * dsin( VER_MODE*PI*GRID_CZ(k)/GRID_CZ(KE) )
       endif

       pott(k,i,j) = temp(k,i,j) * ( Pstd / pres(k,i,j) )**RovCP
       rhot(k,i,j) = dens(k,i,j) * pott(k,i,j)

    enddo
    enddo
    enddo

    call ATMOS_vars_put( dens, momx, momy, momz, rhot, qtrc  )

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*) 

    return
  end subroutine MKEXP_lambwave

end program lambwave
