!-------------------------------------------------------------------------------
!> Program Consistency With Continuity (CWC) Test for SCALE-LES ver.3
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-02-16 (Y.Miyamoto) [new] constructed from the cold bubble test
!!
!<
!-------------------------------------------------------------------------------
program cwc
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
  call MKEXP_cwc

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
  subroutine MKEXP_cwc
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
       Pstd   => CONST_Pstd
    use mod_grid, only : &
       KA => GRID_KA, &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KS => GRID_KS, &
       KE => GRID_KE, &
       IS => GRID_IS, &
       IE => GRID_IE, &
       JS => GRID_JS, &
       JE => GRID_JE, &
       GRID_CZ,       &
       GRID_CX,       &
       GRID_CY
    use mod_atmos_vars, only: &
       QA => A_QA,     &
       ATMOS_vars_get, &
       ATMOS_vars_put
    implicit none

    real(8) :: ENV_THETA = 300.D0 ! Potential Temperature of environment [K]
    real(8) :: ENV_QTRC  = 1.D0   ! tracer of environment
    real(8) :: EXT_TBBL  = 5.D0   ! extremum of temperature in bubble [K]
    real(8) :: ZC_BBL = 3.D3      ! center location [m]: z
    real(8) :: XC_BBL = 15.D3     ! center location [m]: x
    real(8) :: YC_BBL = 15.D3     ! center location [m]: y
    real(8) :: ZR_BBL = 2.D3      ! bubble radius   [m]: z
    real(8) :: XR_BBL = 1.D3      ! bubble radius   [m]: x
    real(8) :: YR_BBL = 1.D3      ! bubble radius   [m]: y

    NAMELIST / PARAM_MKEXP_CWC / &
       ENV_THETA, &
       ENV_QTRC,  &
       EXT_TBBL,  &
       ZC_BBL,    &
       XC_BBL,    &
       YC_BBL,    &
       ZR_BBL,    &
       XR_BBL,    &
       YR_BBL

    real(8) :: dens(KA,IA,JA)      ! density     [kg/m3]
    real(8) :: momx(KA,IA,JA)      ! momentum(x) [kg/m3 * m/s]
    real(8) :: momy(KA,IA,JA)      ! momentum(y) [kg/m3 * m/s]
    real(8) :: momz(KA,IA,JA)      ! momentum(z) [kg/m3 * m/s]
    real(8) :: rhot(KA,IA,JA)      ! rho * theta [kg/m3 * K]
    real(8) :: qtrc(KA,IA,JA,QA)   ! tracer mixing ratio [kg/kg],[1/m3]

    real(8) :: pres(KA,IA,JA)    ! pressure [Pa]
    real(8) :: temp(KA,IA,JA)    ! temperature [K]
    real(8) :: pott(KA,IA,JA)    ! potential temperature [K]

    real(8) :: dist, RovP
    real(8) :: DENS_Z0, dhyd, dgrd, dz, tt, pp, dd, d1, d2

    integer :: k, i, j, n
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CWC]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_cwc,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_cwc. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_cwc)


    call ATMOS_vars_get( dens, momx, momy, momz, rhot, qtrc )

    pott(:,:,:)   = ENV_THETA
    momx(:,:,:)   = 0.D0
    momy(:,:,:)   = 0.D0
    momz(:,:,:)   = 0.D0
    qtrc(:,:,:,:) = ENV_QTRC

    RovP = Rdry / (Pstd)**CPovR
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

       pres(k,i,j) = ( dens(k,i,j) * Rdry * pott(k,i,j) )**( CPdry/CVdry ) * ( Pstd )**( -Rdry/CVdry )
       temp(k,i,j) = pres(k,i,j) / dens(k,i,j) * Rdry

    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       dist = ( (GRID_CZ(k)-ZC_BBL)/ZR_BBL )**2.D0 &
            + ( (GRID_CX(i)-XC_BBL)/XR_BBL )**2.D0 &
            + ( (GRID_CY(j)-YC_BBL)/YR_BBL )**2.D0

       if ( dist <= 1.D0 ) then
          pott(k,i,j) = ENV_THETA &
                      - EXT_TBBL * dcos( 0.5D0*PI*sqrt(dist) )**2.D0 &
                      * ( Pstd/pres(k,i,j) )**RovCP
       endif

       rhot(k,i,j) = dens(k,i,j) * pott(k,i,j)

    enddo
    enddo
    enddo

    call ATMOS_vars_put( dens, momx, momy, momz, rhot, qtrc  )

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine MKEXP_cwc

end program cwc
