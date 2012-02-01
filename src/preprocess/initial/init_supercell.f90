!-------------------------------------------------------------------------------
!> Program Supecell Test for SCALE-LES ver.3
!!
!! @par Description
!!          Initial data for Quarter circle shear Super-Cell Test
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
program supercell
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
  use mod_atmos_refstate, only: &
     ATMOS_REFSTATE_setup
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
  call MKEXP_supercell

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
  subroutine MKEXP_supercell
    use mod_stdio, only: &
       IO_get_available_fid, &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L, &
       IO_FILECHR
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       PI     => CONST_PI,     &
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

    character(len=IO_FILECHR) :: ENV_IN_SOUNDING_file
    real(8) :: XC_BBL = 18.D3     ! center location [m]: x
    real(8) :: YC_BBL = 18.D3     ! center location [m]: y
    real(8) :: ZC_BBL = 3.D3      ! center location [m]: z
    real(8) :: XR_BBL = 4.D3      ! bubble radius   [m]: x
    real(8) :: YR_BBL = 4.D3      ! bubble radius   [m]: y
    real(8) :: ZR_BBL = 2.D3      ! bubble radius   [m]: z

    NAMELIST / PARAM_MKEXP_SUPERCELL / &
       ENV_IN_SOUNDING_file, &
       XC_BBL,    &
       YC_BBL,    &
       ZC_BBL,    &
       XR_BBL,    &
       YR_BBL,    &
       ZR_BBL

    integer, parameter :: EXP_klim = 100
    integer            :: EXP_kmax

    real(8) :: EXP_z    (EXP_klim) ! height [m]
    real(8) :: EXP_pres (EXP_klim) ! pressure [hPa]
    real(8) :: EXP_theta(EXP_klim) ! potential temperature [K]
    real(8) :: EXP_qv   (EXP_klim) ! water vapor [g/kg]
    real(8) :: EXP_u    (EXP_klim) ! velocity u [m/s]
    real(8) :: EXP_v    (EXP_klim) ! velocity v [m/s]

    real(8) :: dens(KA,IA,JA)      ! density     [kg/m3]
    real(8) :: momx(KA,IA,JA)      ! momentum(x) [kg/m3 * m/s]
    real(8) :: momy(KA,IA,JA)      ! momentum(y) [kg/m3 * m/s]
    real(8) :: momz(KA,IA,JA)      ! momentum(z) [kg/m3 * m/s]
    real(8) :: rhot(KA,IA,JA)      ! rho * theta [kg/m3 * K]
    real(8) :: qtrc(KA,IA,JA,QA)   ! tracer mixing ratio [kg/kg],[1/m3]

    real(8) :: pott(KA,IA,JA)      ! potential temperature [K]

    real(8) :: pres (KA)
    real(8) :: velx (KA)
    real(8) :: vely (KA)
    real(8) :: theta(KA)
    real(8) :: qv   (KA)

    real(8) :: dist
    real(8) :: gmr, fact1, fact2

    integer :: i, j, k, kref
    integer :: fid, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[supercell]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_SUPERCELL,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_SUPERCELL. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_SUPERCELL)

    call ATMOS_vars_get( dens, momx, momy, momz, rhot, qtrc )

    if( IO_L ) write(IO_FID_LOG,*) '+++ Input sounding file:', trim(ENV_IN_SOUNDING_file)

    !--- Open config file till end
    fid = IO_get_available_fid()
    open( fid,                                 &
          file   = trim(ENV_IN_SOUNDING_file), &
          form   = 'formatted',                &
          status = 'old',                      &
          iostat = ierr                        )

       if ( ierr /= 0 ) then
          if( IO_L ) write(*,*) 'xxx Input file not found!'
       endif

       read(fid,*) EXP_pres(1), EXP_theta(1), EXP_qv(1)

       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface pressure [hPa]',     EXP_pres(1)
       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface pot. temp  [K]',     EXP_theta(1)
       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface water vapor [g/kg]', EXP_qv(1)

       do k = 2, EXP_klim
          read(fid,*,iostat=ierr) EXP_z(k), EXP_theta(k), EXP_qv(k), EXP_u(k), EXP_v(k)
          if ( ierr /= 0 ) exit
       enddo

       EXP_kmax = k - 1

    close(fid)

    gmr      = GRAV / Rdry

    EXP_z(1)    = 0.D0
    EXP_pres(1) = EXP_pres(1) * 1.D+2
    EXP_u(1)    = EXP_u(2)
    EXP_v(1)    = EXP_v(2)

    do kref = 2, EXP_kmax
	    EXP_pres(kref) = EXP_pres(kref-1) &
                      * exp( -gmr / EXP_theta(kref-1) * ( EXP_z(kref)-EXP_z(kref-1) ) )
    enddo

    EXP_qv(:) = EXP_qv(:) * 1.D-3

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Input sounding profiles'
    do k = 1, EXP_klim
       if( IO_L ) write(IO_FID_LOG,'(1x,i3,6(1x,F10.3))') &
       k, EXP_z(k), EXP_pres(k), EXP_theta(k), EXP_qv(k), EXP_u(k), EXP_v(k)
    enddo

    !--- make reference state
    do k = KS, KE
       do kref = 2, EXP_kmax
          if(       GRID_CZ(k) >  EXP_z(kref-1) &
              .AND. GRID_CZ(k) <= EXP_z(kref)   ) then

             fact1 = ( GRID_CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref) - EXP_z(kref-1) )
             fact2 = ( EXP_z(kref) - GRID_CZ(k)   ) / ( EXP_z(kref) - EXP_z(kref-1) )

             theta(k) = EXP_theta(kref-1) * fact1 &
                      + EXP_theta(kref)   * fact2

             pres(k) = EXP_pres(kref-1) &
                     * exp( -gmr / EXP_theta(kref-1) * ( GRID_CZ(k)-EXP_z(kref-1) ) )

             velx(k) = EXP_u(kref)   * fact1 &
                     + EXP_u(kref+1) * fact2
             vely(k) = EXP_v(kref)   * fact1 &
                     + EXP_v(kref+1) * fact2

             qv(k) = EXP_qv(kref)   * fact1 &
                   + EXP_qv(kref+1) * fact2
          endif
       enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Interpolated data'
    do k = 1, KA
       if( IO_L ) write(IO_FID_LOG,'(1x,i3,6(1x,F10.3))') k, GRID_CZ(k), pres(k), theta(k), velx(k), vely(k), qv(k)
    enddo

    momz(:,:,:)   = 0.D0
    qtrc(:,:,:,:) = 0.D0

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       dist = ( (GRID_CX(i)-XC_BBL)/XR_BBL )**2.D0 &
            + ( (GRID_CY(j)-YC_BBL)/YR_BBL )**2.D0 &
            + ( (GRID_CZ(k)-ZC_BBL)/ZR_BBL )**2.D0

       if ( dist > 1.D0 ) then ! out of cold bubble
          pott(k,i,j) = theta(k)
       else
          pott(k,i,j) = theta(k) + 5.D0 * dcos( 0.5D0*PI*sqrt(dist) )**2 &
                      * ( Pstd/pres(k) )**RovCP
       endif
       dens(k,i,j)   = Pstd / Rdry / pott(k,i,j) * ( pres(k)/Pstd )**CVovCP

       rhot(k,i,j)   = dens(k,i,j) * pott(k,i,j)
       momx(k,i,j)   = dens(k,i,j) * velx(k)
       momy(k,i,j)   = dens(k,i,j) * vely(k)
       qtrc(k,i,j,1) = qv(k)
    enddo
    enddo
    enddo

    call ATMOS_vars_put( dens, momx, momy, momz, rhot, qtrc  )

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine MKEXP_supercell

end program supercell
