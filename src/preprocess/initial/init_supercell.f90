!-------------------------------------------------------------------------------
!> Program Supercell Test for SCALE-LES ver.3
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)  [new] follow the supercell test of WRF
!! @li      2012-02-16 (Y.Miyamoto) [mod] added hydrostatic balance calculation
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
  use mod_fileio, only: &
     FIO_setup, &
     FIO_finalize
  use mod_grid, only: &
     GRID_setup
  use mod_comm, only: &
     COMM_setup
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
  call TIME_rapstart('Initialize')

  ! setup file I/O
  call FIO_setup

  ! setup horisontal/veritical grid system
  call GRID_setup

  ! setup mpi communication
  call COMM_setup

  ! setup atmosphere
  call ATMOS_vars_setup

  call TIME_rapend('Initialize')


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
  !> Make initial state for supercell experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_supercell
    use mod_stdio, only: &
       IO_get_available_fid, &
       IO_FILECHR,  &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       PI     => CONST_PI,     &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       RovCP  => CONST_RovCP,  &
       RovCV  => CONST_RovCV,  &
       CVovCP => CONST_CVovCP, &
       CPovCV => CONST_CPovCV, &
       EPSvap => CONST_EPSvap, &
       P00    => CONST_PRE00
    use mod_grid, only : &
       KA  => GRID_KA, &
       IA  => GRID_IA, &
       JA  => GRID_JA, &
       KS  => GRID_KS, &
       KE  => GRID_KE, &
       IS  => GRID_IS, &
       IE  => GRID_IE, &
       JS  => GRID_JS, &
       JE  => GRID_JE, &
       CZ  => GRID_CZ, &
       CX  => GRID_CX, &
       CY  => GRID_CY, &
       FDZ => GRID_FDZ
    use mod_atmos_vars, only: &
       var => atmos_var, &
       QA  => A_QA, &
       I_DENS,      &
       I_MOMX,      &
       I_MOMY,      &
       I_MOMZ,      &
       I_RHOT,      &
       I_QV
    use mod_atmos_hydrostatic, only: &
       hydro_buildrho
    implicit none

    character(len=IO_FILECHR) :: ENV_IN_SOUNDING_file = ''
    real(8) :: EXT_TBBL   =   5.D0  ! extremum of temperature in bubble [K]
    real(8) :: ZC_BBL     =   3.D3  ! center location [m]: z
    real(8) :: XC_BBL     =  15.D3  ! center location [m]: x
    real(8) :: YC_BBL     =  15.D3  ! center location [m]: y
    real(8) :: ZR_BBL     =   2.D3  ! bubble radius   [m]: z
    real(8) :: XR_BBL     =   4.D3  ! bubble radius   [m]: x
    real(8) :: YR_BBL     =   4.D3  ! bubble radius   [m]: y

    NAMELIST / PARAM_MKEXP_SUPERCELL / &
       ENV_IN_SOUNDING_file, &
       EXT_TBBL,  &
       ZC_BBL,    &
       XC_BBL,    &
       YC_BBL,    &
       ZR_BBL,    &
       XR_BBL,    &
       YR_BBL

    integer, parameter :: EXP_klim = 100
    integer            :: EXP_kmax

    real(8) :: EXP_z   (EXP_klim) ! height      [m]
    real(8) :: EXP_dens(EXP_klim) ! density     [kg/m3]
    real(8) :: EXP_pott(EXP_klim) ! potential temperature [K]
    real(8) :: EXP_u   (EXP_klim) ! velocity u  [m/s]
    real(8) :: EXP_v   (EXP_klim) ! velocity v  [m/s]
    real(8) :: EXP_qv  (EXP_klim) ! water vapor [g/kg]
    real(8) :: EXP_pres(EXP_klim) ! pressure    [Pa]

    real(8) :: dens(KA) ! density  [kg/m3]
    real(8) :: pres(KA) ! pressure [Pa]
    real(8) :: pott(KA) ! potential temperature [K]
    real(8) :: temp(KA) ! temperature [K]
    real(8) :: velx(KA) ! velocity u [m/s]
    real(8) :: vely(KA) ! velocity v [m/s]
    real(8) :: qv  (KA) ! water vapor mixing ratio [kg/kg]

    real(8) :: pres_sfc
    real(8) :: pott_sfc
    real(8) :: dist

    real(8) :: fact1, fact2, rdz
    real(8) :: RovP, EXP_dens_s, dens_s, dhyd, dgrd

    real(8), parameter :: criteria = 1.D-10
    integer, parameter :: itelim = 100
    integer            :: kref, ite

    integer :: ierr, fid
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[WARMBUBBLE]/Categ[INIT]'

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

    !--- prepare sounding profile
    if( IO_L ) write(IO_FID_LOG,*) '+++ Input sounding file:', trim(ENV_IN_SOUNDING_file)
    fid = IO_get_available_fid()
    open( fid,                                 &
          file   = trim(ENV_IN_SOUNDING_file), &
          form   = 'formatted',                &
          status = 'old',                      &
          iostat = ierr                        )

       if ( ierr /= 0 ) then
          if( IO_L ) write(*,*) 'xxx Input file not found!'
       endif

       !--- read sounding file till end
       read(fid,*) EXP_pres(1), EXP_pott(1), EXP_qv(1)

       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface pressure [hPa]',     EXP_pres(1)
       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface pot. temp  [K]',     EXP_pott(1)
       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface water vapor [g/kg]', EXP_qv(1)

       do k = 2, EXP_klim
          read(fid,*,iostat=ierr) EXP_z(k), EXP_pott(k), EXP_qv(k), EXP_u(k), EXP_v(k)
          if ( ierr /= 0 ) exit
       enddo

       EXP_kmax = k - 1
    close(fid)

    EXP_pres(1) = EXP_pres(1) * 1.D+2
    EXP_dens(1) = P00 / Rdry / EXP_pott(1) * ( EXP_pres(1)/P00 )**CVovCP

    EXP_z(1)    = 0.D0
    EXP_u(1)    = EXP_u(2)
    EXP_v(1)    = EXP_v(2)
    do k = 1, EXP_klim
       EXP_qv(k) = EXP_qv(k) * 1.D-3
    enddo

    do kref = 2, EXP_kmax
       rdz = 1.D0 / ( EXP_z(kref)-EXP_z(kref-1) )

       EXP_dens_s     = 0.D0
       EXP_dens(kref) = EXP_dens(kref-1) ! first guess

       do ite = 1, itelim
          if ( abs(EXP_dens(kref) - EXP_dens_s) <= criteria ) exit

          EXP_dens_s = EXP_dens(kref) 

          dhyd = + ( P00 * ( EXP_dens(kref-1) * Rdry * EXP_pott(kref-1) / P00 )**CPovCV &
                   - P00 * ( EXP_dens_s       * Rdry * EXP_pott(kref)   / P00 )**CPovCV ) * rdz & ! dp/dz
                 - GRAV * 0.5D0 * ( EXP_dens(kref-1) + EXP_dens_s )                               ! rho*g

          dgrd = - P00 * ( Rdry * EXP_pott(kref) / P00 )**CPovCV * rdz &
                 * CPovCV * EXP_dens_s**RovCV                          &
                 - 0.5D0 * GRAV

          EXP_dens(kref) = EXP_dens_s - dhyd/dgrd
       enddo

       EXP_pres(kref) = P00 * ( EXP_dens(kref) * Rdry * EXP_pott(kref) / P00 )**CPovCV
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Input sounding profiles'
    if( IO_L ) write(IO_FID_LOG,*) '  K       Z[m] rho[kg/m3]      P[Pa]   theta[K]     U[m/s]     V[m/s]  Qv[kg/kg]'
    do k = 1, EXP_kmax
       if( IO_L ) write(IO_FID_LOG,'(1x,i3,7(1x,F10.3))') &
       k, EXP_z(k), EXP_dens(k), EXP_pres(k), EXP_pott(k), EXP_u(k), EXP_v(k), EXP_qv(k)
    enddo

    !--- linear interpolate to model grid
    do k    = KS, KE
    do kref = 2, EXP_kmax

       if (       CZ(k) >  EXP_z(kref-1) &
            .AND. CZ(k) <= EXP_z(kref)   ) then

          fact1 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )
          fact2 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
          rdz   = 1.D0 / ( CZ(k)-EXP_z(kref-1) )

          pott(k) = EXP_pott(kref-1) * fact1 &
                  + EXP_pott(kref)   * fact2
          velx(k) = EXP_u   (kref-1) * fact1 &
                  + EXP_u   (kref)   * fact2
          vely(k) = EXP_v   (kref-1) * fact1 &
                  + EXP_v   (kref)   * fact2
          qv(k)   = EXP_qv  (kref-1) * fact1 &
                  + EXP_qv  (kref)   * fact2
       endif
    enddo
    enddo
    pres_sfc = EXP_pres(1)
    pott_sfc = EXP_pott(1)

    ! make density profile
    call hydro_buildrho( dens(:), pres(:), pott(:), pres_sfc, pott_sfc )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Interpolated data'
    if( IO_L ) write(IO_FID_LOG,*) '  K       Z[m] rho[kg/m3]      P[Pa]   theta[K]     U[m/s]     V[m/s]  Qv[kg/kg]'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(1x,i3,7(1x,F10.3))') k, CZ(k), dens(k), pres(k), pott(k), velx(k), vely(k), qv(k)
    enddo

    ! make warm bubble
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dist = ( (CZ(k)-ZC_BBL)/ZR_BBL )**2.D0 &
            + ( (CX(i)-XC_BBL)/XR_BBL )**2.D0 &
            + ( (CY(j)-YC_BBL)/YR_BBL )**2.D0

       if ( dist <= 1.D0 ) then
          var(k,i,j,I_RHOT) = dens(k) &
                            * ( pott(k) &
                              + EXT_TBBL * cos( 0.5D0*PI*sqrt(dist) )**2 * ( P00/pres(k) )**RovCP )
       else
          var(k,i,j,I_RHOT) = dens(k) * pott(k)
       endif

       var(k,i,j,I_DENS) = dens(k)
       var(k,i,j,I_MOMZ) = 0.D0
       var(k,i,j,I_MOMX) = dens(k) * velx(k)
       var(k,i,j,I_MOMY) = dens(k) * vely(k)


       if ( QA > 0 ) then
          do iq = 1,  QA
             var(k,i,j,5+iq) = 0.D0
          enddo
          var(k,i,j,5+I_QV) = qv(k)
       endif

    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine MKEXP_supercell

end program supercell
