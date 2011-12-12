!-------------------------------------------------------------------------------
!> Program Reference state generator SCALE-LES ver.3
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
program refstate
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_setup
  use mod_process, only: &
     PRC_setup,      &
     PRC_NOMPIstart, &
     PRC_MPIstop
  use mod_const, only: &
     CONST_setup
  use mod_time, only: &
     TIME_setup,           &
     TIME_rapstart,        &
     TIME_rapend,          &
     TIME_rapreport
  use mod_grid, only: &
     GRID_setup
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

  ! start WITHOUT MPI
  call PRC_NOMPIstart

  ! setup process
  call PRC_setup

  ! setup constants
  call CONST_setup

  ! setup time
  call TIME_setup

  ! setup horisontal/veritical grid system
  call GRID_setup


  !########## main ##########

  call TIME_rapstart('Main')

  ! make initial state (restart)
  call MKEXP_refstate

  call TIME_rapend('Main')


  !########## Finalize ##########
  call TIME_rapreport

  ! stop MPI
  call PRC_MPIstop

  stop
  !=============================================================================
contains

  !-----------------------------------------------------------------------------
  !> Make initial state for cold bubble experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_refstate
    use mod_stdio, only: &
       IO_get_available_fid, &
       IO_FILECHR,  &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CONST_UNDEF8,          &
       GRAV   => CONST_GRAV,  &
       Rair   => CONST_Rair,  &
       RovCP  => CONST_RovCP, &
       Pstd   => CONST_Pstd,  &
       Tstd   => CONST_Tstd
    use mod_grid, only : &
       KA   => GRID_KA,   &
       KMAX => GRID_KMAX, &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       CZ   => GRID_CZ
    implicit none

    character(LEN=IO_FILECHR) :: REF_FILENAME = 'refstate.txt' ! output file
    real(8)                   :: TEMP_SFC     = Tstd           ! surface temperature

    NAMELIST / PARAM_MKEXP_refstate / &
       REF_FILENAME, &
       TEMP_SFC

    integer, parameter :: nref = 8
    real(8), parameter :: CZ_ref  (nref) = (/     0.0D0,  &
                                              11000.0D0,  &
                                              20000.0D0,  &
                                              32000.0D0,  &
                                              47000.0D0,  &
                                              51000.0D0,  &
                                              71000.0D0,  &
                                              84852.0D0   /)
    real(8), parameter :: GAMMA(nref) = (/       -6.5D-3, &
                                                  0.0D0 , &
                                                  1.0D-3, &
                                                  2.8D-3, &
                                                  0.0D-3, &
                                                 -2.8D-3, &
                                                 -2.0D-3, &
                                                  0.0D0   /)
    real(8) :: temp_ref(nref)
    real(8) :: pres_ref(nref)

    real(8), allocatable :: temp(:) !! reference temperature
    real(8), allocatable :: pres(:) !! reference pressure
    real(8), allocatable :: dens(:) !! reference density
    real(8), allocatable :: pott(:) !! reference potential temperature

    real(8) :: gmr !! grav / Rair
    integer :: i, k
    integer :: fid, ierr
    !---------------------------------------------------------------------------

    gmr = GRAV / Rair

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING REFERENCE DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[REFSTATE]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_refstate,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_refstate. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_refstate)


    !--- ISA profile
    temp_ref(1) = Tstd
    pres_ref(1) = Pstd

    do i = 2, nref
       temp_ref(i) = temp_ref(i-1) + GAMMA(i-1) * ( CZ_ref(i)-CZ_ref(i-1) )

       if ( GAMMA(i-1) == 0.D0 ) then
	      pres_ref(i) = pres_ref(i-1) * exp( -gmr / temp_ref(i) * ( CZ_ref(i)-CZ_ref(i-1) ) )
       else
	      pres_ref(i) = pres_ref(i-1) * ( temp_ref(i)/temp_ref(i-1) ) ** ( -gmr/GAMMA(i-1) )
       endif
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '###### ICAO International Standard Atmosphere ######'
    if( IO_L ) write(IO_FID_LOG,*) 'height      : lapse rate  : pressure   : temperature'
    do i = 1, nref
       if( IO_L ) write(IO_FID_LOG,'(1x,4(f13.5))') CZ_ref(i), GAMMA(i), pres_ref(i), temp_ref(i)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '####################################################'


    allocate( pres(KA) ); pres = CONST_UNDEF8
    allocate( temp(KA) ); temp = CONST_UNDEF8

    !--- make reference state
    do k = KS, KE
    do i = 2, nref
       if ( CZ(k) > CZ_ref(i-1) .AND. CZ(k) <= CZ_ref(i) ) then

          temp(k) = temp_ref(i-1) + GAMMA(i-1) * ( CZ(k)-CZ_ref(i-1) )

          if ( GAMMA(i-1) == 0.D0 ) then
             pres(k) = pres_ref(i-1) * exp( -gmr/temp_ref(i-1) * ( CZ(k)-CZ_ref(i-1) ) )
          else
             pres(k) = pres_ref(i-1) * ( temp(k)/temp_ref(i-1) ) ** ( -gmr/GAMMA(i-1) )
          endif
       elseif ( CZ(k) <= CZ_ref(1)    ) then

          temp(k) = temp_ref(1) + GAMMA(1) * ( CZ(k)-CZ_ref(1) )
          pres(k) = pres_ref(1) * ( temp(k)/temp_ref(1) ) ** ( -gmr/GAMMA(1) )

       elseif ( CZ(k)  > CZ_ref(nref) ) then

          temp(k) = temp(k-1)
          pres(k) = pres_ref(i-1) * exp( -gmr/temp_ref(i-1) * ( CZ(k)-CZ_ref(i-1) ) )

       endif
    enddo
    enddo


    allocate( pott(KA) ); pott = CONST_UNDEF8
    allocate( dens(KA) ); dens = CONST_UNDEF8

    pott(KS:KE) = temp(KS:KE) * ( Pstd/pres(KS:KE) )**RovCP
    dens(KS:KE) = pres(KS:KE) / ( temp(KS:KE) * Rair )


    if( IO_L ) write(IO_FID_LOG,*) ' height      : pressure   : temperature: density   : pot.temp.'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(5(f13.5))') CZ(k), pres(k), temp(k), dens(k), pott(k)
    enddo

    fid = IO_get_available_fid()
    open (fid, file=trim(REF_FILENAME), form='formatted', iostat=ierr)
       do k = KS, KE
          write(fid,*) CZ(k), dens(k), pott(k)
       enddo
    close(fid)

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine MKEXP_refstate

end program refstate
