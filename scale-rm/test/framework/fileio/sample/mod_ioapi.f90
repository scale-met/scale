!-------------------------------------------------------------------------------
!> module fileio API
!!
!! @par Description
!!          fileio sample
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_ioapi
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: IOAPI
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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> IO API
  subroutine IOAPI()
    use dc_log, only: &
       LogInit
    use gtool_file, only: &
       FileCloseAll
    use scale_precision
    use scale_stdio
    use scale_prof
    use scale_grid_index

    use scale_process, only: &
       PRC_setup,    &
       PRC_MPIstart, &
       PRC_MPIfinish
    use scale_const, only: &
       CONST_setup
    use scale_calendar, only: &
       CALENDAR_setup
    use scale_random, only: &
       RANDOM_setup
    use scale_time, only: &
       TIME_setup
    use scale_grid_index, only: &
       GRID_INDEX_setup
    use scale_grid, only: &
       GRID_setup
    use scale_grid_nest, only: &
       NEST_setup
    use scale_land_grid_index, only: &
       LAND_GRID_INDEX_setup
    use scale_land_grid, only: &
       LAND_GRID_setup
    use scale_urban_grid_index, only: &
       URBAN_GRID_INDEX_setup
    use scale_urban_grid, only: &
       URBAN_GRID_setup
    use scale_tracer, only: &
       TRACER_setup
    use scale_fileio, only: &
       FILEIO_setup, &
       FILEIO_write, &
       FILEIO_read
    implicit none

    character(len=H_MID), parameter :: MODELNAME = "SCALE-RM"

    character(len=H_MID) :: DATATYPE = 'DEFAULT' !< REAL4 or REAL8

    real(RP), allocatable :: DATAforOUTPUT(:,:,:)
    real(RP), allocatable :: DATAforINPUT (:,:,:)

    integer :: k, i, j
    !-----------------------------------------------------------------------------

    ! setup standard I/O
    call IO_setup( MODELNAME )

    ! start MPI
    call PRC_MPIstart

    ! setup process
    call PRC_setup

    ! setup Log
    call LogInit(IO_FID_CONF, IO_FID_LOG, IO_L)

    ! setup constants
    call CONST_setup

    ! setup time
    call TIME_setup( setup_TimeIntegration = .false. )

    call PROF_rapstart('Initialize')

    ! setup horizontal/vertical grid coordinates (cartesian,idealized)
    call GRID_INDEX_setup
    call GRID_setup

    call LAND_GRID_INDEX_setup
    call LAND_GRID_setup

    call URBAN_GRID_INDEX_setup
    call URBAN_GRID_setup

    ! setup file I/O
    call FILEIO_setup

    call PROF_rapend('Initialize')

    !########## main ##########

    call PROF_rapstart('Main')
    if( IO_L ) write(IO_FID_LOG,*) '*** IOAPI test ***'



    ! Write file

    allocate( DATAforOUTPUT(KA,IA,JA) )

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DATAforOUTPUT(k,i,j) = (j-1)*IA*KA + (i-1)*KA + k
       if( IO_L ) write(IO_FID_LOG,*) 'Output(k,i,j,value): ', k, i, j, int(DATAforOUTPUT(k,i,j))
    enddo
    enddo
    enddo

    call FILEIO_write( DATAforOUTPUT(:,:,:), "filename_test", "test outputfile",      & ! [IN]
                       'varname_test', 'description of var', 'unit', 'ZXY', DATATYPE  ) ! [IN]

    deallocate( DATAforOUTPUT )

    call FileCloseAll


    ! Read file

    allocate( DATAforINPUT(KA,IA,JA) )
    DATAforINPUT(:,:,:) = 0.0_RP

    call FILEIO_read( DATAforINPUT(:,:,:),                           & ! [OUT]
                      "filename_test", 'varname_test', 'ZXY', step=1 ) ! [IN]

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       if( IO_L ) write(IO_FID_LOG,*) 'Intput(k,i,j,value): ', k, i, j, int(DATAforINPUT(k,i,j))
    enddo
    enddo
    enddo



    if( IO_L ) write(IO_FID_LOG,*) '*** IOAPI test end ***'
    call PROF_rapend('Main')

    !########## Finalize ##########

    call PROF_rapreport

    call FileCloseAll

    ! stop MPI
    call PRC_MPIfinish

    return
  end subroutine IOAPI

end module mod_ioapi
