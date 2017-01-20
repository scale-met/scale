!-------------------------------------------------------------------------------
!> module LANDUSE
!!
!! @par Description
!!          Land use category module
!!          Manage land/lake/urban/PFT fraction and PFT index
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_landuse
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LANDUSE_setup
  public :: LANDUSE_calc_fact
  public :: LANDUSE_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: LANDUSE_fact_ocean(:,:) !< ocean factor
  real(RP), public, allocatable :: LANDUSE_fact_land (:,:) !< land  factor
  real(RP), public, allocatable :: LANDUSE_fact_urban(:,:) !< urban factor

  real(RP), public, allocatable :: LANDUSE_frac_land (:,:) !< land  fraction
  real(RP), public, allocatable :: LANDUSE_frac_lake (:,:) !< lake  fraction
  real(RP), public, allocatable :: LANDUSE_frac_urban(:,:) !< urban fraction

  integer,  public              :: LANDUSE_PFT_mosaic = 2   !< number of PFT mosaic
  integer,  public              :: LANDUSE_PFT_nmax   = 15  !< number of plant functional type(PFT)

  real(RP), public, allocatable :: LANDUSE_frac_PFT (:,:,:) !< fraction of PFT for each mosaic
  integer,  public, allocatable :: LANDUSE_index_PFT(:,:,:) !< index of PFT for each mosaic

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LANDUSE_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: LANDUSE_IN_BASENAME  = ''                  !< basename of the input  file
  logical,                private :: LANDUSE_IN_CHECK_COORDINATES = .true.      !< switch for check of coordinates
  character(len=H_LONG),  private :: LANDUSE_OUT_BASENAME = ''                  !< basename of the output file
  character(len=H_MID),   private :: LANDUSE_OUT_TITLE    = 'SCALE-RM LANDUSE'  !< title    of the output file
  character(len=H_SHORT), private :: LANDUSE_OUT_DTYPE    = 'DEFAULT'           !< REAL4 or REAL8
  logical,                private :: LANDUSE_AllOcean     = .false.
  logical,                private :: LANDUSE_AllLand      = .false.
  logical,                private :: LANDUSE_AllUrban     = .false.
  logical,                private :: LANDUSE_MosaicWorld  = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LANDUSE_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_LANDUSE / &
       LANDUSE_IN_BASENAME,          &
       LANDUSE_IN_CHECK_COORDINATES, &
       LANDUSE_OUT_BASENAME,         &
       LANDUSE_OUT_DTYPE,            &
       LANDUSE_PFT_mosaic,           &
       LANDUSE_PFT_nmax,             &
       LANDUSE_AllOcean,             &
       LANDUSE_AllLand,              &
       LANDUSE_AllUrban,             &
       LANDUSE_MosaicWorld

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[LANDUSE] / Categ[COUPLER] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LANDUSE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LANDUSE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LANDUSE)

    allocate( LANDUSE_frac_land (IA,JA) )
    allocate( LANDUSE_frac_lake (IA,JA) )
    allocate( LANDUSE_frac_urban(IA,JA) )
    LANDUSE_frac_land (:,:) = 0.0_RP
    LANDUSE_frac_lake (:,:) = 0.0_RP
    LANDUSE_frac_urban(:,:) = 0.0_RP

    allocate( LANDUSE_index_PFT(IA,JA,LANDUSE_PFT_mosaic) )
    allocate( LANDUSE_frac_PFT (IA,JA,LANDUSE_PFT_mosaic) )
    LANDUSE_frac_PFT (:,:,:) = 0.0_RP
    LANDUSE_frac_PFT (:,:,1) = 1.0_RP ! tentative, mosaic is off
    LANDUSE_index_PFT(:,:,:) = 1      ! default

    allocate( LANDUSE_fact_ocean(IA,JA) )
    allocate( LANDUSE_fact_land (IA,JA) )
    allocate( LANDUSE_fact_urban(IA,JA) )
    LANDUSE_fact_ocean(:,:) = 0.0_RP
    LANDUSE_fact_land (:,:) = 0.0_RP
    LANDUSE_fact_urban(:,:) = 0.0_RP


    if    ( LANDUSE_AllOcean ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Assume all grids are ocean'
       call LANDUSE_calc_fact
    elseif( LANDUSE_AllLand ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Assume all grids are land'
       LANDUSE_frac_land (:,:) = 1.0_RP
       call LANDUSE_calc_fact
    elseif( LANDUSE_AllUrban ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Assume all grids are land'
       LANDUSE_frac_land (:,:) = 1.0_RP
       if( IO_L ) write(IO_FID_LOG,*) '*** Assume all lands are urban'
       LANDUSE_frac_urban(:,:) = 1.0_RP
       call LANDUSE_calc_fact
    elseif( LANDUSE_MosaicWorld ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Assume all grids have ocean, land, and urban'
       LANDUSE_frac_land (:,:) = 0.5_RP
       LANDUSE_frac_urban(:,:) = 0.5_RP
       call LANDUSE_calc_fact
    else
       ! read from file
       call LANDUSE_read
    endif

    return
  end subroutine LANDUSE_setup

  !-----------------------------------------------------------------------------
  subroutine LANDUSE_calc_fact
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ calculate landuse factor'

    ! tentative treatment for lake fraction
    LANDUSE_frac_land (:,:) = LANDUSE_frac_land(:,:) * ( 1.0_RP - LANDUSE_frac_lake(:,:) )

    ! make factors
    LANDUSE_fact_ocean(:,:) = ( 1.0_RP - LANDUSE_frac_land(:,:) )
    LANDUSE_fact_land (:,:) = (          LANDUSE_frac_land(:,:) ) * ( 1.0_RP - LANDUSE_frac_urban(:,:) )
    LANDUSE_fact_urban(:,:) = (          LANDUSE_frac_land(:,:) ) * (          LANDUSE_frac_urban(:,:) )

    return
  end subroutine LANDUSE_calc_fact

  !-----------------------------------------------------------------------------
  !> Read landuse data
  subroutine LANDUSE_read
    use scale_fileio, only: &
       FILEIO_open, &
       FILEIO_read, &
       FILEIO_flush, &
       FILEIO_check_coordinates, &
       FILEIO_close
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: temp(IA,JA)

    character(len=H_SHORT) :: varname
    integer :: fid
    integer :: p, tag
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input landuse file ***'

    if ( LANDUSE_IN_BASENAME /= '' ) then

       call FILEIO_open( fid, LANDUSE_IN_BASENAME )

       call FILEIO_read( LANDUSE_frac_land(:,:),         & ! [OUT]
                         fid, 'FRAC_LAND',  'XY', step=1 ) ! [IN]
       call FILEIO_read( LANDUSE_frac_lake(:,:),         & ! [OUT]
                         fid, 'FRAC_LAKE',  'XY', step=1 ) ! [IN]
       call FILEIO_read( LANDUSE_frac_urban(:,:),        & ! [OUT]
                         fid, 'FRAC_URBAN', 'XY', step=1 ) ! [IN]

       call FILEIO_read( LANDUSE_fact_ocean(:,:),            & ! [OUT]
                         fid, 'FRAC_OCEAN_abs', 'XY', step=1 ) ! [IN]
       call FILEIO_read( LANDUSE_fact_land(:,:),             & ! [OUT]
                         fid, 'FRAC_LAND_abs',  'XY', step=1 ) ! [IN]
       call FILEIO_read( LANDUSE_fact_urban(:,:),            & ! [OUT]
                         fid, 'FRAC_URBAN_abs', 'XY', step=1 ) ! [IN]

       call FILEIO_flush( fid )

       if ( LANDUSE_IN_CHECK_COORDINATES ) then
          call FILEIO_check_coordinates( fid )
       end if

       call COMM_vars8( LANDUSE_frac_land (:,:), 1 )
       call COMM_vars8( LANDUSE_frac_lake (:,:), 2 )
       call COMM_vars8( LANDUSE_frac_urban(:,:), 3 )
       call COMM_vars8( LANDUSE_fact_ocean(:,:), 4 )
       call COMM_vars8( LANDUSE_fact_land (:,:), 5 )
       call COMM_vars8( LANDUSE_fact_urban(:,:), 6 )
       call COMM_wait ( LANDUSE_frac_land (:,:), 1 )
       call COMM_wait ( LANDUSE_frac_lake (:,:), 2 )
       call COMM_wait ( LANDUSE_frac_urban(:,:), 3 )
       call COMM_wait ( LANDUSE_fact_ocean(:,:), 4 )
       call COMM_wait ( LANDUSE_fact_land (:,:), 5 )
       call COMM_wait ( LANDUSE_fact_urban(:,:), 6 )

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p

          call FILEIO_read( LANDUSE_frac_PFT(:,:,p),   & ! [OUT]
                            fid, varname, 'XY', step=1 ) ! [IN]
          call FILEIO_flush( fid )

          tag = 7 + (p-1)*LANDUSE_PFT_mosaic
          call COMM_vars8( LANDUSE_frac_PFT (:,:,p), tag )
          call COMM_wait ( LANDUSE_frac_PFT (:,:,p), tag )

          write(varname,'(A9,I1.1)') 'INDEX_PFT', p

          call FILEIO_read( temp(:,:),                 & ! [OUT]
                            fid, varname, 'XY', step=1 ) ! [IN]
          call FILEIO_flush( fid )

          tag = 8 + (p-1)*LANDUSE_PFT_mosaic
          call COMM_vars8( temp(:,:), tag )
          call COMM_wait ( temp(:,:), tag )

          LANDUSE_index_PFT(:,:,p) = int(temp(:,:))
       enddo

       call FILEIO_close( fid )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** landuse file is not specified.'
       if( IO_L ) write(IO_FID_LOG,*) '*** Assume all grids are ocean'
    endif

    return
  end subroutine LANDUSE_read

  !-----------------------------------------------------------------------------
  !> Write landuse data
  subroutine LANDUSE_write
    use scale_fileio, only: &
       FILEIO_create, &
       FILEIO_def_var, &
       FILEIO_enddef, &
       FILEIO_write_var, &
       FILEIO_close
    implicit none

    real(RP) :: temp(IA,JA)

    integer :: fid, vid(6+LANDUSE_PFT_mosaic*2)
    character(len=H_SHORT) :: varname
    integer :: p
    !---------------------------------------------------------------------------

    if ( LANDUSE_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output landuse file ***'

       call FILEIO_create( fid, LANDUSE_OUT_BASENAME, LANDUSE_OUT_TITLE, &
                           LANDUSE_OUT_DTYPE, nozcoord=.true. )

       call FILEIO_def_var( fid, vid(1), 'FRAC_LAND',  'LAND fraction',  '1', 'XY', LANDUSE_OUT_DTYPE )
       call FILEIO_def_var( fid, vid(2), 'FRAC_LAKE',  'LAKE fraction',  '1', 'XY', LANDUSE_OUT_DTYPE )
       call FILEIO_def_var( fid, vid(3), 'FRAC_URBAN', 'URBAN fraction', '1', 'XY', LANDUSE_OUT_DTYPE )
       call FILEIO_def_var( fid, vid(4), 'FRAC_OCEAN_abs', 'absolute OCEAN fraction', '1', 'XY', LANDUSE_OUT_DTYPE )
       call FILEIO_def_var( fid, vid(5), 'FRAC_LAND_abs',  'absolute LAND fraction',  '1', 'XY', LANDUSE_OUT_DTYPE )
       call FILEIO_def_var( fid, vid(6), 'FRAC_URBAN_abs', 'absolute URBAN fraction', '1', 'XY', LANDUSE_OUT_DTYPE )

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p
          call FILEIO_def_var( fid, vid(7+2*(p-1)), varname, 'PFT fraction', '1', 'XY', LANDUSE_OUT_DTYPE )
          write(varname,'(A9,I1.1)') 'INDEX_PFT', p
          call FILEIO_def_var( fid, vid(8+2*(p-1)), varname, 'PFT index',    '1', 'XY', LANDUSE_OUT_DTYPE )
       end do

       call FILEIO_enddef( fid )

       call FILEIO_write_var( fid, vid(1), LANDUSE_frac_land (:,:), 'FRAC_LAND',  'XY' )
       call FILEIO_write_var( fid, vid(2), LANDUSE_frac_lake (:,:), 'FRAC_LAKE',  'XY' )
       call FILEIO_write_var( fid, vid(3), LANDUSE_frac_urban(:,:), 'FRAC_URBAN', 'XY' )
       call FILEIO_write_var( fid, vid(4), LANDUSE_fact_ocean(:,:), 'FRAC_OCEAN_abs', 'XY' )
       call FILEIO_write_var( fid, vid(5), LANDUSE_fact_land (:,:), 'FRAC_LAND_abs',  'XY' )
       call FILEIO_write_var( fid, vid(6), LANDUSE_fact_urban(:,:), 'FRAC_URBAN_abs', 'XY' )

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p
          call FILEIO_write_var( fid, vid(7+2*(p-1)), LANDUSE_frac_PFT(:,:,p), varname, 'XY' )
          write(varname,'(A9,I1.1)') 'INDEX_PFT', p
          temp(:,:) = real(LANDUSE_index_PFT(:,:,p),kind=RP)
          call FILEIO_write_var( fid, vid(8+2*(p-1)), temp(:,:),               varname, 'XY' )
       end do

       call FILEIO_close( fid )

    end if

    return
  end subroutine LANDUSE_write

end module scale_landuse
