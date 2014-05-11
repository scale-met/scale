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
  public :: LANDUSE_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: LANDUSE_frac_land (:,:) !< land  fraction
  real(RP), public, allocatable :: LANDUSE_frac_lake (:,:) !< lake  fraction
  real(RP), public, allocatable :: LANDUSE_frac_urban(:,:) !< urban fraction

  integer,  public, save        :: LANDUSE_PFT_mosaic = 2   !< number of PFT mosaic
  integer,  public, save        :: LANDUSE_PFT_nmax   = 4   !< number of plant functional type(PFT)

  real(RP), public, allocatable :: LANDUSE_frac_PFT (:,:,:) !< fraction of PFT for each mosaic
  integer,  public, allocatable :: LANDUSE_index_PFT(:,:,:) !< index    of PFT for each mosaic

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LANDUSE_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private :: LANDUSE_IN_BASENAME  = ''                  !< basename of the input  file
  character(len=H_LONG), private :: LANDUSE_OUT_BASENAME = ''                  !< basename of the output file
  character(len=H_MID),  private :: LANDUSE_OUT_TITLE    = 'SCALE-LES LANDUSE' !< title    of the output file
  character(len=H_MID),  private :: LANDUSE_OUT_DTYPE    = 'DEFAULT'           !< REAL4 or REAL8
  logical,               private :: LANDUSE_AllLand      = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LANDUSE_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_LANDUSE / &
       LANDUSE_IN_BASENAME,  &
       LANDUSE_OUT_BASENAME, &
       LANDUSE_OUT_DTYPE,    &
       LANDUSE_PFT_mosaic,   &
       LANDUSE_PFT_nmax,     &
       LANDUSE_AllLand

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[LANDUSE]/Categ[GRID]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LANDUSE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LANDUSE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_LANDUSE)

    allocate( LANDUSE_frac_land (IA,JA) )
    allocate( LANDUSE_frac_lake (IA,JA) )
    allocate( LANDUSE_frac_urban(IA,JA) )
    LANDUSE_frac_land (:,:) = 0.0_RP
    LANDUSE_frac_lake (:,:) = 0.0_RP
    LANDUSE_frac_urban(:,:) = 0.0_RP

    allocate( LANDUSE_index_PFT(IA,JA,LANDUSE_PFT_mosaic) )
    allocate( LANDUSE_frac_PFT (IA,JA,LANDUSE_PFT_mosaic) )
    LANDUSE_frac_PFT (:,:,:) = 0.0_RP
    LANDUSE_frac_PFT (:,:,1) = 1.0_RP
    LANDUSE_index_PFT(:,:,:) = 1

    ! read from file
    call LANDUSE_read

    if ( LANDUSE_AllLand ) then
       LANDUSE_frac_land (:,:) = 1.0_RP
    endif

    return
  end subroutine LANDUSE_setup

  !-----------------------------------------------------------------------------
  !> Read landuse data
  subroutine LANDUSE_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: temp(IA,JA)

    character(len=H_SHORT) :: varname
    integer :: p
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input landuse file ***'

    if ( LANDUSE_IN_BASENAME /= '' ) then

       call FILEIO_read( LANDUSE_frac_land(:,:),                        & ! [OUT]
                         LANDUSE_IN_BASENAME, 'FRAC_LAND', 'XY', step=1 ) ! [IN]
       call FILEIO_read( LANDUSE_frac_lake(:,:),                        & ! [OUT]
                         LANDUSE_IN_BASENAME, 'FRAC_LAKE', 'XY', step=1 ) ! [IN]
       call FILEIO_read( LANDUSE_frac_urban(:,:),                        & ! [OUT]
                         LANDUSE_IN_BASENAME, 'FRAC_URBAN', 'XY', step=1 ) ! [IN]

       ! fill IHALO & JHALO
       call COMM_vars8( LANDUSE_frac_land (:,:), 1 )
       call COMM_vars8( LANDUSE_frac_lake (:,:), 2 )
       call COMM_vars8( LANDUSE_frac_urban(:,:), 3 )
       call COMM_wait ( LANDUSE_frac_land (:,:), 1 )
       call COMM_wait ( LANDUSE_frac_lake (:,:), 2 )
       call COMM_wait ( LANDUSE_frac_urban(:,:), 3 )

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p

          call FILEIO_read( LANDUSE_frac_PFT(:,:,p),                   & ! [OUT]
                            LANDUSE_IN_BASENAME, varname, 'XY', step=1 ) ! [IN]

          call COMM_vars8( LANDUSE_frac_PFT (:,:,p), 4 )
          call COMM_wait ( LANDUSE_frac_PFT (:,:,p), 4 )

          write(varname,'(A9,I1.1)') 'INDEX_PFT', p

          call FILEIO_read( temp(:,:),                                 & ! [OUT]
                            LANDUSE_IN_BASENAME, varname, 'XY', step=1 ) ! [IN]

          call COMM_vars8( temp(:,:), 5 )
          call COMM_wait ( temp(:,:), 5 )

          LANDUSE_index_PFT(:,:,p) = int(temp(:,:))
       enddo

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
       FILEIO_write
    implicit none

    real(RP) :: temp(IA,JA)

    character(len=H_SHORT) :: varname
    integer :: p
    !---------------------------------------------------------------------------

    if ( LANDUSE_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output landuse file ***'

       call FILEIO_write( LANDUSE_frac_land (:,:), LANDUSE_OUT_BASENAME, LANDUSE_OUT_TITLE, & ! [IN]
                          'FRAC_LAND',  'LAND fraction',  '0-1', 'XY',   LANDUSE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( LANDUSE_frac_lake (:,:), LANDUSE_OUT_BASENAME, LANDUSE_OUT_TITLE, & ! [IN]
                          'FRAC_LAKE',  'LAKE fraction',  '0-1', 'XY',   LANDUSE_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( LANDUSE_frac_urban(:,:), LANDUSE_OUT_BASENAME, LANDUSE_OUT_TITLE, & ! [IN]
                          'FRAC_URBAN', 'URBAN fraction', '0-1', 'XY',   LANDUSE_OUT_DTYPE  ) ! [IN]

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p

          call FILEIO_write( LANDUSE_frac_PFT(:,:,p), LANDUSE_OUT_BASENAME, LANDUSE_OUT_TITLE, & ! [IN]
                             varname, 'PFT fraction', '0-1', 'XY',          LANDUSE_OUT_DTYPE  ) ! [IN]

          write(varname,'(A9,I1.1)') 'INDEX_PFT', p
          temp(:,:) = real(LANDUSE_index_PFT(:,:,p),kind=RP)

          call FILEIO_write( temp(:,:), LANDUSE_OUT_BASENAME,   LANDUSE_OUT_TITLE, & ! [IN]
                             varname, 'PFT index', 'IDX', 'XY', LANDUSE_OUT_DTYPE  ) ! [IN]
       enddo

    endif

    return
  end subroutine LANDUSE_write

end module scale_landuse
