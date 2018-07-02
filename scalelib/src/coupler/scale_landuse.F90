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
#include "scalelib.h"
module scale_landuse
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LANDUSE_setup
  public :: LANDUSE_calc_fact
  public :: LANDUSE_fillhalo
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
  integer,  public              :: LANDUSE_PFT_nmax   = 17  !< number of plant functional type(PFT)

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
  logical,                private :: LANDUSE_IN_AGGREGATE                       !< switch to use aggregated file
  logical,                private :: LANDUSE_IN_CHECK_COORDINATES = .true.      !< switch for check of coordinates
  character(len=H_LONG),  private :: LANDUSE_OUT_BASENAME = ''                  !< basename of the output file
  logical,                private :: LANDUSE_OUT_AGGREGATE                      !< switch to use aggregated file
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
    use scale_file, only: &
       FILE_AGGREGATE
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_LANDUSE / &
       LANDUSE_IN_BASENAME,          &
       LANDUSE_IN_AGGREGATE,         &
       LANDUSE_IN_CHECK_COORDINATES, &
       LANDUSE_OUT_BASENAME,         &
       LANDUSE_OUT_AGGREGATE,        &
       LANDUSE_OUT_DTYPE,            &
       LANDUSE_PFT_mosaic,           &
       LANDUSE_PFT_nmax,             &
       LANDUSE_AllOcean,             &
       LANDUSE_AllLand,              &
       LANDUSE_AllUrban,             &
       LANDUSE_MosaicWorld

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LANDUSE_setup",*) 'Setup'

    LANDUSE_IN_AGGREGATE  = FILE_AGGREGATE
    LANDUSE_OUT_AGGREGATE = FILE_AGGREGATE

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LANDUSE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LANDUSE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LANDUSE_setup",*) 'Not appropriate names in namelist PARAM_LANDUSE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_LANDUSE)

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
       LOG_INFO("LANDUSE_setup",*) 'Assume all grids are ocean'
       call LANDUSE_calc_fact
    elseif( LANDUSE_AllLand ) then
       LOG_INFO("LANDUSE_setup",*) 'Assume all grids are land'
       LOG_INFO("LANDUSE_setup",*) 'Assume land PFT is 1 (bare ground)'
       LANDUSE_frac_land (:,:) = 1.0_RP
       call LANDUSE_calc_fact
    elseif( LANDUSE_AllUrban ) then
       LOG_INFO("LANDUSE_setup",*) 'Assume all grids are land'
       LOG_INFO("LANDUSE_setup",*) 'Assume land PFT is 1 (bare ground)'
       LANDUSE_frac_land (:,:) = 1.0_RP
       LOG_INFO("LANDUSE_setup",*) 'Assume all lands are urban'
       LANDUSE_frac_urban(:,:) = 1.0_RP
       call LANDUSE_calc_fact
    elseif( LANDUSE_MosaicWorld ) then
       LOG_INFO("LANDUSE_setup",*) 'Assume all grids have ocean, land, and urban'
       LOG_INFO("LANDUSE_setup",*) 'Assume land PFT is 1 (bare ground)'
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

    LOG_NEWLINE
    LOG_INFO("LANDUSE_calc_fact",*) 'calculate landuse factor'

    ! tentative treatment: The area of the lake is treated as the ocean
    LANDUSE_frac_land (:,:) = LANDUSE_frac_land(:,:) * ( 1.0_RP - LANDUSE_frac_lake(:,:) )

    ! make factors
    LANDUSE_fact_ocean(:,:) = ( 1.0_RP - LANDUSE_frac_land(:,:) )
    LANDUSE_fact_land (:,:) = (          LANDUSE_frac_land(:,:) ) * ( 1.0_RP - LANDUSE_frac_urban(:,:) )
    LANDUSE_fact_urban(:,:) = (          LANDUSE_frac_land(:,:) ) * (          LANDUSE_frac_urban(:,:) )

    return
  end subroutine LANDUSE_calc_fact

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine LANDUSE_fillhalo( FILL_BND )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    logical, intent(in), optional :: FILL_BND

    real(RP) :: temp(IA,JA)

    logical  :: FILL_BND_
    integer  :: p
    !---------------------------------------------------------------------------

    FILL_BND_ = .true.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    call COMM_vars8( LANDUSE_frac_land (:,:), 1 )
    call COMM_vars8( LANDUSE_frac_lake (:,:), 2 )
    call COMM_vars8( LANDUSE_frac_urban(:,:), 3 )
    call COMM_vars8( LANDUSE_fact_ocean(:,:), 4 )
    call COMM_vars8( LANDUSE_fact_land (:,:), 5 )
    call COMM_vars8( LANDUSE_fact_urban(:,:), 6 )

    call COMM_wait ( LANDUSE_frac_land (:,:), 1, FILL_BND_ )
    call COMM_wait ( LANDUSE_frac_lake (:,:), 2, FILL_BND_ )
    call COMM_wait ( LANDUSE_frac_urban(:,:), 3, FILL_BND_ )
    call COMM_wait ( LANDUSE_fact_ocean(:,:), 4, FILL_BND_ )
    call COMM_wait ( LANDUSE_fact_land (:,:), 5, FILL_BND_ )
    call COMM_wait ( LANDUSE_fact_urban(:,:), 6, FILL_BND_ )

    do p = 1, LANDUSE_PFT_mosaic
       temp(:,:) = real(LANDUSE_index_PFT(:,:,p),kind=RP)

       call COMM_vars8( LANDUSE_frac_PFT(:,:,p), 7+2*(p-1) )
       call COMM_vars8( temp            (:,:)  , 8+2*(p-1) )

       call COMM_wait ( LANDUSE_frac_PFT(:,:,p), 7+2*(p-1), FILL_BND_ )
       call COMM_wait ( temp            (:,:)  , 8+2*(p-1), FILL_BND_ )

       LANDUSE_index_PFT(:,:,p) = int(temp(:,:)+1.E-3_RP,kind=4)
    enddo

    return
  end subroutine LANDUSE_fillhalo

  !-----------------------------------------------------------------------------
  !> Read landuse data
  subroutine LANDUSE_read
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush, &
       FILE_CARTESC_check_coordinates, &
       FILE_CARTESC_close
    implicit none

    real(RP) :: temp(IA,JA)

    character(len=H_SHORT) :: varname

    integer  :: fid
    integer  :: p, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LANDUSE_read",*) 'Input landuse file '

    if ( LANDUSE_IN_BASENAME /= '' ) then

       call FILE_CARTESC_open( LANDUSE_IN_BASENAME, fid, aggregate=LANDUSE_IN_AGGREGATE )

       call FILE_CARTESC_read( fid, 'FRAC_LAND',  'XY', LANDUSE_frac_land(:,:) )
       call FILE_CARTESC_read( fid, 'FRAC_LAKE',  'XY', LANDUSE_frac_lake(:,:) )
       call FILE_CARTESC_read( fid, 'FRAC_URBAN', 'XY', LANDUSE_frac_urban(:,:) )

       call FILE_CARTESC_read( fid, 'FRAC_OCEAN_abs', 'XY', LANDUSE_fact_ocean(:,:) )
       call FILE_CARTESC_read( fid, 'FRAC_LAND_abs',  'XY', LANDUSE_fact_land(:,:) )
       call FILE_CARTESC_read( fid, 'FRAC_URBAN_abs', 'XY', LANDUSE_fact_urban(:,:) )

       call FILE_CARTESC_flush( fid )

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p

          call FILE_CARTESC_read( fid, varname, 'XY', LANDUSE_frac_PFT(:,:,p) )

          write(varname,'(A9,I1.1)') 'INDEX_PFT', p

          call FILE_CARTESC_read( fid, varname, 'XY', temp(:,:) )
          call FILE_CARTESC_flush( fid ) ! for non-blocking I/O

          do j = JS, JE
          do i = IS, IE
             LANDUSE_index_PFT(i,j,p) = int(temp(i,j)+1.E-3_RP,kind=4)
          enddo
          enddo
       enddo

       if ( LANDUSE_IN_CHECK_COORDINATES ) then
          call FILE_CARTESC_check_coordinates( fid )
       end if

       call FILE_CARTESC_close( fid )

       call LANDUSE_fillhalo( FILL_BND=.false. )

    else
       LOG_INFO_CONT(*) 'landuse file is not specified.'
       LOG_INFO_CONT(*) 'Assume all grids are ocean'
    endif

    return
  end subroutine LANDUSE_read

  !-----------------------------------------------------------------------------
  !> Write landuse data
  subroutine LANDUSE_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_create, &
       FILE_CARTESC_def_var, &
       FILE_CARTESC_enddef, &
       FILE_CARTESC_write_var, &
       FILE_CARTESC_close
    implicit none

    real(RP) :: temp(IA,JA)

    integer                :: vid(6+LANDUSE_PFT_mosaic*2)
    character(len=H_SHORT) :: varname

    integer :: fid
    integer :: p
    !---------------------------------------------------------------------------

    if ( LANDUSE_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("LANDUSE_write",*) 'Output landuse file '

       call LANDUSE_fillhalo( FILL_BND=.false. )

       call FILE_CARTESC_create( &
            LANDUSE_OUT_BASENAME, LANDUSE_OUT_TITLE, LANDUSE_OUT_DTYPE, & ! [IN]
            fid,                                                        & ! [OUT]
            haszcoord=.false., aggregate=LANDUSE_OUT_AGGREGATE          ) ! [IN]

       call FILE_CARTESC_def_var( fid, 'FRAC_LAND'     , 'LAND fraction'          , '1', 'XY', LANDUSE_OUT_DTYPE, vid(1), standard_name="land_area_fraction" )
       call FILE_CARTESC_def_var( fid, 'FRAC_LAKE'     , 'LAKE fraction'          , '1', 'XY', LANDUSE_OUT_DTYPE, vid(2) )
       call FILE_CARTESC_def_var( fid, 'FRAC_URBAN'    , 'URBAN fraction'         , '1', 'XY', LANDUSE_OUT_DTYPE, vid(3) )
       call FILE_CARTESC_def_var( fid, 'FRAC_OCEAN_abs', 'absolute OCEAN fraction', '1', 'XY', LANDUSE_OUT_DTYPE, vid(4) )
       call FILE_CARTESC_def_var( fid, 'FRAC_LAND_abs' , 'absolute LAND fraction' , '1', 'XY', LANDUSE_OUT_DTYPE, vid(5) )
       call FILE_CARTESC_def_var( fid, 'FRAC_URBAN_abs', 'absolute URBAN fraction', '1', 'XY', LANDUSE_OUT_DTYPE, vid(6) )

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p
          call FILE_CARTESC_def_var( fid, varname, 'PFT fraction', '1', 'XY', LANDUSE_OUT_DTYPE, vid(7+2*(p-1)) )
          write(varname,'(A9,I1.1)') 'INDEX_PFT', p
          call FILE_CARTESC_def_var( fid, varname, 'PFT index',    '1', 'XY', LANDUSE_OUT_DTYPE, vid(8+2*(p-1)) )
       end do

       call FILE_CARTESC_enddef( fid )

       call FILE_CARTESC_write_var( fid, vid(1), LANDUSE_frac_land (:,:), 'FRAC_LAND'     , 'XY' )
       call FILE_CARTESC_write_var( fid, vid(2), LANDUSE_frac_lake (:,:), 'FRAC_LAKE'     , 'XY' )
       call FILE_CARTESC_write_var( fid, vid(3), LANDUSE_frac_urban(:,:), 'FRAC_URBAN'    , 'XY' )
       call FILE_CARTESC_write_var( fid, vid(4), LANDUSE_fact_ocean(:,:), 'FRAC_OCEAN_abs', 'XY' )
       call FILE_CARTESC_write_var( fid, vid(5), LANDUSE_fact_land (:,:), 'FRAC_LAND_abs' , 'XY' )
       call FILE_CARTESC_write_var( fid, vid(6), LANDUSE_fact_urban(:,:), 'FRAC_URBAN_abs', 'XY' )

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p
          call FILE_CARTESC_write_var( fid, vid(7+2*(p-1)), LANDUSE_frac_PFT(:,:,p), varname, 'XY' )
          write(varname,'(A9,I1.1)') 'INDEX_PFT', p
          temp(:,:) = real(LANDUSE_index_PFT(:,:,p),kind=RP)
          call FILE_CARTESC_write_var( fid, vid(8+2*(p-1)), temp(:,:),               varname, 'XY' )
       end do

       call FILE_CARTESC_close( fid )

    end if

    return
  end subroutine LANDUSE_write

end module scale_landuse
