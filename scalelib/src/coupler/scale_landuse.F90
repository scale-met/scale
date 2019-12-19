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
  real(RP), public, allocatable :: LANDUSE_fact_lake (:,:) !< lake  factor

  logical,  public, allocatable :: LANDUSE_exists_ocean(:,:) !< ocean calculation flag
  logical,  public, allocatable :: LANDUSE_exists_land (:,:) !< land  calculation flag
  logical,  public, allocatable :: LANDUSE_exists_urban(:,:) !< urban calculation flag
  logical,  public, allocatable :: LANDUSE_exists_lake (:,:) !< lake  calculation flag

  real(RP), public, allocatable :: LANDUSE_frac_land (:,:) !< land  fraction
  real(RP), public, allocatable :: LANDUSE_frac_urban(:,:) !< urban fraction
  real(RP), public, allocatable :: LANDUSE_frac_lake (:,:) !< lake  fraction

  integer,  public, parameter   :: LANDUSE_index_OCEAN  =  0 !< ocean index
  integer,  public, parameter   :: LANDUSE_index_URBAN  = -1 !< urban index
  integer,  public, parameter   :: LANDUSE_index_LAKE   = -2 !< lake  index

  integer,  public, parameter   :: LANDUSE_PFT_nmin   = -2  !< minimum number of PFT type
  integer,  public              :: LANDUSE_PFT_nmax   = 17  !< number of plant functional type(PFT)
  integer,  public              :: LANDUSE_PFT_mosaic = 2   !< number of PFT mosaic

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
  logical,                private :: LANDUSE_AllLake      = .false.
  logical,                private :: LANDUSE_Ignore_Lake  = .false.
  logical,                private :: LANDUSE_MosaicWorld  = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LANDUSE_setup( &
       OCEAN_do, &
       URBAN_do, &
       LAKE_do   )
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_AGGREGATE
    implicit none

    logical, intent(in) :: OCEAN_do
    logical, intent(in) :: URBAN_do
    logical, intent(in) :: LAKE_do

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
       LANDUSE_AllLake,              &
       LANDUSE_Ignore_Lake,          &
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
    allocate( LANDUSE_frac_urban(IA,JA) )
    allocate( LANDUSE_frac_lake (IA,JA) )
    LANDUSE_frac_land (:,:) = 0.0_RP
    LANDUSE_frac_urban(:,:) = 0.0_RP
    LANDUSE_frac_lake (:,:) = 0.0_RP

    allocate( LANDUSE_index_PFT(IA,JA,LANDUSE_PFT_mosaic) )
    allocate( LANDUSE_frac_PFT (IA,JA,LANDUSE_PFT_mosaic) )
    LANDUSE_frac_PFT (:,:,:) = 0.0_RP
    LANDUSE_frac_PFT (:,:,1) = 1.0_RP ! tentative, mosaic is off
    LANDUSE_index_PFT(:,:,:) = 1      ! default

    allocate( LANDUSE_fact_ocean(IA,JA) )
    allocate( LANDUSE_fact_land (IA,JA) )
    allocate( LANDUSE_fact_urban(IA,JA) )
    allocate( LANDUSE_fact_lake (IA,JA) )
    LANDUSE_fact_ocean(:,:) = 0.0_RP
    LANDUSE_fact_land (:,:) = 0.0_RP
    LANDUSE_fact_urban(:,:) = 0.0_RP
    LANDUSE_fact_lake (:,:) = 0.0_RP

    allocate( LANDUSE_exists_ocean(IA,JA) )
    allocate( LANDUSE_exists_land (IA,JA) )
    allocate( LANDUSE_exists_urban(IA,JA) )
    allocate( LANDUSE_exists_lake (IA,JA) )
    LANDUSE_exists_ocean(:,:) = .false.
    LANDUSE_exists_land (:,:) = .false.
    LANDUSE_exists_urban(:,:) = .false.
    LANDUSE_exists_lake (:,:) = .false.



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
       LOG_INFO("LANDUSE_setup",*) 'Assume all lands are urban'
       LOG_INFO("LANDUSE_setup",*) 'Assume land PFT is urban: ', LANDUSE_index_URBAN
       LANDUSE_frac_land (:,:)   = 1.0_RP
       LANDUSE_frac_urban(:,:)   = 1.0_RP
       LANDUSE_index_PFT (:,:,:) = LANDUSE_index_URBAN

       call LANDUSE_calc_fact

    elseif( LANDUSE_AllLake ) then

       LOG_INFO("LANDUSE_setup",*) 'Assume all grids are land'
       LOG_INFO("LANDUSE_setup",*) 'Assume all lands are lake'
       LOG_INFO("LANDUSE_setup",*) 'Assume land PFT is lake: ', LANDUSE_index_LAKE
       LANDUSE_frac_land (:,:)   = 1.0_RP
       LANDUSE_frac_lake (:,:)   = 1.0_RP
       LANDUSE_index_PFT (:,:,:) = LANDUSE_index_LAKE

       call LANDUSE_calc_fact

    elseif( LANDUSE_MosaicWorld ) then

       LOG_INFO("LANDUSE_setup",*) 'Assume all grids have ocean, land, and urban'
!       LOG_INFO("LANDUSE_setup",*) 'Assume all grids have ocean, land, urban, and lake'
       LOG_INFO("LANDUSE_setup",*) 'Assume land PFT is 1 (bare ground)'
       LANDUSE_frac_land (:,:) = 0.5_RP
       LANDUSE_frac_urban(:,:) = 0.5_RP
!       LANDUSE_frac_lake (:,:) = 0.25_RP

       call LANDUSE_calc_fact

    else ! default: read from file

       call LANDUSE_read( OCEAN_do, URBAN_do, LAKE_do )
       call LANDUSE_calc_fact

    endif

    return
  end subroutine LANDUSE_setup

  !-----------------------------------------------------------------------------
  subroutine LANDUSE_calc_fact
    implicit none

    real(RP) :: fact_soil
    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LANDUSE_calc_fact",*) 'calculate landuse factor'

    ! make factors
    !$omp parallel do private(fact_soil)
    do j = 1, JA
    do i = 1, IA
       LANDUSE_fact_ocean(i,j) = max( 1.0_RP - LANDUSE_frac_land(i,j), 0.0_RP )
       LANDUSE_fact_lake (i,j) = LANDUSE_frac_land(i,j) * LANDUSE_frac_lake(i,j)
       fact_soil = max( LANDUSE_frac_land(i,j) - LANDUSE_fact_lake(i,j), 0.0_RP )
       LANDUSE_fact_urban(i,j) = fact_soil * LANDUSE_frac_urban(i,j)
       LANDUSE_fact_land (i,j) = max( fact_soil - LANDUSE_fact_urban(i,j), 0.0_RP )

       if( LANDUSE_fact_ocean(i,j) > 0.0_RP ) LANDUSE_exists_ocean(i,j) = .true.
       if( LANDUSE_fact_land (i,j) > 0.0_RP ) LANDUSE_exists_land (i,j) = .true.
       if( LANDUSE_fact_urban(i,j) > 0.0_RP ) LANDUSE_exists_urban(i,j) = .true.
       if( LANDUSE_fact_lake (i,j) > 0.0_RP ) LANDUSE_exists_lake (i,j) = .true.
    enddo
    enddo

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
    call COMM_vars8( LANDUSE_fact_lake (:,:), 7 )

    call COMM_wait ( LANDUSE_frac_land (:,:), 1, FILL_BND_ )
    call COMM_wait ( LANDUSE_frac_lake (:,:), 2, FILL_BND_ )
    call COMM_wait ( LANDUSE_frac_urban(:,:), 3, FILL_BND_ )
    call COMM_wait ( LANDUSE_fact_ocean(:,:), 4, FILL_BND_ )
    call COMM_wait ( LANDUSE_fact_land (:,:), 5, FILL_BND_ )
    call COMM_wait ( LANDUSE_fact_urban(:,:), 6, FILL_BND_ )
    call COMM_wait ( LANDUSE_fact_lake (:,:), 7, FILL_BND_ )

    do p = 1, LANDUSE_PFT_mosaic
       temp(:,:) = real(LANDUSE_index_PFT(:,:,p),kind=RP)

       call COMM_vars8( LANDUSE_frac_PFT(:,:,p), 8+2*(p-1) )
       call COMM_vars8( temp            (:,:)  , 9+2*(p-1) )

       call COMM_wait ( LANDUSE_frac_PFT(:,:,p), 8+2*(p-1), FILL_BND_ )
       call COMM_wait ( temp            (:,:)  , 9+2*(p-1), FILL_BND_ )

       LANDUSE_index_PFT(:,:,p) = aint(temp(:,:),kind=4)
    enddo

    return
  end subroutine LANDUSE_fillhalo

  !-----------------------------------------------------------------------------
  !> Read landuse data
  subroutine LANDUSE_read( &
       OCEAN_do, &
       URBAN_do, &
       LAKE_do   )
    use scale_file_cartesC, only: &
       FILE_CARTESC_open,              &
       FILE_CARTESC_read,              &
       FILE_CARTESC_flush,             &
       FILE_CARTESC_check_coordinates, &
       FILE_CARTESC_close
    use scale_prc, only: &
       PRC_abort
    implicit none

    logical, intent(in) :: OCEAN_do
    logical, intent(in) :: URBAN_do
    logical, intent(in) :: LAKE_do

    real(RP) :: temp(IA,JA)

    character(len=H_SHORT) :: varname

    integer  :: fid
    integer  :: p, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LANDUSE_read",*) 'Input landuse file '

    if ( LANDUSE_IN_BASENAME /= '' ) then

       call FILE_CARTESC_open( LANDUSE_IN_BASENAME, fid, aggregate=LANDUSE_IN_AGGREGATE )

       call FILE_CARTESC_read( fid, 'FRAC_LAND',  'XY', LANDUSE_frac_land (:,:) )
       call FILE_CARTESC_read( fid, 'FRAC_LAKE',  'XY', LANDUSE_frac_lake (:,:) )
       call FILE_CARTESC_read( fid, 'FRAC_URBAN', 'XY', LANDUSE_frac_urban(:,:) )

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

!!$       if ( .not. OCEAN_do ) then
!!$          !$omp parallel do private(frac_ocean)
!!$          do j = 1, JA
!!$          do i = 1, IA
!!$             do p = LANDUSE_PFT_mosaic, 1, -1
!!$                frac_ocean = 1.0_RP - LANDUSE_frac_land(i,j)
!!$                if ( frac_ocean > LANDUSE_frac_PFT(i,j,p) ) then
!!$                   if ( p < LANDUSE_PFT_mosaic ) then
!!$                      LANDUSE_index_PFT(i,j,p+1) = LANDUSE_index_PFT(i,j,p)
!!$                      LANDUSE_frac_PFT (i,j,p+1) = LANDUSE_frac_PFT (i,j,p)
!!$                   end if
!!$                   LANDUSE_index_PFT(i,j,p) = LANDUSE_index_OCEAN
!!$                   LANDUSE_frac_PFT (i,j,p) = frac_ocean
!!$                else
!!$                   exit
!!$                end if
!!$             end do
!!$             LANDUSE_frac_PFT(i,j,:) = LANDUSE_frac_PFT(i,j,:) / sum( LANDUSE_frac_PFT(i,j,:) )
!!$          end do
!!$          end do
!!$          LANDUSE_frac_land(:,:) = 1.0_RP
!!$       end if

       if ( .not. URBAN_do ) then
          !$omp parallel do
          do j = 1, JA
          do i = 1, IA
             do p = LANDUSE_PFT_mosaic, 1, -1
                if ( LANDUSE_frac_urban(i,j) > LANDUSE_frac_PFT(i,j,p) ) then
                   if ( p < LANDUSE_PFT_mosaic ) then
                      LANDUSE_index_PFT(i,j,p+1) = LANDUSE_index_PFT(i,j,p)
                      LANDUSE_frac_PFT (i,j,p+1) = LANDUSE_frac_PFT (i,j,p)
                   end if
                   LANDUSE_index_PFT(i,j,p) = LANDUSE_index_URBAN
                   LANDUSE_frac_PFT (i,j,p) = LANDUSE_frac_urban(i,j)
                end if
             end do
             LANDUSE_frac_PFT(i,j,:) = LANDUSE_frac_PFT(i,j,:) / sum( LANDUSE_frac_PFT(i,j,:) )
          end do
          end do
          LANDUSE_frac_urban(:,:) = 0.0_RP
       end if

       if ( .not. LAKE_do ) then
          if ( LANDUSE_Ignore_Lake .or. (.not. OCEAN_do) ) then
             !$omp parallel do
             do j = 1, JA
             do i = 1, IA
                do p = LANDUSE_PFT_mosaic, 1, -1
                   if ( LANDUSE_frac_lake(i,j) > LANDUSE_frac_PFT(i,j,p) ) then
                      if ( p < LANDUSE_PFT_mosaic ) then
                         LANDUSE_index_PFT(i,j,p+1) = LANDUSE_index_PFT(i,j,p)
                         LANDUSE_frac_PFT (i,j,p+1) = LANDUSE_frac_PFT (i,j,p)
                      end if
                      LANDUSE_index_PFT(i,j,p) = LANDUSE_index_LAKE
                      LANDUSE_frac_PFT (i,j,p) = LANDUSE_frac_lake(i,j)
                   else
                      exit
                   end if
                end do
                LANDUSE_frac_PFT(i,j,:) = LANDUSE_frac_PFT(i,j,:) / sum( LANDUSE_frac_PFT(i,j,:) )
             end do
             end do
          else
             ! lake is assumed to be ocean
             !$omp parallel do
             do j = 1, JA
             do i = 1, IA
                if( LANDUSE_frac_land(i,j) == 0.0_RP .and. LANDUSE_frac_urban(i,j) > 0.0_RP)then
                   LOG_ERROR("LANDUSE_read",*) 'Not appropriate original landuse (w/ Lake). Bug!',i,j
                   call PRC_abort
                endif

                LANDUSE_frac_land(i,j)  = min( max( LANDUSE_frac_land(i,j) * (1.0_RP-LANDUSE_frac_lake(i,j)), 0.0_RP), 1.0_RP)

                if( LANDUSE_frac_land(i,j) == 0.0_RP .and. LANDUSE_frac_urban(i,j) > 0.0_RP)then
                   LOG_ERROR("LANDUSE_read",*) 'Not appropriate replaced landuse (w/o Lake). Bug!',i,j
                   call PRC_abort
                endif

             end do
             end do
          end if
          LANDUSE_frac_lake(:,:) = 0.0_RP
       end if

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
       FILE_CARTESC_create,    &
       FILE_CARTESC_def_var,   &
       FILE_CARTESC_enddef,    &
       FILE_CARTESC_write_var, &
       FILE_CARTESC_close
    implicit none

    real(RP) :: temp(IA,JA)

    integer                :: vid(7+LANDUSE_PFT_mosaic*2)
    character(len=H_SHORT) :: varname

    integer :: fid
    integer :: p
    !---------------------------------------------------------------------------

    if ( LANDUSE_OUT_BASENAME /= '' .and. LANDUSE_OUT_BASENAME /= LANDUSE_IN_BASENAME ) then

       LOG_NEWLINE
       LOG_INFO("LANDUSE_write",*) 'Output landuse file '

       call LANDUSE_fillhalo( FILL_BND=.false. )

       call FILE_CARTESC_create( LANDUSE_OUT_BASENAME, LANDUSE_OUT_TITLE, LANDUSE_OUT_DTYPE, & ! [IN]
                                 fid,                                                        & ! [OUT]
                                 haszcoord=.false., aggregate=LANDUSE_OUT_AGGREGATE          ) ! [IN]

       call FILE_CARTESC_def_var( fid, 'FRAC_LAND'     , 'LAND fraction'          , '1', 'XY', LANDUSE_OUT_DTYPE, vid(1), standard_name="land_area_fraction" )
       call FILE_CARTESC_def_var( fid, 'FRAC_LAKE'     , 'LAKE fraction'          , '1', 'XY', LANDUSE_OUT_DTYPE, vid(2) )
       call FILE_CARTESC_def_var( fid, 'FRAC_URBAN'    , 'URBAN fraction'         , '1', 'XY', LANDUSE_OUT_DTYPE, vid(3) )
       call FILE_CARTESC_def_var( fid, 'FRAC_OCEAN_abs', 'absolute OCEAN fraction', '1', 'XY', LANDUSE_OUT_DTYPE, vid(4) )
       call FILE_CARTESC_def_var( fid, 'FRAC_LAND_abs' , 'absolute LAND fraction' , '1', 'XY', LANDUSE_OUT_DTYPE, vid(5) )
       call FILE_CARTESC_def_var( fid, 'FRAC_URBAN_abs', 'absolute URBAN fraction', '1', 'XY', LANDUSE_OUT_DTYPE, vid(6) )
       call FILE_CARTESC_def_var( fid, 'FRAC_LAKE_abs' , 'absolute LAKE fraction' , '1', 'XY', LANDUSE_OUT_DTYPE, vid(7) )

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p
          call FILE_CARTESC_def_var( fid, varname, 'PFT fraction', '1', 'XY', LANDUSE_OUT_DTYPE, vid(8+2*(p-1)) )
          write(varname,'(A9,I1.1)') 'INDEX_PFT', p
          call FILE_CARTESC_def_var( fid, varname, 'PFT index',    '1', 'XY', LANDUSE_OUT_DTYPE, vid(9+2*(p-1)) )
       end do

       call FILE_CARTESC_enddef( fid )

       call FILE_CARTESC_write_var( fid, vid(1), LANDUSE_frac_land (:,:), 'FRAC_LAND'     , 'XY' )
       call FILE_CARTESC_write_var( fid, vid(2), LANDUSE_frac_lake (:,:), 'FRAC_LAKE'     , 'XY' )
       call FILE_CARTESC_write_var( fid, vid(3), LANDUSE_frac_urban(:,:), 'FRAC_URBAN'    , 'XY' )
       call FILE_CARTESC_write_var( fid, vid(4), LANDUSE_fact_ocean(:,:), 'FRAC_OCEAN_abs', 'XY' )
       call FILE_CARTESC_write_var( fid, vid(5), LANDUSE_fact_land (:,:), 'FRAC_LAND_abs' , 'XY' )
       call FILE_CARTESC_write_var( fid, vid(6), LANDUSE_fact_urban(:,:), 'FRAC_URBAN_abs', 'XY' )
       call FILE_CARTESC_write_var( fid, vid(7), LANDUSE_fact_lake (:,:), 'FRAC_LAKE_abs' , 'XY' )

       do p = 1, LANDUSE_PFT_mosaic
          write(varname,'(A8,I1.1)') 'FRAC_PFT', p
          call FILE_CARTESC_write_var( fid, vid(8+2*(p-1)), LANDUSE_frac_PFT(:,:,p), varname, 'XY' )
          write(varname,'(A9,I1.1)') 'INDEX_PFT', p
          temp(:,:) = real(LANDUSE_index_PFT(:,:,p),kind=RP)
          call FILE_CARTESC_write_var( fid, vid(9+2*(p-1)), temp(:,:),               varname, 'XY' )
       end do

       call FILE_CARTESC_close( fid )

    end if

    return
  end subroutine LANDUSE_write

end module scale_landuse
