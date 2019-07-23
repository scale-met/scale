!-------------------------------------------------------------------------------
!> Module SNO (RM) PLUGIN
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          SCALE NetCDF Operator (SNO)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_snoplugin_hgridope
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use mod_sno_h, only: &
     axisinfo, &
     iteminfo
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SNOPLGIN_hgridope_setup
  public :: SNOPLGIN_hgridope_setcoef
  public :: SNOPLGIN_hgridope_alloc
  public :: SNOPLGIN_hgridope_dealloc
  public :: SNOPLGIN_hgridope_llinterp

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
  character(len=H_SHORT), private              :: SNOPLGIN_hgridope_type        = 'OFF'   ! type of average
                                                                                ! 'OFF'    : disable
                                                                                ! 'LATLON' : remap to latitude-longitude grid
  real(RP),               private              :: SNOPLGIN_hgridope_lat_start   = -1.0_RP
  real(RP),               private              :: SNOPLGIN_hgridope_lat_end     = -1.0_RP
  real(RP),               private              :: SNOPLGIN_hgridope_dlat        = -1.0_RP
  real(RP),               private              :: SNOPLGIN_hgridope_lon_start   = -1.0_RP
  real(RP),               private              :: SNOPLGIN_hgridope_lon_end     = -1.0_RP
  real(RP),               private              :: SNOPLGIN_hgridope_dlon        = -1.0_RP

  integer,                private              :: SNOPLGIN_hgridope_nintrp      = 5       ! number of interpolation point
  integer,                private              :: SNOPLGIN_hgridope_weight      = 2       ! weighting factor for interpolation

  integer,                private              :: naxis_ll
  type(axisinfo),         private              :: ainfo_ll(8)
  type(iteminfo),         private              :: dinfo_ll

  real(RP),               private              :: lon_rng
  real(RP),               private              :: lat_rng
  integer,                private              :: lon_num
  integer,                private              :: lat_num

  integer,                private              :: imax_ref
  integer,                private              :: jmax_ref
  integer,                private              :: imax_new
  integer,                private              :: jmax_new
  integer,                private, allocatable :: idx_i(:,:,:)
  integer,                private, allocatable :: idx_j(:,:,:)
  real(RP),               private, allocatable :: hfact(:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_hgridope_setup( &
       nprocs_x_out,    &
       nprocs_y_out,    &
       output_single,   &
       output_grads,    &
       output_gradsctl, &
       enable_plugin,   &
       do_output        )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in)    :: nprocs_x_out             ! x length of 2D processor topology (output)
    integer, intent(in)    :: nprocs_y_out             ! y length of 2D processor topology (output)
    logical, intent(in)    :: output_single            ! output single file when using MPI?
    logical, intent(in)    :: output_grads
    logical, intent(in)    :: output_gradsctl
    logical, intent(out)   :: enable_plugin
    logical, intent(inout) :: do_output

    namelist / PARAM_SNOPLGIN_HGRIDOPE / &
       SNOPLGIN_hgridope_type,      &
       SNOPLGIN_hgridope_lat_start, &
       SNOPLGIN_hgridope_lat_end,   &
       SNOPLGIN_hgridope_dlat,      &
       SNOPLGIN_hgridope_lon_start, &
       SNOPLGIN_hgridope_lon_end,   &
       SNOPLGIN_hgridope_dlon,      &
       SNOPLGIN_hgridope_nintrp,    &
       SNOPLGIN_hgridope_weight

    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SNOPLGIN_hgridope_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SNOPLGIN_HGRIDOPE,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
       LOG_INFO("SNOPLGIN_hgridope_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'Not appropriate names in namelist PARAM_SNOPLGIN_HGRIDOPE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_SNOPLGIN_HGRIDOPE)

    LOG_NEWLINE
    select case(SNOPLGIN_hgridope_type)
    case('OFF')

       LOG_INFO("SNOPLGIN_hgridope_setup",*) 'SNOPLGIN_hgridope_type     : OFF'
       enable_plugin = .false.

    case('LATLON')

       LOG_INFO("SNOPLGIN_hgridope_setup",*) 'SNOPLGIN_hgridope_type     : remap to lat-lon'
       enable_plugin = .true.

       if (      SNOPLGIN_hgridope_lat_start < -90.0_RP &
            .OR. SNOPLGIN_hgridope_lat_start >  90.0_RP &
            .OR. SNOPLGIN_hgridope_lat_end   < -90.0_RP &
            .OR. SNOPLGIN_hgridope_lat_end   >  90.0_RP ) then
          LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'latitude should be within the range between -90 [deg] anf 90 [deg].'
          call PRC_abort
       endif

       if ( SNOPLGIN_hgridope_lat_start > SNOPLGIN_hgridope_lat_end ) then
          LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'SNOPLGIN_hgridope_lat_start should be smaller than SNOPLGIN_hgridope_lat_end.'
          call PRC_abort
       endif

       if (      SNOPLGIN_hgridope_lon_start < -180.0_RP &
            .OR. SNOPLGIN_hgridope_lon_start >  180.0_RP &
            .OR. SNOPLGIN_hgridope_lon_end   < -180.0_RP &
            .OR. SNOPLGIN_hgridope_lon_end   >  180.0_RP ) then
          LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'Longitude should be within the range between -180 [deg] anf 180 [deg].'
          call PRC_abort
       endif

       if ( SNOPLGIN_hgridope_lon_start > SNOPLGIN_hgridope_lon_end ) then
          LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'SNOPLGIN_hgridope_lon_start should be smaller than SNOPLGIN_hgridope_lon_end.'
          call PRC_abort
       endif

       if (      SNOPLGIN_hgridope_dlat <= 0.0_RP &
            .OR. SNOPLGIN_hgridope_dlon <= 0.0_RP ) then
          LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'delta(latitude) and delat(longitude) should be positive.'
          call PRC_abort
       endif

       lon_rng = SNOPLGIN_hgridope_lon_end - SNOPLGIN_hgridope_lon_start
       lon_num = int( lon_rng / SNOPLGIN_hgridope_dlon / nprocs_x_out )

       if ( mod( int( lon_rng / SNOPLGIN_hgridope_dlon ), nprocs_x_out ) /= 0 ) then
          LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'the number of x-process should be able to divide longitudinal size.', &
                                                 nprocs_x_out, int( lon_rng / SNOPLGIN_hgridope_dlon )
          call PRC_abort
       endif

       lat_rng = SNOPLGIN_hgridope_lat_end - SNOPLGIN_hgridope_lat_start 
       lat_num = int( lat_rng / SNOPLGIN_hgridope_dlat / nprocs_y_out )

       if ( mod( int( lat_rng / SNOPLGIN_hgridope_dlat ), nprocs_y_out ) /= 0 ) then
          LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'the number of y-process should be able to divide latitudinal size.', &
                                                 nprocs_y_out, int( lat_rng / SNOPLGIN_hgridope_dlat )
          call PRC_abort
       endif

    case default
       LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'the name of SNOPLGIN_hgridope_type is not appropriate : ', trim(SNOPLGIN_hgridope_type)
       LOG_ERROR_CONT(*) 'you can choose OFF,NUMBER,DAILY,MONTHLY,ANNUAL'
       call PRC_abort
    end select

    if ( enable_plugin ) then
       if    ( output_grads ) then
          LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'This plugin only supports NetCDF format output.'
          call PRC_abort
       elseif( output_gradsctl ) then
          LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'This plugin ignores control file. please turn off output_gradsctl.'
          call PRC_abort
       else
          if ( .NOT. output_single .AND. nprocs_x_out * nprocs_y_out /= 1 ) then
             LOG_ERROR("SNOPLGIN_hgridope_setup",*) 'To use this plugin, the number of output file must be 1.'
             call PRC_abort
          endif
       endif

       ! do not output original data
       do_output = .false.
    endif

    return
  end subroutine SNOPLGIN_hgridope_setup

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_hgridope_setcoef( &
       ismaster,      &
       output_single, &
       nprocs_x_out,  &
       nprocs_y_out,  &
       ngrids_x_out,  &
       ngrids_y_out,  &
       nowrank_x,     &
       nowrank_y,     &
       hinfo,         &
       naxis,         &
       ainfo,         &
       debug          )
    use scale_file_h, only: &
       FILE_REAL8
    use scale_const, only: &
       CONST_D2R,   &
       CONST_RADIUS
    use scale_interp, only: &
       INTERP_setup,   &
       INTERP_factor2d
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo
    use mod_sno_comm, only: &
       SNO_comm_globalaxis
    implicit none

    logical,          intent(in)  :: ismaster                              ! master process?                    (execution)
    logical,          intent(in)  :: output_single                         ! output single file when using MPI?
    integer,          intent(in)  :: nprocs_x_out                          ! x length of 2D processor topology (output)
    integer,          intent(in)  :: nprocs_y_out                          ! y length of 2D processor topology (output)
    integer,          intent(in)  :: ngrids_x_out                          ! number of x-axis grids per process (output,sometimes including halo)
    integer,          intent(in)  :: ngrids_y_out                          ! number of y-axis grids per process (output,sometimes including halo)
    integer,          intent(in)  :: nowrank_x                             ! current process number at x-axis
    integer,          intent(in)  :: nowrank_y                             ! current process number at y-axis
    type(commoninfo), intent(in)  :: hinfo                                 ! common information                 (input)
    integer,          intent(in)  :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(in)  :: ainfo   (naxis)                       ! axis information                   (input)
    logical,          intent(in)  :: debug

    real(RP), allocatable :: lon_ref(:,:) ! [rad]
    real(RP), allocatable :: lat_ref(:,:) ! [rad]
    real(RP), allocatable :: lon_new(:,:) ! [rad]
    real(RP), allocatable :: lat_new(:,:) ! [rad]

    real(RP) :: dxy, dx, dy ! [m]
    real(RP) :: clat        ! [rad]

    type(axisinfo), allocatable :: ainfo_all(:)

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SNOPLGIN_hgridope_setcoef",*) 'Setup remapping coefficient'

    ! set new axis set

    ainfo_ll(1)%varname     = 'lon'
    ainfo_ll(1)%description = 'longitude'
    ainfo_ll(1)%units       = 'degree'
    ainfo_ll(1)%datatype    = FILE_REAL8
    ainfo_ll(1)%dim_rank    = 1
    ainfo_ll(1)%dim_name(1) = 'lon'
    ainfo_ll(1)%transpose   = .false.
    ainfo_ll(1)%regrid      = .false.
    ainfo_ll(1)%has_bounds  = .false.
    ainfo_ll(1)%is_bounds   = .false.

    if ( nowrank_x == 1 ) then
       ainfo_ll(1)%dim_size(1) = lon_num + 1
       allocate( ainfo_ll(1)%AXIS_1d(ainfo_ll(1)%dim_size(1)) )
       ainfo_ll(1)%AXIS_1d(1) = SNOPLGIN_hgridope_lon_start &
                              + lon_rng / nprocs_x_out * ( nowrank_x - 1 )
    else
       ainfo_ll(1)%dim_size(1) = lon_num
       allocate( ainfo_ll(1)%AXIS_1d(ainfo_ll(1)%dim_size(1)) )
       ainfo_ll(1)%AXIS_1d(1) = SNOPLGIN_hgridope_lon_start + SNOPLGIN_hgridope_dlon &
                              + lon_rng / nprocs_x_out * ( nowrank_x - 1 )
    endif
    do i = 2, ainfo_ll(1)%dim_size(1)
       ainfo_ll(1)%AXIS_1d(i) = ainfo_ll(1)%AXIS_1d(i-1) + SNOPLGIN_hgridope_dlon
    enddo

    ainfo_ll(2)%varname     = 'lat'
    ainfo_ll(2)%description = 'latitude'
    ainfo_ll(2)%units       = 'degree'
    ainfo_ll(2)%datatype    = FILE_REAL8
    ainfo_ll(2)%dim_rank    = 1
    ainfo_ll(2)%dim_name(1) = 'lat'
    ainfo_ll(2)%transpose   = .false.
    ainfo_ll(2)%regrid      = .false.
    ainfo_ll(2)%has_bounds  = .false.
    ainfo_ll(2)%is_bounds   = .false.

    if ( nowrank_y == 1 ) then
       ainfo_ll(2)%dim_size(1) = lat_num + 1
       allocate( ainfo_ll(2)%AXIS_1d(ainfo_ll(2)%dim_size(1)) )
       ainfo_ll(2)%AXIS_1d(1) = SNOPLGIN_hgridope_lat_start &
                              + lat_rng / nprocs_y_out * ( nowrank_y - 1 )
    else
       ainfo_ll(2)%dim_size(1) = lat_num
       allocate( ainfo_ll(2)%AXIS_1d(ainfo_ll(2)%dim_size(1)) )
       ainfo_ll(2)%AXIS_1d(1) = SNOPLGIN_hgridope_lat_start + SNOPLGIN_hgridope_dlat &
                              + lat_rng / nprocs_y_out * ( nowrank_y - 1 )
    endif
    do j = 2, ainfo_ll(2)%dim_size(1)
       ainfo_ll(2)%AXIS_1d(j) = ainfo_ll(2)%AXIS_1d(j-1) + SNOPLGIN_hgridope_dlat
    enddo

    naxis_ll = 2

    do n = 1, naxis
       select case(ainfo(n)%varname)
       case('z','zh','pressure','oz','lz','uz')
          naxis_ll = naxis_ll + 1

          ainfo_ll(naxis_ll)%varname     = ainfo(n)%varname
          ainfo_ll(naxis_ll)%description = ainfo(n)%description
          ainfo_ll(naxis_ll)%units       = ainfo(n)%units
          ainfo_ll(naxis_ll)%datatype    = ainfo(n)%datatype
          ainfo_ll(naxis_ll)%dim_rank    = ainfo(n)%dim_rank
          ainfo_ll(naxis_ll)%dim_name(1) = ainfo(n)%dim_name(1)
          ainfo_ll(naxis_ll)%dim_size(1) = ainfo(n)%dim_size(1)
          ainfo_ll(naxis_ll)%transpose   = ainfo(n)%transpose
          ainfo_ll(naxis_ll)%regrid      = ainfo(n)%regrid
          allocate( ainfo_ll(naxis_ll)%AXIS_1d(ainfo_ll(naxis_ll)%dim_size(1)) )
          ainfo_ll(naxis_ll)%AXIS_1d(:)  = ainfo(n)%AXIS_1d(:)
       end select
    enddo

    ! set remapping coefficient
    clat = ainfo_ll(2)%AXIS_1d(ainfo_ll(2)%dim_size(1)/2) * CONST_D2R

    dx  = SNOPLGIN_hgridope_dlon * CONST_D2R * CONST_RADIUS * cos(clat)
    dy  = SNOPLGIN_hgridope_dlat * CONST_D2R * CONST_RADIUS * cos(clat)
    dxy = sqrt( dx*dx + dy*dy ) * 10.0_RP

    call INTERP_setup( SNOPLGIN_hgridope_weight, & ! [IN]
                       search_limit = dxy        ) ! [IN]

    imax_ref = ngrids_x_out * nprocs_x_out
    jmax_ref = ngrids_y_out * nprocs_y_out
    imax_new = ainfo_ll(1)%dim_size(1)
    jmax_new = ainfo_ll(2)%dim_size(1)

    allocate( lon_ref(imax_ref,jmax_ref) )
    allocate( lat_ref(imax_ref,jmax_ref) )
    allocate( lon_new(imax_new,jmax_new) )
    allocate( lat_new(imax_new,jmax_new) )

    allocate( idx_i(imax_new,jmax_new,SNOPLGIN_hgridope_nintrp) )
    allocate( idx_j(imax_new,jmax_new,SNOPLGIN_hgridope_nintrp) )
    allocate( hfact(imax_new,jmax_new,SNOPLGIN_hgridope_nintrp) )

    allocate( ainfo_all(naxis) )

    call SNO_comm_globalaxis( ismaster,      & ! [IN]
                              output_single, & ! [IN]
                              nprocs_x_out,  & ! [IN]
                              nprocs_y_out,  & ! [IN]
                              hinfo,         & ! [IN]
                              naxis,         & ! [IN]
                              ainfo    (:),  & ! [IN]
                              ainfo_all(:),  & ! [OUT]
                              bcast = .true. ) ! [IN]

    do n = 1, naxis
       if    ( ainfo_all(n)%varname == 'lon' ) then
          lon_ref(:,:) = ainfo_all(n)%AXIS_2d(:,:) * CONST_D2R
       elseif( ainfo_all(n)%varname == 'lat' ) then
          lat_ref(:,:) = ainfo_all(n)%AXIS_2d(:,:) * CONST_D2R
       endif
    enddo

    do i = 1, imax_new
    do j = 1, jmax_new
       lon_new(i,j) = ainfo_ll(1)%AXIS_1d(i) * CONST_D2R
       lat_new(i,j) = ainfo_ll(2)%AXIS_1d(j) * CONST_D2R
    enddo
    enddo

    call INTERP_factor2d( SNOPLGIN_hgridope_nintrp,   & ! [IN]
                          imax_ref, jmax_ref,         & ! [IN]
                          imax_new, jmax_new,         & ! [IN]
                          lon_ref(:,:), lat_ref(:,:), & ! [IN]
                          lon_new(:,:), lat_new(:,:), & ! [IN]
                          idx_i(:,:,:), idx_j(:,:,:), & ! [OUT]
                          hfact(:,:,:)                ) ! [OUT]

    deallocate( lon_ref )
    deallocate( lat_ref )
    deallocate( lon_new )
    deallocate( lat_new )

    return
  end subroutine SNOPLGIN_hgridope_setcoef

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_hgridope_alloc( &
       dinfo, &
       debug  )
    use mod_sno_h, only: &
       iteminfo
    implicit none

    type(iteminfo), intent(in)  :: dinfo ! variable information               (input)
    logical,        intent(in)  :: debug

    integer  :: kmax
    !---------------------------------------------------------------------------

    if ( debug ) then
       LOG_INFO("SNOPLGIN_hgridope_alloc",*) 'Allocate temporal array'
    endif

    if ( dinfo%dim_rank == 1 ) then

       ! do nothing

    elseif( dinfo%dim_rank == 2 ) then

       allocate( dinfo_ll%VAR_2d(imax_new,jmax_new) )
       dinfo_ll%VAR_2d(:,:) = 0.0_RP

    elseif( dinfo%dim_rank == 3 ) then

       if ( dinfo%transpose ) then
          kmax = dinfo%dim_size(3)
       else
          kmax = dinfo%dim_size(1)
       endif
       allocate( dinfo_ll%VAR_3d(kmax,imax_new,jmax_new) )

       dinfo_ll%VAR_3d(:,:,:) = 0.0_RP

    endif

    return
  end subroutine SNOPLGIN_hgridope_alloc

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_hgridope_dealloc( &
       debug  )
    use mod_sno_h, only: &
       iteminfo
    implicit none

    logical,        intent(in)  :: debug
    !---------------------------------------------------------------------------

    if ( debug ) then
       LOG_INFO("SNOPLGIN_hgridope_dealloc",*) 'Deallocate temporal array'
    endif

    if( allocated(dinfo_ll%VAR_1d) ) deallocate( dinfo_ll%VAR_1d )
    if( allocated(dinfo_ll%VAR_2d) ) deallocate( dinfo_ll%VAR_2d )
    if( allocated(dinfo_ll%VAR_3d) ) deallocate( dinfo_ll%VAR_3d )

    return
  end subroutine SNOPLGIN_hgridope_dealloc

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_hgridope_llinterp( &
       ismaster,      &
       dirpath,       &
       basename,      &
       output_single, &
       output_grads,  &
       nowrank,       &
       nowstep,       &
       nprocs_x_out,  &
       nprocs_y_out,  &
       nhalos_x,      &
       nhalos_y,      &
       hinfo,         &
       dinfo,         &
       debug          )
    use scale_interp, only: &
       INTERP_interp2d
    use mod_sno_h, only: &
       commoninfo, &
       iteminfo
    use mod_sno_vars, only: &
       SNO_vars_write
    use mod_sno_comm, only: &
       SNO_comm_globalvars
    implicit none

    logical,          intent(in)    :: ismaster                              ! master process?                    (execution)
    character(len=*), intent(in)    :: dirpath                               ! directory path                     (output)
    character(len=*), intent(in)    :: basename                              ! basename of file                   (output)
    logical,          intent(in)    :: output_single                         ! output single file when using MPI?
    logical,          intent(in)    :: output_grads
    integer,          intent(in)    :: nowrank                               ! current rank                       (output)
    integer,          intent(in)    :: nowstep                               ! current step                       (output)
    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology  (output)
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology  (output)
    integer,          intent(in)    :: nhalos_x                              ! number of x-axis halo grids        (global domain)
    integer,          intent(in)    :: nhalos_y                              ! number of y-axis halo grids        (global domain)
    type(commoninfo), intent(in)    :: hinfo                                 ! common information                 (input)
    type(iteminfo),   intent(in)    :: dinfo                                 ! variable information               (input)
    logical,          intent(in)    :: debug

    type(iteminfo) :: dinfo_all

    logical  :: do_output, finalize, add_rm_attr
    integer  :: k, t
    !---------------------------------------------------------------------------

    do_output = .false.

    ! set variable information
    dinfo_ll%varname     = dinfo%varname
    dinfo_ll%description = dinfo%description
    dinfo_ll%units       = dinfo%units
    dinfo_ll%datatype    = dinfo%datatype
    dinfo_ll%dim_rank    = dinfo%dim_rank
    dinfo_ll%transpose   = .true.
    dinfo_ll%step_nmax   = dinfo%step_nmax
    do t = 1, dinfo%step_nmax
       dinfo_ll%time_start(t) = dinfo%time_start(t)
       dinfo_ll%time_end  (t) = dinfo%time_end  (t)
    enddo
    dinfo_ll%dt          = dinfo%dt
    dinfo_ll%time_units  = dinfo%time_units

    if ( dinfo%dim_rank == 1 ) then

       ! do nothing

    elseif( dinfo%dim_rank == 2 ) then

       dinfo_ll%dim_name(1) = "lon"
       dinfo_ll%dim_name(2) = "lat"
       dinfo_ll%dim_size(1) = imax_new
       dinfo_ll%dim_size(2) = jmax_new

       call SNO_comm_globalvars( ismaster,      & ! [IN]
                                 output_single, & ! [IN]
                                 nprocs_x_out,  & ! [IN]
                                 nprocs_y_out,  & ! [IN]
                                 dinfo,         & ! [IN]
                                 dinfo_all,     & ! [OUT]
                                 bcast = .true. ) ! [IN]

       call INTERP_interp2d( SNOPLGIN_hgridope_nintrp,  & ! [IN]
                             imax_ref, jmax_ref,        & ! [IN]
                             imax_new, jmax_new,        & ! [IN]
                             idx_i           (:,:,:),   & ! [IN]
                             idx_j           (:,:,:),   & ! [IN]
                             hfact           (:,:,:),   & ! [IN]
                             dinfo_all%VAR_2d(:,:),     & ! [IN]
                             dinfo_ll%VAR_2d (:,:)      ) ! [OUT]

       do_output = .true.

    elseif( dinfo%dim_rank == 3 ) then

       dinfo_ll%dim_name(1) = "lon"
       dinfo_ll%dim_size(1) = imax_new
       dinfo_ll%dim_name(2) = "lat"
       dinfo_ll%dim_size(2) = jmax_new

       if ( dinfo%transpose ) then
          dinfo_ll%dim_name(3) = dinfo%dim_name(3)
          dinfo_ll%dim_size(3) = dinfo%dim_size(3)
       else
          dinfo_ll%dim_name(3) = dinfo%dim_name(1)
          dinfo_ll%dim_size(3) = dinfo%dim_size(1)
       endif

       call SNO_comm_globalvars( ismaster,      & ! [IN]
                                 output_single, & ! [IN]
                                 nprocs_x_out,  & ! [IN]
                                 nprocs_y_out,  & ! [IN]
                                 dinfo,         & ! [IN]
                                 dinfo_all,     & ! [OUT]
                                 bcast = .true. ) ! [IN]

       do k = 1, dinfo_ll%dim_size(3)
          call INTERP_interp2d( SNOPLGIN_hgridope_nintrp,  & ! [IN]
                                imax_ref, jmax_ref,        & ! [IN]
                                imax_new, jmax_new,        & ! [IN]
                                idx_i           (:,:,:),   & ! [IN]
                                idx_j           (:,:,:),   & ! [IN]
                                hfact           (:,:,:),   & ! [IN]
                                dinfo_all%VAR_3d(k,:,:),   & ! [IN]
                                dinfo_ll%VAR_3d (k,:,:)    ) ! [OUT]
       enddo

       do_output = .true.
    endif

    ! output

    if ( do_output ) then
       finalize    = ( nowstep == dinfo_ll%step_nmax )
       add_rm_attr = .false.

       call SNO_vars_write( ismaster,                   & ! [IN] from MPI
                            dirpath,                    & ! [IN] from namelist
                            basename,                   & ! [IN] from namelist
                            output_single,              & ! [IN] from namelist
                            output_grads,               & ! [IN] from namelist
                            nowrank,                    & ! [IN]
                            nowstep,                    & ! [IN]
                            finalize,                   & ! [IN]
                            add_rm_attr,                & ! [IN]
                            nprocs_x_out, nprocs_y_out, & ! [IN] from namelist
                            nhalos_x,     nhalos_y,     & ! [IN] from SNO_file_getinfo
                            hinfo,                      & ! [IN] from SNO_file_getinfo
                            naxis_ll,                   & ! [IN] from SNO_file_getinfo
                            ainfo_ll(1:naxis_ll),       & ! [IN] from SNO_axis_getinfo
                            dinfo_ll,                   & ! [IN] from SNO_vars_getinfo
                            debug                       ) ! [IN]
    endif

    return
  end subroutine SNOPLGIN_hgridope_llinterp

end module mod_snoplugin_hgridope
