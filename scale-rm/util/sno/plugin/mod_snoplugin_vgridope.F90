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
module mod_snoplugin_vgridope
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use mod_sno_h, only: &
     lev_limit, &
     axisinfo,  &
     iteminfo
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SNOPLGIN_vgridope_setup
  public :: SNOPLGIN_vgridope_setcoef
  public :: SNOPLGIN_vgridope_alloc
  public :: SNOPLGIN_vgridope_dealloc
  public :: SNOPLGIN_vgridope_updatecoef
  public :: SNOPLGIN_vgridope_vinterp

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
  character(len=H_SHORT), private              :: SNOPLGIN_vgridope_type        = 'OFF'   ! type of average
                                                                                ! 'OFF'    : disable
                                                                                ! 'ZLEV'   : remap to z-level grid
                                                                                ! 'PLEV'   : remap to pressure-level grid
  integer,                private              :: SNOPLGIN_vgridope_lev_num             = -1
  real(RP),               private              :: SNOPLGIN_vgridope_lev_data(lev_limit) = -1.0_RP

  type(axisinfo),         private, allocatable :: ainfo_v(:)
  type(iteminfo),         private              :: dinfo_v

  integer,                private              :: znum
  integer,                private              :: kmax_ref
  integer,                private              :: imax_ref
  integer,                private              :: jmax_ref
  integer,                private              :: kmax_new

  integer,                private, allocatable :: idx_Z (:,:)
  integer,                private, allocatable :: idx_Zh(:,:)
  real(RP),               private, allocatable :: Zfact (:)
  real(RP),               private, allocatable :: Zhfact(:)

  real(RP),               private, allocatable :: Z_ref (:)
  real(RP),               private, allocatable :: Zh_ref(:)

  integer,                private, allocatable :: idx_P (:,:,:,:)
  integer,                private, allocatable :: idx_Ph(:,:,:,:)
  real(RP),               private, allocatable :: Pfact (:,:,:)
  real(RP),               private, allocatable :: Phfact(:,:,:)

  real(RP),               private, allocatable :: PRES_ref  (:,:,:)
  real(RP),               private, allocatable :: PRESh_ref (:,:,:)
  real(RP),               private, allocatable :: height_ref(:,:,:)

  real(RP),               private, allocatable :: SFC_PRES_ref(:,:)

  real(RP),               private, allocatable :: LnPaxis(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_vgridope_setup( &
       nprocs_x_out,    &
       nprocs_y_out,    &
       output_grads,    &
       output_gradsctl, &
       enable_plugin,   &
       do_output        )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in)    :: nprocs_x_out             ! x length of 2D processor topology (output)
    integer, intent(in)    :: nprocs_y_out             ! y length of 2D processor topology (output)
    logical, intent(in)    :: output_grads
    logical, intent(in)    :: output_gradsctl
    logical, intent(out)   :: enable_plugin
    logical, intent(inout) :: do_output

    namelist / PARAM_SNOPLGIN_VGRIDOPE / &
       SNOPLGIN_vgridope_type,        &
       SNOPLGIN_vgridope_lev_num,     &
       SNOPLGIN_vgridope_lev_data

    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SNOPLGIN_vgridope_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SNOPLGIN_VGRIDOPE,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
       LOG_INFO("SNOPLGIN_vgridope_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SNOPLGIN_vgridope_setup",*) 'Not appropriate names in namelist PARAM_SNOPLGIN_VGRIDOPE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_SNOPLGIN_VGRIDOPE)

    LOG_NEWLINE
    select case(SNOPLGIN_vgridope_type)
    case('OFF')

       LOG_INFO("SNOPLGIN_vgridope_setup",*) 'SNOPLGIN_vgridope_type     : OFF'
       enable_plugin = .false.

    case('ZLEV')

       LOG_INFO("SNOPLGIN_vgridope_setup",*) 'SNOPLGIN_vgridope_type     : remap to z-level'
       enable_plugin = .true.

       if ( SNOPLGIN_vgridope_lev_num <= 0 ) then
          LOG_ERROR("SNOPLGIN_vgridope_setup",*) 'The number of vertical layers for interpolation should be positive: SNOPLGIN_vgridope_lev_num'
          call PRC_abort
       endif
       if ( all( SNOPLGIN_vgridope_lev_data(:) < 0.0_RP ) ) then
          LOG_ERROR("SNOPLGIN_vgridope_setup",*) 'At least one vertical layer for interpolation should be specified: SNOPLGIN_vgridope_lev_data'
          call PRC_abort
       endif

    case('PLEV')

       LOG_INFO("SNOPLGIN_vgridope_setup",*) 'SNOPLGIN_vgridope_type     : remap to pressure-level'
       enable_plugin = .true.

       if ( SNOPLGIN_vgridope_lev_num <= 0 ) then
          LOG_ERROR("SNOPLGIN_vgridope_setup",*) 'The number of vertical layers for interpolation should be positive.'
          call PRC_abort
       endif
       if ( all( SNOPLGIN_vgridope_lev_data(:) < 0.0_RP ) ) then
          LOG_ERROR("SNOPLGIN_vgridope_setup",*) 'At least one vertical layer for interpolation should be specified: SNOPLGIN_vgridope_lev_data'
          call PRC_abort
       endif

    case default
       LOG_ERROR("SNOPLGIN_vgridope_setup",*) 'the name of SNOPLGIN_vgridope_type is not appropriate : ', trim(SNOPLGIN_vgridope_type)
       LOG_ERROR_CONT(*) 'you can choose OFF,ZLEV,PLEV'
       call PRC_abort
    end select

    if ( enable_plugin ) then
       if    ( output_grads ) then
          LOG_ERROR("SNOPLGIN_vgridope_setup",*) 'This plugin only supports NetCDF format output.'
          call PRC_abort
       elseif( output_gradsctl ) then
          LOG_ERROR("SNOPLGIN_vgridope_setup",*) 'This plugin ignores control file. please turn off output_gradsctl.'
          call PRC_abort
       endif

       ! do not output original data
       do_output = .false.
    endif

    return
  end subroutine SNOPLGIN_vgridope_setup

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_vgridope_setcoef( &
       ngrids_z_out, &
       ngrids_x_out, &
       ngrids_y_out, &
       naxis,        &
       ainfo,        &
       debug         )
    use scale_file_h, only: &
       FILE_REAL8
    use scale_interp, only: &
       INTERP_setup,    &
       INTERP_factor1d
    use mod_sno_h, only: &
       axisinfo
    implicit none

    integer,        intent(in)  :: ngrids_z_out                          ! number of z-axis grids per process (output,sometimes including halo)
    integer,        intent(in)  :: ngrids_x_out                          ! number of x-axis grids per process (output,sometimes including halo)
    integer,        intent(in)  :: ngrids_y_out                          ! number of y-axis grids per process (output,sometimes including halo)
    integer,        intent(in)  :: naxis                                 ! number of axis variables           (input)
    type(axisinfo), intent(in)  :: ainfo   (naxis)                       ! axis information                   (input)
    logical,        intent(in)  :: debug

    integer  :: k, n
    !---------------------------------------------------------------------------

    ! set region size
    kmax_ref = ngrids_z_out
    imax_ref = ngrids_x_out
    jmax_ref = ngrids_y_out

    kmax_new = SNOPLGIN_vgridope_lev_num

    select case( trim(SNOPLGIN_vgridope_type ) )
    case('ZLEV')

       LOG_NEWLINE
       LOG_INFO("SNOPLGIN_vgridope_setcoef",*) 'Setup remapping coefficient (height)'

       allocate( Z_ref (  kmax_ref) )
       allocate( Zh_ref(0:kmax_ref) )

       allocate( idx_Z (kmax_new,2) )
       allocate( idx_Zh(kmax_new,2) )
       allocate( Zfact (kmax_new  ) )
       allocate( Zhfact(kmax_new  ) )

       ! set basic axis
       allocate( ainfo_v( naxis ) )
       ainfo_v(:) = ainfo(:)

       do n = 1, naxis
          select case( trim(ainfo(n)%varname) )
          case('z')
             Z_ref (:) = ainfo(n)%AXIS_1d(:)

             znum = n

             ! rewrite axis
             if( allocated( ainfo_v(znum)%AXIS_1d ) ) deallocate( ainfo_v(znum)%AXIS_1d )
             ainfo_v(znum)%dim_size(1) = kmax_new
             allocate( ainfo_v(znum)%AXIS_1d(kmax_new) )

             do k = 1, kmax_new
                ainfo_v(znum)%AXIS_1d(k) = SNOPLGIN_vgridope_lev_data(k)
             enddo
          case('zh')
             Zh_ref(:) = ainfo(n)%AXIS_1d(:)
          endselect
       enddo

       ! set remapping coefficient
       call INTERP_setup( 2 ) ! [IN] not used

       call INTERP_factor1d( kmax_ref, 1, kmax_ref,    & ! [IN]
                             kmax_new, 1, kmax_new,    & ! [IN]
                             Z_ref(:),                 & ! [IN]
                             ainfo_v(znum)%AXIS_1d(:), & ! [IN]
                             idx_Z(:,:),               & ! [OUT]
                             Zfact(:),                 & ! [OUT]
                             flag_extrap = .false.     ) ! [IN]

       call INTERP_factor1d( kmax_ref+1, 1, kmax_ref+1, & ! [IN]
                             kmax_new,   1, kmax_new,   & ! [IN]
                             Zh_ref(:),                 & ! [IN]
                             ainfo_v(znum)%AXIS_1d(:),  & ! [IN]
                             idx_Zh(:,:),               & ! [OUT]
                             Zhfact(:),                 & ! [OUT]
                             flag_extrap = .false.      ) ! [IN]

    case('PLEV')

       LOG_NEWLINE
       LOG_INFO("SNOPLGIN_vgridope_setcoef",*) 'Setup remapping coefficient (pressure)'

       allocate( PRES_ref  (  kmax_ref,imax_ref,jmax_ref) )
       allocate( PRESh_ref (0:kmax_ref,imax_ref,jmax_ref) )
       allocate( height_ref(  kmax_ref,imax_ref,jmax_ref) )

       PRES_ref  (:,:,:) = -1.0_RP
       PRESh_ref (:,:,:) = -1.0_RP
       height_ref(:,:,:) = -1.0_RP

       allocate( SFC_PRES_ref(imax_ref,jmax_ref) )

       SFC_PRES_ref(:,:) = -1.0_RP

       allocate( idx_P (kmax_new,2,imax_ref,jmax_ref) )
       allocate( idx_Ph(kmax_new,2,imax_ref,jmax_ref) )
       allocate( Pfact (kmax_new,  imax_ref,jmax_ref) )
       allocate( Phfact(kmax_new,  imax_ref,jmax_ref) )

       ! set basic axis
       allocate( ainfo_v( naxis+1 ) )
       ainfo_v(2:naxis+1) = ainfo(:)

       znum = 1

       ! add pressure axis
       ainfo_v(znum)%varname     = 'pressure'
       ainfo_v(znum)%description = 'Pressure Level'
       ainfo_v(znum)%units       = 'Pa'
       ainfo_v(znum)%datatype    = FILE_REAL8
       ainfo_v(znum)%dim_rank    = 1
       ainfo_v(znum)%dim_name(1) = 'pressure'
       ainfo_v(znum)%transpose   = .false.
       ainfo_v(znum)%regrid      = .false.

       if( allocated( ainfo_v(znum)%AXIS_1d ) ) deallocate( ainfo_v(znum)%AXIS_1d )
       ainfo_v(znum)%dim_size(1) = kmax_new
       allocate( ainfo_v(znum)%AXIS_1d(kmax_new) )

       do k = 1, kmax_new
          ainfo_v(znum)%AXIS_1d(k) = SNOPLGIN_vgridope_lev_data(k)
       enddo

       allocate( LnPaxis(kmax_new) )

       ! logarithmic
       LnPaxis(:) = - log( ainfo_v(znum)%AXIS_1d(:) )

       ! geopotential height
       do n = 1, naxis
          select case( trim(ainfo(n)%varname) )
          case('height')
             height_ref(:,:,:) = ainfo(n)%AXIS_3d(:,:,:)
          endselect
       enddo

       ! set remapping coefficient
       call INTERP_setup( 2 ) ! [IN] not used

    end select

    return
  end subroutine SNOPLGIN_vgridope_setcoef

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_vgridope_alloc( &
       dinfo, &
       debug  )
    use mod_sno_h, only: &
       iteminfo
    implicit none

    type(iteminfo), intent(in)  :: dinfo ! variable information               (input)
    logical,        intent(in)  :: debug

    character(len=H_SHORT) :: zaxis_orgname

    integer  :: zaxis_orgsize
    !---------------------------------------------------------------------------

    if ( debug ) then
       LOG_INFO("SNOPLGIN_vgridope_alloc",*) 'Allocate temporal array'
    endif

    if ( dinfo%dim_rank == 1 ) then

       ! do nothing

    elseif( dinfo%dim_rank == 2 ) then

       allocate( dinfo_v%VAR_2d(imax_ref,jmax_ref) )
       dinfo_v%VAR_2d(:,:) = 0.0_RP

    elseif( dinfo%dim_rank == 3 ) then

       if ( dinfo%transpose ) then
          zaxis_orgname = dinfo%dim_name(3)
          zaxis_orgsize = dinfo%dim_size(3)
       else
          zaxis_orgname = dinfo%dim_name(1)
          zaxis_orgsize = dinfo%dim_size(1)
       endif

       select case( trim(zaxis_orgname) )
       case('z')
          allocate( dinfo_v%VAR_3d(kmax_new,     imax_ref,jmax_ref) )
       case('zh')
          allocate( dinfo_v%VAR_3d(kmax_new,     imax_ref,jmax_ref) )
       case default
          allocate( dinfo_v%VAR_3d(zaxis_orgsize,imax_ref,jmax_ref) )
       endselect

       dinfo_v%VAR_3d(:,:,:) = 0.0_RP

    endif

    return
  end subroutine SNOPLGIN_vgridope_alloc

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_vgridope_dealloc( &
       debug  )
    use mod_sno_h, only: &
       iteminfo
    implicit none

    logical,        intent(in)  :: debug
    !---------------------------------------------------------------------------

    if ( debug ) then
       LOG_INFO("SNOPLGIN_vgridope_dealloc",*) 'Deallocate temporal array'
    endif

    if( allocated(dinfo_v%VAR_1d) ) deallocate( dinfo_v%VAR_1d )
    if( allocated(dinfo_v%VAR_2d) ) deallocate( dinfo_v%VAR_2d )
    if( allocated(dinfo_v%VAR_3d) ) deallocate( dinfo_v%VAR_3d )

    return
  end subroutine SNOPLGIN_vgridope_dealloc

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_vgridope_updatecoef( &
       basename_in,   &
       nowvar,        &
       nowstep,       &
       nprocs_x_in,   &
       nprocs_y_in,   &
       ngrids_x,      &
       ngrids_y,      &
       nhalos_x,      &
       nhalos_y,      &
       hinfo,         &
       ngrids_x_out,  &
       ngrids_y_out,  &
       ngrids_xh_out, &
       ngrids_yh_out, &
       nvars,         &
       dinfo,         &
       localmap,      &
       readflag,      &
       debug          )
    use scale_prc, only: &
       PRC_abort
    use scale_interp, only: &
       INTERP_factor1d
    use mod_sno_h, only: &
       commoninfo, &
       iteminfo
    use mod_sno_vars, only: &
       SNO_vars_alloc,   &
       SNO_vars_dealloc, &
       SNO_vars_read
    implicit none

    character(len=H_LONG), intent(in) :: basename_in   ! Basename of the input  file

    integer,    intent(in) :: nowvar        ! target variable number
    integer,    intent(in) :: nowstep       ! current step
    integer,    intent(in) :: nprocs_x_in   ! x length of 2D processor topology (input)
    integer,    intent(in) :: nprocs_y_in   ! y length of 2D processor topology (input)
    integer,    intent(in) :: ngrids_x      ! size of x-axis grids              (global,sometimes including halo)
    integer,    intent(in) :: ngrids_y      ! size of y-axis grids              (global,sometimes including halo)
    integer,    intent(in) :: nhalos_x      ! size of x-axis halo grids         (global,sometimes have a size)
    integer,    intent(in) :: nhalos_y      ! size of y-axis halo grids         (global,sometimes have a size)

    type(commoninfo), intent(in) :: hinfo

    integer,    intent(in) :: ngrids_x_out  ! size of x-axis grids              (output,sometimes including halo)
    integer,    intent(in) :: ngrids_y_out  ! size of y-axis grids              (output,sometimes including halo)
    integer,    intent(in) :: ngrids_xh_out ! size of x-axis grids, staggard    (output,sometimes including halo)
    integer,    intent(in) :: ngrids_yh_out ! size of y-axis grids, staggard    (output,sometimes including halo)

    integer,        intent(in)  :: nvars          ! number of item variables
    type(iteminfo), intent(in)  :: dinfo(nvars)   ! item information

    integer(2), intent(in) :: localmap (ngrids_x_out,ngrids_y_out,3)  ! mapping table
    logical,    intent(in) :: readflag (nprocs_x_in,nprocs_y_in)      ! flag to read each input file
    logical,    intent(in) :: debug

    type(iteminfo) :: pinfo

    integer :: i, j, k, v
    !---------------------------------------------------------------------------

    select case( trim(SNOPLGIN_vgridope_type) )
    case('PLEV')

       ! update PRES and SFC_PRES
       do v = 1, nvars

          select case( trim(dinfo(v)%varname) )
          case('PRES')

             if( v == nowvar ) then
                PRES_ref(:,:,:) = dinfo(v)%VAR_3d(:,:,:)
             else
                pinfo = dinfo(v)

                call SNO_vars_alloc( ngrids_x_out,  ngrids_y_out,  & ! [IN]    from SNO_map_getsize_local
                                     ngrids_xh_out, ngrids_yh_out, & ! [IN]    from SNO_map_getsize_local
                                     pinfo,                        & ! [INOUT] from SNO_vars_getinfo
                                     debug                         ) ! [IN]

                call SNO_vars_read( basename_in,                  & ! [IN]    from namelist
                                    nowstep,                      & ! [IN]
                                    nprocs_x_in,   nprocs_y_in,   & ! [IN]    from SNO_file_getinfo
                                    ngrids_x,      ngrids_y,      & ! [IN]    from SNO_file_getinfo
                                    nhalos_x,      nhalos_y,      & ! [IN]    from SNO_file_getinfo
                                    hinfo,                        & ! [IN]    from SNO_file_getinfo
                                    ngrids_x_out,  ngrids_y_out,  & ! [IN]    from SNO_map_getsize_local
                                    ngrids_xh_out, ngrids_yh_out, & ! [IN]    from SNO_map_getsize_local
                                    pinfo,                        & ! [INOUT] from SNO_vars_getinfo
                                    localmap(:,:,:),              & ! [IN]    from SNO_map_settable_local
                                    readflag(:,:),                & ! [IN]    from SNO_map_settable_local
                                    debug                         ) ! [IN]

                PRES_ref(:,:,:) = pinfo%VAR_3d(:,:,:)

                call SNO_vars_dealloc( pinfo, & ! [INOUT] from SNO_vars_getinfo
                                       debug  ) ! [IN]
             endif

          case('SFC_PRES')

             if( v == nowvar ) then
                SFC_PRES_ref(:,:) = dinfo(v)%VAR_2d(:,:)
             else
                pinfo = dinfo(v)

                call SNO_vars_alloc( ngrids_x_out,  ngrids_y_out,  & ! [IN]    from SNO_map_getsize_local
                                     ngrids_xh_out, ngrids_yh_out, & ! [IN]    from SNO_map_getsize_local
                                     pinfo,                        & ! [INOUT] from SNO_vars_getinfo
                                     debug                         ) ! [IN]

                call SNO_vars_read( basename_in,                  & ! [IN]    from namelist
                                    nowstep,                      & ! [IN]
                                    nprocs_x_in,   nprocs_y_in,   & ! [IN]    from SNO_file_getinfo
                                    ngrids_x,      ngrids_y,      & ! [IN]    from SNO_file_getinfo
                                    nhalos_x,      nhalos_y,      & ! [IN]    from SNO_file_getinfo
                                    hinfo,                        & ! [IN]    from SNO_file_getinfo
                                    ngrids_x_out,  ngrids_y_out,  & ! [IN]    from SNO_map_getsize_local
                                    ngrids_xh_out, ngrids_yh_out, & ! [IN]    from SNO_map_getsize_local
                                    pinfo,                        & ! [INOUT] from SNO_vars_getinfo
                                    localmap(:,:,:),              & ! [IN]    from SNO_map_settable_local
                                    readflag(:,:),                & ! [IN]    from SNO_map_settable_local
                                    debug                         ) ! [IN]

                SFC_PRES_ref(:,:) = pinfo%VAR_2d(:,:)

                call SNO_vars_dealloc( pinfo, & ! [INOUT] from SNO_vars_getinfo
                                       debug  ) ! [IN]
             endif

          endselect

       enddo

       ! check PRES and SFC_PRES
       if( all( PRES_ref(:,:,:) < 0.0_RP ) .or. all( SFC_PRES_ref(:,:) < 0.0_RP ) ) then
          LOG_ERROR("SNOPLGIN_vgridope_setcoef",*) 'Not found required pressure data. Check "PRES" and "SFC_PRES"!'
          call PRC_abort
       endif

       ! update pressure at half-level
       PRESh_ref(0,:,:) = SFC_PRES_ref(:,:)
       do k = 1, kmax_ref-1
          PRESh_ref(k,:,:) = PRESh_ref(k-1,:,:) + ( PRES_ref(k,:,:) - PRES_ref(k+1,:,:) )
       enddo
       PRESh_ref(kmax_ref,:,:) = PRESh_ref(kmax_ref-1,:,:) + ( PRES_ref(kmax_ref-1,:,:) - PRES_ref(kmax_ref,:,:) )

       ! logarithmic
       PRES_ref (:,:,:) = - log( PRES_ref (:,:,:) )
       PRESh_ref(:,:,:) = - log( PRESh_ref(:,:,:) )

       SFC_PRES_ref(:,:) = - log( SFC_PRES_ref(:,:) )

       ! update remapping coefficient
       do j = 1, jmax_ref
       do i = 1, imax_ref
          call INTERP_factor1d( kmax_ref, 1, kmax_ref, & ! [IN]
                                kmax_new, 1, kmax_new, & ! [IN]
                                PRES_ref(:,i,j),       & ! [IN]
                                LnPaxis(:),            & ! [IN]
                                idx_P(:,:,i,j),        & ! [OUT]
                                Pfact(:,i,j),          & ! [OUT]
                                flag_extrap = .false.  ) ! [IN]

          call INTERP_factor1d( kmax_ref+1, 1, kmax_ref+1, & ! [IN]
                                kmax_new,   1, kmax_new,   & ! [IN]
                                PRESh_ref(:,i,j),          & ! [IN]
                                LnPaxis(:),                & ! [IN]
                                idx_Ph(:,:,i,j),           & ! [OUT]
                                Phfact(:,i,j),             & ! [OUT]
                                flag_extrap = .false.      ) ! [IN]
       enddo
       enddo

    endselect

    return
  end subroutine SNOPLGIN_vgridope_updatecoef

  !-----------------------------------------------------------------------------
  subroutine SNOPLGIN_vgridope_vinterp( &
       dirpath,       &
       basename,      &
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
       INTERP_interp1d
    use mod_sno_h, only: &
       commoninfo, &
       iteminfo
    use mod_sno_vars, only: &
       SNO_vars_write
    implicit none

    character(len=*), intent(in)    :: dirpath                               ! directory path                     (output)
    character(len=*), intent(in)    :: basename                              ! basename of file                   (output)
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

    character(len=H_SHORT) :: zaxis_orgname

    logical  :: do_output, finalize, add_rm_attr
    integer  :: zaxis_orgsize
    integer  :: i, j, t
    !---------------------------------------------------------------------------

    do_output = .false.

    ! set variable information
    dinfo_v%varname     = dinfo%varname
    dinfo_v%description = dinfo%description
    dinfo_v%units       = dinfo%units
    dinfo_v%datatype    = dinfo%datatype
    dinfo_v%dim_rank    = dinfo%dim_rank
    dinfo_v%transpose   = .true.
    dinfo_v%step_nmax   = dinfo%step_nmax
    do t = 1, dinfo%step_nmax
       dinfo_v%time_start(t) = dinfo%time_start(t)
       dinfo_v%time_end  (t) = dinfo%time_end  (t)
    enddo
    dinfo_v%dt          = dinfo%dt
    dinfo_v%time_units  = dinfo%time_units

    if ( dinfo%dim_rank == 1 ) then

       ! do nothing

    elseif( dinfo%dim_rank == 2 ) then

       ! do output file but no vertical interpolation
       dinfo_v%dim_name(1) = dinfo%dim_name(1)
       dinfo_v%dim_name(2) = dinfo%dim_name(2)
       dinfo_v%dim_size(1) = imax_ref
       dinfo_v%dim_size(2) = jmax_ref

       dinfo_v%VAR_2d(:,:) = dinfo%VAR_2d(:,:)

       do_output = .true.

    elseif( dinfo%dim_rank == 3 ) then

       if ( dinfo%transpose ) then
          dinfo_v%dim_name(1) = dinfo%dim_name(1)
          dinfo_v%dim_name(2) = dinfo%dim_name(2)
          zaxis_orgname = dinfo%dim_name(3)
          zaxis_orgsize = dinfo%dim_size(3)
       else
          dinfo_v%dim_name(1) = dinfo%dim_name(2)
          dinfo_v%dim_name(2) = dinfo%dim_name(3)
          zaxis_orgname = dinfo%dim_name(1)
          zaxis_orgsize = dinfo%dim_size(1)
       endif
       dinfo_v%dim_size(1) = imax_ref
       dinfo_v%dim_size(2) = jmax_ref

       select case( trim(SNOPLGIN_vgridope_type) )

       case('ZLEV')

          select case( trim(zaxis_orgname) )
          case('z')

             dinfo_v%dim_name(3) = 'z'
             dinfo_v%dim_size(3) = kmax_new

             do j = 1, jmax_ref
             do i = 1, imax_ref
                call INTERP_interp1d( kmax_ref, 1, kmax_ref,    & ! [IN]
                                      kmax_new, 1, kmax_new,    & ! [IN]
                                      idx_Z(:,:),               & ! [IN]
                                      Zfact(:),                 & ! [IN]
                                      Z_ref(:),                 & ! [IN]
                                      ainfo_v(znum)%AXIS_1d(:), & ! [IN]
                                      dinfo%VAR_3d  (:,i,j),    & ! [IN]
                                      dinfo_v%VAR_3d(:,i,j),    & ! [OUT]
                                      logwgt = .false.          ) ! [IN]
             end do
             end do

          case('zh')

             dinfo_v%dim_name(3) = 'z'
             dinfo_v%dim_size(3) = kmax_new

             do j = 1, jmax_ref
             do i = 1, imax_ref
                call INTERP_interp1d( kmax_ref+1, 1, kmax_ref+1, & ! [IN]
                                      kmax_new,   1, kmax_new,   & ! [IN]
                                      idx_Zh(:,:),               & ! [IN]
                                      Zhfact(:),                 & ! [IN]
                                      Zh_ref(:),                 & ! [IN]
                                      ainfo_v(znum)%AXIS_1d(:),  & ! [IN]
                                      dinfo%VAR_3d  (:,i,j),     & ! [IN]
                                      dinfo_v%VAR_3d(:,i,j),     & ! [OUT]
                                      logwgt = .false.           ) ! [IN]
             end do
             end do

          case default

             dinfo_v%dim_name(3) = zaxis_orgname
             dinfo_v%dim_size(3) = zaxis_orgsize

             dinfo_v%VAR_3d(:,:,:) = dinfo%VAR_3d(:,:,:)

          end select

       case('PLEV')

          select case( trim(zaxis_orgname) )
          case('z')

             ! make geopotential height insted of PRES
             if( trim(dinfo_v%varname) == 'PRES' ) then

                dinfo_v%varname     = 'GPH'
                dinfo_v%description = 'geopotential height'
                dinfo_v%units       = 'm'
                dinfo_v%datatype    = dinfo%datatype
                dinfo_v%dim_rank    = dinfo%dim_rank
                dinfo_v%transpose   = .true.
                dinfo_v%step_nmax   = dinfo%step_nmax
                do t = 1, dinfo%step_nmax
                   dinfo_v%time_start(t) = dinfo%time_start(t)
                   dinfo_v%time_end  (t) = dinfo%time_end  (t)
                enddo
                dinfo_v%dt          = dinfo%dt
                dinfo_v%time_units  = dinfo%time_units
                dinfo_v%dim_name(3) = 'pressure'
                dinfo_v%dim_size(3) = kmax_new

                do j = 1, jmax_ref
                do i = 1, imax_ref
                   call INTERP_interp1d( kmax_ref, 1, kmax_ref, & ! [IN]
                                         kmax_new, 1, kmax_new, & ! [IN]
                                         idx_P(:,:,i,j),        & ! [IN]
                                         Pfact(:,i,j),          & ! [IN]
                                         PRES_ref(:,i,j),       & ! [IN]
                                         LnPaxis(:),            & ! [IN]
                                         height_ref(:,i,j),     & ! [IN]
                                         dinfo_v%VAR_3d(:,i,j), & ! [OUT]
                                         logwgt = .false.       ) ! [IN]
                end do
                end do

             else

                dinfo_v%dim_name(3) = 'pressure'
                dinfo_v%dim_size(3) = kmax_new

                do j = 1, jmax_ref
                do i = 1, imax_ref
                   call INTERP_interp1d( kmax_ref, 1, kmax_ref, & ! [IN]
                                         kmax_new, 1, kmax_new, & ! [IN]
                                         idx_P(:,:,i,j),        & ! [IN]
                                         Pfact(:,i,j),          & ! [IN]
                                         PRES_ref(:,i,j),       & ! [IN]
                                         LnPaxis(:),            & ! [IN]
                                         dinfo%VAR_3d  (:,i,j), & ! [IN]
                                         dinfo_v%VAR_3d(:,i,j), & ! [OUT]
                                         logwgt = .false.       ) ! [IN]
                end do
                end do

             endif

          case('zh')

             dinfo_v%dim_name(3) = 'pressure'
             dinfo_v%dim_size(3) = kmax_new

             do j = 1, jmax_ref
             do i = 1, imax_ref
                call INTERP_interp1d( kmax_ref+1, 1, kmax_ref+1, & ! [IN]
                                      kmax_new,   1, kmax_new,   & ! [IN]
                                      idx_Ph(:,:,i,j),           & ! [IN]
                                      Phfact(:,i,j),             & ! [IN]
                                      PRESh_ref(:,i,j),          & ! [IN]
                                      LnPaxis(:),                & ! [IN]
                                      dinfo%VAR_3d  (:,i,j),     & ! [IN]
                                      dinfo_v%VAR_3d(:,i,j),     & ! [OUT]
                                      logwgt = .false.           ) ! [IN]
             end do
             end do

          case default

             dinfo_v%dim_name(3) = zaxis_orgname
             dinfo_v%dim_size(3) = zaxis_orgsize

             dinfo_v%VAR_3d(:,:,:) = dinfo%VAR_3d(:,:,:)

          end select

       endselect

       do_output = .true.

    endif

    ! output

    if ( do_output ) then
       finalize    = ( nowstep == dinfo_v%step_nmax )
       add_rm_attr = .true.

       call SNO_vars_write( dirpath,                    & ! [IN] from namelist
                            basename,                   & ! [IN] from namelist
                            output_grads,               & ! [IN] from namelist
                            nowrank,                    & ! [IN]
                            nowstep,                    & ! [IN]
                            finalize,                   & ! [IN]
                            add_rm_attr,                & ! [IN]
                            nprocs_x_out, nprocs_y_out, & ! [IN] from namelist
                            nhalos_x,     nhalos_y,     & ! [IN] from SNO_file_getinfo
                            hinfo,                      & ! [IN] from SNO_file_getinfo
                            size(ainfo_v(:)),           & ! [IN] from SNO_file_getinfo
                            ainfo_v(:),                 & ! [IN] from SNO_axis_getinfo
                            dinfo_v,                    & ! [IN] from SNO_vars_getinfo
                            debug                       ) ! [IN]
    endif

    return
  end subroutine SNOPLGIN_vgridope_vinterp

end module mod_snoplugin_vgridope
