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
  integer,                private, allocatable :: idx_k (:,:)
  integer,                private, allocatable :: idx_kh(:,:)
  real(RP),               private, allocatable :: vfact (:)
  real(RP),               private, allocatable :: vhfact(:)

  real(RP),               private, allocatable :: Z_ref (:)
  real(RP),               private, allocatable :: Zh_ref(:)

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
       nvars,        &
       ainfo,        &
       dinfo,        &
       debug         )
    use scale_file_h, only: &
       FILE_REAL8
    use scale_prc, only: &
       PRC_abort
    use scale_interp, only: &
       INTERP_setup,    &
       INTERP_factor1d
    use scale_interp_vert, only: &
       INTERP_VERT_alloc_pres,   &
       INTERP_VERT_setcoef_pres
    use mod_sno_h, only: &
       axisinfo
    implicit none

    integer,        intent(in)  :: ngrids_z_out                          ! number of z-axis grids per process (output,sometimes including halo)
    integer,        intent(in)  :: ngrids_x_out                          ! number of x-axis grids per process (output,sometimes including halo)
    integer,        intent(in)  :: ngrids_y_out                          ! number of y-axis grids per process (output,sometimes including halo)
    integer,        intent(in)  :: naxis                                 ! number of axis variables           (input)
    integer,        intent(in)  :: nvars                                 ! number of item variables           (input)
    type(axisinfo), intent(in)  :: ainfo   (naxis)                       ! axis information                   (input)
    type(iteminfo), intent(in)  :: dinfo   (nvars)                       ! item information                   (input)
    logical,        intent(in)  :: debug

    real(RP), allocatable :: PRES_ref    (:,:,:)
    real(RP), allocatable :: PRESh_ref   (:,:,:)
    real(RP), allocatable :: SFC_PRES_ref(:,:)

    integer  :: i, j, k, n
    !---------------------------------------------------------------------------

    ! set region size
    kmax_ref = ngrids_z_out
    imax_ref = ngrids_x_out
    jmax_ref = ngrids_y_out

    kmax_new = SNOPLGIN_vgridope_lev_num

    ! set basic axis
    allocate( ainfo_v( naxis ) )
    ainfo_v(:) = ainfo(:)

    select case(SNOPLGIN_vgridope_type)
    case('ZLEV')

       LOG_NEWLINE
       LOG_INFO("SNOPLGIN_vgridope_setcoef",*) 'Setup remapping coefficient (height)'

       allocate( Z_ref (  kmax_ref) )
       allocate( Zh_ref(0:kmax_ref) )

       allocate( idx_k (kmax_new,2) )
       allocate( idx_kh(kmax_new,2) )
       allocate( vfact (kmax_new  ) )
       allocate( vhfact(kmax_new  ) )

       do n = 1, naxis
          if    ( ainfo(n)%varname == 'z'  ) then
             Z_ref (:) = ainfo(n)%AXIS_1d(:)

             znum = n
             ! rewrite axis
             deallocate( ainfo_v(znum)%AXIS_1d )
             ainfo_v(znum)%dim_size(1) = kmax_new
             allocate( ainfo_v(znum)%AXIS_1d(kmax_new) )

             do i = 1, kmax_new
                ainfo_v(znum)%AXIS_1d(i) = SNOPLGIN_vgridope_lev_data(i)
             enddo
          elseif( ainfo(n)%varname == 'zh' ) then
             Zh_ref(:) = ainfo(n)%AXIS_1d(:)
          endif
       enddo

       ! set remapping coefficient
       call INTERP_setup( 2 ) ! [IN] not used

       call INTERP_factor1d( kmax_ref, 1, kmax_ref,    & ! [IN]
                             kmax_new, 1, kmax_new,    & ! [IN]
                             Z_ref(:),                 & ! [IN]
                             ainfo_v(znum)%AXIS_1d(:), & ! [IN]
                             idx_k(:,:),               & ! [OUT]
                             vfact(:),                 & ! [OUT]
                             flag_extrap = .false.     ) ! [IN]

       call INTERP_factor1d( kmax_ref+1, 1, kmax_ref+1, & ! [IN]
                             kmax_new,   1, kmax_new,   & ! [IN]
                             Zh_ref(:),                 & ! [IN]
                             ainfo_v(znum)%AXIS_1d(:),  & ! [IN]
                             idx_kh(:,:),               & ! [IN]
                             vhfact(:),                 & ! [IN]
                             flag_extrap = .false.      ) ! [IN]

    case('PLEV')

       LOG_NEWLINE
       LOG_INFO("SNOPLGIN_vgridope_setcoef",*) 'Setup remapping coefficient (pressure)'

       ! set new axis set
       ainfo_v(1)%varname     = 'p'
       ainfo_v(1)%description = 'pressure'
       ainfo_v(1)%units       = 'Pa'
       ainfo_v(1)%datatype    = FILE_REAL8
       ainfo_v(1)%dim_rank    = 1
       ainfo_v(1)%dim_name(1) = 'pressure'
       ainfo_v(1)%transpose   = .false.
       ainfo_v(1)%regrid      = .false.

       ainfo_v(1)%dim_size(1) = kmax_new
       allocate( ainfo_v(1)%AXIS_1d(kmax_new) )

       do i = 1, kmax_new
          ainfo_v(1)%AXIS_1d(i) = SNOPLGIN_vgridope_lev_data(i)
       enddo

       allocate( PRES_ref (  kmax_ref,imax_ref,jmax_ref) )
       allocate( PRESh_ref(0:kmax_ref,imax_ref,jmax_ref) )

       allocate( SFC_PRES_ref(imax_ref,jmax_ref) )

       do n = 1, nvars
          if    ( dinfo(n)%varname == 'PRES' ) then
             PRES_ref(:,:,:) = ainfo(n)%AXIS_3d(:,:,:)
          elseif( dinfo(n)%varname == 'SFC_PRES' ) then
             SFC_PRES_ref(:,:) = ainfo(n)%AXIS_2d(:,:)
          else
             LOG_ERROR("SNOPLGIN_vgridope_setcoef",*) 'Not found pressure data at vertical level and the surface. Check "PRES" and "SFC_PRES"!'
             call PRC_abort
          endif
       enddo

       ! set pressure at half-level
       PRESh_ref(0,:,:) = SFC_PRES_ref(:,:)
       do k = 1, kmax_ref-1
          PRESh_ref(k,:,:) = PRESh_ref(k-1,:,:) + ( PRES_ref(k,:,:) - PRES_ref(k+1,:,:) )
       end do
       PRESh_ref(kmax_ref,:,:) = PRESh_ref(kmax_ref-1,:,:) + ( PRES_ref(kmax_ref-1,:,:) - PRES_ref(kmax_ref,:,:) )

       ! set remapping coefficient
       call INTERP_VERT_alloc_pres( kmax_new, & ! [IN]
                                    kmax_ref, & ! [IN]
                                    imax_ref, & ! [IN]
                                    jmax_ref  ) ! [IN]

       call INTERP_VERT_setcoef_pres( kmax_new,              & ! [IN]
                                      kmax_ref, 1, kmax_ref, & ! [IN]
                                      imax_ref, 1, imax_ref, & ! [IN]
                                      jmax_ref, 1, jmax_ref, & ! [IN]
                                      PRES_ref (:,:,:),      & ! [IN]
                                      PRESh_ref(:,:,:),      & ! [IN]
                                      SFC_PRES_ref(:,:),     & ! [IN]
                                      ainfo_v(1)%AXIS_1d(:)  ) ! [IN]

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

       select case(trim(zaxis_orgname))
       case('z')
          allocate( dinfo_v%VAR_3d(kmax_new,     imax_ref,jmax_ref) )
       case('zh')
          allocate( dinfo_v%VAR_3d(kmax_new,     imax_ref,jmax_ref) )
       case default
          allocate( dinfo_v%VAR_3d(zaxis_orgsize,imax_ref,jmax_ref) )
       end select

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

       select case(trim(zaxis_orgname))
       case('z')
          dinfo_v%dim_name(3) = "z"
          dinfo_v%dim_size(3) = kmax_new

          do j = 1, jmax_ref
          do i = 1, imax_ref
             call INTERP_interp1d( kmax_ref, 1, kmax_ref,    & ! [IN]
                                   kmax_new, 1, kmax_new,    & ! [IN]
                                   idx_k(:,:),               & ! [IN]
                                   vfact(:),                 & ! [IN]
                                   Z_ref(:),                 & ! [IN]
                                   ainfo_v(znum)%AXIS_1d(:), & ! [IN]
                                   dinfo%VAR_3d  (:,i,j),    & ! [IN]
                                   dinfo_v%VAR_3d(:,i,j),    & ! [OUT]
                                   logwgt = .false.          ) ! [OUT]
          end do
          end do
       case('zh')
          dinfo_v%dim_name(3) = "z"
          dinfo_v%dim_size(3) = kmax_new

          do j = 1, jmax_ref
          do i = 1, imax_ref
             call INTERP_interp1d( kmax_ref+1, 1, kmax_ref+1, & ! [IN]
                                   kmax_new,   1, kmax_new,   & ! [IN]
                                   idx_kh(:,:),               & ! [IN]
                                   vhfact(:),                 & ! [IN]
                                   Zh_ref(:),                 & ! [IN]
                                   ainfo_v(znum)%AXIS_1d(:),  & ! [IN]
                                   dinfo%VAR_3d  (:,i,j),     & ! [IN]
                                   dinfo_v%VAR_3d(:,i,j),     & ! [OUT]
                                   logwgt = .false.           ) ! [OUT]
          end do
          end do
       case default
          dinfo_v%dim_name(3) = zaxis_orgname
          dinfo_v%dim_size(3) = zaxis_orgsize

          dinfo_v%VAR_3d(:,:,:) = dinfo%VAR_3d(:,:,:)
       end select

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
