!-------------------------------------------------------------------------------------------
!> module NET2G netcdf
!!
!! @par Description
!!          netCDF I/O module for post-process of scale
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-02-03 (R.Yoshida)  original
!!
!<
!-------------------------------------------------------------------------------------------
module mod_net2g_netcdf
  !-----------------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use netcdf

  use mod_net2g_vars
  use mod_net2g_error
  use mod_net2g_anal
  use mod_net2g_setup, only:  &
      set_flag_bnd,           &
      set_index_readbuf,      &
      set_index_readbuf_grid, &
      set_index_netcdf

  !-----------------------------------------------------------------------------------------
  implicit none
  private
  !++ included parameters
#include "inc_net2g.h"
  !-----------------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: netcdf_setup
  public :: netcdf_var_dim
  public :: netcdf_retrieve_dims
  public :: netcdf_read_grid
  public :: netcdf_read_var

  !-----------------------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: irank2ixjy
  !-----------------------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(DP),allocatable, private :: p_3d(:,:,:,:)         ! partial grid data
  real(DP),allocatable, private :: p_3d_urban(:,:,:,:)   ! partial grid data
  real(DP),allocatable, private :: p_3d_land(:,:,:,:)    ! partial grid data
  real(DP),allocatable, private :: p_2d(:,:,:)           ! partial grid data
  real(DP),allocatable, private :: p_2dt(:,:)            ! partial grid data

  !-----------------------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------------------

  !> read from netcdf: retrieve dimension size
  !-----------------------------------------------------------------------------------------
  subroutine netcdf_retrieve_dims( &
      ncfile,          & ! [in ]
      nxp, nxgp, nxh,  & ! [out]
      nyp, nygp, nyh,  & ! [out]
      nz,  nzg,  nzh,  & ! [out]
      uz, lz, ks, ke   ) ! [out]
    implicit none

    character(CLNG), intent(in) :: ncfile
    integer, intent(out) :: nxp, nxgp, nxh
    integer, intent(out) :: nyp, nygp, nyh
    integer, intent(out) :: nz,  nzg,  nzh
    integer, intent(out) :: uz,  lz
    integer, intent(out) :: ks,  ke

    integer :: ncid, dimid
    integer :: istat
    !---------------------------------------------------------------------------

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
 
    istat = nf90_inq_dimid ( ncid,  'x',  dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nxp )
    istat = nf90_inq_dimid ( ncid,  'CX', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nxgp )
    nxh = ( nxgp - nxp ) / 2

    istat = nf90_inq_dimid ( ncid,  'y',  dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nyp )
    istat = nf90_inq_dimid ( ncid,  'CY', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nygp )
    nyh = ( nygp - nyp ) / 2

    istat = nf90_inq_dimid ( ncid,  'z',  dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nz )
    istat = nf90_inq_dimid ( ncid,  'CZ', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nzg )
    nzh = ( nzg - nz ) / 2
    ks = 1
    ke = nz

    istat = nf90_inq_dimid ( ncid,  'uz', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=uz )

    istat = nf90_inq_dimid ( ncid,  'lz', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=lz )

    istat = nf90_close(ncid)
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    return
  end subroutine netcdf_retrieve_dims

  !> setting of flag boundary
  !-----------------------------------------------------------------------------------------
  subroutine netcdf_var_dim( &
      ncfile,    & ! [in ]
      varname,   & ! [in ]
      ndim       ) ! [out]
    implicit none

    character(CLNG), intent(in) :: ncfile
    character(CMID), intent(in) :: varname
    integer, intent(out)        :: ndim


    integer :: ncid, varid
    integer :: istat
    !---------------------------------------------------------------------------

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
 
    istat = nf90_inq_varid( ncid, trim(varname), varid )
    istat = nf90_inquire_variable(ncid, varid, ndims=ndim )
    if (istat .ne. nf90_noerr) then
       if ( LOUT ) write( FID_LOG, '(1x,A,A)' ) "ERROR: Such a variable was not found: ", trim(varname)
       call handle_err(istat, __LINE__)
    endif

    istat = nf90_close(ncid)
    return
  end subroutine netcdf_var_dim

  !> read from netcdf: retrieve dimension size
  !-----------------------------------------------------------------------------------------
  subroutine netcdf_setup( &
      mnxp,   & ! [in]
      mnyp,   & ! [in]
      nz      ) ! [in]
    implicit none

    integer, intent(in) :: mnxp, mnyp
    integer, intent(in) :: nz(3)
    !---------------------------------------------------------------------------

    allocate( p_2dt       (mnxp, mnyp          ) )
    allocate( p_2d        (mnxp, mnyp, 1       ) )
    allocate( p_3d        (mnxp, mnyp, nz(1), 1) )
    allocate( p_3d_urban  (mnxp, mnyp, nz(2), 1) )
    allocate( p_3d_land   (mnxp, mnyp, nz(3), 1) )

    return
  end subroutine netcdf_setup

  !> read from netcdf: grid data
  !-----------------------------------------------------------------------------------------
  subroutine netcdf_read_grid( &
      imnge,        & ! [in ]
      nm,           & ! [in ]
      nxgp, nygp,   & ! [in ]
      zlev,         & ! [out]
      cz, cdz,      & ! [out]
      p_cx, p_cdx,  & ! [out]
      p_cy, p_cdy   ) ! [out]
    implicit none

    integer, intent(in) :: imnge
    integer, intent(in) :: nm
    integer, intent(in) :: nxgp, nygp
    real(SP), intent(out) :: zlev(:)
    real(SP), intent(out) :: cz(:),   cdz(:)
    real(SP), intent(out) :: p_cx(:), p_cdx(:)
    real(SP), intent(out) :: p_cy(:), p_cdy(:)

    integer :: ncid, varid
    integer :: istat
    integer :: is, ie, js, je
    character(CLNG) :: ncfile
    character(6)    :: num
    !---------------------------------------------------------------------------

    write ( num,'(I6.6)' ) imnge
    ncfile = trim(IDIR)//"/"//trim(HISTORY_DEFAULT_BASENAME)//".pe"//num//".nc"
    if ( LOUT ) write( FID_LOG, '(1x,A,A)' ) "+++ Target File (grd): ", trim(ncfile)

    call set_index_readbuf_grid( nm, nxgp, nygp, is, ie, js, je )

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
 
    istat = nf90_inq_varid( ncid, 'z', varid )
    istat = nf90_get_var( ncid, varid, zlev )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CZ', varid )
    istat = nf90_get_var( ncid, varid, cz )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CDZ', varid )
    istat = nf90_get_var( ncid, varid, cdz )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CX', varid )
    istat = nf90_get_var( ncid, varid, p_cx(is:ie) )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CDX', varid )
    istat = nf90_get_var( ncid, varid, p_cdx(is:ie) )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CY', varid )
    istat = nf90_get_var( ncid, varid, p_cy(js:je) )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CDY', varid )
    istat = nf90_get_var( ncid, varid, p_cdy(js:je) )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_close(ncid)
    return
  end subroutine netcdf_read_grid


  !> read from netcdf: variables data
  !-----------------------------------------------------------------------------------------
  subroutine netcdf_read_var( &
      imnge,      & ! [in ]
      nm,         & ! [in ]
      nxp,  nyp,  & ! [in ]
      mnxp, mnyp, & ! [in ]
      it,         & ! [in ]
      zz,         & ! [in ]
      nz,         & ! [in ]
      varname,    & ! [in ]
      atype,      & ! [in ]
      vtype,      & ! [in ]
      p_var       ) ! [out]
    implicit none

    integer, intent(in) :: imnge, nm, it
    integer, intent(in) :: nxp, nyp, mnxp, mnyp
    integer, intent(in) :: zz, nz(3)
    character(CMID), intent(in)  :: varname
    integer,         intent(in)  :: atype, vtype
    real(SP),        intent(out) :: p_var(:,:)

    integer :: start_3d(4)
    integer :: start_2d(3)
    integer :: start_2dt(2)
    integer :: count_3d(4)
    integer :: count_2d(3)
    integer :: count_urban(4)
    integer :: count_land(4)
    integer :: count_height(3)
    integer :: count_tpmsk(2)

    integer :: is, ie, js, je
    integer :: isn, ien, jsn, jen, nzn
    integer :: ix, jy
    character(CLNG) :: ncfile
    character(6)    :: num

    integer :: ncid, varid
    integer :: istat
    logical :: flag_bnd = .false.
    !---------------------------------------------------------------------------

    write ( num,'(I6.6)' ) imnge
    ncfile = trim(IDIR)//"/"//trim(HISTORY_DEFAULT_BASENAME)//".pe"//num//".nc"
    if ( LOUT ) write( FID_LOG, '(1x,A,A)' ) "+++ Target File (var): ", trim(ncfile)

    call irank2ixjy( imnge, ix, jy )
    call set_flag_bnd( ix, jy, flag_bnd )
    call set_index_readbuf( nm, nxp, nyp, mnxp, mnyp,       &
                            is, ie, js, je, im_bnd=flag_bnd )
    call set_index_netcdf( atype, vtype, ix, jy, nxp, nyp,       & ! [in]
                           mnxp, mnyp, it, nz, zz, varname,      & ! [in]
                           isn, ien, jsn, jen, nzn,              & ! [out]
                           start_3d, start_2d, start_2dt,        & ! [out]
                           count_3d, count_2d, count_urban,      & ! [out]
                           count_land, count_height, count_tpmsk ) ! [out]

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
 
    istat = nf90_inq_varid( ncid, trim(varname), varid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    select case( vtype )
    case ( vt_urban )
       istat = nf90_get_var( ncid, varid, p_3d_urban(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_urban )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       call simple_analysis( atype, is,  ie,  js,  je,  &
                             isn, jsn, nzn, p_3d_urban, p_var )

    case ( vt_land )
       istat = nf90_get_var( ncid, varid, p_3d_land(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_land )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       call simple_analysis( atype, is,  ie,  js,  je,  &
                             isn, jsn, nzn, p_3d_land, p_var )

    case ( vt_height )
       istat = nf90_get_var( ncid, varid, p_2d(isn:ien,jsn:jen,1), start=start_2d, count=count_height )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(is:ie,js:je) = real( p_2d(isn:ien,jsn:jen,1) )

    case ( vt_tpmsk )
       istat = nf90_get_var( ncid, varid, p_2dt(isn:ien,jsn:jen), start=start_2dt, count=count_tpmsk )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(is:ie,js:je) = real( p_2dt(isn:ien,jsn:jen) )

    case ( vt_3d )
       istat = nf90_get_var( ncid, varid, p_3d(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_3d )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       call simple_analysis( atype, is,  ie,  js,  je,  &
                             isn, jsn, nzn, p_3d, p_var )

    case ( vt_2d )
       istat = nf90_get_var( ncid, varid, p_2d(isn:ien,jsn:jen,1), start=start_2d, count=count_2d )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(is:ie,js:je) = real( p_2d(isn:ien,jsn:jen,1) )

    case default
       call err_abort( 0, __LINE__ )

    end select

    istat = nf90_close(ncid)
    return
  end subroutine netcdf_read_var


  !> assignment of tile position corresponding to the rank
  !-----------------------------------------------------------------------------------------
  subroutine irank2ixjy( &
      irank,    & ! [in ]
      ix,       & ! [out]
      jy        ) ! [out]
    implicit none

    integer, intent(in)  :: irank
    integer, intent(out) :: ix, jy

    integer :: xproc, yproc
    !---------------------------------------------------------------------------

    xproc = PRC_NUM_X
    yproc = PRC_NUM_Y

    ix = mod( irank, xproc )  + 1
    jy = int( irank / xproc ) + 1

    return
  end subroutine irank2ixjy

end module mod_net2g_netcdf
