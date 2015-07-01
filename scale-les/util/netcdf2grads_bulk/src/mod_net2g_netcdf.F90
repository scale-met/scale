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
  public :: netcdf_read_ref

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

    if ( LOUT .and. LOG_DBUG ) write( FID_LOG, '(1x,A,A)' ) "+++ Target File (dim): ", trim(ncfile)
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
       write( FID_LOG, '(1x,A,A)' ) "ERROR: Such a variable was not found: ", trim(varname)
       write( *, '(1x,A,A)' ) "ERROR: Such a variable was not found: ", trim(varname)
       call handle_err(istat, __LINE__)
    endif

    istat = nf90_close(ncid)
    return
  end subroutine netcdf_var_dim

  !> allocate temporaly arraies
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
      ctype,      & ! [in ]
      vtype,      & ! [in ]
      p_var       ) ! [out]
    implicit none

    integer, intent(in) :: imnge, nm, it
    integer, intent(in) :: nxp, nyp, mnxp, mnyp
    integer, intent(in) :: zz, nz(3)
    character(CMID), intent(in)  :: varname
    integer,         intent(in)  :: atype, ctype, vtype
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
    real(SP) :: lev = 0.0
    character(CLNG) :: ncfile
    character(6)    :: num
    logical         :: logwgt = .false.

    integer :: ncid, varid
    integer :: istat
    logical :: flag_bnd = .false.
    !---------------------------------------------------------------------------

    write ( num,'(I6.6)' ) imnge
    ncfile = trim(IDIR)//"/"//trim(HISTORY_DEFAULT_BASENAME)//".pe"//num//".nc"
    if ( LOUT .and. LOG_DBUG ) write( FID_LOG, '(1x,A,A)' ) "+++ Target File (var): ", trim(ncfile)

    if ( atype == a_conv ) then
       if ( ctype == c_height ) then
          lev = real( zz )
       elseif ( ctype == c_pres ) then
          lev = real( zz ) * 100.D0
       endif

       select case ( trim(varname) )
       case("PRES","pres", "DENS","dens")
          logwgt = .true.
       case default
          logwgt = .false.
       end select
    endif

    call irank2ixjy( imnge, ix, jy )
    call set_flag_bnd( ix, jy, flag_bnd )
    call set_index_readbuf( nm, nxp, nyp, mnxp, mnyp,       &
                            is, ie, js, je, im_bnd=flag_bnd )
    call set_index_netcdf( atype, vtype, ix, jy, nxp, nyp,       &
                           mnxp, mnyp, it, nz, zz, varname,      &
                           isn, ien, jsn, jen, nzn,              &
                           start_3d, start_2d, start_2dt,        &
                           count_3d, count_2d, count_urban,      &
                           count_land, count_height, count_tpmsk )

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, trim(varname), varid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    select case( vtype )
    case ( vt_urban )
       istat = nf90_get_var( ncid, varid, p_3d_urban(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_urban )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       call anal_simple( atype, is,  ie,  js,  je,  &
                         isn, jsn, nzn, p_3d_urban, p_var )

    case ( vt_land )
       istat = nf90_get_var( ncid, varid, p_3d_land(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_land )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       call anal_simple( atype, is,  ie,  js,  je,  &
                         isn, jsn, nzn, p_3d_land, p_var )

!    case ( vt_height )
!       istat = nf90_get_var( ncid, varid, p_2d(isn:ien,jsn:jen,1), start=start_2d, count=count_height )
!       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
!       p_var(is:ie,js:je) = real( p_2d(isn:ien,jsn:jen,1) )

    case ( vt_tpmsk )
       istat = nf90_get_var( ncid, varid, p_2dt(isn:ien,jsn:jen), start=start_2dt, count=count_tpmsk )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(is:ie,js:je) = real( p_2dt(isn:ien,jsn:jen) )

    case ( vt_3d, vt_height )
       if ( vtype == vt_3d ) then
          istat = nf90_get_var( ncid, varid, p_3d(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_3d )
       elseif ( vtype == vt_height ) then
          istat = nf90_get_var( ncid, varid, p_3d(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_height )
       endif
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

       if ( atype == a_conv ) then
          call anal_ref_interp( ctype, lev, p_3d, logwgt, is, ie, js, je, &
                                isn, ien, jsn, jen, nzn, nm, p_var )
       else
          call anal_simple( atype, is,  ie,  js,  je,  &
                            isn, jsn, nzn, p_3d, p_var )
       endif

    case ( vt_2d )
       istat = nf90_get_var( ncid, varid, p_2d(isn:ien,jsn:jen,1), start=start_2d, count=count_2d )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(is:ie,js:je) = real( p_2d(isn:ien,jsn:jen,1) )

    case default
       call err_abort( 0, __LINE__, loc_netcdf )

    end select

    istat = nf90_close(ncid)
    return
  end subroutine netcdf_read_var

  !> read from netcdf: reference data for convert
  !-----------------------------------------------------------------------------------------
  subroutine netcdf_read_ref( &
      imnge,      & ! [in ]
      nxp,  nyp,  & ! [in ]
      mnxp, mnyp, & ! [in ]
      it,         & ! [in ]
      nz,         & ! [in ]
      ctype,      & ! [in ]
      p_var       ) ! [out]
    implicit none

    integer, intent(in) :: imnge
    integer, intent(in) :: nxp, nyp, mnxp, mnyp
    integer, intent(in) :: it
    integer, intent(in) :: nz(3)
    integer,         intent(in)  :: ctype
    real(SP),        intent(out) :: p_var(:,:,:)

    integer :: start_3d(4)
    integer :: start_2d(3)
    integer :: count_3d(4)
    integer :: count_2d(3)

    integer :: isn, ien, jsn, jen, nzn
    integer :: ix, jy
    character(CLNG) :: ncfile
    character(6)    :: num

    integer :: xproc, yproc
    integer :: ncid, varid
    integer :: istat
    logical :: flag_bnd = .false.
    !---------------------------------------------------------------------------

    xproc = PRC_NUM_X
    yproc = PRC_NUM_Y

    write ( num,'(I6.6)' ) imnge
    ncfile = trim(IDIR)//"/"//trim(HISTORY_DEFAULT_BASENAME)//".pe"//num//".nc"
    if ( LOUT .and. LOG_DBUG ) write( FID_LOG, '(1x,A,A)' ) "+++ Target File (ref): ", trim(ncfile)

    call irank2ixjy( imnge, ix, jy )

    nzn = nz(1)
    isn = 1
    jsn = 1
    if ( flag_bnd ) then
       if ( ix == 1 .or. ix == xproc ) then
          ien = mnxp
       else
          ien = mnxp - 2
       endif
       if ( jy == 1 .or. jy == yproc ) then
          jen = mnyp
       else
          jen = mnyp - 2
       endif
    else
       ien = nxp
       jen = nyp
    endif

    start_2d(1:3) = (/   1,   1,   1      /)
    start_3d(1:4) = (/   1,   1,   1,  it /)
    count_2d(1:3) = (/ ien, jen, nzn      /)
    count_3d(1:4) = (/ ien, jen, nzn,  1  /)

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    select case( ctype )
    case ( c_height )
       istat = nf90_inq_varid( ncid, 'height', varid )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       istat = nf90_get_var( ncid, varid, p_3d(isn:ien,jsn:jen,1:nzn,1), start=start_2d, count=count_2d )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(:,:,:) = real( p_3d(isn:ien,jsn:jen,1:nzn,1) )

    case ( c_pres )
       istat = nf90_inq_varid( ncid, 'PRES', varid )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       istat = nf90_get_var( ncid, varid, p_3d(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_3d )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(:,:,:) = real( p_3d(isn:ien,jsn:jen,1:nzn,1) )

    case default
       call err_abort( 0, __LINE__, loc_netcdf )

    end select

    istat = nf90_close(ncid)
    return
  end subroutine netcdf_read_ref

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
