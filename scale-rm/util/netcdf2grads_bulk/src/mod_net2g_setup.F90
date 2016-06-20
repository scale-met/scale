!-------------------------------------------------------------------------------------------
!> module NET2G setup
!!
!! @par Description
!!          Procedure setting module for post-process of scale
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-02-03 (R.Yoshida)  original
!!
!<
!-------------------------------------------------------------------------------------------
module mod_net2g_setup
  !-----------------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_net2g_vars
  use mod_net2g_error

  !-----------------------------------------------------------------------------------------
  implicit none
  private
  !++ included parameters
#include "inc_net2g.h"
  !-----------------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: set_vtype
  public :: set_atype
  public :: set_ctype
  public :: set_flag_bnd
  public :: set_index
  public :: set_index_grid
  public :: set_index_gathered
  public :: set_index_readbuf
  public :: set_index_readbuf_grid
  public :: set_index_netcdf
  public :: set_rank_manage
  public :: set_array_size
  public :: set_calender

  !-----------------------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------------------

  !> setting of flag boundary
  !-----------------------------------------------------------------------------------------
  subroutine set_vtype( &
      ndim,      & ! [in ]
      varname,   & ! [in ]
      vtype,     & ! [out]
      atype      ) ! [inout]
    implicit none

    integer, intent(in)         :: ndim
    character(CMID), intent(in) :: varname
    integer, intent(out)        :: vtype
    integer, intent(inout)      :: atype

    character(CMID)             :: vname
    !---------------------------------------------------------------------------

    vname=trim(varname)
    if(vname(1:6)=="height") vname="height"
    if(vname(1:3)=="lon")    vname="lon"
    if(vname(1:3)=="lat")    vname="lat"

    select case( trim(vname) )
    case ( "TRL_URB", "TBL_URB", "TGL_URB" )
       vtype = vt_urban
    case ( "LAND_TEMP", "LAND_WATER" )
       vtype = vt_land
    case ( "height" )
       vtype = vt_height
    case ( "topo", "lsmask", "lon", "lat" )
       vtype = vt_tpmsk
       Z_MERGE_OUT = .false.
    case default
       if( ndim == 4 ) then
          vtype = vt_3d
       elseif( ndim == 3 ) then
          vtype = vt_2d
          Z_MERGE_OUT = .false.
          ANALYSIS    = "SLICE"
          atype       = a_slice
       else
          call err_abort( 0, __LINE__, loc_setup )
       end if
    end select

    return
  end subroutine set_vtype

  !> setting of flag analysis
  !-----------------------------------------------------------------------------------------
  subroutine set_atype( &
      atype    ) ! [out]
    implicit none

    integer, intent(out) :: atype
    !---------------------------------------------------------------------------

    select case( trim(Z_LEV_TYPE) )
    case ( "ZLEV", "zlev" )
       atype = a_conv

    case ( "PLEV", "plev" )
       atype = a_conv

    case ( "ORIGINAL", "original" )
       atype = a_slice
       if ( LOUT ) write( FID_LOG, '(1x,A)' ) "+++ ANALYSYS TYPE: slice"

    case ( "ANAL", "anal" )
       select case( trim(ANALYSIS) )
       case ( "MAX", "max", "MAXIMUM", "maximum" )
          if ( LOUT ) write( FID_LOG, '(1x,A)' ) "+++ ANALYSYS TYPE: column max"
          atype = a_max
          ZCOUNT  = 1
          Z_MERGE_OUT = .false.
       case ( "MIN", "min", "MINIMUM", "minimum" )
          if ( LOUT ) write( FID_LOG, '(1x,A)' ) "+++ ANALYSYS TYPE: column min"
          atype = a_min
          ZCOUNT  = 1
          Z_MERGE_OUT = .false.
       case ( "SUM", "sum", "SUMMATION", "summation" )
          if ( LOUT ) write( FID_LOG, '(1x,A)' ) "+++ ANALYSYS TYPE: column sum"
          atype = a_sum
          ZCOUNT  = 1
          Z_MERGE_OUT = .false.
       case ( "AVE", "ave", "AVERAGE", "average" )
          if ( LOUT ) write( FID_LOG, '(1x,A)' ) "+++ ANALYSYS TYPE: column ave"
          atype = a_ave
          ZCOUNT  = 1
          Z_MERGE_OUT = .false.
       case default
          write (*, *) "ERROR: specified analysis type is not appropiate"
          write (*, *) "***** ", trim(ANALYSIS)
          call err_abort( 1, __LINE__, loc_setup )
       end select

    case default
       write (*, *) "ERROR: specified Z_LEV_TYPE is not appropiate"
       write (*, *) "***** ", trim(Z_LEV_TYPE)
       call err_abort( 1, __LINE__, loc_setup )

    end select

    return
  end subroutine set_atype

  !> setting of flag level convert
  !-----------------------------------------------------------------------------------------
  subroutine set_ctype( &
      ctype    ) ! [out]
    implicit none

    integer, intent(out) :: ctype
    !---------------------------------------------------------------------------

    select case( trim(Z_LEV_TYPE) )
    case ( "ZLEV", "zlev" )
       if ( LOUT ) write( FID_LOG, '(1x,A)' ) "+++ reference var: height"
       ctype = c_height
    case ( "PLEV", "plev" )
       if ( LOUT ) write( FID_LOG, '(1x,A)' ) "+++ reference var: PRES"
       ctype = c_pres
    case default
       ctype = -1
       if ( LOUT ) write( FID_LOG, '(1x,A)' ) "+++ No Vertical Interpolation: "
    end select

    return
  end subroutine set_ctype

  !> setting of flag boundary
  !-----------------------------------------------------------------------------------------
  subroutine set_flag_bnd( &
      ix, jy,    & ! [in ]
      flag_bnd   ) ! [out]
    implicit none

    integer, intent(in)  :: ix, jy
    logical, intent(out) :: flag_bnd

    integer :: xproc, yproc
    !---------------------------------------------------------------------------

    xproc = PRC_NUM_X
    yproc = PRC_NUM_Y

    flag_bnd = .false.

    if ( HIST_BND ) then
       if ( ix == 1 .or. ix == xproc ) then
          flag_bnd = .true.
       endif
       if ( jy == 1 .or. jy == yproc ) then
          flag_bnd = .true.
       endif
    endif

    return
  end subroutine set_flag_bnd

  !> setting of indices
  !-----------------------------------------------------------------------------------------
  subroutine set_index( &
      ix,  jy,      & ! [in ]
      nxp, nyp,     & ! [in ]
      is,  ie,      & ! [out]
      js,  je       ) ! [out]
    implicit none

    integer, intent(in)  :: ix, jy
    integer, intent(in)  :: nxp,  nyp        ! partial data size
    integer, intent(out) :: is, ie
    integer, intent(out) :: js, je

    integer :: xproc, yproc
    !---------------------------------------------------------------------------

    xproc = PRC_NUM_X
    yproc = PRC_NUM_Y

    if ( HIST_BND ) then
       if ( ix == 1 ) then
          is = 1
          ie = nxp
       elseif ( ix == xproc ) then
          is = (ix-1)*(nxp-2)+2 + 1
          ie = (ix-1)*(nxp-2)+2 + nxp
       else
          is = (ix-1)*nxp+2 + 1
          ie = (ix-1)*nxp+2 + nxp
       endif

       if ( jy == 1 ) then
          js = 1
          je = nyp
       elseif ( jy == yproc ) then
          js = (jy-1)*(nyp-2)+2 + 1
          je = (jy-1)*(nyp-2)+2 + nyp
       else
          js = (jy-1)*nyp+2 + 1
          je = (jy-1)*nyp+2 + nyp
       endif
    else
       is = (ix-1)*nxp + 1
       ie = (ix-1)*nxp + nxp
       js = (jy-1)*nyp + 1
       je = (jy-1)*nyp + nyp
    endif

    return
  end subroutine set_index

  !> setting of indices for grid data
  !-----------------------------------------------------------------------------------------
  subroutine set_index_grid( &
      ix,   jy,      & ! [in ]
      nxgp, nygp,    & ! [in ]
      is,   ie,      & ! [out]
      js,   je       ) ! [out]
    implicit none

    integer, intent(in)  :: ix, jy
    integer, intent(in)  :: nxgp, nygp       ! partial data size for grids
    integer, intent(out) :: is, ie
    integer, intent(out) :: js, je

    integer :: xproc, yproc
    !---------------------------------------------------------------------------

    xproc = PRC_NUM_X
    yproc = PRC_NUM_Y

    is = (ix-1)*nxgp + 1
    ie = (ix-1)*nxgp + nxgp
    js = (jy-1)*nygp*xproc + 1
    je = (jy-1)*nygp*xproc + nygp

    return
  end subroutine set_index_grid

  !> setting of indices for grid data
  !-----------------------------------------------------------------------------------------
  subroutine set_index_gathered( &
      ix,   jy,      & ! [in ]
      mnxp, mnyp,    & ! [in ]
      is,   ie,      & ! [out]
      js,   je       ) ! [out]
    implicit none

    integer, intent(in)  :: ix, jy
    integer, intent(in)  :: mnxp, mnyp       ! maximum partial data size
    integer, intent(out) :: is, ie
    integer, intent(out) :: js, je

    integer :: xproc, yproc
    !---------------------------------------------------------------------------

    xproc = PRC_NUM_X
    yproc = PRC_NUM_Y

    if ( HIST_BND ) then
       if ( ix == 1 .or. ix == xproc ) then
          is = 1
          ie = mnxp
       else
          is = 1
          ie = mnxp - 2
       endif

       if ( jy == 1 ) then
          js = 1
          je = mnyp
       elseif ( jy == yproc ) then
          js = (ix-1)*mnyp + (jy-1)*mnyp*xproc + 1
          je = (ix-1)*mnyp + (jy-1)*mnyp*xproc + mnxp
       else
          js = (ix-1)*mnyp + (jy-1)*mnyp*xproc + 1
          je = (ix-1)*mnyp + (jy-1)*mnyp*xproc + (mnxp-2)
       endif
    else
       is = 1
       ie = (mnxp-2)

       js = (ix-1)*mnyp + (jy-1)*mnyp*xproc + 1
       je = (ix-1)*mnyp + (jy-1)*mnyp*xproc + (mnyp-2)
    endif

    return
  end subroutine set_index_gathered

  !> setting of indices for reading buffer
  !-----------------------------------------------------------------------------------------
  subroutine set_index_readbuf( &
      nm,          & ! [in ]
      nxp,  nyp,   & ! [in ]
      mnxp, mnyp,  & ! [in ]
      is,   ie,    & ! [out]
      js,   je,    & ! [out]
      im_bnd       ) ! [in ] optional
    implicit none

    integer, intent(in)  :: nm               ! num of loop for manage ranks
    integer, intent(in)  :: nxp,  nyp        ! partial data size
    integer, intent(in)  :: mnxp, mnyp       ! maximum partial data size
    integer, intent(out) :: is, ie           ! start index, end index
    integer, intent(out) :: js, je           ! start index, end index
    logical, intent(in), optional :: im_bnd  ! flag of boundary (edge tile)

    logical :: not_bnd
    !---------------------------------------------------------------------------

    not_bnd = .true.
    if ( present(im_bnd) ) then
       if ( im_bnd ) not_bnd = .false.
    endif

    if ( HIST_BND ) then
       if ( not_bnd ) then
          is = 1
          ie = nxp
          js = (nm-1)*mnyp + 1
          je = (nm-1)*mnyp + nyp
       else
          is = 1
          ie = mnxp
          js = (nm-1)*mnyp + 1
          je = (nm-1)*mnyp + mnyp
       endif
    else
       is = 1
       ie = nxp
       js = (nm-1)*nyp + 1
       je = (nm-1)*nyp + nyp
    endif

    return
  end subroutine set_index_readbuf


  !> setting of indices for reading buffer for grid
  !-----------------------------------------------------------------------------------------
  subroutine set_index_readbuf_grid( &
      nm,          & ! [in ]
      nxgp, nygp,  & ! [in ]
      is,   ie,    & ! [out]
      js,   je     ) ! [out]
    implicit none

    integer, intent(in)  :: nm               ! num of loop for manage ranks
    integer, intent(in)  :: nxgp, nygp       ! partial data size for grids
    integer, intent(out) :: is, ie           ! start index, end index
    integer, intent(out) :: js, je           ! start index, end index
    !---------------------------------------------------------------------------

    is = (nm-1)*nxgp + 1
    ie = (nm-1)*nxgp + nxgp
    js = (nm-1)*nygp + 1
    je = (nm-1)*nygp + nygp

    return
  end subroutine set_index_readbuf_grid

  !> setting of indices for reading buffer (netcdf) with counts
  !-----------------------------------------------------------------------------------------
  subroutine set_index_netcdf( &
      atype, vtype,  & ! [in ]
      ix,    jy,     & ! [in ]
      nxp,   nyp,    & ! [in ]
      mnxp,  mnyp,   & ! [in ]
      it,            & ! [in ]
      nz, zz,        & ! [in ]
      varname,       & ! [in ]
      is, ie,        & ! [out]
      js, je,        & ! [out]
      nzn,           & ! [out]
      start_3d,      & ! [out]
      start_2d,      & ! [out]
      start_2dt,     & ! [out]
      count_3d,      & ! [out]
      count_2d,      & ! [out]
      count_urban,   & ! [out]
      count_land,    & ! [out]
      count_height,  & ! [out]
      count_tpmsk    ) ! [out]
    implicit none

    integer, intent(in)  :: atype, vtype
    integer, intent(in)  :: ix, jy
    integer, intent(in)  :: nxp, nyp
    integer, intent(in)  :: mnxp, mnyp
    integer, intent(in)  :: it
    integer, intent(in)  :: nz(3), zz        !nz: 1=atom, 2=urban, 3=land
    character(CMID), intent(in) :: varname
    integer, intent(out) :: is, ie           ! start index, end index
    integer, intent(out) :: js, je           ! start index, end index
    integer, intent(out) :: nzn              ! number of z-dimension in netcdf
    integer, intent(out) :: start_3d(4)      ! start index for reading
    integer, intent(out) :: start_2d(3)      ! start index for reading
    integer, intent(out) :: start_2dt(2)     ! start index for reading
    integer, intent(out) :: count_3d(4)      ! data count for reading
    integer, intent(out) :: count_2d(3)      ! data count for reading
    integer, intent(out) :: count_urban(4)   ! data count for reading
    integer, intent(out) :: count_land(4)    ! data count for reading
    integer, intent(out) :: count_height(3)  ! data count for reading
    integer, intent(out) :: count_tpmsk(2)   ! data count for reading

    integer :: xproc, yproc
    !---------------------------------------------------------------------------

    xproc = PRC_NUM_X
    yproc = PRC_NUM_Y
    is = 1
    js = 1

    if ( HIST_BND ) then
       if ( ix == 1 .or. ix == xproc ) then
          ie = mnxp
       else
          ie = mnxp - 2
       endif
       if ( jy == 1 .or. jy == yproc ) then
          je = mnyp
       else
          je = mnyp - 2
       endif
    else
       ie = nxp
       je = nyp
    endif

    count_tpmsk(1:2) = (/ ie, je    /)
    start_2d   (1:3) = (/ 1,  1, it /)
    start_2dt  (1:2) = (/ 1,  1     /)

    select case( atype )
    case ( a_slice )
       count_3d    (1:4) = (/ ie, je, 1,  1  /)
       count_2d    (1:3) = (/ ie, je,     1  /)
       count_urban (1:4) = (/ ie, je, 1,  1  /)
       count_land  (1:4) = (/ ie, je, 1,  1  /)
       count_height(1:3) = (/ ie, je, 1      /)
       start_3d    (1:4) = (/ 1,  1,  zz, it /)
       nzn = 1
    case ( a_max, a_min, a_sum, a_ave, a_conv )
       count_3d    (1:4) = (/ ie, je, nz(1), 1  /)
       count_2d    (1:3) = (/ ie, je,        1  /)
       count_urban (1:4) = (/ ie, je, nz(2), 1  /)
       count_land  (1:4) = (/ ie, je, nz(3), 1  /)
       count_height(1:3) = (/ ie, je, nz(1)     /)
       start_3d    (1:4) = (/ 1,  1,  1,     it /)

       select case( vtype )
       case ( vt_urban )
          nzn = nz(2)
       case ( vt_land )
          nzn = nz(3)
       case ( vt_3d, vt_height )
          nzn = nz(1)
       case ( vt_2d, vt_tpmsk )
          write (*, *) "ERROR: specified anal-type is not appropiate for the var"
          write (*, *) "***** ", trim(ANALYSIS), " --> ", trim(varname)
          call err_abort( 1, __LINE__, loc_setup )
       case default
          call err_abort( 0, __LINE__, loc_setup )
       end select

    case default
       call err_abort( 0, __LINE__, loc_setup )
    end select

    return
  end subroutine set_index_netcdf

  !> setting of ranks should be managed in the process
  !-----------------------------------------------------------------------------------------
  subroutine set_rank_manage( &
      irank,     & ! [in ]
      nmnge,     & ! [in ]
      rk_mnge    ) ! [out]
    implicit none

    integer, intent(in)  :: irank
    integer, intent(in)  :: nmnge
    integer, intent(out) :: rk_mnge(:)

    integer :: n
    !---------------------------------------------------------------------------

    if ( LOUT ) write( FID_LOG, '(1X,A,I5)') &
                "+++ number of ranks to manage: ", nmnge
    do n = 1, nmnge
       rk_mnge(n) = irank * nmnge + (n-1)
       if ( LOUT ) write( FID_LOG, '(1X,A,I5,A,I5)') &
                   "+++ myrank: ", irank, " -->  manage: ", rk_mnge(n)
    enddo

    return
  end subroutine set_rank_manage

  !> setting of indices
  !-----------------------------------------------------------------------------------------
  subroutine set_array_size( &
      nxp,  nyp,     & ! [in ]
      nxgp, nygp,    & ! [in ]
      nx,   ny,      & ! [out]
      mnx,  mny,     & ! [out]
      mnxp, mnyp,    & ! [out]
      nxg_tproc,     & ! [out]
      nyg_tproc      ) ! [out]
    implicit none

    integer, intent(in)  :: nxp, nxgp
    integer, intent(in)  :: nyp, nygp
    integer, intent(out) :: nx, mnx, mnxp, nxg_tproc
    integer, intent(out) :: ny, mny, mnyp, nyg_tproc

    integer :: xproc, yproc, tproc
    !---------------------------------------------------------------------------

    xproc = PRC_NUM_X
    yproc = PRC_NUM_Y
    tproc = xproc * yproc

    if ( HIST_BND ) then
       mnxp = nxgp - 2
       mnyp = nygp - 2
       nx = (mnxp-2) * xproc + 4
       ny = (mnyp-2) * yproc + 4
    else
       mnxp = nxp
       mnyp = nyp
       nx = mnxp * xproc
       ny = mnyp * yproc
    endif

    mnx = mnxp * xproc
    mny = mnyp * yproc
    nxg_tproc = nxgp * tproc
    nyg_tproc = nygp * tproc

    if ( LOUT ) write( FID_LOG, '(1X,"+++ nx:",I7,2X,"ny:",I7)') nx, ny
    if ( LOUT ) write( FID_LOG, '(1X,A)') ""

    return
  end subroutine set_array_size

  !> setting of calender indices
  !-----------------------------------------------------------------------------------------
  subroutine set_calender( &
      yy,       & ! [in ]
      mm,       & ! [in ]
      dd,       & ! [in ]
      hh,       & ! [in ]
      mn,       & ! [in ]
      sc,       & ! [in ]
      STIME,    & ! [out]
      FTIME     ) ! [out]
    implicit none

    integer, intent(in) :: yy, mm, dd
    integer, intent(in) :: hh, mn, sc
    character(*), intent(out) :: STIME
    character(*), intent(out) :: FTIME

    character(2) :: csc, cmn, chh, cdd, cmm2
    character(4) :: cyy
    !---------------------------------------------------------------------------

    write(csc, '(I2.2)') sc
    write(cmn, '(I2.2)') mn
    write(chh, '(I2.2)') hh
    write(cdd, '(I2.2)') dd
    write(cmm2,'(I2.2)') mm
    write(cyy, '(I4.4)') yy
    STIME = chh//':'//cmn//'Z'//cdd//cmm(mm)//cyy
    FTIME = cyy//cmm2//cdd//chh//cmn//csc

    return
  end subroutine set_calender

end module mod_net2g_setup
