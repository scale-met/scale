!-------------------------------------------------------------------------------
!> Module SNO (RM)
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
module mod_sno_grads
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SNO_grads_setup
  public :: SNO_grads_write
  public :: SNO_grads_write_ctl
  public :: SNO_grads_netcdfctl

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: SNO_grads_calc_timechar

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private              :: GRADS_grd_fid = -1

  real(SP), private, allocatable :: TMPDATA(:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine SNO_grads_setup( &
       nprocs_x_out,   &
       nprocs_y_out,   &
       output_grads,   &
       output_gradsctl )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in)  :: nprocs_x_out             ! x length of 2D processor topology (output)
    integer, intent(in)  :: nprocs_y_out             ! y length of 2D processor topology (output)
    logical, intent(in)  :: output_grads
    logical, intent(in)  :: output_gradsctl
    !---------------------------------------------------------------------------

    if ( output_grads .OR. output_gradsctl ) then
       LOG_NEWLINE
       LOG_INFO("SNO_grads_setup",*) 'Output file as GrADS format, instead of NeCDF'

       if ( nprocs_x_out * nprocs_y_out /= 1 ) then
          LOG_ERROR("SNO_grads_setup",*) 'To output file as GrADS format, the number of output file must be 1.'
          call PRC_abort
       endif
    endif

    return
  end subroutine SNO_grads_setup

  !-----------------------------------------------------------------------------
  subroutine SNO_grads_write( &
       dirpath,  &
       nowstep,  &
       finalize, &
       hinfo,    &
       naxis,    &
       ainfo,    &
       dinfo,    &
       debug     )
    use scale_prc, only: &
       PRC_abort
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo,   &
       iteminfo
    implicit none

    character(len=*), intent(in)  :: dirpath                               ! directory path                     (output)
    integer,          intent(in)  :: nowstep                               ! current step                       (output)
    logical,          intent(in)  :: finalize                              ! finalize in this step?
    type(commoninfo), intent(in)  :: hinfo                                 ! common information                 (input)
    integer,          intent(in)  :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(in)  :: ainfo(naxis)                          ! axis information                   (input)
    type(iteminfo),   intent(in)  :: dinfo                                 ! variable information               (input)
    logical,          intent(in)  :: debug

    character(len=H_LONG) :: grdname
    integer(8)            :: recsize

    integer  :: imax, jmax, kmax
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if ( dinfo%dim_rank == 2 ) then
       imax = size(dinfo%VAR_2d(:,:),1)
       jmax = size(dinfo%VAR_2d(:,:),2)
       kmax = 1
    elseif( dinfo%dim_rank == 3 ) then
       imax = size(dinfo%VAR_3d(:,:,:),2)
       jmax = size(dinfo%VAR_3d(:,:,:),3)
       kmax = size(dinfo%VAR_3d(:,:,:),1)
    endif

    ! Open data file
    if ( nowstep == 1 ) then
       if ( dirpath == '' ) then
          grdname = trim(dinfo%varname)//'.grd'
       else
          grdname = trim(dirpath)//'/'//trim(dinfo%varname)//'.grd'
       endif

       recsize = int(imax,kind=8) * int(jmax,kind=8) * int(kmax,kind=8) * 4_8

       LOG_INFO("SNO_grads_write",*) 'filename : ', trim(grdname)

       GRADS_grd_fid = IO_get_available_fid()
       open( unit   = GRADS_grd_fid, &
             file   = trim(grdname), &
             form   = 'unformatted', &
             access = 'direct',      &
             recl   = recsize,       &
             status = 'unknown'      )

       allocate( TMPDATA(imax,jmax,kmax) )
    endif



    if ( dinfo%dim_rank == 2 ) then
       do j = 1, jmax
       do i = 1, imax
          TMPDATA(i,j,1) = real(dinfo%VAR_2d(i,j),kind=SP)
       enddo
       enddo
    elseif( dinfo%dim_rank == 3 ) then
       do k = 1, kmax
       do j = 1, jmax
       do i = 1, imax
          TMPDATA(i,j,k) = real(dinfo%VAR_3d(k,i,j),kind=SP)
       enddo
       enddo
       enddo
    endif

    write(GRADS_grd_fid,rec=nowstep) TMPDATA(:,:,:)



    ! Close data file
    if ( finalize ) then
       deallocate( TMPDATA )

       close(GRADS_grd_fid)

       call SNO_grads_write_ctl( dirpath, & ! [IN]
                                 hinfo,   & ! [IN]
                                 naxis,   & ! [IN]
                                 ainfo,   & ! [IN]
                                 dinfo,   & ! [IN]
                                 debug    ) ! [IN]
    endif

    return
  end subroutine SNO_grads_write

  !-----------------------------------------------------------------------------
  subroutine SNO_grads_write_ctl( &
       dirpath, &
       hinfo,   &
       naxis,   &
       ainfo,   &
       dinfo,   &
       debug    )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       CONST_UNDEF, &
       CONST_D2R,   &
       CONST_RADIUS
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo,   &
       iteminfo
    implicit none

    character(len=*), intent(in)  :: dirpath                               ! directory path                     (output)
    type(commoninfo), intent(in)  :: hinfo                                 ! common information                 (input)
    integer,          intent(in)  :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(in)  :: ainfo(naxis)                          ! axis information                   (input)
    type(iteminfo),   intent(in)  :: dinfo                                 ! variable information               (input)
    logical,          intent(in)  :: debug

    character(len=H_LONG)  :: grdname
    character(len=H_LONG)  :: ctlname

    character(len=H_SHORT) :: idim, jdim, kdim
    integer                :: imax, jmax, kmax
    integer                :: imax_, jmax_

    real(RP)               :: latstart, latend
    real(RP)               :: lonstart, lonend
    real(RP)               :: clat, clon
    real(RP)               :: dx, dy
    real(RP)               :: dlat, dlon

    character(len=20)      :: cdate, dhour

    logical  :: written
    integer  :: fid
    integer  :: k, n
    !---------------------------------------------------------------------------

    if ( dinfo%dim_rank == 2 ) then
       imax = size(dinfo%VAR_2d(:,:),1)
       jmax = size(dinfo%VAR_2d(:,:),2)
       kmax = 1

       idim = trim(dinfo%dim_name(1))
       jdim = trim(dinfo%dim_name(2))
       kdim = 'none'
    elseif( dinfo%dim_rank == 3 ) then
       imax = size(dinfo%VAR_3d(:,:,:),2)
       jmax = size(dinfo%VAR_3d(:,:,:),3)
       kmax = size(dinfo%VAR_3d(:,:,:),1)

       if ( dinfo%transpose ) then
          idim = trim(dinfo%dim_name(1))
          jdim = trim(dinfo%dim_name(2))
          kdim = trim(dinfo%dim_name(3))
       else
          idim = trim(dinfo%dim_name(2))
          jdim = trim(dinfo%dim_name(3))
          kdim = trim(dinfo%dim_name(1))
       endif
    endif

    grdname = trim(dinfo%varname)//'.grd'
    if ( dirpath == '' ) then
       ctlname = trim(dinfo%varname)//'.ctl'
    else
       ctlname = trim(dirpath)//'/'//trim(dinfo%varname)//'.ctl'
    endif

    do n = 1, naxis
       if    ( ainfo(n)%varname == 'lat' ) then
          latstart = minval( ainfo(n)%AXIS_2d(:,:) )
          latend   = maxval( ainfo(n)%AXIS_2d(:,:) )
       elseif( ainfo(n)%varname == 'lon' ) then
          lonstart = minval( ainfo(n)%AXIS_2d(:,:) )
          lonend   = maxval( ainfo(n)%AXIS_2d(:,:) )
       endif
    enddo
    clat = 0.5_RP * ( latstart + latend ) * CONST_D2R
    clon = 0.5_RP * ( lonstart + lonend ) * CONST_D2R


    !##### write control file #####

    LOG_INFO("SNO_grads_write_ctl",*) 'filename : ', trim(ctlname)

    fid = IO_get_available_fid()
    open( unit   = fid,                    &
          file   = trim(ctlname),          &
          form   = 'formatted',            &
          status = 'replace'               )

       write(fid,'(A)')       'DSET ^'//trim(grdname)
       write(fid,'(A)')       'TITLE SCALE-RM data output'
       write(fid,'(A)')       'OPTIONS BIG_ENDIAN'
       write(fid,'(A,E12.5)') 'UNDEF ', CONST_UNDEF

       !--- XDEF

       written = .false.
       do n = 1, naxis
          if ( ainfo(n)%varname == idim ) then
             dx   = ainfo(n)%AXIS_1d(imax/2+1) - ainfo(n)%AXIS_1d(imax/2)
             dlon = dx / ( CONST_RADIUS * cos(clat) ) / CONST_D2R

             if    ( hinfo%minfo_mapping_name == 'lambert_conformal_conic' ) then
                imax_ = int(imax*1.1_RP)
                write(fid,'(A,I5,A,1x,F9.2,1x,F9.3)') 'XDEF ', imax_, ' LINEAR', lonstart, dlon
             elseif( hinfo%minfo_mapping_name == 'polar_stereographic' ) then
                write(fid,'(A)') 'XDEF   720 LINEAR -179.5 0.5'
             else
                write(fid,'(A,I5,A,1x,F9.2,1x,F9.3)') 'XDEF ', imax, ' LINEAR', lonstart, dlon
             endif
             written = .true.
          endif
       enddo

       if ( .NOT. written ) then
          LOG_ERROR("SNO_grads_write_ctl",*) 'AXIS data for XDEF not found. ', trim(idim)
          call PRC_abort
       endif

       !--- YDEF

       written = .false.
       do n = 1, naxis
          if ( ainfo(n)%varname == jdim ) then
             dy    = ainfo(n)%AXIS_1d(jmax/2+1) - ainfo(n)%AXIS_1d(jmax/2)
             dlat  = dy / ( CONST_RADIUS * cos(clat) ) / CONST_D2R

             if    ( hinfo%minfo_mapping_name == 'lambert_conformal_conic' ) then
                jmax_ = int(jmax*0.9_RP)
                write(fid,'(A,I5,A,1x,F9.2,1x,F9.3)') 'YDEF ', jmax_, ' LINEAR', latstart, dlat
             elseif( hinfo%minfo_mapping_name == 'polar_stereographic' ) then
                latstart = real(int(latstart),kind=RP)
                jmax_ = int( ( 90.0_RP - abs(latstart) ) / 0.5_RP )
                write(fid,'(A,I5,A,1x,F9.2,1x,F9.3)') 'YDEF ', jmax_, ' LINEAR', latstart, 0.5_RP
             else
                write(fid,'(A,I5,A,1x,F9.2,1x,F9.3)') 'YDEF ', jmax, ' LINEAR', latstart, dlat
             endif
             written = .true.
          endif
       enddo

       if ( .NOT. written ) then
          LOG_ERROR("SNO_grads_write_ctl",*) 'AXIS data for YDEF not found. ', trim(idim)
          call PRC_abort
       endif

       !--- ZDEF

       written = .false.
       do n = 1, naxis
          if ( ainfo(n)%varname == kdim ) then
             write(fid,'(A,I5,A)') 'ZDEF ', kmax, ' LEVELS'
             write(fid,'(10(1x,F9.3))') (ainfo(n)%AXIS_1d(k),k=1,kmax)
             written = .true.
          endif
       enddo

       if ( .NOT. written ) then
          write(fid,'(A,I5,A,2I5)') 'ZDEF ', 1, ' LINEAR', 1, 1
       endif

       !--- TDEF

       call SNO_grads_calc_timechar( dinfo%time_units, dinfo%dt, & ! [IN]
                                     cdate, dhour                ) ! [OUT]

       write(fid,'(A,I5,3(1x,A))') 'TDEF ', dinfo%step_nmax, ' LINEAR ', trim(cdate), trim(dhour)

       !--- PDEF

       if ( hinfo%minfo_mapping_name == 'lambert_conformal_conic' ) then
          write(fid,'(A,2(1x,I5),A,2(1x,F9.2),2(1x,I5),5(1x,F9.2))') 'PDEF',                                       &
                                                                     imax,                                         &
                                                                     jmax,                                         &
                                                                     ' LCC',                                       &
                                                                     hinfo%minfo_latitude_of_projection_origin(1), &
                                                                     hinfo%minfo_longitude_of_central_meridian(1), &
                                                                     imax/2,                                       &
                                                                     jmax/2,                                       &
                                                                     hinfo%minfo_standard_parallel            (1), &
                                                                     hinfo%minfo_standard_parallel            (2), &
                                                                     hinfo%minfo_longitude_of_central_meridian(1), &
                                                                     dx,                                           &
                                                                     dy
       elseif( hinfo%minfo_mapping_name == 'polar_stereographic' ) then
          write(fid,'(A,2(1x,I5),A,2(1x,F9.2),2(1x,I5),3(1x,F9.2))') 'PDEF',                                               &
                                                                     imax,                                                 &
                                                                     jmax,                                                 &
                                                                     ' PSE',                                               &
                                                                     hinfo%minfo_latitude_of_projection_origin        (1), &
                                                                     hinfo%minfo_straight_vertical_longitude_from_pole(1), &
                                                                     imax/2,                                               &
                                                                     jmax/2,                                               &
                                                                     dx/1000.0_RP,                                         &
                                                                     dy/1000.0_RP,                                         &
                                                                     sign(1.0_DP,hinfo%minfo_standard_parallel(1))
       endif

       !--- VARS

       write(fid,'(A,I5)') 'VARS ', 1
       if ( kmax == 1 ) then
          write(fid,'(A,2I5,1X,A)') trim(dinfo%varname), 0   , 99, trim(dinfo%description)
       else
          write(fid,'(A,2I5,1X,A)') trim(dinfo%varname), kmax, 99, trim(dinfo%description)
       endif
       write(fid,'(A)') 'ENDVARS '

    close(fid)

    return
  end subroutine SNO_grads_write_ctl

  !-----------------------------------------------------------------------------
  subroutine SNO_grads_netcdfctl( &
       dirpath,  &
       basename, &
       hinfo,    &
       naxis,    &
       ainfo,    &
       nvars,    &
       dinfo,    &
       debug     )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       CONST_UNDEF, &
       CONST_D2R,   &
       CONST_RADIUS
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo,   &
       iteminfo
    implicit none

    character(len=*), intent(in)  :: dirpath                               ! directory path                     (output)
    character(len=*), intent(in)  :: basename                              ! basename of file                   (output)
    type(commoninfo), intent(in)  :: hinfo                                 ! common information                 (input)
    integer,          intent(in)  :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(in)  :: ainfo(naxis)                          ! axis information                   (input)
    integer,          intent(in)  :: nvars                                 ! number of variables                (input)
    type(iteminfo),   intent(in)  :: dinfo(nvars)                          ! variable information               (input)
    logical,          intent(in)  :: debug

    character(len=H_LONG)  :: grdname
    character(len=H_LONG)  :: ctlname

    character(len=H_SHORT) :: idim, jdim, kdim
    integer                :: imax, jmax, kmax

    real(RP)               :: latstart, latend
    real(RP)               :: lonstart, lonend
    real(RP)               :: clat, clon
    real(RP)               :: dx, dy
    real(RP)               :: dlat, dlon

    character(len=20)      :: cdate, dhour

    character(len=20)      :: dimorder
    character(len=H_SHORT) :: kdim_
    integer                :: kmax_
    integer                :: count

    logical  :: written
    integer  :: fid
    integer  :: k, n, v
    !---------------------------------------------------------------------------

    kmax = 1
    kdim = 'none'
    do v = 1, nvars
       if ( dinfo(v)%dim_rank == 2 ) then
          imax = size(dinfo(v)%VAR_2d(:,:),1)
          jmax = size(dinfo(v)%VAR_2d(:,:),2)

          idim = trim(dinfo(v)%dim_name(1))
          jdim = trim(dinfo(v)%dim_name(2))
       elseif( dinfo(v)%dim_rank == 3 ) then
          imax  = size(dinfo(v)%VAR_3d(:,:,:),2)
          jmax  = size(dinfo(v)%VAR_3d(:,:,:),3)
          kmax_ = size(dinfo(v)%VAR_3d(:,:,:),1)

          if ( dinfo(v)%transpose ) then
             idim  = trim(dinfo(v)%dim_name(1))
             jdim  = trim(dinfo(v)%dim_name(2))
             kdim_ = trim(dinfo(v)%dim_name(3))
          else
             idim  = trim(dinfo(v)%dim_name(2))
             jdim  = trim(dinfo(v)%dim_name(3))
             kdim_ = trim(dinfo(v)%dim_name(1))
          endif

          if ( kmax_ > kmax ) then
             kmax = kmax_
             kdim = kdim_
          endif
       endif
    enddo

    if ( basename == '' ) then
       LOG_ERROR("SNO_grads_netcdfctl",*) 'Namelist parameter basename_out in PARAM_SNO is empty. Check!'
       call PRC_abort
    endif

    grdname = trim(basename)//'.pe000000.nc'
    if ( dirpath == '' ) then
       ctlname = trim(basename)//'.ctl'
    else
       ctlname = trim(dirpath)//'/'//trim(basename)//'.ctl'
    endif

    do n = 1, naxis
       if    ( ainfo(n)%varname == 'lat' ) then
          latstart = minval( ainfo(n)%AXIS_2d(:,:) )
          latend   = maxval( ainfo(n)%AXIS_2d(:,:) )
       elseif( ainfo(n)%varname == 'lon' ) then
          lonstart = minval( ainfo(n)%AXIS_2d(:,:) )
          lonend   = maxval( ainfo(n)%AXIS_2d(:,:) )
       endif
    enddo
    clat = 0.5_RP * ( latstart + latend ) * CONST_D2R
    clon = 0.5_RP * ( lonstart + lonend ) * CONST_D2R


    !##### write control file #####

    fid = IO_get_available_fid()
    open( unit   = fid,                    &
          file   = trim(ctlname),          &
          form   = 'formatted',            &
          status = 'replace'               )

       write(fid,'(A)')       'DSET ^'//trim(grdname)
       write(fid,'(A)')       'TITLE SCALE-RM data output'
       write(fid,'(A)')       'DTYPE netcdf'
       write(fid,'(A,E12.5)') 'UNDEF ', CONST_UNDEF

       !--- XDEF

       written = .false.
       do n = 1, naxis
          if ( ainfo(n)%varname == idim ) then
             dx   = ainfo(n)%AXIS_1d(imax/2+1) - ainfo(n)%AXIS_1d(imax/2)
             dlon = dx / ( CONST_RADIUS * cos(clat) ) / CONST_D2R

             write(fid,'(A,I5,A,1x,F9.2,1x,F9.3)') 'XDEF ', int(imax*1.1_RP), ' LINEAR', lonstart, dlon
             written = .true.
          endif
       enddo

       if ( .NOT. written ) then
          LOG_ERROR("SNO_grads_netcdfctl",*) '[SNO_grads_write_ctl] AXIS data for XDEF not found. ', trim(idim)
          call PRC_abort
       endif

       !--- YDEF

       written = .false.
       do n = 1, naxis
          if ( ainfo(n)%varname == jdim ) then
             dy    = ainfo(n)%AXIS_1d(jmax/2+1) - ainfo(n)%AXIS_1d(jmax/2)
             dlat  = dy / ( CONST_RADIUS * cos(clat) ) / CONST_D2R

             write(fid,'(A,I5,A,1x,F9.2,1x,F9.3)') 'YDEF ', jmax, ' LINEAR', latstart, dlat
             written = .true.
          endif
       enddo

       if ( .NOT. written ) then
          LOG_ERROR("SNO_grads_netcdfctl",*) '[SNO_grads_write_ctl] AXIS data for YDEF not found. ', trim(idim)
          call PRC_abort
       endif

       !--- ZDEF

       written = .false.
       do n = 1, naxis
          if ( ainfo(n)%varname == kdim ) then
             write(fid,'(A,I5,A)') 'ZDEF ', kmax, ' LEVELS'
             write(fid,'(10(1x,F9.3))') (ainfo(n)%AXIS_1d(k),k=1,kmax)
             written = .true.
          endif
       enddo

       if ( .NOT. written ) then
          write(fid,'(A,I5,A,2I5)') 'ZDEF ', 1, ' LINEAR', 1, 1
       endif

       !--- TDEF

       call SNO_grads_calc_timechar( dinfo(1)%time_units, dinfo(1)%dt, & ! [IN]
                                     cdate, dhour                      ) ! [OUT]

       write(fid,'(A,I5,3(1x,A))') 'TDEF ', dinfo(1)%step_nmax, ' LINEAR ', trim(cdate), trim(dhour)

       !--- PDEF

       if ( hinfo%minfo_mapping_name == 'lambert_conformal_conic' ) then
          write(fid,'(A,2(1x,I5),A,2(1x,F9.2),2(1x,I5),5(1x,F9.2))') 'PDEF',                                       &
                                                                     imax,                                         &
                                                                     jmax,                                         &
                                                                     ' LCC',                                       &
                                                                     hinfo%minfo_latitude_of_projection_origin(1), &
                                                                     hinfo%minfo_longitude_of_central_meridian(1), &
                                                                     imax/2,                                       &
                                                                     jmax/2,                                       &
                                                                     hinfo%minfo_standard_parallel            (1), &
                                                                     hinfo%minfo_standard_parallel            (2), &
                                                                     hinfo%minfo_longitude_of_central_meridian(1), &
                                                                     dx,                                           &
                                                                     dy
       elseif( hinfo%minfo_mapping_name == 'polar_stereographic' ) then
          write(fid,'(A,2(1x,I5),A,2(1x,F9.2),2(1x,I5),5(1x,F9.2))') 'PDEF',                                               &
                                                                     imax,                                                 &
                                                                     jmax,                                                 &
                                                                     ' LCC',                                               &
                                                                     hinfo%minfo_latitude_of_projection_origin        (1), &
                                                                     hinfo%minfo_straight_vertical_longitude_from_pole(1), &
                                                                     imax/2,                                               &
                                                                     jmax/2,                                               &
                                                                     hinfo%minfo_standard_parallel                    (1), &
                                                                     hinfo%minfo_standard_parallel                    (1), &
                                                                     hinfo%minfo_straight_vertical_longitude_from_pole(1), &
                                                                     dx,                                                   &
                                                                     dy
       endif

       !--- VARS

       count = 0
       do v = 1, nvars
          if    ( dinfo(v)%dim_rank == 1 ) then
             cycle ! skip
          elseif( dinfo(v)%dim_rank == 2 ) then
             if( size(dinfo(v)%VAR_2d(:,:),1) /= imax ) cycle ! skip
             if( size(dinfo(v)%VAR_2d(:,:),2) /= jmax ) cycle ! skip
             if ( dinfo(v)%step_nmax > 1 ) then
                if( abs(dinfo(v)%dt-dinfo(1)%dt) > 1.E-5_DP ) cycle ! skip
             endif
          elseif( dinfo(v)%dim_rank == 3 ) then
             if( size(dinfo(v)%VAR_3d(:,:,:),2) /= imax ) cycle ! skip
             if( size(dinfo(v)%VAR_3d(:,:,:),3) /= jmax ) cycle ! skip
             if( size(dinfo(v)%VAR_3d(:,:,:),1) /= kmax ) cycle ! skip
             if ( dinfo(v)%transpose ) then
                if ( dinfo(v)%step_nmax > 1 ) then
                   if( abs(dinfo(v)%dt-dinfo(1)%dt) > 1.E-5_DP ) cycle ! skip
                endif
             else
                if ( dinfo(v)%step_nmax > 1 ) then
                   if( abs(dinfo(v)%dt-dinfo(1)%dt) > 1.E-5_DP ) cycle ! skip
                endif
             endif
          endif

          count = count + 1
       enddo

       write(fid,'(A,I5)') 'VARS ', count
       do v = 1, nvars
          if    ( dinfo(v)%dim_rank == 1 ) then
             cycle ! skip
          elseif( dinfo(v)%dim_rank == 2 ) then
             kmax_ = 0

             if( size(dinfo(v)%VAR_2d(:,:),1) /= imax ) cycle ! skip
             if( size(dinfo(v)%VAR_2d(:,:),2) /= jmax ) cycle ! skip

             if ( dinfo(v)%step_nmax > 1 ) then
                if( abs(dinfo(v)%dt-dinfo(1)%dt) > 1.E-5_DP ) cycle ! skip

                dimorder = 't,y,x'
             else
                dimorder = 'y,x'
             endif
          elseif( dinfo(v)%dim_rank == 3 ) then
             kmax_ = kmax

             if( size(dinfo(v)%VAR_3d(:,:,:),2) /= imax ) cycle ! skip
             if( size(dinfo(v)%VAR_3d(:,:,:),3) /= jmax ) cycle ! skip
             if( size(dinfo(v)%VAR_3d(:,:,:),1) /= kmax ) cycle ! skip

             if ( dinfo(v)%transpose ) then
                if ( dinfo(v)%step_nmax > 1 ) then
                   if( abs(dinfo(v)%dt-dinfo(1)%dt) > 1.E-5_DP ) cycle ! skip

                   dimorder = 't,z,y,x'
                else
                   dimorder = 'z,y,x'
                endif
             else
                if ( dinfo(v)%step_nmax > 1 ) then
                   if( abs(dinfo(v)%dt-dinfo(1)%dt) > 1.E-5_DP ) cycle ! skip

                   dimorder = 't,y,x,z'
                else
                   dimorder = 'y,x,z'
                endif
             endif
          endif

          write(fid,'(A,I5,2(1X,A))') trim(dinfo(v)%varname)//'=>'//trim(dinfo(v)%varname), &
                                      kmax_, trim(dimorder), trim(dinfo(v)%description)
       enddo
       write(fid,'(A)') 'ENDVARS '

    close(fid)

    return
  end subroutine SNO_grads_netcdfctl

  !-----------------------------------------------------------------------------
  subroutine SNO_grads_calc_timechar( &
       time_units, &
       dt,         &
       cdate,      &
       dhour       )
    use scale_calendar, only: &
       CALENDAR_daysec2date,   &
       CALENDAR_adjust_daysec, &
       CALENDAR_CFunits2sec
    use scale_const, only: &
       CONST_EPS
    implicit none

    character(len=*),  intent(in)  :: time_units
    real(DP),          intent(in)  :: dt
    character(len=20), intent(out) :: cdate
    character(len=20), intent(out) :: dhour

    integer  :: start_day     !< absolute day    (time=0)
    real(DP) :: start_sec     !< absolute second (time=0)
    integer  :: start_date(6) !< date            (time=t)
    real(DP) :: start_ms      !< subsecond       (time=t)

    character(len=3) :: mlist(12)
    data mlist / 'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC' /
    !---------------------------------------------------------------------------

    start_day = 0
    start_sec = CALENDAR_CFunits2sec( cftime=0.0_DP, cfunits=time_units, offset_year=0 )

    call CALENDAR_adjust_daysec( start_day, start_sec ) ! [INOUT]

    call CALENDAR_daysec2date  ( start_date, start_ms,                & ! [OUT]
                                 start_day,  start_sec, offset_year=0 ) ! [IN]

    write(cdate,'(I2.2,A1,I2.2,A1,I2.2,A3,I4.4)') start_date(4),        &
                                                  ':',                  &
                                                  start_date(5),        &
                                                  'Z',                  &
                                                  start_date(3),        &
                                                  mlist(start_date(2)), &
                                                  start_date(1)

    if(      int(dt/3600.0_DP/24.0_DP) >= 1 .and. mod(dt,86400.0_DP) < CONST_EPS ) then
       write(dhour,'(I3)') int(dt/3600.0_DP/24.0_DP) ! day
       dhour = trim(dhour)//'dy'
    else if( int(dt/3600.0_DP) >= 1         .and. mod(dt,3600.0_DP) < CONST_EPS ) then
       write(dhour,'(I3)') int(dt/3600.0_DP) ! hour
       dhour = trim(dhour)//'hr'
    else if( int(dt/60.0_DP) >= 1           .and. mod(dt, 60.0_DP) < CONST_EPS ) then
       write(dhour,'(I3)') int(dt/60.0_DP)   ! minute
       dhour = trim(dhour)//'mn'
    else
       LOG_WARN("SNO_grads_calc_timechar",*) 'Output interval is not supported by GrADS.'
       !write(dhour,'(I3)') max(int(dt/60.0_DP),1) ! less than 1 minute
       write(dhour,'(I8)') int(dt) ! less than 1 minute
       dhour = trim(dhour)//'sec'
    endif

    return
  end subroutine SNO_grads_calc_timechar

end module mod_sno_grads
