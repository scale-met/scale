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
module mod_sno_axis
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SNO_axis_getinfo
  public :: SNO_axis_alloc
  public :: SNO_axis_dealloc
  public :: SNO_axis_copy
  public :: SNO_axis_read
  public :: SNO_axis_define
  public :: SNO_axis_write

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine SNO_axis_getinfo( &
       ismaster,  &
       basename,  &
       naxis,     &
       axisname,  &
       ainfo,     &
       debug      )
    use mpi
    use scale_file_h, only: &
       FILE_dtypelist
    use scale_file, only: &
       FILE_Get_Datainfo
    use scale_process, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD, &
       PRC_MPIstop
    use mod_sno_h, only: &
       dim_limit, &
       axisinfo
    use mod_sno, only: &
       SNO_read_bcast_1d
    implicit none

    logical,          intent(in)  :: ismaster                 ! master process?                    (execution)
    character(len=*), intent(in)  :: basename                 ! basename of file                   (input)
    integer,          intent(in)  :: naxis                    ! number of axis variables           (input)
    character(len=*), intent(in)  :: axisname(naxis)          ! name   of axis variables           (input)
    type(axisinfo),   intent(out) :: ainfo   (naxis)          ! axis information                   (input)
    logical,          intent(in)  :: debug

    integer :: nowrank, nowstep
    integer :: ierr
    integer :: n, d
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** [SNO_axis_getinfo] Read information of axis'

    if ( ismaster ) then
       nowrank = 0 ! first file
       nowstep = 0

       do n = 1, naxis
          if ( debug ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** read info : ', trim(axisname(n))
          endif

          call FILE_Get_Datainfo( basename,                           & ! [IN]
                                  axisname(n),                        & ! [IN]
                                  nowrank,                            & ! [IN]
                                  nowstep,                            & ! [IN]
                                  description = ainfo(n)%description, & ! [OUT]
                                  units       = ainfo(n)%units,       & ! [OUT]
                                  datatype    = ainfo(n)%datatype,    & ! [OUT]
                                  dim_rank    = ainfo(n)%dim_rank,    & ! [OUT]
                                  dim_name    = ainfo(n)%dim_name(:), & ! [OUT]
                                  dim_size    = ainfo(n)%dim_size(:)  ) ! [OUT]

          ainfo(n)%varname = axisname(n)

          ainfo(n)%transpose = .false.
          if ( ainfo(n)%dim_rank > 2 ) then
             if (       ainfo(n)%dim_name(1) /= 'z'   &
                  .AND. ainfo(n)%dim_name(1) /= 'zh'  &
                  .AND. ainfo(n)%dim_name(1) /= 'oz'  &
                  .AND. ainfo(n)%dim_name(1) /= 'ozh' &
                  .AND. ainfo(n)%dim_name(1) /= 'lz'  &
                  .AND. ainfo(n)%dim_name(1) /= 'lzh' &
                  .AND. ainfo(n)%dim_name(1) /= 'uz'  &
                  .AND. ainfo(n)%dim_name(1) /= 'uzh' ) then
                ainfo(n)%transpose = .true.
             endif
          endif

          ainfo(n)%regrid = .false.
          do d = 1, ainfo(n)%dim_rank
             select case(ainfo(n)%dim_name(d))
             case('x','xh','y','yh','CX','CY','FX','FY','FDX','FDY')
                ainfo(n)%regrid = .true.
             end select
          enddo

       enddo
    endif

    do n = 1, naxis
       call MPI_BCAST( ainfo(n)%varname    , H_SHORT          , MPI_CHARACTER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( ainfo(n)%description, H_MID            , MPI_CHARACTER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( ainfo(n)%units      , H_SHORT          , MPI_CHARACTER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( ainfo(n)%datatype   , 1                , MPI_INTEGER  , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( ainfo(n)%dim_rank   , 1                , MPI_INTEGER  , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( ainfo(n)%dim_name(:), H_SHORT*dim_limit, MPI_CHARACTER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( ainfo(n)%dim_size(:), dim_limit        , MPI_INTEGER  , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( ainfo(n)%transpose  , 1                , MPI_LOGICAL  , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( ainfo(n)%regrid     , 1                , MPI_LOGICAL  , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

       if ( debug ) then
          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '*** Axis No.', n
          if( IO_L ) write(IO_FID_LOG,*) '*** varname     : ', trim(ainfo(n)%varname)
          if( IO_L ) write(IO_FID_LOG,*) '*** description : ', trim(ainfo(n)%description)
          if( IO_L ) write(IO_FID_LOG,*) '*** units       : ', trim(ainfo(n)%units)
          if( IO_L ) write(IO_FID_LOG,*) '*** datatype    : ', trim(FILE_dtypelist(ainfo(n)%datatype))
          if( IO_L ) write(IO_FID_LOG,*) '*** dim_rank    : ', ainfo(n)%dim_rank
          do d = 1, ainfo(n)%dim_rank
             if( IO_L ) write(IO_FID_LOG,*) '*** dim No.', d
             if( IO_L ) write(IO_FID_LOG,*) '*** + dim_name  : ', trim(ainfo(n)%dim_name(d))
             if( IO_L ) write(IO_FID_LOG,*) '*** + dim_size  : ', ainfo(n)%dim_size(d)
          enddo
          if( IO_L ) write(IO_FID_LOG,*) '*** transpose   : ', ainfo(n)%transpose
          if( IO_L ) write(IO_FID_LOG,*) '*** regrid      : ', ainfo(n)%regrid
       endif

       if ( ainfo(n)%dim_rank == 1 .AND. .NOT. ainfo(n)%regrid ) then

          allocate( ainfo(n)%AXIS_1d( ainfo(n)%dim_size(1) ) )

          call SNO_read_bcast_1d( ismaster,             & ! [IN]
                                  basename,             & ! [IN]
                                  ainfo(n)%varname,     & ! [IN]
                                  ainfo(n)%datatype,    & ! [IN]
                                  ainfo(n)%dim_size(1), & ! [IN]
                                  ainfo(n)%AXIS_1d(:)   ) ! [OUT]

          if ( debug ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** value : ', ainfo(n)%AXIS_1d(:)
          endif
       endif
    enddo

    return
  end subroutine SNO_axis_getinfo

  !-----------------------------------------------------------------------------
  subroutine SNO_axis_alloc( &
       nprocs_x_out,  &
       nprocs_y_out,  &
       ngrids_x_out,  &
       ngrids_y_out,  &
       ngrids_xh_out, &
       ngrids_yh_out, &
       hinfo,         &
       naxis,         &
       ainfo,         &
       debug          )
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo
    implicit none

    integer,          intent(in)    :: nprocs_x_out             ! x length of 2D processor topology  (output)
    integer,          intent(in)    :: nprocs_y_out             ! y length of 2D processor topology  (output)
    integer,          intent(in)    :: ngrids_x_out             ! number of x-axis grids per process (output,sometimes including halo)
    integer,          intent(in)    :: ngrids_y_out             ! number of y-axis grids per process (output,sometimes including halo)
    integer,          intent(in)    :: ngrids_xh_out            ! number of x-axis grids per process (output,sometimes including halo,staggered)
    integer,          intent(in)    :: ngrids_yh_out            ! number of y-axis grids per process (output,sometimes including halo,staggered)
    type(commoninfo), intent(in)    :: hinfo                    ! common information                 (input)
    integer,          intent(in)    :: naxis                    ! number of axis variables           (input)
    type(axisinfo),   intent(inout) :: ainfo(naxis)             ! axis information                   (input)
    logical,          intent(in)    :: debug

    integer  :: IA_out, JA_out
    integer  :: gout1, gout2, gout3
    integer  :: n
    !---------------------------------------------------------------------------

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** [SNO_axis_alloc] Allocate axis array'
    endif

    IA_out = ( hinfo%gridsize(2) ) / nprocs_x_out + 2 * hinfo%halosize(2)
    JA_out = ( hinfo%gridsize(3) ) / nprocs_y_out + 2 * hinfo%halosize(3)

    do n = 1, naxis
       if ( ainfo(n)%regrid ) then
          if ( ainfo(n)%dim_rank == 1 ) then

             if    ( ainfo(n)%dim_name(1) == 'x'   ) then
                gout1 = ngrids_x_out
             elseif( ainfo(n)%dim_name(1) == 'xh'  ) then
                gout1 = ngrids_xh_out
             elseif( ainfo(n)%dim_name(1) == 'y'   ) then
                gout1 = ngrids_y_out
             elseif( ainfo(n)%dim_name(1) == 'yh'  ) then
                gout1 = ngrids_yh_out
             elseif( ainfo(n)%dim_name(1) == 'CX'  ) then
                gout1 = IA_out
             elseif( ainfo(n)%dim_name(1) == 'CY'  ) then
                gout1 = JA_out
             elseif( ainfo(n)%dim_name(1) == 'FX'  ) then
                gout1 = IA_out+1
             elseif( ainfo(n)%dim_name(1) == 'FY'  ) then
                gout1 = JA_out+1
             elseif( ainfo(n)%dim_name(1) == 'FDX' ) then
                gout1 = IA_out-1
             elseif( ainfo(n)%dim_name(1) == 'FDY' ) then
                gout1 = JA_out-1
             endif

             allocate( ainfo(n)%AXIS_1d(gout1) )
             ainfo(n)%AXIS_1d(:) = 0.0_RP

          elseif( ainfo(n)%dim_rank == 2 ) then

             if    ( ainfo(n)%dim_name(1) == 'x'  ) then
                gout1 = ngrids_x_out
             elseif( ainfo(n)%dim_name(1) == 'xh' ) then
                gout1 = ngrids_xh_out
             else
                gout1 = ainfo(n)%dim_size(1)
             endif

             if    ( ainfo(n)%dim_name(2) == 'y'  ) then
                gout2 = ngrids_y_out
             elseif( ainfo(n)%dim_name(2) == 'yh' ) then
                gout2 = ngrids_yh_out
             else
                gout2 = ainfo(n)%dim_size(2)
             endif

             allocate( ainfo(n)%AXIS_2d(gout1,gout2) )
             ainfo(n)%AXIS_2d(:,:) = 0.0_RP

          elseif( ainfo(n)%dim_rank == 3 ) then

             if ( ainfo(n)%transpose ) then
                gout1 = ainfo(n)%dim_size(3)

                if    ( ainfo(n)%dim_name(1) == 'x'  ) then
                   gout2 = ngrids_x_out
                elseif( ainfo(n)%dim_name(1) == 'xh' ) then
                   gout2 = ngrids_xh_out
                else
                   gout2 = ainfo(n)%dim_size(1)
                endif

                if    ( ainfo(n)%dim_name(2) == 'y'  ) then
                   gout3 = ngrids_y_out
                elseif( ainfo(n)%dim_name(2) == 'yh' ) then
                   gout3 = ngrids_yh_out
                else
                   gout3 = ainfo(n)%dim_size(2)
                endif
             else
                gout1 = ainfo(n)%dim_size(1)

                if    ( ainfo(n)%dim_name(2) == 'x'  ) then
                   gout2 = ngrids_x_out
                elseif( ainfo(n)%dim_name(2) == 'xh' ) then
                   gout2 = ngrids_xh_out
                else
                   gout2 = ainfo(n)%dim_size(2)
                endif

                if    ( ainfo(n)%dim_name(3) == 'y'  ) then
                   gout3 = ngrids_y_out
                elseif( ainfo(n)%dim_name(3) == 'yh' ) then
                   gout3 = ngrids_yh_out
                else
                   gout3 = ainfo(n)%dim_size(3)
                endif
             endif

             allocate( ainfo(n)%AXIS_3d(gout1,gout2,gout3) )
             ainfo(n)%AXIS_3d(:,:,:) = 0.0_RP

          endif
       endif
    enddo

    return
  end subroutine SNO_axis_alloc

  !-----------------------------------------------------------------------------
  subroutine SNO_axis_dealloc( &
       naxis, &
       ainfo, &
       debug  )
    use mod_sno_h, only: &
       axisinfo
    implicit none

    integer,          intent(in)    :: naxis                    ! number of axis variables           (input)
    type(axisinfo),   intent(inout) :: ainfo(naxis)             ! axis information                   (input)
    logical,          intent(in)    :: debug

    integer  :: n
    !---------------------------------------------------------------------------

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** [SNO_axis_dealloc] Deallocate axis array'
    endif

    do n = 1, naxis
       if ( ainfo(n)%regrid ) then
          if( allocated(ainfo(n)%AXIS_1d) ) deallocate( ainfo(n)%AXIS_1d )
          if( allocated(ainfo(n)%AXIS_2d) ) deallocate( ainfo(n)%AXIS_2d )
          if( allocated(ainfo(n)%AXIS_3d) ) deallocate( ainfo(n)%AXIS_3d )
       endif
    enddo

    return
  end subroutine SNO_axis_dealloc

  !-----------------------------------------------------------------------------
  subroutine SNO_axis_copy( &
       nprocs_x_out,  &
       nprocs_y_out,  &
       px,            &
       py,            &
       hinfo,         &
       naxis,         &
       ainfo,         &
       debug          )
    use scale_process, only: &
       PRC_MPIstop
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo
    implicit none

    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology  (output)
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology  (output)
    integer,          intent(in)    :: px                                    ! x index  in 2D processor topology
    integer,          intent(in)    :: py                                    ! y index  in 2D processor topology
    type(commoninfo), intent(in)    :: hinfo                                 ! common information                 (input)
    integer,          intent(in)    :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(inout) :: ainfo   (naxis)                       ! axis information                   (input)
    logical,          intent(in)    :: debug

    integer  :: IA_out, JA_out
    integer  :: IMAX_out, JMAX_out
    integer  :: IS, IE
    integer  :: JS, JE
    logical  :: exist

    integer  :: n, nn
    !---------------------------------------------------------------------------

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** [SNO_axis_copy] Copy axis array (local)'
    endif

    IMAX_out = ( hinfo%gridsize(2) ) / nprocs_x_out
    JMAX_out = ( hinfo%gridsize(3) ) / nprocs_y_out
    IA_out   = ( hinfo%gridsize(2) ) / nprocs_x_out + 2 * hinfo%halosize(2)
    JA_out   = ( hinfo%gridsize(3) ) / nprocs_y_out + 2 * hinfo%halosize(3)

    do n = 1, naxis

       exist = .false.

       if ( ainfo(n)%varname == 'CX' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'CXG' ) then
                IS = (px-1) * IMAX_out + 1
                IE = (px-1) * IMAX_out + IA_out

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(IS:IE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS CXG not found! : necessary for CX'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'CY' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'CYG' ) then
                JS = (py-1) * JMAX_out + 1
                JE = (py-1) * JMAX_out + JA_out

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(JS:JE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS CYG not found! : necessary for CY'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'FX' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'FXG' ) then
                IS = (px-1) * IMAX_out + 1
                IE = (px-1) * IMAX_out + IA_out + 1

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(IS:IE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS FXG not found! : necessary for FX'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'FY' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'FYG' ) then
                JS = (py-1) * JMAX_out + 1
                JE = (py-1) * JMAX_out + JA_out + 1

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(JS:JE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS FYG not found! : necessary for FY'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'CDX' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'CDXG' ) then
                IS = (px-1) * IMAX_out + 1
                IE = (px-1) * IMAX_out + IA_out

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(IS:IE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS CDXG not found! : necessary for CDX'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'CDY' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'CDYG' ) then
                JS = (py-1) * JMAX_out + 1
                JE = (py-1) * JMAX_out + JA_out

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(JS:JE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS CDYG not found! : necessary for CDY'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'FDX' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'FDXG' ) then
                IS = (px-1) * IMAX_out + 1
                IE = (px-1) * IMAX_out + IA_out - 1

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(IS:IE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS FDXG not found! : necessary for FDX'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'FDY' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'FDYG' ) then
                JS = (py-1) * JMAX_out + 1
                JE = (py-1) * JMAX_out + JA_out - 1

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(JS:JE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS FDYG not found! : necessary for FDY'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'CBFX' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'CBFXG' ) then
                IS = (px-1) * IMAX_out + 1
                IE = (px-1) * IMAX_out + IA_out

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(IS:IE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS CBFXG not found! : necessary for CBFX'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'CBFY' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'CBFYG' ) then
                JS = (py-1) * JMAX_out + 1
                JE = (py-1) * JMAX_out + JA_out

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(JS:JE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS CBFYG not found! : necessary for CBFY'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'FBFX' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'FBFXG' ) then
                IS = (px-1) * IMAX_out + 1
                IE = (px-1) * IMAX_out + IA_out + 1

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(IS:IE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS FBFXG not found! : necessary for FBFX'
             call PRC_MPIstop
          endif

       elseif( ainfo(n)%varname == 'FBFY' ) then

          do nn = 1, naxis
             if ( ainfo(nn)%varname == 'FBFYG' ) then
                JS = (py-1) * JMAX_out + 1
                JE = (py-1) * JMAX_out + JA_out + 1

                ainfo(n)%AXIS_1d(:) = ainfo(nn)%AXIS_1d(JS:JE)
                exist = .true.
             endif
          enddo

          if ( .NOT. exist ) then
             write(*,*) 'xxx [SNO_axis_read] AXIS FBFYG not found! : necessary for FBFY'
             call PRC_MPIstop
          endif

       endif
    enddo

    return
  end subroutine SNO_axis_copy

  !-----------------------------------------------------------------------------
  subroutine SNO_axis_read( &
       basename,      &
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
       naxis,         &
       ainfo,         &
       localmap,      &
       readflag,      &
       debug          )
    use scale_process, only: &
       PRC_MPIstop
    use mod_sno_h, only: &
       I_map_p,    &
       I_map_i,    &
       I_map_j,    &
       commoninfo, &
       axisinfo
    use mod_sno, only: &
       SNO_calc_localsize, &
       SNO_read_map_1d,    &
       SNO_read_map_2d,    &
       SNO_read_map_3d
    implicit none

    character(len=*), intent(in)    :: basename                              ! basename of file                   (input)
    integer,          intent(in)    :: nprocs_x_in                           ! x length of 2D processor topology  (input)
    integer,          intent(in)    :: nprocs_y_in                           ! y length of 2D processor topology  (input)
    integer,          intent(in)    :: ngrids_x                              ! number of x-axis grids             (global domain,sometimes including halo)
    integer,          intent(in)    :: ngrids_y                              ! number of y-axis grids             (global domain,sometimes including halo)
    integer,          intent(in)    :: nhalos_x                              ! number of x-axis halo grids        (global domain)
    integer,          intent(in)    :: nhalos_y                              ! number of y-axis halo grids        (global domain)
    type(commoninfo), intent(in)    :: hinfo                                 ! common information                 (input)
    integer,          intent(in)    :: ngrids_x_out                          ! number of x-axis grids per process (output,sometimes including halo)
    integer,          intent(in)    :: ngrids_y_out                          ! number of y-axis grids per process (output,sometimes including halo)
    integer,          intent(in)    :: ngrids_xh_out                         ! number of x-axis grids per process (output,sometimes including halo)
    integer,          intent(in)    :: ngrids_yh_out                         ! number of y-axis grids per process (output,sometimes including halo)
    integer,          intent(in)    :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(inout) :: ainfo   (naxis)                       ! axis information                   (input)
    integer(2),       intent(in)    :: localmap(ngrids_x_out,ngrids_y_out,3) ! mapping table from input to output (local domain)
    logical,          intent(in)    :: readflag(nprocs_x_in,nprocs_y_in)     ! local domain requires the input file?
    logical,          intent(in)    :: debug

    integer(2), allocatable :: localmap_1d(:,:)
    integer(2), allocatable :: localmap_2d(:,:,:)

    integer  :: ngrids_x_in       ! number of x-axis grids per process (input,without halo)
    integer  :: ngrids_y_in       ! number of y-axis grids per process (input,without halo)
    integer  :: ngrids_xh_in      ! number of x-axis grids per process (input,without halo)
    integer  :: ngrids_yh_in      ! number of y-axis grids per process (input,without halo)
    integer  :: staggered_x_in
    integer  :: staggered_y_in
    integer  :: staggered_x_out
    integer  :: staggered_y_out

    integer  :: gin1, gin2, gin3
    integer  :: gout1, gout2, gout3
    integer  :: stgin1, stgin2, stgin3
    integer  :: stgout1, stgout2, stgout3
    logical  :: readflag_1d

    integer  :: p, px, py
    integer  :: n
    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** [SNO_axis_read] Read axis array (local)'
    endif

    do py = 1, nprocs_y_in
       do px = 1, nprocs_x_in
          p = (py-1) * nprocs_x_in + px - 1

          call SNO_calc_localsize( nprocs_x_in,  nprocs_y_in, & ! [IN] from namelist
                                   px,           py,          & ! [IN]
                                   ngrids_x,     ngrids_y,    & ! [IN] from SNO_file_getinfo
                                   nhalos_x,     nhalos_y,    & ! [IN] from SNO_file_getinfo
                                   hinfo,                     & ! [IN] from SNO_file_getinfo
                                   ngrids_x_in,  ngrids_y_in, & ! [OUT]
                                   ngrids_xh_in, ngrids_yh_in ) ! [OUT]

          staggered_x_in  = 0
          staggered_y_in  = 0
          staggered_x_out = 0
          staggered_y_out = 0
          if ( ngrids_x_in  /= ngrids_xh_in  ) staggered_x_in  = 1
          if ( ngrids_y_in  /= ngrids_yh_in  ) staggered_y_in  = 1
          if ( ngrids_x_out /= ngrids_xh_out ) staggered_x_out = 1
          if ( ngrids_y_out /= ngrids_yh_out ) staggered_y_out = 1

          if ( readflag(px,py) ) then

             do n = 1, naxis
                if ( ainfo(n)%regrid ) then
                   if ( ainfo(n)%dim_rank == 1 ) then

                      readflag_1d = .false.

                      if    ( ainfo(n)%varname == 'x' ) then
                         gin1    = ngrids_x_in
                         gout1   = ngrids_x_out
                         stgin1  = 0
                         stgout1 = 0

                         allocate( localmap_1d(gout1,2) )

                         do i = 1, ngrids_x_out
                            localmap_1d(i,I_map_p) = localmap(i,1,I_map_p)
                            localmap_1d(i,I_map_i) = localmap(i,1,I_map_i)
                         enddo

                         readflag_1d = .true.
                      elseif( ainfo(n)%varname == 'xh' ) then
                         gin1    = ngrids_xh_in
                         gout1   = ngrids_xh_out
                         stgin1  = staggered_x_in
                         stgout1 = staggered_x_out

                         allocate( localmap_1d(gout1,2) )

                         do i = 1, ngrids_x_out
                            localmap_1d(i+stgout1,I_map_p) = localmap(i,1,I_map_p)
                            localmap_1d(i+stgout1,I_map_i) = localmap(i,1,I_map_i) + stgin1
                         enddo

                         if ( stgout1 > 0 ) then
                            if ( stgin1 > 0 ) then
                               do i = 1, stgout1
                                  localmap_1d(i,I_map_p) = localmap(1,1,I_map_p)
                                  localmap_1d(i,I_map_i) = i
                               enddo
                            else
                               do i = 1, stgout1
                                  localmap_1d(i,I_map_p) = -1
                                  localmap_1d(i,I_map_i) = -1
                               enddo
                            endif
                         endif

                         readflag_1d = .true.
                      elseif( ainfo(n)%varname == 'y' ) then
                         gin1    = ngrids_y_in
                         gout1   = ngrids_y_out
                         stgin1  = 0
                         stgout1 = 0

                         allocate( localmap_1d(gout1,2) )

                         do j = 1, ngrids_y_out
                            localmap_1d(j,I_map_p) = localmap(1,j,I_map_p)
                            localmap_1d(j,I_map_i) = localmap(1,j,I_map_j)
                         enddo

                         readflag_1d = .true.
                      elseif( ainfo(n)%varname == 'yh' ) then
                         gin1    = ngrids_yh_in
                         gout1   = ngrids_yh_out
                         stgin1  = staggered_y_in
                         stgout1 = staggered_y_out

                         allocate( localmap_1d(gout1,2) )

                         do j = 1, ngrids_y_out
                            localmap_1d(j+stgout1,I_map_p) = localmap(1,j,I_map_p)
                            localmap_1d(j+stgout1,I_map_i) = localmap(1,j,I_map_j) + stgin1
                         enddo

                         if ( stgout1 > 0 ) then
                            if ( stgin1 > 0 ) then
                               do i = 1, stgout1
                                  localmap_1d(i,I_map_p) = localmap(1,1,I_map_p)
                                  localmap_1d(i,I_map_i) = i
                               enddo
                            else
                               do i = 1, stgout1
                                  localmap_1d(i,I_map_p) = -1
                                  localmap_1d(i,I_map_i) = -1
                               enddo
                            endif
                         endif

                         readflag_1d = .true.
                      elseif( ainfo(n)%varname == 'CX' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'CY' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'FX' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'FY' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'CDX' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'CDY' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'FDX' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'FDY' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'CBFX' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'CBFY' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'FBFX' ) then
                         ! tentatively copy from global axis, without mapping
                      elseif( ainfo(n)%varname == 'FBFY' ) then
                         ! tentatively copy from global axis, without mapping
                      endif

                      if ( readflag_1d ) then
!                          if ( debug ) then
!                             if( IO_L ) write(IO_FID_LOG,*)
!                             if( IO_L ) write(IO_FID_LOG,*) '*** Axis No.', n, ' : ', trim(ainfo(n)%varname)
!                             if( IO_L ) write(IO_FID_LOG,*) '*** staggered_x_in  = ', staggered_x_in
!                             if( IO_L ) write(IO_FID_LOG,*) '*** staggered_y_in  = ', staggered_y_in
!                             if( IO_L ) write(IO_FID_LOG,*) '*** staggered_x_out = ', staggered_x_out
!                             if( IO_L ) write(IO_FID_LOG,*) '*** staggered_y_out = ', staggered_y_out
!                             if( IO_L ) write(IO_FID_LOG,*) '*** stgin1          = ', stgin1
!                             if( IO_L ) write(IO_FID_LOG,*) '*** stgout1         = ', stgout1
!
!                             if( IO_L ) write(IO_FID_LOG,*) 'localmap(rank)'
!                             do i = 1, gout1
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                             do i = 1, gout1
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap_1d(i,I_map_p)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!
!                             if( IO_L ) write(IO_FID_LOG,*) 'localmap(i-index)'
!                             do i = 1, gout1
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                             do i = 1, gout1
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap_1d(i,I_map_i)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                          endif

                         call SNO_read_map_1d( basename, p, 1,        & ! [IN]
                                               ainfo(n)%varname,      & ! [IN]
                                               ainfo(n)%datatype,     & ! [IN]
                                               gin1,                  & ! [IN]
                                               gout1,                 & ! [IN]
                                               localmap_1d     (:,:), & ! [IN]
                                               ainfo(n)%AXIS_1d(:)    ) ! [INOUT]

!                          if ( debug ) then
!                             if( IO_L ) write(IO_FID_LOG,*) 'AXIS_1d'
!                             do i = 1, gout1
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3.3,F10.1)') i, ainfo(n)%AXIS_1d(i)
!                             enddo
!                          endif

                         deallocate( localmap_1d )
                      endif

                   elseif( ainfo(n)%dim_rank == 2 ) then

                      if    ( ainfo(n)%dim_name(1) == 'x'  ) then
                         gin1    = ngrids_x_in
                         gout1   = ngrids_x_out
                         stgin1  = 0
                         stgout1 = 0
                      elseif( ainfo(n)%dim_name(1) == 'xh' ) then
                         gin1    = ngrids_xh_in
                         gout1   = ngrids_xh_out
                         stgin1  = staggered_x_in
                         stgout1 = staggered_x_out
                      else
                         gin1    = ainfo(n)%dim_size(1)
                         gout1   = ainfo(n)%dim_size(1)
                         stgin1  = 0
                         stgout1 = 0
                      endif

                      if    ( ainfo(n)%dim_name(2) == 'y'  ) then
                         gin2    = ngrids_y_in
                         gout2   = ngrids_y_out
                         stgin2  = 0
                         stgout2 = 0
                      elseif( ainfo(n)%dim_name(2) == 'yh' ) then
                         gin2    = ngrids_yh_in
                         gout2   = ngrids_yh_out
                         stgin2  = staggered_y_in
                         stgout2 = staggered_y_out
                      else
                         gin2    = ainfo(n)%dim_size(2)
                         gout2   = ainfo(n)%dim_size(2)
                         stgin2  = 0
                         stgout2 = 0
                      endif

                      allocate( localmap_2d(gout1,gout2,3) )

                      do j = 1, ngrids_y_out
                      do i = 1, ngrids_x_out
                         localmap_2d(i+stgout1,j+stgout2,I_map_p) = localmap(i,j,I_map_p)
                         localmap_2d(i+stgout1,j+stgout2,I_map_i) = localmap(i,j,I_map_i) + stgin1
                         localmap_2d(i+stgout1,j+stgout2,I_map_j) = localmap(i,j,I_map_j) + stgin2
                      enddo
                      enddo

                      if ( stgout1 > 0 ) then
                         if ( stgin1 > 0 ) then
                            do j = 1, ngrids_y_out
                            do i = 1, stgout1
                               localmap_2d(i,j+stgout2,I_map_p) = localmap(1,j,I_map_p)
                               localmap_2d(i,j+stgout2,I_map_i) = i
                               localmap_2d(i,j+stgout2,I_map_j) = localmap(1,j,I_map_j) + stgin2
                            enddo
                            enddo
                         else
                            do j = 1, ngrids_y_out
                            do i = 1, stgout1
                               localmap_2d(i,j+stgout2,I_map_p) = -1
                               localmap_2d(i,j+stgout2,I_map_i) = -1
                               localmap_2d(i,j+stgout2,I_map_j) = -1
                            enddo
                            enddo
                         endif
                      endif

                      if ( stgout2 > 0 ) then
                         if ( stgin2 > 0 ) then
                            do j = 1, stgout2
                            do i = 1, ngrids_x_out
                               localmap_2d(i+stgout1,j,I_map_p) = localmap(i,1,I_map_p)
                               localmap_2d(i+stgout1,j,I_map_i) = localmap(i,1,I_map_i) + stgin1
                               localmap_2d(i+stgout1,j,I_map_j) = j
                            enddo
                            enddo
                         else
                            do j = 1, stgout2
                            do i = 1, ngrids_x_out
                               localmap_2d(i+stgout1,j,I_map_p) = -1
                               localmap_2d(i+stgout1,j,I_map_i) = -1
                               localmap_2d(i+stgout1,j,I_map_j) = -1
                            enddo
                            enddo
                         endif
                      endif

                      if ( stgout1 > 0 .AND. stgout2 > 0 ) then
                         if ( stgin1 > 0 .AND. stgin2 > 0 ) then
                            do j = 1, stgout2
                            do i = 1, stgout1
                               localmap_2d(i,j,I_map_p) = localmap(1,1,I_map_p)
                               localmap_2d(i,j,I_map_i) = i
                               localmap_2d(i,j,I_map_j) = j
                            enddo
                            enddo
                         else
                            do j = 1, stgout2
                            do i = 1, stgout1
                               localmap_2d(i,j,I_map_p) = -1
                               localmap_2d(i,j,I_map_i) = -1
                               localmap_2d(i,j,I_map_j) = -1
                            enddo
                            enddo
                         endif
                      endif

!                       if ( debug ) then
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          if( IO_L ) write(IO_FID_LOG,*) '*** Axis No.', n, ' : ', trim(ainfo(n)%varname)
!                          if( IO_L ) write(IO_FID_LOG,*) '*** staggered_x_in  = ', staggered_x_in
!                          if( IO_L ) write(IO_FID_LOG,*) '*** staggered_y_in  = ', staggered_y_in
!                          if( IO_L ) write(IO_FID_LOG,*) '*** staggered_x_out = ', staggered_x_out
!                          if( IO_L ) write(IO_FID_LOG,*) '*** staggered_y_out = ', staggered_y_out
!                          if( IO_L ) write(IO_FID_LOG,*) '*** stgin1          = ', stgin1
!                          if( IO_L ) write(IO_FID_LOG,*) '*** stgout1         = ', stgout1
!                          if( IO_L ) write(IO_FID_LOG,*) '*** stgin2          = ', stgin2
!                          if( IO_L ) write(IO_FID_LOG,*) '*** stgout2         = ', stgout2
!
!                          if( IO_L ) write(IO_FID_LOG,*) 'localmap(rank)'
!                          if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
!                          do i = 1, gout1
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                          enddo
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          do j = 1, gout2
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
!                             do i = 1, gout1
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap_2d(i,j,I_map_p)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                          enddo
!
!                          if( IO_L ) write(IO_FID_LOG,*) 'localmap(i-index)'
!                          if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
!                          do i = 1, gout1
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                          enddo
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          do j = 1, gout2
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
!                             do i = 1, gout1
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap_2d(i,j,I_map_i)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                          enddo
!
!                          if( IO_L ) write(IO_FID_LOG,*) 'localmap(j-index)'
!                          if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
!                          do i = 1, gout1
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                          enddo
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          do j = 1, gout2
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
!                             do i = 1, gout1
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap_2d(i,j,I_map_j)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                          enddo
!                       endif

                      call SNO_read_map_2d( basename, p, 1,          & ! [IN]
                                            ainfo(n)%varname,        & ! [IN]
                                            ainfo(n)%datatype,       & ! [IN]
                                            gin1,  gin2,             & ! [IN]
                                            gout1, gout2,            & ! [IN]
                                            localmap_2d     (:,:,:), & ! [IN]
                                            ainfo(n)%AXIS_2d(:,:)    ) ! [INOUT]

!                       if ( debug ) then
!                          if( IO_L ) write(IO_FID_LOG,*) 'AXIS_2d'
!                          if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
!                          do i = 1, gout1
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                          enddo
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          do j = 1, gout2
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
!                             do i = 1, gout1
!                                if( IO_L ) write(IO_FID_LOG,'(1x,F8.3)',advance='no') ainfo(n)%AXIS_2d(i,j)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                          enddo
!                       endif

                      deallocate( localmap_2d )

                   elseif( ainfo(n)%dim_rank == 3 ) then

                      if ( ainfo(n)%transpose ) then
                         gin1  = ainfo(n)%dim_size(3)
                         gout1 = ainfo(n)%dim_size(3)

                         if    ( ainfo(n)%dim_name(1) == 'x'  ) then
                            gin2    = ngrids_x_in
                            gout2   = ngrids_x_out
                            stgin2  = 0
                            stgout2 = 0
                         elseif( ainfo(n)%dim_name(1) == 'xh' ) then
                            gin2    = ngrids_xh_in
                            gout2   = ngrids_xh_out
                            stgin2  = staggered_x_in
                            stgout2 = staggered_x_out
                         else
                            gin2    = ainfo(n)%dim_size(1)
                            gout2   = ainfo(n)%dim_size(1)
                            stgin2  = 0
                            stgout2 = 0
                         endif

                         if    ( ainfo(n)%dim_name(2) == 'y'  ) then
                            gin3    = ngrids_y_in
                            gout3   = ngrids_y_out
                            stgin3  = 0
                            stgout3 = 0
                         elseif( ainfo(n)%dim_name(2) == 'yh' ) then
                            gin3    = ngrids_yh_in
                            gout3   = ngrids_yh_out
                            stgin3  = staggered_y_in
                            stgout3 = staggered_y_out
                         else
                            gin3    = ainfo(n)%dim_size(2)
                            gout3   = ainfo(n)%dim_size(2)
                            stgin3  = 0
                            stgout3 = 0
                         endif
                      else
                         gin1  = ainfo(n)%dim_size(1)
                         gout1 = ainfo(n)%dim_size(1)

                         if    ( ainfo(n)%dim_name(2) == 'x'  ) then
                            gin2    = ngrids_x_in
                            gout2   = ngrids_x_out
                            stgin2  = 0
                            stgout2 = 0
                         elseif( ainfo(n)%dim_name(2) == 'xh' ) then
                            gin2    = ngrids_xh_in
                            gout2   = ngrids_xh_out
                            stgin2  = staggered_x_in
                            stgout2 = staggered_x_out
                         else
                            gin2    = ainfo(n)%dim_size(2)
                            gout2   = ainfo(n)%dim_size(2)
                            stgin2  = 0
                            stgout2 = 0
                         endif

                         if    ( ainfo(n)%dim_name(3) == 'y'  ) then
                            gin3    = ngrids_y_in
                            gout3   = ngrids_y_out
                            stgin3  = 0
                            stgout3 = 0
                         elseif( ainfo(n)%dim_name(3) == 'yh' ) then
                            gin3    = ngrids_yh_in
                            gout3   = ngrids_yh_out
                            stgin3  = staggered_y_in
                            stgout3 = staggered_y_out
                         else
                            gin3    = ainfo(n)%dim_size(3)
                            gout3   = ainfo(n)%dim_size(3)
                            stgin3  = 0
                            stgout3 = 0
                         endif
                      endif

                      allocate( localmap_2d(gout2,gout3,3) )

                      do j = 1, ngrids_y_out
                      do i = 1, ngrids_x_out
                         localmap_2d(i+stgout2,j+stgout3,I_map_p) = localmap(i,j,I_map_p)
                         localmap_2d(i+stgout2,j+stgout3,I_map_i) = localmap(i,j,I_map_i) + stgin2
                         localmap_2d(i+stgout2,j+stgout3,I_map_j) = localmap(i,j,I_map_j) + stgin3
                      enddo
                      enddo

                      if ( stgout2 > 0 ) then
                         if ( stgin2 > 0 ) then
                            do j = 1, ngrids_y_out
                            do i = 1, stgout2
                               localmap_2d(i,j+stgout3,I_map_p) = localmap(1,j,I_map_p)
                               localmap_2d(i,j+stgout3,I_map_i) = i
                               localmap_2d(i,j+stgout3,I_map_j) = localmap(1,j,I_map_j) + stgin3
                            enddo
                            enddo
                         else
                            do j = 1, ngrids_y_out
                            do i = 1, stgout2
                               localmap_2d(i,j+stgout3,I_map_p) = -1
                               localmap_2d(i,j+stgout3,I_map_i) = -1
                               localmap_2d(i,j+stgout3,I_map_j) = -1
                            enddo
                            enddo
                         endif
                      endif

                      if ( stgout3 > 0 ) then
                         if ( stgin3 > 0 ) then
                            do j = 1, stgout3
                            do i = 1, ngrids_x_out
                               localmap_2d(i+stgout2,j,I_map_p) = localmap(i,1,I_map_p)
                               localmap_2d(i+stgout2,j,I_map_i) = localmap(i,1,I_map_i) + stgin2
                               localmap_2d(i+stgout2,j,I_map_j) = j
                            enddo
                            enddo
                         else
                            do j = 1, stgout3
                            do i = 1, ngrids_x_out
                               localmap_2d(i+stgout2,j,I_map_p) = -1
                               localmap_2d(i+stgout2,j,I_map_i) = -1
                               localmap_2d(i+stgout2,j,I_map_j) = -1
                            enddo
                            enddo
                         endif
                      endif

                      if ( stgout2 > 0 .AND. stgout3 > 0 ) then
                         if ( stgin2 > 0 .AND. stgin3 > 0 ) then
                            do j = 1, stgout3
                            do i = 1, stgout2
                               localmap_2d(i,j,I_map_p) = localmap(1,1,I_map_p)
                               localmap_2d(i,j,I_map_i) = i
                               localmap_2d(i,j,I_map_j) = j
                            enddo
                            enddo
                         else
                            do j = 1, stgout3
                            do i = 1, stgout2
                               localmap_2d(i,j,I_map_p) = -1
                               localmap_2d(i,j,I_map_i) = -1
                               localmap_2d(i,j,I_map_j) = -1
                            enddo
                            enddo
                         endif
                      endif

                      call SNO_read_map_3d( basename, p, 1,          & ! [IN]
                                            ainfo(n)%varname,        & ! [IN]
                                            ainfo(n)%datatype,       & ! [IN]
                                            ainfo(n)%transpose,      & ! [IN]
                                            gin1,  gin2,  gin3,      & ! [IN]
                                            gout1, gout2, gout3,     & ! [IN]
                                            localmap_2d     (:,:,:), & ! [IN]
                                            ainfo(n)%AXIS_3d(:,:,:)  ) ! [INOUT]

!                       if ( debug ) then
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          if( IO_L ) write(IO_FID_LOG,*) '*** Axis No.', n, ' : ', trim(ainfo(n)%varname)
!
!                          if( IO_L ) write(IO_FID_LOG,*) 'AXIS_3d'
!                          if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
!                          do i = 1, gout2
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                          enddo
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          do j = 1, gout3
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
!                             do i = 1, gout2
!                                if( IO_L ) write(IO_FID_LOG,'(1x,F8.3)',advance='no') ainfo(n)%AXIS_3d(1,i,j)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                          enddo
!
!                          if( IO_L ) write(IO_FID_LOG,*) 'localmap(rank)'
!                          if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
!                          do i = 1, gout2
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                          enddo
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          do j = 1, gout3
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
!                             do i = 1, gout2
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap_2d(i,j,I_map_p)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                          enddo
!
!                          if( IO_L ) write(IO_FID_LOG,*) 'localmap(i-index)'
!                          if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
!                          do i = 1, gout2
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                          enddo
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          do j = 1, gout3
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
!                             do i = 1, gout2
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap_2d(i,j,I_map_i)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                          enddo
!
!                          if( IO_L ) write(IO_FID_LOG,*) 'localmap(j-index)'
!                          if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
!                          do i = 1, gout2
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
!                          enddo
!                          if( IO_L ) write(IO_FID_LOG,*)
!                          do j = 1, gout3
!                             if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
!                             do i = 1, gout2
!                                if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap_2d(i,j,I_map_j)
!                             enddo
!                             if( IO_L ) write(IO_FID_LOG,*)
!                          enddo
!                       endif

                      deallocate( localmap_2d )

                   endif ! dim rank?
                endif ! regrid?
             enddo ! axis loop

          endif ! readflag?
       enddo ! px
    enddo ! py

    return
  end subroutine SNO_axis_read

  !-----------------------------------------------------------------------------
  subroutine SNO_axis_define( &
       fid,   &
       naxis, &
       ainfo, &
       debug  )
    use scale_file, only: &
       FILE_Def_Axis,                &
       FILE_Def_AssociatedCoordinate
    use mod_sno_h, only: &
       axisinfo
    implicit none

    integer,          intent(in)    :: fid                      ! number of axis variables           (input)
    integer,          intent(in)    :: naxis                    ! number of axis variables           (input)
    type(axisinfo),   intent(in)    :: ainfo(naxis)             ! axis information                   (input)
    logical,          intent(in)    :: debug

    integer  :: dsize, drank
    integer  :: n

    intrinsic size
    !---------------------------------------------------------------------------

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** [SNO_axis_define] define axis'
    endif

    do n = 1, naxis
       if ( ainfo(n)%dim_rank == 1 ) then
          dsize = size(ainfo(n)%AXIS_1d(:))

          call FILE_Def_Axis( fid,                  & ! [IN]
                              ainfo(n)%varname,     & ! [IN]
                              ainfo(n)%description, & ! [IN]
                              ainfo(n)%units,       & ! [IN]
                              ainfo(n)%dim_name(1), & ! [IN]
                              ainfo(n)%datatype,    & ! [IN]
                              dsize                 ) ! [IN]
       else
          drank = ainfo(n)%dim_rank

          call FILE_Def_AssociatedCoordinate( fid,                        & ! [IN]
                                              ainfo(n)%varname,           & ! [IN]
                                              ainfo(n)%description,       & ! [IN]
                                              ainfo(n)%units,             & ! [IN]
                                              ainfo(n)%dim_name(1:drank), & ! [IN]
                                              ainfo(n)%datatype           ) ! [IN]
       endif
    enddo

    return
  end subroutine SNO_axis_define

  !-----------------------------------------------------------------------------
  subroutine SNO_axis_write( &
       fid,   &
       naxis, &
       ainfo, &
       debug  )
    use scale_file, only: &
       FILE_Write_Axis,                &
       FILE_Write_AssociatedCoordinate
    use mod_sno_h, only: &
       axisinfo
    implicit none

    integer,          intent(in)    :: fid                      ! number of axis variables           (input)
    integer,          intent(in)    :: naxis                    ! number of axis variables           (input)
    type(axisinfo),   intent(in)    :: ainfo(naxis)             ! axis information                   (input)
    logical,          intent(in)    :: debug

    integer, parameter :: start(3) = 1
    integer  :: n
    !---------------------------------------------------------------------------

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** [SNO_axis_write] write axis'
    endif

    do n = 1, naxis
       if    ( ainfo(n)%dim_rank == 1 ) then

          call FILE_Write_Axis( fid,                & ! [IN]
                                ainfo(n)%varname,   & ! [IN]
                                ainfo(n)%AXIS_1d(:) ) ! [IN]

       elseif( ainfo(n)%dim_rank == 2 ) then

          call FILE_Write_AssociatedCoordinate( fid,                   & ! [IN]
                                                ainfo(n)%varname,      & ! [IN]
                                                ainfo(n)%AXIS_2d(:,:), & ! [IN]
                                                start(2:3)             ) ! [IN]

       elseif( ainfo(n)%dim_rank == 3 ) then

          call FILE_Write_AssociatedCoordinate( fid,                     & ! [IN]
                                                ainfo(n)%varname,        & ! [IN]
                                                ainfo(n)%AXIS_3d(:,:,:), & ! [IN]
                                                start(1:3)               ) ! [IN]

       endif
    enddo

    return
  end subroutine SNO_axis_write

end module mod_sno_axis