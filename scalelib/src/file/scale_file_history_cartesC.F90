!-------------------------------------------------------------------------------
!> module file / history_cartesC
!!
!! @par Description
!!          History output module for the cartesianC grid
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_file_history_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_ocean_grid_cartesC_index
  use scale_land_grid_cartesC_index
  use scale_urban_grid_cartesC_index
  use scale_process, only: &
     PRC_abort
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: FILE_HISTORY_CARTESC_setup
  public :: FILE_HISTORY_CARTESC_Set_Pres

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: FILE_HISTORY_CARTESC_set_dims
  private :: FILE_HISTORY_CARTESC_set_axes
  private :: FILE_HISTORY_CARTESC_set_axes_attributes

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,          parameter :: nzs = 3
  character(len=8), parameter :: zs(nzs) = (/ "model   ", &
                                              "z       ", &
                                              "pressure"  /)

  integer               :: FILE_HISTORY_CARTESC_PRES_nlayer = 0
  real(RP), allocatable :: FILE_HISTORY_CARTESC_PRES_val(:)

  integer  :: im,   jm,  km
  integer  :: ims,  ime
  integer  :: jms,  jme
  integer  :: imh,  jmh
  integer  :: imsh, jmsh

  integer  :: FILE_HISTORY_CARTESCORY_STARTDATE(6) !< start time [YYYY MM DD HH MM SS]
  real(DP) :: FILE_HISTORY_CARTESCORY_STARTMS      !< subsecond part of start time [millisec]

  logical  :: FILE_HISTORY_CARTESC_BOUNDARY = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FILE_HISTORY_CARTESC_setup
    use scale_file_h, only: &
       FILE_HSHORT
    use scale_file_history, only: &
       FILE_HISTORY_Setup, &
       FILE_HISTORY_Set_NowDate, &
       FILE_HISTORY_truncate_1D, &
       FILE_HISTORY_truncate_2D, &
       FILE_HISTORY_truncate_3D
    use scale_process, only: &
       PRC_masterrank, &
       PRC_myrank,     &
       PRC_LOCAL_COMM_WORLD
    use scale_rm_process, only: &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y, &
       PRC_HAS_W,      &
       PRC_HAS_S
    use scale_time, only: &
       TIME_NOWDATE,     &
       TIME_NOWMS,       &
       TIME_NOWSTEP,     &
       TIME_DTSEC,       &
       TIME_STARTDAYSEC
    use scale_calendar, only: &
       CALENDAR_get_name
    use scale_interp_vert, only: &
       INTERP_VERT_alloc_pres
    implicit none

    integer, parameter :: nlayer_max = 300
    real(RP)           :: FILE_HISTORY_CARTESC_PRES(nlayer_max) !> pressure level to output [hPa]

    NAMELIST / PARAM_FILE_HISTORY_CARTESC / &
       FILE_HISTORY_CARTESC_PRES_nlayer, &
       FILE_HISTORY_CARTESC_PRES,        &
       FILE_HISTORY_CARTESC_BOUNDARY

    character(len=H_MID) :: FILE_HISTORY_CARTESCORY_H_TITLE = 'SCALE-RM FILE_HISTORY_CARTESC OUTPUT' !< title of the output file
    character(len=H_MID) :: FILE_HISTORY_CARTESCORY_T_SINCE

    character(len=FILE_HSHORT) :: calendar
    real(DP) :: start_daysec
    integer  :: ierr
    integer  :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[FILE_HISTORY_CARTESCORY] / Categ[ATMOS-RM IO] / Origin[SCALElib]'

    FILE_HISTORY_CARTESC_PRES(:) = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_FILE_HISTORY_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_FILE_HISTORY_CARTESC. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_FILE_HISTORY_CARTESC)


    ! check pressure coordinate
    if ( FILE_HISTORY_CARTESC_PRES_nlayer > 0 ) then
       if ( FILE_HISTORY_CARTESC_PRES_nlayer > nlayer_max ) then
          write(*,'(a,i3)') 'xxx FILE_HISTORY_CARTESC_PRES_nlayer must be <= ', nlayer_max
          call PRC_abort
       end if
       allocate( FILE_HISTORY_CARTESC_PRES_val(FILE_HISTORY_CARTESC_PRES_nlayer) )

       do k = 1, FILE_HISTORY_CARTESC_PRES_nlayer
          if ( FILE_HISTORY_CARTESC_PRES(k) <= 0.0_RP ) then
             write(*,'(a,i3,f7.1)') 'xxx Invalid value found in pressure coordinate! (k,value)=', k, FILE_HISTORY_CARTESC_PRES(k)
             call PRC_abort
          elseif ( FILE_HISTORY_CARTESC_PRES(k+1) >= FILE_HISTORY_CARTESC_PRES(k) ) then
             write(*,'(a,i3,2f7.1)') 'xxx The value of pressure coordinate must be descending order! ', &
                  '(k,value[k],value[k+1])=', k, FILE_HISTORY_CARTESC_PRES(k), FILE_HISTORY_CARTESC_PRES(k+1)
             call PRC_abort
          endif
          FILE_HISTORY_CARTESC_PRES_val(k) = FILE_HISTORY_CARTESC_PRES(k) * 100.0_RP ! [hPa->Pa]
       enddo

       call INTERP_VERT_alloc_pres( FILE_HISTORY_CARTESC_PRES_nlayer, IA, JA ) ! [IN]
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** FILE_HISTORY_CARTESC_PRES_nlayer is not set.'
       if( IO_L ) write(IO_FID_LOG,*) '*** Output with pressure coordinate is disabled'
    endif



    FILE_HISTORY_CARTESCORY_STARTDATE(:) = TIME_NOWDATE
    FILE_HISTORY_CARTESCORY_STARTMS      = TIME_NOWMS

    start_daysec = TIME_STARTDAYSEC
    if ( TIME_NOWDATE(1) > 0 ) then
       write(FILE_HISTORY_CARTESCORY_T_SINCE,'(I4.4,5(A1,I2.2))') TIME_NOWDATE(1), &
                                                             '-', TIME_NOWDATE(2), &
                                                             '-', TIME_NOWDATE(3), &
                                                             ' ', TIME_NOWDATE(4), &
                                                             ':', TIME_NOWDATE(5), &
                                                             ':', TIME_NOWDATE(6)
       start_daysec = TIME_NOWMS
    else
       FILE_HISTORY_CARTESCORY_T_SINCE = ''
    endif

    if ( FILE_HISTORY_CARTESC_BOUNDARY ) then
       ims  = ISB
       ime  = IEB
       jms  = JSB
       jme  = JEB

       imsh = ims
       jmsh = jms

       im   = IMAXB
       jm   = JMAXB
       imh  = im
       jmh  = jm
    else
       ims  = IS
       ime  = IE
       jms  = JS
       jme  = JE

       if ( PRC_HAS_W .OR. PRC_PERIODIC_X ) then
          imsh = ims
       else
          imsh = ims - 1 ! including i = IS-1
       endif
       if ( PRC_HAS_S .OR. PRC_PERIODIC_Y ) then
          jmsh = jms
       else
          jmsh = jms - 1 ! include j = JS-1
       endif

       im   = ime - ims  + 1
       jm   = jme - jms  + 1
       imh  = ime - imsh + 1
       jmh  = jme - jmsh + 1
    endif

    ! get calendar name
    call CALENDAR_get_name( calendar )

    call FILE_HISTORY_Setup( &
         FILE_HISTORY_CARTESCORY_H_TITLE,              & ! [IN]
         H_SOURCE, H_INSTITUTE,                        & ! [IN]
         start_daysec, TIME_DTSEC,                     & ! [IN]
         time_since = FILE_HISTORY_CARTESCORY_T_SINCE, & ! [IN]
         calendar = calendar,                          & ! [IN]
         default_zcoord = 'model',                     & ! [IN]
         myrank = PRC_myrank                           ) ! [IN]

    call FILE_HISTORY_Set_NowDate( TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP )

    call FILE_HISTORY_CARTESC_set_dims

    call FILE_HISTORY_CARTESC_set_axes

    call FILE_HISTORY_CARTESC_set_axes_attributes

    FILE_HISTORY_truncate_1D => FILE_HISTORY_CARTESC_truncate_1D
    FILE_HISTORY_truncate_2D => FILE_HISTORY_CARTESC_truncate_2D
    FILE_HISTORY_truncate_3D => FILE_HISTORY_CARTESC_truncate_3D

    return
  end subroutine FILE_HISTORY_CARTESC_setup

  !-----------------------------------------------------------------------------
  !> set hydrostatic pressure for pressure coordinate
  !-----------------------------------------------------------------------------
  subroutine FILE_HISTORY_CARTESC_Set_Pres( &
       PRES,    &
       PRESH,   &
       SFC_PRES )
    use scale_interp_vert, only: &
       INTERP_VERT_setcoef_pres
    implicit none

    real(RP), intent(in) :: PRES    (:,:,:) ! pressure at the full level [Pa]
    real(RP), intent(in) :: PRESH   (:,:,:) ! pressure at the half level [Pa]
    real(RP), intent(in) :: SFC_PRES(  :,:) ! surface pressure           [Pa]
    !---------------------------------------------------------------------------

    if ( FILE_HISTORY_CARTESC_PRES_nlayer > 0 ) then
       call INTERP_VERT_setcoef_pres( FILE_HISTORY_CARTESC_PRES_nlayer, & ! [IN]
                                      KA, KS, KE,                       & ! [IN]
                                      IA, IS, IE,                       & ! [IN]
                                      JA, JS, JE,                       & ! [IN]
                                      PRES         (:,:,:),             & ! [IN]
                                      PRESH        (:,:,:),             & ! [IN]
                                      SFC_PRES     (:,:)  ,             & ! [IN]
                                      FILE_HISTORY_CARTESC_PRES_val(:)  ) ! [IN]
    endif

    return
  end subroutine FILE_HISTORY_CARTESC_Set_Pres

  ! private routines

  !-----------------------------------------------------------------------------
  !> set dimension information
  !-----------------------------------------------------------------------------
  subroutine FILE_HISTORY_CARTESC_set_dims
    use scale_process, only: &
       PRC_myrank
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_HAS_W, &
       PRC_HAS_S
    use scale_file_history, only: &
       FILE_HIStoRY_set_dim
    use scale_mapprojection, only: &
       MAPPROJECTION_get_attributes
    implicit none

    character(len=H_SHORT) :: mapping

    character(len=H_SHORT) :: dims(3,3)

    integer :: start(3,3), count(3,3)
    integer :: xs, xc, ys, yc
    !---------------------------------------------------------------------------

    ! get start and count for x and y
    if ( FILE_HISTORY_CARTESC_BOUNDARY ) then
       xs = ISGB
       ys = JSGB
       xc = IMAXB
       yc = JMAXB
    else
       ! for the case the shared-file contains no halos
       xs = PRC_2Drank(PRC_myrank,1) * IMAX + 1 ! no IHALO
       xc = IMAX
       ys = PRC_2Drank(PRC_myrank,2) * JMAX + 1 ! no JHALO
       yc = JMAX
    end if

    ! get mapping name
    call MAPPROJECTION_get_attributes( mapping )


    !  Vertical 1D
    start(1,1) = 1
    dims (1,1) = "z"
    count(1,1) = KMAX
    call FILE_HISTORY_Set_Dim( "Z",   1, 1, dims(:,:), zs(:), start(:,:), count(:,:) ) ! [IN]
    dims (1,1) = "zh"
    count(1,1) = KMAX + 1
    call FILE_HISTORY_Set_Dim( "ZH",  1, 1, dims(:,:), zs(:), start(:,:), count(:,:) ) ! [IN]

    dims (1,1) = "oz"
    count(1,1) = OKMAX
    call FILE_HISTORY_Set_Dim( "OZ",  1, 1, dims(:,:), zs(:), start(:,:), count(:,:) ) ! [IN]
    dims (1,1) = "ozh"
    count(1,1) = OKMAX + 1
    call FILE_HISTORY_Set_Dim( "OZH", 1, 1, dims(:,:), zs(:), start(:,:), count(:,:) ) ! [IN]

    dims (1,1) = "lz"
    count(1,1) = LKMAX
    call FILE_HISTORY_Set_Dim( "LZ",  1, 1, dims(:,:), zs(:), start(:,:), count(:,:) ) ! [IN]
    dims (1,1) = "lzh"
    count(1,1) = LKMAX + 1
    call FILE_HISTORY_Set_Dim( "LZH", 1, 1, dims(:,:), zs(:), start(:,:), count(:,:) ) ! [IN]

    dims (1,1) = "uz"
    count(1,1) = UKMAX
    call FILE_HISTORY_Set_Dim( "UZ",  1, 1, dims(:,:), zs(:), start(:,:), count(:,:) ) ! [IN]
    dims (1,1) = "uzh"
    count(1,1) = UKMAX + 1
    call FILE_HISTORY_Set_Dim( "UZH", 1, 1, dims(:,:), zs(:), start(:,:), count(:,:) ) ! [IN]

    ! X, Y
    start(1,:) = xs
    start(2,:) = ys
    dims (1,:) = 'lon'
    dims (2,:) = 'lat'
    count(1,:) = xc
    count(2,:) = yc
    call FILE_HISTORY_Set_Dim( "XY", 2, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area", location="face" ) ! [IN]

    start(3,:) = 1
    dims (3,:) = (/ "height    ", "z         ", "pressure  " /)
    count(3,:) = (/ KMAX,   KMAX,   FILE_HISTORY_CARTESC_PRES_nlayer /)
    call FILE_HISTORY_Set_Dim( "ZXY",  3, nzs, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, & ! [IN]
                               area="cell_area", area_x="cell_area_xyz_x", area_y="cell_area_xyz_y", volume="cell_volume", location="face" ) ! [IN]
    dims (3,:) = (/ "height_xyw", "zh        ", "pressure  " /)
    count(3,:) = (/ KMAX+1, KMAX+1, FILE_HISTORY_CARTESC_PRES_nlayer /)
    call FILE_HISTORY_Set_Dim( "ZHXY", 3, nzs, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area", volume="cell_volume_xyw", location="face" )

    dims (3,1) = "oz"
    count(3,1) = OKMAX
    call FILE_HISTORY_Set_Dim( "OXY",  3, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area", volume="cell_volume_xyo", location="face", grid="ocean" ) ! [IN]
    dims (3,1) = "ozh"
    count(3,1) = OKMAX + 1
    call FILE_HISTORY_Set_Dim( "OHXY", 3, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area", location="face", grid="ocean" ) ! [IN]
    dims (3,1) = "lz"
    count(3,1) = LKMAX
    call FILE_HISTORY_Set_Dim( "LXY",  3, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area", volume="cell_volume_xyl", location="face", grid="land" ) ! [IN]
    dims (3,1) = "lzh"
    count(3,1) = LKMAX + 1
    call FILE_HISTORY_Set_Dim( "LHXY", 3, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area", location="face", grid="land" ) ! [IN]
    dims (3,1) = "uz"
    count(3,1) = UKMAX
    call FILE_HISTORY_Set_Dim( "UXY",  3, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area", volume="cell_volume_xyu", grid="urban" ) ! [IN]
    dims (3,1) = "uzh"
    count(3,1) = UKMAX + 1
    call FILE_HISTORY_Set_Dim( "UHXY", 3, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area", location="face", grid="urban" ) ! [IN]

    ! XH, Y
    dims(1,:) = 'lon_uy'
    dims(2,:) = 'lat_uy'
    if ( PRC_HAS_W ) then
       start(1,:) = xs+1
       count(1,:) = xc
    else
       start(1,:) = xs
       count(1,:) = xc+1
    endif
    call FILE_HISTORY_Set_Dim( "XHY", 2, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area_uy", location="edge1" ) ! [IN]

    dims (3,:) = (/ "height_uyz", "z         ", "pressure  " /)
    count(3,:) = (/ KMAX,   KMAX,   FILE_HISTORY_CARTESC_PRES_nlayer /)
    call FILE_HISTORY_Set_Dim( "ZXHY",  3, nzs, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area_uy", volume="cell_volume_uyz", location="edge1" ) ! [IN]
    dims (3,:) = (/ "height_uyw", "zh        ", "pressure  " /)
    count(3,:) = (/ KMAX+1, KMAX+1, FILE_HISTORY_CARTESC_PRES_nlayer /)
    call FILE_HISTORY_Set_Dim( "ZHXHY", 3, nzs, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area_uy", location="edge1" ) ! [IN]

    ! X, YH
    dims (1,:) = 'lon_xv'
    dims (2,:) = 'lat_xv'
    start(1,:) = xs
    count(1,:) = xc
    if ( PRC_HAS_S ) then
       start(2,:) = ys+1
       count(2,:) = yc
    else
       start(2,:) = ys
       count(2,:) = yc+1
    endif
    call FILE_HISTORY_Set_Dim( "XYH", 2, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area_xv", location="edge2" ) ! [IN]

    dims (3,:) = (/ "height_xvz", "z         ", "pressure  " /)
    count(3,:) = (/ KMAX,   KMAX,   FILE_HISTORY_CARTESC_PRES_nlayer /)
    call FILE_HISTORY_Set_Dim( "ZXYH",  3, nzs, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area_xv", volume="cell_volume_xvz", location="edge2" ) ! [IN]
    dims (3,:) = (/ "height_xvw", "zh        ", "pressure  " /)
    count(3,:) = (/ KMAX+1, KMAX+1, FILE_HISTORY_CARTESC_PRES_nlayer /)
    call FILE_HISTORY_Set_Dim( "ZHXYH", 3, nzs, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, area="cell_area_xv", location="edge2" ) ! [IN]

    ! XH, YH
    dims(1,:) = 'lon_uv'
    dims(2,:) = 'lat_uv'
    if ( PRC_HAS_W ) then
       start(1,:) = xs+1
       count(1,:) = xc
    else
       start(1,:) = xs
       count(1,:) = xc+1
    endif
    call FILE_HISTORY_Set_Dim( "XHYH", 2, 1, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, location="face" ) ! [IN]

    dims (3,:) = (/ "height_uvz", "z         ", "pressure  " /)
    count(3,:) = (/ KMAX,   KMAX,   FILE_HISTORY_CARTESC_PRES_nlayer /)
    call FILE_HISTORY_Set_Dim( "ZXHYH",  3, nzs, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, & ![IN]
                               area="cell_area_uv", area_x="cell_area_uvz_x", area_y="cell_area_uvz_y", location="node" ) ! [IN]
    dims (3,:) = (/ "height_uvw", "zh        ", "pressure  " /)
    count(3,:) = (/ KMAX+1, KMAX+1, FILE_HISTORY_CARTESC_PRES_nlayer /)
    call FILE_HISTORY_Set_Dim( "ZHXHYH", 3, nzs, dims(:,:), zs(:), start(:,:), count(:,:), mapping=mapping, location="node" ) ! [IN]

    return
  end subroutine FILE_HISTORY_CARTESC_set_dims

  !-----------------------------------------------------------------------------
  !> truncate 1D data to history buffer
  subroutine FILE_HISTORY_CARTESC_truncate_1D( &
       src, &
       dim_type, zcoord, fill_halo, &
       dst )
    implicit none

    real(RP),         intent(in) :: src(:)
    character(len=*), intent(in) :: dim_type
    character(len=*), intent(in) :: zcoord
    logical,          intent(in) :: fill_halo ! ignored

    real(DP), intent(out) :: dst(:)

    integer  :: ksize
    integer  :: kstart
    integer  :: k
    !---------------------------------------------------------------------------

    ! select dimension
    select case ( dim_type )
    case ('Z')
        ksize  = KMAX
        kstart = KS
    case ('ZH')
       ksize  = KMAX+1
       kstart = KS-1
    case ('OZ')
       ksize  = OKMAX
       kstart = OKS
    case ('OZH')
       ksize  = OKMAX+1
       kstart = OKS-1
    case ('LZ')
       ksize  = LKMAX
       kstart = LKS
    case ('LZH')
       ksize  = LKMAX+1
       kstart = LKS-1
    case ('UZ')
       ksize  = UKMAX
       kstart = UKS
    case ('UZH')
       ksize  = UKMAX+1
       kstart = UKS-1
    case default
       write(*,*) 'xxx [FILE_HISTORY_CARTESC_truncate_1D] dim_type is invalid: ', trim(dim_type)
       call PRC_abort
    end select

    do k = 1, ksize
       dst(k) = src(kstart+k-1)
    enddo

    return
  end subroutine FILE_HISTORY_CARTESC_truncate_1D

  !-----------------------------------------------------------------------------
  !> truncate 2D data to history buffer
  subroutine FILE_HISTORY_CARTESC_truncate_2D( &
       src, &
       dim_type, zcoord, fill_halo, &
       dst )
    use scale_file_h, only: &
       RMISS => FILE_RMISS
    implicit none

    real(RP),         intent(in) :: src(:,:)
    character(len=*), intent(in) :: dim_type
    character(len=*), intent(in) :: zcoord    ! ignored
    logical,          intent(in) :: fill_halo

    real(DP), intent(out) :: dst(:)

    integer :: isize, jsize
    integer :: istart, jstart
    integer :: i, j
    !---------------------------------------------------------------------------

    ! select dimension
    select case( dim_type )
    case ( 'XY', 'XYH' )
       isize  = im
       istart = ims
    case ( 'XHY', 'XHYH' )
       isize  = imh
       istart = imsh
    case default
       write(*,*) 'xxx [FILE_HISTORY_CARTESC_truncate_2D] dim_type is invalid: ', trim(dim_type)
       call PRC_abort
    end select

    select case ( dim_type )
    case ( 'XY', 'XHY' )
       jsize  = jm
       jstart = jms
    case ( 'XYH', 'XHYH' )
       jsize  = jmh
       jstart = jmsh
    case default
       write(*,*) 'xxx [FILE_HISTORY_CARTESC_truncate_2D] dim_type is invalid: ', trim(dim_type)
       call PRC_abort
    end select

    !$omp parallel do
    do j = 1, jsize
    do i = 1, isize
       dst((j-1)*isize+i) = src(istart+i-1,jstart+j-1)
    enddo
    enddo

    if ( fill_halo ) then
       ! W halo
       do j = 1, jsize
       do i = 1, IS-istart
          dst((j-1)*isize+i) = RMISS
       enddo
       enddo
       ! E halo
       do j = 1, jsize
       do i = IE-istart+2, ime-istart+1
          dst((j-1)*isize+i) = RMISS
       enddo
       enddo
       ! S halo
       do j = 1, JS-jstart
       do i = 1, isize
          dst((j-1)*isize+i) = RMISS
       enddo
       enddo
       ! N halo
       do j = JE-jstart+2, jme-jstart+1
       do i = 1, isize
          dst((j-1)*isize+i) = RMISS
       enddo
       enddo

    end if

    return
  end subroutine FILE_HISTORY_CARTESC_truncate_2D

  !-----------------------------------------------------------------------------
  !> truncate 3D data to history buffer
  subroutine FILE_HISTORY_CARTESC_truncate_3D( &
       src, &
       dim_type, zcoord, fill_halo, &
       dst )
    use scale_file_h, only: &
       RMISS => FILE_RMISS
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CZ, &
       ATMOS_GRID_CARTESC_FZ
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_CZ, &
       ATMOS_GRID_CARTESC_REAL_FZ
    use scale_interp_vert, only: &
       INTERP_VERT_xi2z,   &
       INTERP_VERT_xi2p,   &
       INTERP_VERT_xih2zh, &
       INTERP_VERT_xih2p,  &
       INTERP_available
    implicit none

    real(RP),         intent(in)  :: src(:,:,:)
    character(len=*), intent(in)  :: dim_type
    character(len=*), intent(in)  :: zcoord
    logical,          intent(in)  :: fill_halo
    real(DP),         intent(out) :: dst(:)

    real(RP) :: src_Z(KA,IA,JA)
    real(RP) :: src_P(FILE_HISTORY_CARTESC_PRES_nlayer,IA,JA)

    integer  :: isize,  jsize,  ksize
    integer  :: istart, jstart, kstart
    integer  :: i, j, k
    !---------------------------------------------------------------------------

    ! select dimension
    if ( index( dim_type, 'XH' ) == 0 ) then
       isize  = im
       istart = ims
    else
       isize  = imh
       istart = imsh
    end if

    if ( index( dim_type, 'YH' ) == 0 ) then
       jsize  = jm
       jstart = jms
    else
       jsize  = jmh
       jstart = jmsh
    end if

    select case( dim_type(1:1) )
    case ( 'Z' )
       ksize  = KMAX
       kstart = KS
    case('O')
       ksize  = OKMAX
       kstart = OKS
    case('L')
       ksize  = LKMAX
       kstart = LKS
    case('U')
       ksize  = UKMAX
       kstart = UKS
    case default
       write(*,*) 'xxx [FILE_HISTORY_CARTESC_truncate_3D] dim_type is invalid: ', trim(dim_type)
       call PRC_abort
    end select
    if ( dim_type(2:2) == 'H' ) then
       ksize = ksize + 1
       kstart = kstart - 1
    end if


    if ( ksize == KMAX .and. zcoord == "z" .and. INTERP_available ) then ! z*->z interpolation (full level)

       call PROF_rapstart('FILE_O_interp', 2)
       call INTERP_VERT_xi2z( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                              ATMOS_GRID_CARTESC_CZ(:),          & ! [IN]
                              ATMOS_GRID_CARTESC_REAL_CZ(:,:,:), & ! [IN]
                              src    (:,:,:), & ! [IN]
                              src_Z  (:,:,:)  ) ! [OUT]
       call PROF_rapend  ('FILE_O_interp', 2)

       !$omp parallel do
       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = src_Z(kstart+k-1,istart+i-1,jstart+j-1)
       enddo
       enddo
       enddo

    else if( ksize == KMAX+1 .and. zcoord == "z" .and. INTERP_available ) then ! z*->z interpolation (half level)


       call PROF_rapstart('FILE_O_interp', 2)
       call INTERP_VERT_xih2zh( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                ATMOS_GRID_CARTESC_FZ(:),          & ! [IN]
                                ATMOS_GRID_CARTESC_REAL_FZ(:,:,:), & ! [IN]
                                src    (:,:,:), & ! [IN]
                                src_Z  (:,:,:)  ) ! [OUT]
       call PROF_rapend  ('FILE_O_interp', 2)

       !$omp parallel do
       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = src_Z(kstart+k-1,istart+i-1,jstart+j-1)
       enddo
       enddo
       enddo

    elseif( ksize == KMAX .and. zcoord == "pressure" ) then ! z*->p interpolation (full level)
       ksize = FILE_HISTORY_CARTESC_PRES_nlayer
       if ( ksize == 0 ) then
          write(*,*) 'xxx FILE_HISTORY_CARTESC_PRES_nlayer must be set to output variable with the pressure coordinate'
          call PRC_abort
       end if

       call PROF_rapstart('FILE_O_interp', 2)
       call INTERP_VERT_xi2p( FILE_HISTORY_CARTESC_PRES_nlayer, & ! [IN]
                              KA,                               & ! [IN]
                              IA, ISB, IEB,                     & ! [IN]
                              JA, JSB, JEB,                     & ! [IN]
                              src  (:,:,:),                     & ! [IN]
                              src_P(:,:,:)                      ) ! [OUT]
       call PROF_rapend  ('FILE_O_interp', 2)

       !$omp parallel do
       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = src_P(k,istart+i-1,jstart+j-1)
       enddo
       enddo
       enddo

    elseif( ksize == KMAX+1 .and. zcoord == "pressure" ) then ! z*->p interpolation (half level)
       ksize = FILE_HISTORY_CARTESC_PRES_nlayer
       if ( ksize == 0 ) then
          write(*,*) 'xxx FILE_HISTORY_CARTESC_PRES_nlayer must be set to output variable with the pressure coordinate'
          call PRC_abort
       end if

       call PROF_rapstart('FILE_O_interp', 2)
       call INTERP_VERT_xih2p( FILE_HISTORY_CARTESC_PRES_nlayer, & ! [IN]
                               KA,                               & ! [IN]
                               IA, ISB, IEB,                     & ! [IN]
                               JA, JSB, JEB,                     & ! [IN]
                               src  (:,:,:),                     & ! [IN]
                               src_P(:,:,:)                      ) ! [OUT]
       call PROF_rapend  ('FILE_O_interp', 2)

       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = src_P(k,istart+i-1,jstart+j-1)
       enddo
       enddo
       enddo

    else ! no interpolation

       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = src(kstart+k-1,istart+i-1,jstart+j-1)
       enddo
       enddo
       enddo

    endif

    if ( fill_halo ) then
       ! W halo
       do k = 1, ksize
       do j = 1, jsize
       do i = 1, IS-istart
          dst((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
       enddo
       enddo
       enddo
       ! E halo
       do k = 1, ksize
       do j = 1, jsize
       do i = IE-istart+2, ime-istart+1
          dst((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
       enddo
       enddo
       enddo
       ! S halo
       do k = 1, ksize
       do j = 1, JS-jstart
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
       enddo
       enddo
       enddo
       ! N halo
       do k = 1, ksize
       do j = JE-jstart+2, jme-jstart+1
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
       enddo
       enddo
       enddo
    endif

    return
  end subroutine FILE_HISTORY_CARTESC_truncate_3D

  !-----------------------------------------------------------------------------
  !> get information of axis coordinate to history file
  !  only register the axis and coordinate variables into internal buffers
  !  The actual write happens later when calling FILE_HISTORY_CARTESC_write
  subroutine FILE_HISTORY_CARTESC_set_axes
    use scale_file_history, only: &
       FILE_HISTORY_AGGREGATE, &
       FILE_HISTORY_Set_Axis, &
       FILE_HISTORY_Set_AssociatedCoordinate
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       D2R => CONST_D2R
    use scale_process, only: &
       PRC_myrank
    use scale_rm_process, only: &
       PRC_2Drank,     &
       PRC_NUM_X,      &
       PRC_NUM_Y,      &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y, &
       PRC_HAS_W,      &
       PRC_HAS_S
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CZ,    &
       ATMOS_GRID_CARTESC_CX,    &
       ATMOS_GRID_CARTESC_CY,    &
       ATMOS_GRID_CARTESC_FZ,    &
       ATMOS_GRID_CARTESC_FX,    &
       ATMOS_GRID_CARTESC_FY,    &
       ATMOS_GRID_CARTESC_CDZ,   &
       ATMOS_GRID_CARTESC_CDX,   &
       ATMOS_GRID_CARTESC_CDY,   &
       ATMOS_GRID_CARTESC_FDZ,   &
       ATMOS_GRID_CARTESC_FDX,   &
       ATMOS_GRID_CARTESC_FDY,   &
       ATMOS_GRID_CARTESC_CBFZ,  &
       ATMOS_GRID_CARTESC_CBFX,  &
       ATMOS_GRID_CARTESC_CBFY,  &
       ATMOS_GRID_CARTESC_FBFZ,  &
       ATMOS_GRID_CARTESC_FBFX,  &
       ATMOS_GRID_CARTESC_FBFY,  &
       ATMOS_GRID_CARTESC_CXG,   &
       ATMOS_GRID_CARTESC_CYG,   &
       ATMOS_GRID_CARTESC_FXG,   &
       ATMOS_GRID_CARTESC_FYG,   &
       ATMOS_GRID_CARTESC_CDXG,  &
       ATMOS_GRID_CARTESC_CDYG,  &
       ATMOS_GRID_CARTESC_FDXG,  &
       ATMOS_GRID_CARTESC_FDYG,  &
       ATMOS_GRID_CARTESC_CBFXG, &
       ATMOS_GRID_CARTESC_CBFYG, &
       ATMOS_GRID_CARTESC_FBFXG, &
       ATMOS_GRID_CARTESC_FBFYG
    use scale_ocean_grid_cartesC, only: &
       OCEAN_GRID_CARTESC_CZ, &
       OCEAN_GRID_CARTESC_FZ, &
       OCEAN_GRID_CARTESC_CDZ
    use scale_land_grid_cartesC, only: &
       LAND_GRID_CARTESC_CZ, &
       LAND_GRID_CARTESC_FZ, &
       LAND_GRID_CARTESC_CDZ
    use scale_urban_grid_cartesC, only: &
       URBAN_GRID_CARTESC_CZ, &
       URBAN_GRID_CARTESC_FZ, &
       URBAN_GRID_CARTESC_CDZ
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_CZ,    &
       ATMOS_GRID_CARTESC_REAL_FZ,    &
       ATMOS_GRID_CARTESC_REAL_LON,   &
       ATMOS_GRID_CARTESC_REAL_LONUY, &
       ATMOS_GRID_CARTESC_REAL_LONXV, &
       ATMOS_GRID_CARTESC_REAL_LONUV, &
       ATMOS_GRID_CARTESC_REAL_LAT,   &
       ATMOS_GRID_CARTESC_REAL_LATUY, &
       ATMOS_GRID_CARTESC_REAL_LATXV, &
       ATMOS_GRID_CARTESC_REAL_LATUV, &
       AREA   => ATMOS_GRID_CARTESC_REAL_AREA,   &
       AREAUY => ATMOS_GRID_CARTESC_REAL_AREAUY, &
       AREAXV => ATMOS_GRID_CARTESC_REAL_AREAXV, &
       AREAZUY_X => ATMOS_GRID_CARTESC_REAL_AREAZUY_X, &
       AREAZXV_Y => ATMOS_GRID_CARTESC_REAL_AREAZXV_Y, &
       AREAWUY_X => ATMOS_GRID_CARTESC_REAL_AREAWUY_X, &
       AREAWXV_Y => ATMOS_GRID_CARTESC_REAL_AREAWXV_Y, &
       AREAZXY_X => ATMOS_GRID_CARTESC_REAL_AREAZXY_X, &
       AREAZUV_Y => ATMOS_GRID_CARTESC_REAL_AREAZUV_Y, &
       AREAZUV_X => ATMOS_GRID_CARTESC_REAL_AREAZUV_X, &
       AREAZXY_Y => ATMOS_GRID_CARTESC_REAL_AREAZXY_Y, &
       VOL    => ATMOS_GRID_CARTESC_REAL_VOL,    &
       VOLWXY => ATMOS_GRID_CARTESC_REAL_VOLWXY, &
       VOLZUY => ATMOS_GRID_CARTESC_REAL_VOLZUY, &
       VOLZXV => ATMOS_GRID_CARTESC_REAL_VOLZXV
    use scale_ocean_grid_cartesC_real, only: &
       VOLO => OCEAN_GRID_CARTESC_REAL_VOL
    use scale_land_grid_cartesC_real, only: &
       VOLL => LAND_GRID_CARTESC_REAL_VOL
    use scale_urban_grid_cartesC_real, only: &
       VOLU => URBAN_GRID_CARTESC_REAL_VOL
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_landuse, only: &
       LANDUSE_frac_land
    implicit none

    real(RP)         :: AXIS (imh,jmh,0:KMAX)
    real(RP)         :: AXISO(im, jm, OKMAX)
    real(RP)         :: AXISL(im, jm, LKMAX)
    real(RP)         :: AXISU(im, jm, UKMAX)
    character(len=2) :: AXIS_name(3)

    integer :: rankidx(2)
    integer :: start(3,4) !> 1: FF, 2: HF, 3: FH, 4: HH (x,y)
    integer :: startX, startY, startZ
    integer :: startXH, startYH
    integer :: XAG, YAG
    integer :: XAGH, YAGH

    real(RP) :: z_bnds(2,KA), zh_bnds(2,0:KA)
    real(RP) :: oz_bnds(2,OKA), ozh_bnds(2,0:OKA)
    real(RP) :: lz_bnds(2,LKA), lzh_bnds(2,0:LKA)
    real(RP) :: uz_bnds(2,UKA), uzh_bnds(2,0:UKA)
    real(RP) :: x_bnds(2,IA), xh_bnds(2,0:IA)
    real(RP) :: y_bnds(2,JA), yh_bnds(2,0:JA)

    real(RP) :: FDXG(0:IAG), FDYG(0:JAG)
    real(RP) :: FDX(0:IA), FDY(0:JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    ! For parallel I/O, some variables are written by a subset of processes.
    ! 1. Only PRC_myrank 0 writes all z axes
    ! 2. Only south-most processes (rankidx(2) == 0) write x axes
    !         rankidx(1) == 0           writes west HALO
    !         rankidx(1) == PRC_NUM_X-1 writes east HALO
    !         others                    writes without HALO
    ! 3. Only west-most processes (rankidx(1) == 0) write y axes
    !         rankidx(1) == 0           writes south HALO
    !         rankidx(1) == PRC_NUM_Y-1 writes north HALO
    !         others                    writes without HALO

    if ( FILE_HISTORY_AGGREGATE ) then

       startZ = 1

       if ( FILE_HISTORY_CARTESC_BOUNDARY ) then
          startX  = ISGB ! global subarray starting index
          startY  = JSGB ! global subarray starting index
          startXH = startX
          startYH = startY
          XAG     = IAGB
          YAG     = JAGB
          XAGH    = XAG
          YAGH    = YAG
       else
          startX  = PRC_2Drank(PRC_myrank,1) * IMAX + 1
          startY  = PRC_2Drank(PRC_myrank,2) * JMAX + 1
          startXH = startX
          startYH = startY
          XAG     = IMAXG
          YAG     = JMAXG
          XAGH    = XAG
          YAGH    = YAG

          if ( .NOT. PRC_PERIODIC_X ) then
             XAGH = XAGH + 1
             if( PRC_HAS_W ) startXH = startXH + 1
          endif

          if ( .NOT. PRC_PERIODIC_Y ) then
             YAGH = YAGH + 1
             if( PRC_HAS_S ) startYH = startYH + 1
          endif
       endif

       ! for shared-file parallel I/O, only a part of rank writes variables
       if ( PRC_myrank > 0 ) then ! only rank 0 writes Z axes
          startZ = -1
       endif
       if ( rankidx(2) > 0 ) then ! only south-most processes write
          startX  = -1
          startXH = -1
       endif
       if ( rankidx(1) > 0 ) then ! only west-most processes write
          startY  = -1
          startYH = -1
       endif
    else
       startZ  = 1
       startX  = 1
       startY  = 1
       startXH = startX
       startYH = startY
       XAG     = im
       YAG     = jm
       XAGH    = imh
       YAGH    = jmh
    endif


    ! bounds
    do k = KS, KE
       z_bnds(1,k) = ATMOS_GRID_CARTESC_FZ(k-1)
       z_bnds(2,k) = ATMOS_GRID_CARTESC_FZ(k  )
    end do
    do k = KS-1, KE
       zh_bnds(1,k) = ATMOS_GRID_CARTESC_CZ(k  )
       zh_bnds(2,k) = ATMOS_GRID_CARTESC_CZ(k+1)
    end do

    do k = OKS, OKE
       oz_bnds(1,k) = OCEAN_GRID_CARTESC_FZ(k-1)
       oz_bnds(2,k) = OCEAN_GRID_CARTESC_FZ(k  )
    end do
    ozh_bnds(1,OKS-1) = OCEAN_GRID_CARTESC_FZ(OKS-1)
    do k = OKS-1, OKE-1
       ozh_bnds(2,k  ) = OCEAN_GRID_CARTESC_CZ(k+1)
       ozh_bnds(1,k+1) = OCEAN_GRID_CARTESC_CZ(k+1)
    end do
    ozh_bnds(2,OKE) = OCEAN_GRID_CARTESC_FZ(OKE)

    do k = LKS, LKE
       lz_bnds(1,k) = LAND_GRID_CARTESC_FZ(k-1)
       lz_bnds(2,k) = LAND_GRID_CARTESC_FZ(k  )
    end do
    lzh_bnds(1,LKS-1) = LAND_GRID_CARTESC_FZ(LKS-1)
    do k = LKS-1, LKE-1
       lzh_bnds(2,k  ) = LAND_GRID_CARTESC_CZ(k+1)
       lzh_bnds(1,k+1) = LAND_GRID_CARTESC_CZ(k+1)
    end do
    lzh_bnds(2,LKE) = LAND_GRID_CARTESC_FZ(LKE)

    do k = UKS, UKE
       uz_bnds(1,k) = URBAN_GRID_CARTESC_FZ(k-1)
       uz_bnds(2,k) = URBAN_GRID_CARTESC_FZ(k  )
    end do
    uzh_bnds(1,UKS-1) = URBAN_GRID_CARTESC_FZ(UKS-1)
    do k = UKS-1, UKE-1
       uzh_bnds(2,k  ) = URBAN_GRID_CARTESC_CZ(k+1)
       uzh_bnds(1,k+1) = URBAN_GRID_CARTESC_CZ(k+1)
    end do
    uzh_bnds(2,UKE) = URBAN_GRID_CARTESC_FZ(UKE)

    do i = ims, ime
       x_bnds(1,i) = ATMOS_GRID_CARTESC_FX(i-1)
       x_bnds(2,i) = ATMOS_GRID_CARTESC_FX(i  )
    end do
    if ( imsh == 0 ) then
       xh_bnds(1,0) = ATMOS_GRID_CARTESC_FX(0)
    else
       xh_bnds(1,imsh) = ATMOS_GRID_CARTESC_CX(imsh)
    end if
    do i = imsh, ime-1
       xh_bnds(2,i  ) = ATMOS_GRID_CARTESC_CX(i+1)
       xh_bnds(1,i+1) = ATMOS_GRID_CARTESC_CX(i+1)
    end do
    if ( ime == IA ) then
       xh_bnds(2,ime) = ATMOS_GRID_CARTESC_FX(IA)
    else
       xh_bnds(2,ime) = ATMOS_GRID_CARTESC_CX(ime+1)
    end if

    do j = jms, jme
       y_bnds(1,j) = ATMOS_GRID_CARTESC_FY(j-1)
       y_bnds(2,j) = ATMOS_GRID_CARTESC_FY(j  )
    end do
    if ( jmsh == 0 ) then
       yh_bnds(1,jmsh) = ATMOS_GRID_CARTESC_FY(jmsh)
    else
       yh_bnds(1,jmsh) = ATMOS_GRID_CARTESC_CY(jmsh)
    end if
    do j = jmsh, jme-1
       yh_bnds(2,j  ) = ATMOS_GRID_CARTESC_CY(j+1)
       yh_bnds(1,j+1) = ATMOS_GRID_CARTESC_CY(j+1)
    end do
    if ( jme == JA ) then
       yh_bnds(2,jme) = ATMOS_GRID_CARTESC_FY(jme)
    else
       yh_bnds(2,jme) = ATMOS_GRID_CARTESC_CY(jme+1)
    end if


    FDXG(1:IAG-1) = ATMOS_GRID_CARTESC_FDXG(:)
    FDXG(0  ) = UNDEF
    FDXG(IAG) = UNDEF
    FDYG(1:JAG-1) = ATMOS_GRID_CARTESC_FDYG(:)
    FDYG(0  ) = UNDEF
    FDYG(JAG) = UNDEF

    FDX(1:IA-1) = ATMOS_GRID_CARTESC_FDX(:)
    FDX(0 ) = FDXG(IS_inG-IHALO-1)
    FDX(IA) = FDXG(IE_inG+IHALO  )
    FDY(1:JA-1) = ATMOS_GRID_CARTESC_FDY(:)
    FDY(0 ) = FDYG(JS_inG-JHALO-1)
    FDY(JA) = FDYG(JE_inG+JHALO  )


    ! for the shared-file I/O method, the axes are global (gsize)
    ! for one-file-per-process I/O method, the axes size is equal to the local buffer size

    call FILE_HISTORY_Set_Axis( 'z',   'Z',               'm', 'z',   ATMOS_GRID_CARTESC_CZ (KS  :KE), &
                                bounds=z_bnds (:,KS  :KE), gsize=KMAX   , start=startZ                 )
    call FILE_HISTORY_Set_Axis( 'zh',  'Z (half level)',  'm', 'zh',  ATMOS_GRID_CARTESC_FZ (KS-1:KE), &
                                bounds=zh_bnds(:,KS-1:KE), gsize=KMAX+1 , start=startZ                 )

    if ( FILE_HISTORY_CARTESC_PRES_nlayer > 0 ) then
       call FILE_HISTORY_Set_Axis( 'pressure', 'Pressure', 'hPa', 'pressure', FILE_HISTORY_CARTESC_PRES_val(:)/100.0_RP, &
                                    gsize=FILE_HISTORY_CARTESC_PRES_nlayer, start=startZ, down=.true.                    )
    endif

    call FILE_HISTORY_Set_Axis( 'oz',  'OZ',              'm', 'oz',  OCEAN_GRID_CARTESC_CZ(OKS  :OKE), &
                                bounds=oz_bnds (:,OKS  :OKE), gsize=OKMAX  , start=startZ, down=.true.  )
    call FILE_HISTORY_Set_Axis( 'ozh', 'OZ (half level)', 'm', 'ozh', OCEAN_GRID_CARTESC_FZ(OKS-1:OKE), &
                                bounds=ozh_bnds(:,OKS-1:OKE), gsize=OKMAX+1, start=startZ, down=.true.  )

    call FILE_HISTORY_Set_Axis( 'lz',  'LZ',              'm', 'lz',  LAND_GRID_CARTESC_CZ(LKS  :LKE),  &
                                bounds=lz_bnds (:,LKS  :LKE), gsize=LKMAX  , start=startZ, down=.true.  )
    call FILE_HISTORY_Set_Axis( 'lzh', 'LZ (half level)', 'm', 'lzh', LAND_GRID_CARTESC_FZ(LKS-1:LKE),  &
                                bounds=lzh_bnds(:,LKS-1:LKE), gsize=LKMAX+1, start=startZ, down=.true.  )

    call FILE_HISTORY_Set_Axis( 'uz',  'UZ',              'm', 'uz',  URBAN_GRID_CARTESC_CZ(UKS  :UKE), &
                                bounds=uz_bnds (:,UKS  :UKE), gsize=UKMAX  , start=startZ, down=.true.  )
    call FILE_HISTORY_Set_Axis( 'uzh', 'UZ (half level)', 'm', 'uzh', URBAN_GRID_CARTESC_FZ(UKS-1:UKE), &
                                bounds=uzh_bnds(:,UKS-1:UKE), gsize=UKMAX+1, start=startZ, down=.true.  )

    call FILE_HISTORY_Set_Axis( 'x',   'X',               'm', 'x',   ATMOS_GRID_CARTESC_CX (ims :ime), &
                                bounds=x_bnds (:,ims :ime), gsize=XAG , start=startX                    )
    call FILE_HISTORY_Set_Axis( 'xh',  'X (half level)',  'm', 'xh',  ATMOS_GRID_CARTESC_FX (imsh:ime), &
                                bounds=xh_bnds(:,imsh:ime), gsize=XAGH, start=startXH                   )

    call FILE_HISTORY_Set_Axis( 'y',   'Y',               'm', 'y',   ATMOS_GRID_CARTESC_CY (jms :jme), &
                                bounds=y_bnds (:,jms :jme), gsize=YAG , start=startY                    )
    call FILE_HISTORY_Set_Axis( 'yh',  'Y (half level)',  'm', 'yh',  ATMOS_GRID_CARTESC_FY (jmsh:jme), &
                                bounds=yh_bnds(:,jmsh:jme), gsize=YAGH, start=startYH                   )

    ! axes below always include halos when written to file regardless of PRC_PERIODIC_X/PRC_PERIODIC_Y
    call FILE_HISTORY_Set_Axis( 'CZ',   'Atmos Grid Center Position Z', 'm', 'CZ',  ATMOS_GRID_CARTESC_CZ,   gsize=KA,      start=startZ )
    call FILE_HISTORY_Set_Axis( 'FZ',   'Atmos Grid Face Position Z',   'm', 'FZ',  ATMOS_GRID_CARTESC_FZ,   gsize=KA+1,    start=startZ )
    call FILE_HISTORY_Set_Axis( 'CDZ',  'Grid Cell length Z',           'm', 'CZ',  ATMOS_GRID_CARTESC_CDZ,  gsize=KA,      start=startZ )
    call FILE_HISTORY_Set_Axis( 'FDZ',  'Grid distance Z',              'm', 'FDZ', ATMOS_GRID_CARTESC_FDZ,  gsize=KA-1,    start=startZ )
    call FILE_HISTORY_Set_Axis( 'CBFZ', 'Boundary factor Center Z',     '1', 'CZ',  ATMOS_GRID_CARTESC_CBFZ, gsize=KA,      start=startZ )
    call FILE_HISTORY_Set_Axis( 'FBFZ', 'Boundary factor Face Z',       '1', 'FZ',  ATMOS_GRID_CARTESC_FBFZ, gsize=KA+1,    start=startZ )

    call FILE_HISTORY_Set_Axis( 'OCZ',  'Ocean Grid Center Position Z', 'm', 'OCZ', OCEAN_GRID_CARTESC_CZ,  gsize=OKMAX,   start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'OFZ',  'Ocean Grid Face Position Z',   'm', 'OFZ', OCEAN_GRID_CARTESC_FZ,  gsize=OKMAX+1, start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'OCDZ', 'Ocean Grid Cell length Z',     'm', 'OCZ', OCEAN_GRID_CARTESC_CDZ, gsize=OKMAX,   start=startZ              )

    call FILE_HISTORY_Set_Axis( 'LCZ',  'Land Grid Center Position Z',  'm', 'LCZ', LAND_GRID_CARTESC_CZ,  gsize=LKMAX,   start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'LFZ',  'Land Grid Face Position Z',    'm', 'LFZ', LAND_GRID_CARTESC_FZ,  gsize=LKMAX+1, start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'LCDZ', 'Land Grid Cell length Z',      'm', 'LCZ', LAND_GRID_CARTESC_CDZ, gsize=LKMAX,   start=startZ              )

    call FILE_HISTORY_Set_Axis( 'UCZ',  'Urban Grid Center Position Z', 'm', 'UCZ', URBAN_GRID_CARTESC_CZ,  gsize=UKMAX,   start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'UFZ',  'Urban Grid Face Position Z',   'm', 'UFZ', URBAN_GRID_CARTESC_FZ,  gsize=UKMAX+1, start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'UCDZ', 'Urban Grid Cell length Z',     'm', 'UCZ', URBAN_GRID_CARTESC_CDZ, gsize=UKMAX,   start=startZ              )

    if ( FILE_HISTORY_AGGREGATE ) then
       call FILE_HISTORY_Set_Axis( 'CX',   'Atmos Grid Center Position X', 'm', 'CX', ATMOS_GRID_CARTESC_CXG,   gsize=IAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'CY',   'Atmos Grid Center Position Y', 'm', 'CY', ATMOS_GRID_CARTESC_CYG,   gsize=JAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'FX',   'Atmos Grid Face Position X',   'm', 'FX', ATMOS_GRID_CARTESC_FXG,   gsize=IAG+1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'FY',   'Atmos Grid Face Position Y',   'm', 'FY', ATMOS_GRID_CARTESC_FYG,   gsize=JAG+1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'CDX',  'Grid Cell length X',           'm', 'CX', ATMOS_GRID_CARTESC_CDXG,  gsize=IAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'CDY',  'Grid Cell length Y',           'm', 'CY', ATMOS_GRID_CARTESC_CDYG,  gsize=JAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'FDX',  'Grid distance X',              'm', 'FX',                    FDXG,  gsize=IAG+1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'FDY',  'Grid distance Y',              'm', 'FY',                    FDYG,  gsize=JAG+1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'CBFX', 'Boundary factor Center X',     '1', 'CX', ATMOS_GRID_CARTESC_CBFXG, gsize=IAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'CBFY', 'Boundary factor Center Y',     '1', 'CY', ATMOS_GRID_CARTESC_CBFYG, gsize=JAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'FBFX', 'Boundary factor Face X',       '1', 'FX', ATMOS_GRID_CARTESC_FBFXG, gsize=IAG+1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'FBFY', 'Boundary factor Face Y',       '1', 'FY', ATMOS_GRID_CARTESC_FBFYG, gsize=JAG+1, start=startZ )
    else
       call FILE_HISTORY_Set_Axis( 'CX',   'Atmos Grid Center Position X', 'm', 'CX', ATMOS_GRID_CARTESC_CX   )
       call FILE_HISTORY_Set_Axis( 'CY',   'Atmos Grid Center Position Y', 'm', 'CY', ATMOS_GRID_CARTESC_CY   )
       call FILE_HISTORY_Set_Axis( 'FX',   'Atmos Grid Face Position X',   'm', 'FX', ATMOS_GRID_CARTESC_FX   )
       call FILE_HISTORY_Set_Axis( 'FY',   'Atmos Grid Face Position Y',   'm', 'FY', ATMOS_GRID_CARTESC_FY   )
       call FILE_HISTORY_Set_Axis( 'CDX',  'Grid Cell length X',           'm', 'CX', ATMOS_GRID_CARTESC_CDX  )
       call FILE_HISTORY_Set_Axis( 'CDY',  'Grid Cell length Y',           'm', 'CY', ATMOS_GRID_CARTESC_CDY  )
       call FILE_HISTORY_Set_Axis( 'FDX',  'Grid distance X',              'm', 'FX',                    FDX  )
       call FILE_HISTORY_Set_Axis( 'FDY',  'Grid distance Y',              'm', 'FY',                    FDY  )
       call FILE_HISTORY_Set_Axis( 'CBFX', 'Boundary factor Center X',     '1', 'CX', ATMOS_GRID_CARTESC_CBFX )
       call FILE_HISTORY_Set_Axis( 'CBFY', 'Boundary factor Center Y',     '1', 'CY', ATMOS_GRID_CARTESC_CBFY )
       call FILE_HISTORY_Set_Axis( 'FBFX', 'Boundary factor Face X',       '1', 'FX', ATMOS_GRID_CARTESC_FBFX )
       call FILE_HISTORY_Set_Axis( 'FBFY', 'Boundary factor Face Y',       '1', 'FY', ATMOS_GRID_CARTESC_FBFY )
    endif

    call FILE_HISTORY_Set_Axis('CXG',   'Grid Center Position X (global)',   'm', 'CXG',  ATMOS_GRID_CARTESC_CXG,   gsize=IAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('CYG',   'Grid Center Position Y (global)',   'm', 'CYG',  ATMOS_GRID_CARTESC_CYG,   gsize=JAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('FXG',   'Grid Face Position X (global)',     'm', 'FXG',  ATMOS_GRID_CARTESC_FXG,   gsize=IAG+1, start=startZ )
    call FILE_HISTORY_Set_Axis('FYG',   'Grid Face Position Y (global)',     'm', 'FYG',  ATMOS_GRID_CARTESC_FYG,   gsize=JAG+1, start=startZ )
    call FILE_HISTORY_Set_Axis('CDXG',  'Grid Cell length X (global)',       'm', 'CXG',  ATMOS_GRID_CARTESC_CDXG,  gsize=IAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('CDYG',  'Grid Cell length Y (global)',       'm', 'CYG',  ATMOS_GRID_CARTESC_CDYG,  gsize=JAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('FDXG',  'Grid distance X (global)',          'm', 'FDXG',                    FDXG,  gsize=IAG+1, start=startZ )
    call FILE_HISTORY_Set_Axis('FDYG',  'Grid distance Y (global)',          'm', 'FDYG',                    FDYG,  gsize=JAG+1, start=startZ )
    call FILE_HISTORY_Set_Axis('CBFXG', 'Boundary factor Center X (global)', '1', 'CXG',  ATMOS_GRID_CARTESC_CBFXG, gsize=IAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG',  ATMOS_GRID_CARTESC_CBFYG, gsize=JAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('FBFXG', 'Boundary factor Face X (global)',   '1', 'FXG',  ATMOS_GRID_CARTESC_FBFXG, gsize=IAG+1, start=startZ )
    call FILE_HISTORY_Set_Axis('FBFYG', 'Boundary factor Face Y (global)',   '1', 'FYG',  ATMOS_GRID_CARTESC_FBFYG, gsize=JAG+1, start=startZ )



    ! associate coordinates
    if ( FILE_HISTORY_AGGREGATE ) then
       if ( FILE_HISTORY_CARTESC_BOUNDARY ) then
          start(1,:) = ISGB   ! global subarray starting index
          start(2,:) = JSGB   ! global subarray starting index
       else
          start(1,:) = PRC_2Drank(PRC_myrank,1) * IMAX + 1 ! no IHALO
          start(2,:) = PRC_2Drank(PRC_myrank,2) * JMAX + 1 ! no JHALO
          if ( (.NOT. PRC_PERIODIC_X) .AND. PRC_HAS_W ) then
             start(1,2) = start(1,2) + 1
             start(1,4) = start(1,4) + 1
          endif
          if ( (.NOT. PRC_PERIODIC_Y) .AND. PRC_HAS_S ) then
             start(2,3) = start(2,3) + 1
             start(2,4) = start(2,4) + 1
          endif
       endif
       start(3,:) = 1
    else
       start(:,:) = 1
    endif

    do k = 1, KMAX
       AXIS(1:im,1:jm,k) = ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,ims:ime,jms:jme)
    enddo
    AXIS_name(1:3) = (/'x ', 'y ', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height', 'height above ground level',                        &
                                                'm', AXIS_name(1:3), AXIS(1:im,1:jm,1:KMAX), start=start(:,1) )

    do k = 0, KMAX
       AXIS(1:im,1:jm,k) = ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,ims:ime,jms:jme)
    enddo
    AXIS_name(1:3) = (/'x ', 'y ', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_xyw', 'height above ground level (half level xyw)',    &
                                                'm' , AXIS_name(1:3), AXIS(1:im,1:jm,0:KMAX), start=start(:,1) )

    do k = 1, KMAX
    do j = 1, jm
    do i = 1, min(imh,IA-imsh)
       AXIS(i,j,k) = ( ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+i-1,jms+j-1) + ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+i,jms+j-1) ) * 0.5_RP
    enddo
    enddo
    enddo
    if ( imh == IA-imsh+1 ) then
       do k = 1, KMAX
       do j = 1, jm
          AXIS(imh,j,k) = ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+imh-1,jms+j-1)
       enddo
       enddo
    endif
    AXIS_name(1:3) = (/'xh', 'y ', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_uyz', 'height above ground level (half level uyz)',    &
                                                'm', AXIS_name(1:3), AXIS(1:imh,1:jm,1:KMAX), start=start(:,2) )

    do k = 1, KMAX
    do j = 1, min(jmh,JA-jmsh)
    do i = 1, im
       AXIS(i,j,k) = ( ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,ims+i-1,jmsh+j-1) + ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,ims+i-1,jmsh+j) ) * 0.5_RP
    enddo
    enddo
    enddo
    if ( jmh == JA-jmsh+1 ) then
       do k = 1, KMAX
       do i = 1, im
          AXIS(i,jmh,k) = ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,ims+i-1,jmsh+jmh-1)
       enddo
       enddo
    endif
    AXIS_name(1:3) = (/'x ', 'yh', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_xvz', 'height above ground level (half level xvz)',    &
                                                'm', AXIS_name(1:3), AXIS(1:im,1:jmh,1:KMAX), start=start(:,3) )

    do k = 1, KMAX
    do j = 1, min(jmh,JA-jmsh)
    do i = 1, min(imh,IA-imsh)
       AXIS(i,j,k) = ( ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+i-1,jmsh+j-1) + ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+i  ,jmsh+j-1) &
                     + ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+i-1,jmsh+j  ) + ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+i  ,jmsh+j  ) ) * 0.25_RP
    enddo
    enddo
    enddo
    if ( jmh == JA-jmsh+1 ) then
       do k = 1, KMAX
       do i = 1, min(imh,IA-imsh)
          AXIS(i,jmh,k) = ( ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+i-1,jmsh+jmh-1) + ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+i,jmsh+jmh-1) ) * 0.5_RP
       enddo
       enddo
    endif
    if ( imh == IA-imsh+1 ) then
       do k = 1, KMAX
       do j = 1, min(jmh,JA-jmsh)
          AXIS(imh,j,k) = ( ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+imh-1,jmsh+j-1) + ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+imh-1,jmsh+j) ) * 0.5_RP
       enddo
       enddo
    endif
    if ( imh == IA-imsh+1 .AND. jmh == JA-jmsh+1 ) then
       do k = 1, KMAX
          AXIS(imh,jmh,k) = ATMOS_GRID_CARTESC_REAL_CZ(k+KS-1,imsh+imh-1,jmsh+jmh-1)
       enddo
    endif
    AXIS_name(1:3) = (/'xh', 'yh', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_uvz', 'height above ground level (half level uvz)',     &
                                                'm', AXIS_name(1:3), AXIS(1:imh,1:jmh,1:KMAX), start=start(:,4) )

    do k = 0, KMAX
    do j = 1, jm
    do i = 1, min(imh,IA-imsh)
       AXIS(i,j,k) = ( ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+i-1,jms+j-1) + ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+i,jms+j-1) ) * 0.5_RP
    enddo
    enddo
    enddo
    if ( imh == IA-imsh+1 ) then
       do k = 0, KMAX
       do j = 1, jm
          AXIS(imh,j,k) = ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+imh-1,jms+j-1)
       enddo
       enddo
    endif
    AXIS_name(1:3) = (/'xh', 'y ', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_uyw', 'height above ground level (half level uyw)',    &
                                                'm', AXIS_name(1:3), AXIS(1:imh,1:jm,0:KMAX), start=start(:,2) )

    do k = 0, KMAX
    do j = 1, min(jmh,JA-jmsh)
    do i = 1, im
       AXIS(i,j,k) = ( ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,ims+i-1,jmsh+j-1) + ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,ims+i-1,jmsh+j) ) * 0.5_RP
    enddo
    enddo
    enddo
    if ( jmh == JA-jmsh+1 ) then
       do k = 0, KMAX
       do i = 1, im
          AXIS(i,jmh,k) = ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,ims+i-1,jmsh+jmh-1)
       enddo
       enddo
    endif
    AXIS_name(1:3) = (/'x ', 'yh', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_xvw', 'height above ground level (half level xvw)',    &
                                                'm', AXIS_name(1:3), AXIS(1:im,1:jmh,0:KMAX), start=start(:,3) )

    do k = 0, KMAX
    do j = 1, min(jmh,JA-jmsh)
    do i = 1, min(imh,IA-imsh)
       AXIS(i,j,k) = ( ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+i-1,jmsh+j-1) + ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+i  ,jmsh+j-1) &
                     + ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+i-1,jmsh+j  ) + ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+i  ,jmsh+j  ) ) * 0.25_RP
    enddo
    enddo
    enddo
    if ( jmh == JA-jmsh+1 ) then
       do k = 0, KMAX
       do i = 1, min(imh,IA-imsh)
          AXIS(i,jmh,k) = ( ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+i-1,jmsh+jmh-1) + ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+i,jmsh+jmh-1) ) * 0.5_RP
       enddo
       enddo
    endif
    if ( imh == IA-imsh+1 ) then
       do k = 0, KMAX
       do j = 1, min(jmh,JA-jmsh)
          AXIS(imh,j,k) = ( ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+imh-1,jmsh+j-1) + ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+imh-1,jmsh+j) ) * 0.5_RP
       enddo
       enddo
    endif
    if ( imh == IA-imsh+1 .AND. jm == JA-jmsh+1 ) then
       do k = 0, KMAX
          AXIS(imh,jmh,k) = ATMOS_GRID_CARTESC_REAL_FZ(k+KS-1,imsh+imh-1,jmsh+jmh-1)
       enddo
    endif
    AXIS_name(1:3) = (/'xh', 'yh', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_uvw', 'height above ground level (half level uvw)', 'm', &
                                                AXIS_name(1:3), AXIS(1:imh,1:jmh,0:KMAX), start=start(:,4)       )

    AXIS(1:im,1:jm,1) = ATMOS_GRID_CARTESC_REAL_LON (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lon', 'longitude', 'degrees_east',                 &
                                                AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1) )

    AXIS(1:imh,1:jm,1) = ATMOS_GRID_CARTESC_REAL_LONUY(imsh:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lon_uy', 'longitude (half level uy)', 'degrees_east', &
                                                AXIS_name(1:2), AXIS(1:imh,1:jm,1), start=start(:,2)   )

    AXIS(1:im,1:jmh,1) = ATMOS_GRID_CARTESC_REAL_LONXV(ims:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'yh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lon_xv', 'longitude (half level xv)', 'degrees_east', &
                                                AXIS_name(1:2), AXIS(1:im,1:jmh,1), start=start(:,3)   )

    AXIS(1:imh,1:jmh,1) = ATMOS_GRID_CARTESC_REAL_LONUV(imsh:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'yh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lon_uv', 'longitude (half level uv)', 'degrees_east', &
                                                AXIS_name(1:2), AXIS(1:imh,1:jmh,1), start=start(:,4)  )

    AXIS(1:im,1:jm,1) = ATMOS_GRID_CARTESC_REAL_LAT (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lat', 'latitude',                     'degrees_north', &
                                                AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1)     )

    AXIS(1:imh,1:jm,1) = ATMOS_GRID_CARTESC_REAL_LATUY(imsh:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lat_uy', 'latitude (half level uy)',  'degrees_north', &
                                                AXIS_name(1:2), AXIS(1:imh,1:jm,1), start=start(:,2)    )

    AXIS(1:im,1:jmh,1) = ATMOS_GRID_CARTESC_REAL_LATXV(ims:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'yh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lat_xv', 'latitude (half level xv)',  'degrees_north', &
                                                AXIS_name(1:2), AXIS(1:im,1:jmh,1), start=start(:,3)    )

    AXIS(1:imh,1:jmh,1) = ATMOS_GRID_CARTESC_REAL_LATUV(imsh:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'yh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lat_uv', 'latitude (half level uv)',  'degrees_north', &
                                                AXIS_name(1:2), AXIS(1:imh,1:jmh,1), start=start(:,4)   )

    AXIS_name(1:2) = (/'x ', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'topo', 'topography', 'm', AXIS_name(1:2),   &
                                                TOPO_Zsfc(ims:ime,jms:jme), start=start(:,1) )

    AXIS_name(1:2) = (/'x ', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lsmask', 'fraction for land-sea mask', '1', AXIS_name(1:2), &
                                                LANDUSE_frac_land(ims:ime,jms:jme), start=start(:,1)         )

    AXIS_name(1:2) = (/'x ','y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area',    'area of grid cell',                  'm2', AXIS_name(1:2), &
                                                AREA(ims:ime,jms:jme),   start=start(:,1)                                   )
    AXIS_name(1:2) = (/'xh','y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_uy', 'area of grid cell (half level uy)',  'm2', AXIS_name(1:2), &
                                                AREAUY(ims:ime,jms:jme), start=start(:,2)                                   )
    AXIS_name(1:2) = (/'x ','yh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_xv', 'area of grid cell (half level xv)',  'm2', AXIS_name(1:2), &
                                                AREAXV(ims:ime,jms:jme), start=start(:,3)                                   )

    do k = 1, KMAX
    do j = 1, jm
    do i = 1, imh
       AXIS(i,j,k) = AREAZUY_X(KS+k-1,imsh+i-1,jms+j-1)
    end do
    end do
    end do
    AXIS_name = (/'xh', 'y ', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_uyz_x', 'area of grid cell face (half level uyz, normal x)', 'm2', &
                                                AXIS_name(1:3), AXIS(1:imh,1:jm,1:KMAX), start=start(:,2)                     )
    do k = 1, KMAX
    do j = 1, jmh
    do i = 1, im
       AXIS(i,j,k) = AREAZXV_Y(KS+k-1,ims+i-1,jmsh+j-1)
    end do
    end do
    end do
    AXIS_name = (/'x ', 'yh', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_xvz_y', 'area of grid cell face (half level xvz, normal y)', 'm2', &
                                                AXIS_name(1:3), AXIS(1:im,1:jmh,1:KMAX), start=start(:,3)                     )
    do k = 0, KMAX
    do j = 1, jmh
    do i = 1, imh
       AXIS(i,j,k) = AREAWUY_X(KS+k-1,imsh+i-1,jmsh+j-1)
    end do
    end do
    end do
    AXIS_name = (/'xh', 'y ', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_uyw_x', 'area of grid cell face (half level uyw, normal x)', 'm2', &
                                                AXIS_name(1:3), AXIS(1:imh,1:jmh,0:KMAX), start=start(:,2)                    )
    do k = 0, KMAX
    do j = 1, jmh
    do i = 1, im
       AXIS(i,j,k) = AREAWXV_Y(KS+k-1,ims+i-1,jmsh+j-1)
    end do
    end do
    end do
    AXIS_name = (/'x ', 'yh', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_xvw_y', 'area of grid cell face (half level xvw, normal y)', 'm2', &
                                                AXIS_name(1:3), AXIS(1:im,1:jmh,0:KMAX), start=start(:,3)                     )
    do k = 1, KMAX
    do j = 1, jm
    do i = 1, im
       AXIS(i,j,k) = AREAZXY_X(KS+k-1,ims+i-1,jms+j-1)
    end do
    end do
    end do
    AXIS_name = (/'x ', 'y ', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_xyz_x', 'area of grid cell face (half level xyz, normal x)', 'm2', &
                                                AXIS_name(1:3), AXIS(1:im,1:jm,1:KMAX), start=start(:,1)                      )
    do k = 1, KMAX
    do j = 1, jmh
    do i = 1, imh
       AXIS(i,j,k) = AREAZUV_Y(KS+k-1,imsh+i-1,jmsh+j-1)
    end do
    end do
    end do
    AXIS_name = (/'xh', 'yh', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_uvz_y', 'area of grid cell face (half level uvz, normal y)', 'm2', &
                                                AXIS_name(1:3), AXIS(1:imh,1:jmh,1:KMAX), start=start(:,4)                    )
    do k = 1, KMAX
    do j = 1, jmh
    do i = 1, imh
       AXIS(i,j,k) = AREAZUV_X(KS+k-1,imsh+i-1,jmsh+j-1)
    end do
    end do
    end do
    AXIS_name = (/'xh', 'yh', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_uvz_x', 'area of grid cell face (half level uvz, normal x)', 'm2', &
                                                AXIS_name(1:3), AXIS(1:imh,1:jmh,1:KMAX), start=start(:,4)                    )
    do k = 1, KMAX
    do j = 1, jm
    do i = 1, im
       AXIS(i,j,k) = AREAZXY_Y(KS+k-1,ims+i-1,jms+j-1)
    end do
    end do
    end do
    AXIS_name = (/'x ', 'y ', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_area_xyz_y', 'area of grid cell face (half level xyz, normal y)', 'm2', &
                                                AXIS_name(1:3), AXIS(1:im,1:jm,1:KMAX), start=start(:,1)                      )

    do k = 1, KMAX
    do j = 1, jm
    do i = 1, im
       AXIS(i,j,k) = VOL(KS+k-1,ims+i-1,jms+j-1)
    end do
    end do
    end do
    AXIS_name = (/ 'x ', 'y ', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_volume',     'volume of grid cell',                  'm3', &
                                                AXIS_name(1:3), AXIS(1:im,1:jm,1:KMAX), start=start(:,1)         )
    do k = 0, KMAX
    do j = 1, jm
    do i = 1, im
       AXIS(i,j,k) = VOLWXY(KS+k-1,ims+i-1,jms+j-1)
    end do
    end do
    end do
    AXIS_name = (/'x ', 'y ', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_volume_xyw', 'volume of grid cell (half level xyw)', 'm3', &
                                                AXIS_name(1:3), AXIS(1:im,1:jm,0:KMAX), start=start(:,1)         )
    do k = 1, KMAX
    do j = 1, jm
    do i = 1, imh
       AXIS(i,j,k) = VOLZUY(KS+k-1,imsh+i-1,jms+j-1)
    end do
    end do
    end do
    AXIS_name = (/'xh', 'y ', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_volume_uyz', 'volume of grid cell (half level uyz)',  'm3', &
                                                AXIS_name(1:3), AXIS(1:imh,1:jm,1:KMAX), start=start(:,2)         )
    do k = 1, KMAX
    do j = 1, jmh
    do i = 1, im
       AXIS(i,j,k) = VOLZXV(KS+k-1,ims+i-1,jmsh+j-1)
    end do
    end do
    end do
    AXIS_name = (/'x ', 'yh', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_volume_xvz', 'volume of grid cell (half level xvz)',  'm3', &
                                                AXIS_name(1:3), AXIS(1:im,1:jmh,1:KMAX), start=start(:,3)         )

    do k = 1, OKMAX
    do j = 1, jm
    do i = 1, im
       AXISO(i,j,k) = VOLO(OKS+k-1,ims+i-1,jms+j-1)
    end do
    end do
    end do
    AXIS_name = (/'x ', 'y ', 'oz'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_volume_xyo', 'volume of grid cell', 'm3', &
                                                AXIS_name(1:3), AXISO(:,:,:), start=start(:,1)  )
    do k = 1, LKMAX
    do j = 1, jm
    do i = 1, im
       AXISL(i,j,k) = VOLL(LKS+k-1,ims+i-1,jms+j-1)
    end do
    end do
    end do
    AXIS_name = (/'x ', 'y ', 'lz'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_volume_xyl', 'volume of grid cell', 'm3', &
                                                AXIS_name(1:3), AXISL(:,:,:), start=start(:,1)  )
    do k = 1, UKMAX
    do j = 1, jm
    do i = 1, im
       AXISU(i,j,k) = VOLU(UKS+k-1,ims+i-1,jms+j-1)
    end do
    end do
    end do
    AXIS_name = (/'x ', 'y ', 'uz'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'cell_volume_xyu', 'volume of grid cell', 'm3', &
                                                AXIS_name(1:3), AXISU(:,:,:), start=start(:,1)  )

    return
  end subroutine FILE_HISTORY_CARTESC_set_axes

  !-----------------------------------------------------------------------------
  subroutine FILE_HISTORY_CARTESC_set_axes_attributes
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_NAME
    use scale_calendar, only: &
       CALENDAR_get_name
    use scale_file_history, only: &
       FILE_HISTORY_AGGREGATE, &
       FILE_HISTORY_Set_Attribute
    use scale_process, only: &
       PRC_myrank
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_rm_process, only: &
       PRC_2Drank,     &
       PRC_NUM_X,      &
       PRC_NUM_Y,      &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y, &
       PRC_HAS_W,      &
       PRC_HAS_E,      &
       PRC_HAS_S,      &
       PRC_HAS_N
    use scale_file, only: &
       FILE_get_CFtunits
    use scale_file_cartesC, only: &
       axisattinfo,        &
       mappinginfo
    use scale_mapprojection, only: &
       MAPPROJECTION_get_attributes
    implicit none

    character(len=34) :: tunits
    character(len=H_SHORT) :: calendar

    type(axisattinfo) :: ainfo(4) ! x, xh, y, yh
    type(mappinginfo) :: minfo
    !---------------------------------------------------------------------------

    call FILE_HISTORY_Set_Attribute( "global", "grid_name", ATMOS_GRID_CARTESC_NAME ) ! [IN]

    if ( FILE_HISTORY_AGGREGATE ) then
       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_rank_x", (/0/) ) ! [IN]
       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_rank_y", (/0/) ) ! [IN]

       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_num_x", (/1/) ) ! [IN]
       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_num_y", (/1/) ) ! [IN]
    else
       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_rank_x", (/PRC_2Drank(PRC_myrank,1)/) ) ! [IN]
       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_rank_y", (/PRC_2Drank(PRC_myrank,2)/) ) ! [IN]

       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_num_x", (/PRC_NUM_X/) ) ! [IN]
       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_num_y", (/PRC_NUM_Y/) ) ! [IN]
    endif

    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_periodic_z", .false.        ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_periodic_x", PRC_PERIODIC_X ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_periodic_y", PRC_PERIODIC_Y ) ! [IN]

    call FILE_HISTORY_Set_Attribute( "global", "scale_atmos_grid_cartesC_index_imaxg", (/IMAXG/) ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_atmos_grid_cartesC_index_jmaxg", (/JMAXG/) ) ! [IN]

                     call FILE_HISTORY_Set_Attribute( "global", "scale_atmos_grid_cartesC_index_kmax", (/KMAX/)  ) ! [IN]
    if ( OKMAX > 0 ) call FILE_HISTORY_Set_Attribute( "global", "scale_ocean_grid_cartesC_index_kmax", (/OKMAX/)  ) ! [IN]
    if ( LKMAX > 0 ) call FILE_HISTORY_Set_Attribute( "global", "scale_land_grid_cartesC_index_kmax",  (/LKMAX/)  ) ! [IN]
    if ( UKMAX > 0 ) call FILE_HISTORY_Set_Attribute( "global", "scale_urban_grid_cartesC_index_kmax", (/UKMAX/)  ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_atmos_grid_cartesC_index_khalo", (/KHALO/) ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_atmos_grid_cartesC_index_ihalo", (/IHALO/) ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_atmos_grid_cartesC_index_jhalo", (/JHALO/) ) ! [IN]

    call FILE_get_CFtunits( FILE_HISTORY_CARTESCORY_STARTDATE(:), tunits )
    call FILE_HISTORY_Set_Attribute( "global", "time_units", tunits )
    call FILE_HISTORY_Set_Attribute( "global", "time_start", (/FILE_HISTORY_CARTESCORY_STARTMS/) )

    call CALENDAR_get_name( calendar )
    if ( calendar /= "" ) call FILE_HISTORY_Set_Attribute( "global", "calendar", calendar )

    if ( PRC_PERIODIC_X ) then
       ainfo(1)%periodic = .true.
       ainfo(2)%periodic = .true.
    else
       ainfo(1)%periodic = .false.
       ainfo(2)%periodic = .false.
    endif

    if ( PRC_PERIODIC_Y ) then
       ainfo(3)%periodic = .true.
       ainfo(4)%periodic = .true.
    else
       ainfo(3)%periodic = .false.
       ainfo(4)%periodic = .false.
    endif

    ! for x
    if ( PRC_PERIODIC_X .OR. .NOT. FILE_HISTORY_CARTESC_BOUNDARY ) then
       ainfo(1)%size_global (1) = IMAX * PRC_NUM_X
       ainfo(1)%start_global(1) = IS_inG - IHALO
       ainfo(1)%halo_global (1) = 0     ! west side
       ainfo(1)%halo_global (2) = 0     ! east side
       ainfo(1)%halo_local  (1) = 0     ! west side
       ainfo(1)%halo_local  (2) = 0     ! east side
    else
       ainfo(1)%size_global (1) = IAG
       ainfo(1)%start_global(1) = ISGA
       ainfo(1)%halo_global (1) = IHALO ! west side
       ainfo(1)%halo_global (2) = IHALO ! east side
       ainfo(1)%halo_local  (1) = IHALO ! west side
       ainfo(1)%halo_local  (2) = IHALO ! east side
       if ( .not. FILE_HISTORY_AGGREGATE ) then
          if( PRC_HAS_W ) ainfo(1)%halo_local(1) = 0
          if( PRC_HAS_E ) ainfo(1)%halo_local(2) = 0
       end if
    endif

    ! for xh
    ainfo(2) = ainfo(1)
    if ( .NOT. PRC_PERIODIC_X .AND. .NOT. FILE_HISTORY_CARTESC_BOUNDARY ) then
       ainfo(2)%size_global (1) = ainfo(2)%size_global (1) + 1
       ainfo(2)%halo_global (1) = ainfo(2)%halo_global (1) + 1
       if ( PRC_HAS_W .and. ( .not. FILE_HISTORY_AGGREGATE ) ) then
          ainfo(2)%start_global(1) = ainfo(2)%start_global(1) + 1
       else
          ainfo(2)%halo_local  (1) = ainfo(2)%halo_local  (1) + 1
       endif
    endif

    ! for y
    if ( PRC_PERIODIC_Y .OR. .NOT. FILE_HISTORY_CARTESC_BOUNDARY ) then
       ainfo(3)%size_global (1) = JMAX * PRC_NUM_Y
       ainfo(3)%start_global(1) = JS_inG - JHALO
       ainfo(3)%halo_global (1) = 0     ! south side
       ainfo(3)%halo_global (2) = 0     ! north side
       ainfo(3)%halo_local  (1) = 0     ! south side
       ainfo(3)%halo_local  (2) = 0     ! north side
    else
       ainfo(3)%size_global (1) = JAG
       ainfo(3)%start_global(1) = JSGA
       ainfo(3)%halo_global (1) = JHALO ! south side
       ainfo(3)%halo_global (2) = JHALO ! north side
       ainfo(3)%halo_local  (1) = JHALO ! south side
       ainfo(3)%halo_local  (2) = JHALO ! north side
       if ( .not. FILE_HISTORY_AGGREGATE ) then
          if( PRC_HAS_S ) ainfo(3)%halo_local(1) = 0
          if( PRC_HAS_N ) ainfo(3)%halo_local(2) = 0
       end if
    endif

    ! for yh
    ainfo(4) = ainfo(3)
    if ( .NOT. PRC_PERIODIC_Y .AND. .NOT. FILE_HISTORY_CARTESC_BOUNDARY ) then
       ainfo(4)%size_global (1) = ainfo(4)%size_global (1) + 1
       ainfo(4)%halo_global (1) = ainfo(4)%halo_global (1) + 1
       if ( PRC_HAS_S .and. ( .not. FILE_HISTORY_AGGREGATE ) ) then
          ainfo(4)%start_global(1) = ainfo(4)%start_global(1) + 1
       else
          ainfo(4)%halo_local  (1) = ainfo(4)%halo_local  (1) + 1
       endif
    endif

    if ( FILE_HISTORY_AGGREGATE ) then
       ainfo(1)%start_global(1) = 1
       ainfo(2)%start_global(1) = 1
       ainfo(3)%start_global(1) = 1
       ainfo(4)%start_global(1) = 1
    end if

    call FILE_HISTORY_Set_Attribute( "x" , "size_global" , ainfo(1)%size_global (:) )
    call FILE_HISTORY_Set_Attribute( "x" , "start_global", ainfo(1)%start_global(:) )
    call FILE_HISTORY_Set_Attribute( "x" , "halo_global" , ainfo(1)%halo_global (:) )
    call FILE_HISTORY_Set_Attribute( "x" , "halo_local"  , ainfo(1)%halo_local  (:) )
    call FILE_HISTORY_Set_Attribute( "x" , "periodic"    , ainfo(1)%periodic        )

    call FILE_HISTORY_Set_Attribute( "xh", "size_global" , ainfo(2)%size_global (:) )
    call FILE_HISTORY_Set_Attribute( "xh", "start_global", ainfo(2)%start_global(:) )
    call FILE_HISTORY_Set_Attribute( "xh", "halo_global" , ainfo(2)%halo_global (:) )
    call FILE_HISTORY_Set_Attribute( "xh", "halo_local"  , ainfo(2)%halo_local  (:) )
    call FILE_HISTORY_Set_Attribute( "xh", "periodic"    , ainfo(2)%periodic        )

    call FILE_HISTORY_Set_Attribute( "y" , "size_global" , ainfo(3)%size_global (:) )
    call FILE_HISTORY_Set_Attribute( "y" , "start_global", ainfo(3)%start_global(:) )
    call FILE_HISTORY_Set_Attribute( "y" , "halo_global" , ainfo(3)%halo_global (:) )
    call FILE_HISTORY_Set_Attribute( "y" , "halo_local"  , ainfo(3)%halo_local  (:) )
    call FILE_HISTORY_Set_Attribute( "y" , "periodic"    , ainfo(3)%periodic        )

    call FILE_HISTORY_Set_Attribute( "yh", "size_global" , ainfo(4)%size_global (:) )
    call FILE_HISTORY_Set_Attribute( "yh", "start_global", ainfo(4)%start_global(:) )
    call FILE_HISTORY_Set_Attribute( "yh", "halo_global" , ainfo(4)%halo_global (:) )
    call FILE_HISTORY_Set_Attribute( "yh", "halo_local"  , ainfo(4)%halo_local  (:) )
    call FILE_HISTORY_Set_Attribute( "yh", "periodic"    , ainfo(4)%periodic        )

    ! map projection info
    call MAPPROJECTION_get_attributes( minfo%mapping_name,                             & ! [OUT]
                              minfo%false_easting                        (1), & ! [OUT]
                              minfo%false_northing                       (1), & ! [OUT]
                              minfo%longitude_of_central_meridian        (1), & ! [OUT]
                              minfo%longitude_of_projection_origin       (1), & ! [OUT]
                              minfo%latitude_of_projection_origin        (1), & ! [OUT]
                              minfo%straight_vertical_longitude_from_pole(1), & ! [OUT]
                              minfo%standard_parallel                    (:)  ) ! [OUT]

    if ( minfo%mapping_name /= "" ) then
       call FILE_HISTORY_Set_Attribute( "x" , "standard_name", "projection_x_coordinate" )
       call FILE_HISTORY_Set_Attribute( "xh", "standard_name", "projection_x_coordinate" )
       call FILE_HISTORY_Set_Attribute( "y" , "standard_name", "projection_y_coordinate" )
       call FILE_HISTORY_Set_Attribute( "yh", "standard_name", "projection_y_coordinate" )

       call FILE_HISTORY_Set_Attribute( minfo%mapping_name, "grid_mapping_name", minfo%mapping_name, add_variable=.true. )

       if ( minfo%false_easting(1) /= UNDEF ) then
          call FILE_HISTORY_Set_Attribute( minfo%mapping_name,    & ! [IN]
                                           "false_easting",       & ! [IN]
                                           minfo%false_easting(:) ) ! [IN]
       endif

       if ( minfo%false_northing(1) /= UNDEF ) then
          call FILE_HISTORY_Set_Attribute( minfo%mapping_name,     & ! [IN]
                                           "false_northing",       & ! [IN]
                                           minfo%false_northing(:) ) ! [IN]
       endif

       if ( minfo%longitude_of_central_meridian(1) /= UNDEF ) then
          call FILE_HISTORY_Set_Attribute( minfo%mapping_name,                    & ! [IN]
                                           "longitude_of_central_meridian",       & ! [IN]
                                           minfo%longitude_of_central_meridian(:) ) ! [IN]
       endif

       if ( minfo%longitude_of_projection_origin(1) /= UNDEF ) then
          call FILE_HISTORY_Set_Attribute( minfo%mapping_name,                     & ! [IN]
                                           "longitude_of_projection_origin",       & ! [IN]
                                           minfo%longitude_of_projection_origin(:) ) ! [IN]
       endif

       if ( minfo%latitude_of_projection_origin(1) /= UNDEF ) then
          call FILE_HISTORY_Set_Attribute( minfo%mapping_name,                    & ! [IN]
                                           "latitude_of_projection_origin",       & ! [IN]
                                           minfo%latitude_of_projection_origin(:) ) ! [IN]
       endif

       if ( minfo%straight_vertical_longitude_from_pole(1) /= UNDEF ) then
          call FILE_HISTORY_Set_Attribute( minfo%mapping_name,                            & ! [IN]
                                           "straight_vertical_longitude_from_pole",       & ! [IN]
                                           minfo%straight_vertical_longitude_from_pole(:) ) ! [IN]
       endif

       if ( minfo%standard_parallel(1) /= UNDEF ) then
          if ( minfo%standard_parallel(2) /= UNDEF ) then
             call FILE_HISTORY_Set_Attribute( minfo%mapping_name,          & ! [IN]
                                              "standard_parallel",         & ! [IN]
                                              minfo%standard_parallel(1:2) ) ! [IN]
          else
             call FILE_HISTORY_Set_Attribute( minfo%mapping_name,          & ! [IN]
                                              "standard_parallel",         & ! [IN]
                                              minfo%standard_parallel(1:1) ) ! [IN]
          endif
       endif
    endif

    ! area and volume
    call FILE_HISTORY_Set_Attribute( "cell_area",    "standard_name", "area" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_area_uy", "standard_name", "area" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_area_xv", "standard_name", "area" ) ! [IN]

    call FILE_HISTORY_Set_Attribute( "cell_area_uyz_x", "standard_name", "area" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_area_xvz_y", "standard_name", "area" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_area_uyw_x", "standard_name", "area" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_area_xvw_y", "standard_name", "area" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_area_xyz_x", "standard_name", "area" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_area_uvz_y", "standard_name", "area" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_area_uvz_x", "standard_name", "area" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_area_xyz_y", "standard_name", "area" ) ! [IN]

    call FILE_HISTORY_Set_Attribute( "cell_volume",     "standard_name", "volume" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_volume_xyw", "standard_name", "volume" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_volume_uyz", "standard_name", "volume" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_volume_xvz", "standard_name", "volume" ) ! [IN]

    call FILE_HISTORY_Set_Attribute( "cell_volume_xyo", "standard_name", "volume" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_volume_xyl", "standard_name", "volume" ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "cell_volume_xyu", "standard_name", "volume" ) ! [IN]

    ! SGRID
    call FILE_HISTORY_Set_Attribute( "grid", "cf_role",             "grid_topology", add_variable=.true. )
    call FILE_HISTORY_Set_Attribute( "grid", "topology_dimension",  (/ 2 /) )
    call FILE_HISTORY_Set_Attribute( "grid", "node_dimensions",     "xh yh" )
    call FILE_HISTORY_Set_Attribute( "grid", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_HISTORY_Set_Attribute( "grid", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_HISTORY_Set_Attribute( "grid", "face_coordinates",    "lon lat" )
    call FILE_HISTORY_Set_Attribute( "grid", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_HISTORY_Set_Attribute( "grid", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_HISTORY_Set_Attribute( "grid", "vertical_dimensions", "z: zh (padding: none)" )

    call FILE_HISTORY_Set_Attribute( "grid_ocean", "cf_role",             "grid_topology", add_variable=.true. )
    call FILE_HISTORY_Set_Attribute( "grid_ocean", "topology_dimension",  (/ 2 /) )
    call FILE_HISTORY_Set_Attribute( "grid_ocean", "node_dimensions",     "xh yh" )
    call FILE_HISTORY_Set_Attribute( "grid_ocean", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_HISTORY_Set_Attribute( "grid_ocean", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_HISTORY_Set_Attribute( "grid_ocean", "face_coordinates",    "lon lat" )
    call FILE_HISTORY_Set_Attribute( "grid_ocean", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_HISTORY_Set_Attribute( "grid_ocean", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_HISTORY_Set_Attribute( "grid_ocean", "vertical_dimensions", "oz: ozh (padding: none)" )

    call FILE_HISTORY_Set_Attribute( "grid_land", "cf_role",             "grid_topology", add_variable=.true. )
    call FILE_HISTORY_Set_Attribute( "grid_land", "topology_dimension",  (/ 2 /) )
    call FILE_HISTORY_Set_Attribute( "grid_land", "node_dimensions",     "xh yh" )
    call FILE_HISTORY_Set_Attribute( "grid_land", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_HISTORY_Set_Attribute( "grid_land", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_HISTORY_Set_Attribute( "grid_land", "face_coordinates",    "lon lat" )
    call FILE_HISTORY_Set_Attribute( "grid_land", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_HISTORY_Set_Attribute( "grid_land", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_HISTORY_Set_Attribute( "grid_land", "vertical_dimensions", "lz: lzh (padding: none)" )

    call FILE_HISTORY_Set_Attribute( "grid_urban", "cf_role",             "grid_topology", add_variable=.true. )
    call FILE_HISTORY_Set_Attribute( "grid_urban", "topology_dimension",  (/ 2 /) )
    call FILE_HISTORY_Set_Attribute( "grid_urban", "node_dimensions",     "xh yh" )
    call FILE_HISTORY_Set_Attribute( "grid_urban", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_HISTORY_Set_Attribute( "grid_urban", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_HISTORY_Set_Attribute( "grid_urban", "face_coordinates",    "lon lat" )
    call FILE_HISTORY_Set_Attribute( "grid_urban", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_HISTORY_Set_Attribute( "grid_urban", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_HISTORY_Set_Attribute( "grid_urban", "vertical_dimensions", "uz: uzh (padding: none)" )

    call FILE_HISTORY_Set_Attribute( "grid_pressure", "cf_role",             "grid_topology", add_variable=.true. )
    call FILE_HISTORY_Set_Attribute( "grid_pressure", "topology_dimension",  (/ 2 /) )
    call FILE_HISTORY_Set_Attribute( "grid_pressure", "node_dimensions",     "xh yh" )
    call FILE_HISTORY_Set_Attribute( "grid_pressure", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_HISTORY_Set_Attribute( "grid_pressure", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_HISTORY_Set_Attribute( "grid_pressure", "face_coordinates",    "lon lat" )
    call FILE_HISTORY_Set_Attribute( "grid_pressure", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_HISTORY_Set_Attribute( "grid_pressure", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_HISTORY_Set_Attribute( "grid_pressure", "vertical_dimensions", "pressure" )

    call FILE_HISTORY_Set_Attribute( "grid_z", "cf_role",             "grid_topology", add_variable=.true. )
    call FILE_HISTORY_Set_Attribute( "grid_z", "topology_dimension",  (/ 2 /) )
    call FILE_HISTORY_Set_Attribute( "grid_z", "node_dimensions",     "xh yh" )
    call FILE_HISTORY_Set_Attribute( "grid_z", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_HISTORY_Set_Attribute( "grid_z", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_HISTORY_Set_Attribute( "grid_z", "face_coordinates",    "lon lat" )
    call FILE_HISTORY_Set_Attribute( "grid_z", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_HISTORY_Set_Attribute( "grid_z", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_HISTORY_Set_Attribute( "grid_z", "vertical_dimensions", "height_xyw: height (padding: none)" )

    call FILE_HISTORY_Set_Attribute( "grid_model", "cf_role",             "grid_topology", add_variable=.true. )
    call FILE_HISTORY_Set_Attribute( "grid_model", "topology_dimension",  (/ 2 /) )
    call FILE_HISTORY_Set_Attribute( "grid_model", "node_dimensions",     "FX FY" )
    call FILE_HISTORY_Set_Attribute( "grid_model", "face_dimensions",     "CX: FY (padding: none) CY: FY (padding: none)" )
    call FILE_HISTORY_Set_Attribute( "grid_model", "vertical_dimensions", "CZ: FZ (padding: none)" )

    call FILE_HISTORY_Set_Attribute( "grid_model_global", "cf_role",             "grid_topology", add_variable=.true. )
    call FILE_HISTORY_Set_Attribute( "grid_model_global", "topology_dimension",  (/ 2 /) )
    call FILE_HISTORY_Set_Attribute( "grid_model_global", "node_dimensions",     "FXG FYG" )
    call FILE_HISTORY_Set_Attribute( "grid_model_global", "face_dimensions",     "CXG: FYG (padding: none) CYG: FYG (padding: none)" )
    call FILE_HISTORY_Set_Attribute( "grid_model_global", "vertical_dimensions", "CZ: FZ (padding: none)" )
    

    return
  end subroutine FILE_HISTORY_CARTESC_set_axes_attributes

end module scale_file_history_cartesC
