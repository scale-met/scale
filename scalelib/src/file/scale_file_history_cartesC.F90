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
  use scale_grid_index
  use scale_ocean_grid_index
  use scale_land_grid_index
  use scale_urban_grid_index
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
                                              "pressure" /)

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
!!$    use scale_grid_cartesC, only: &
!!$       GRID_CARTESC_NAME
    use scale_interpolation, only: &
       INTERP_setup_pres
    use scale_mapproj, only: &
       MPRJ_get_attributes
    implicit none

    character(len=*), parameter :: GRID_CARTESC_NAME = "cartesC"

    integer, parameter :: nlayer_max = 300
    real(RP)           :: FILE_HISTORY_CARTESC_PRES(nlayer_max) !> pressure level to output [hPa]

    NAMELIST / PARAM_FILE_HISTORY_CARTESC / &
       FILE_HISTORY_CARTESC_PRES_nlayer,               &
       FILE_HISTORY_CARTESC_PRES,                      &
       FILE_HISTORY_CARTESC_BOUNDARY

    character(len=H_MID) :: FILE_HISTORY_CARTESCORY_H_TITLE = 'SCALE-RM FILE_HISTORY_CARTESCORY OUTPUT' !< title of the output file
    character(len=H_MID) :: FILE_HISTORY_CARTESCORY_T_SINCE

    character(len=FILE_HSHORT) :: mapping_name
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

       call INTERP_setup_pres( FILE_HISTORY_CARTESC_PRES_nlayer ) ! [IN]
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** FILE_HISTORY_CARTESC_PRES_nlayer is not set.'
       if( IO_L ) write(IO_FID_LOG,*) '*** Output with pressure coordinate is disabled'
    endif



    FILE_HISTORY_CARTESCORY_STARTDATE(:) = TIME_NOWDATE
    FILE_HISTORY_CARTESCORY_STARTMS      = TIME_NOWMS

    start_daysec = TIME_STARTDAYSEC
    if ( TIME_NOWDATE(1) > 0 ) then
       write(FILE_HISTORY_CARTESCORY_T_SINCE,'(I4.4,5(A1,I2.2))')      TIME_NOWDATE(1), &
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

    ! get mapping name
    call MPRJ_get_attributes( mapping_name )

    call FILE_HISTORY_Setup( &
         FILE_HISTORY_CARTESCORY_H_TITLE,              & ! [IN]
         H_SOURCE, H_INSTITUTE,                        & ! [IN]
         GRID_CARTESC_NAME, mapping_name,              & ! [IN]
         start_daysec, TIME_DTSEC,                     & ! [IN]
         time_since = FILE_HISTORY_CARTESCORY_T_SINCE, & ! [IN]
         default_zcoord = 'model',                     & ! [IN]
         myrank = PRC_myrank,                          & ! [IN]
         mpi_comm = PRC_LOCAL_COMM_WORLD               ) ! [IN]

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
       PRES, PRESH, SFC_PRES )
    use scale_interpolation, only: &
       INTERP_update_pres
    implicit none

    real(RP), intent(in) :: PRES    (:,:,:) ! pressure at the full level [Pa]
    real(RP), intent(in) :: PRESH   (:,:,:) ! pressure at the half level [Pa]
    real(RP), intent(in) :: SFC_PRES(  :,:) ! surface pressure           [Pa]
    !---------------------------------------------------------------------------

    if ( FILE_HISTORY_CARTESC_PRES_nlayer > 0 ) then
       call INTERP_update_pres( FILE_HISTORY_CARTESC_PRES_nlayer, & ! [IN]
                                PRES(:,:,:), SFC_PRES(:,:),       & ! [IN]
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
       PRC_IsMaster, &
       PRC_myrank
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_HAS_W, &
       PRC_HAS_S
    use scale_file_history, only: &
       FILE_HIStoRY_set_dim
    implicit none

    character(len=H_SHORT) :: dims(3)

    integer :: start(3), count(3)
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


    !  Virtical 1D
    start(1) = 1
    dims (1) = "z"
    if ( PRC_IsMaster ) then
       count(1) = KMAX
    else
       ! for shared-file parallel I/O, only master rank writes variables with only Z dimension
       count(1) = 0
    end if
    call FILE_HISTORY_Set_Dim( "Z", 1, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]
    dims (1) = "zh"
    if ( PRC_IsMaster ) count(1) = KMAX + 1
    call FILE_HISTORY_Set_Dim( "ZH", 1, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]

    dims (1) = "lz"
    if ( PRC_IsMaster ) count(1) = LKMAX
    call FILE_HISTORY_Set_Dim( "LZ", 1, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]
    dims (1) = "lzh"
    if ( PRC_IsMaster ) count(1) = LKMAX + 1
    call FILE_HISTORY_Set_Dim( "LZH", 1, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]

    dims (1) = "oz"
    if ( PRC_IsMaster ) count(1) = OKMAX
    call FILE_HISTORY_Set_Dim( "OZ", 1, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]
    dims (1) = "ozh"
    if ( PRC_IsMaster ) count(1) = OKMAX + 1
    call FILE_HISTORY_Set_Dim( "OZH", 1, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]

    dims (1) = "uz"
    if ( PRC_IsMaster ) count(1) = UKMAX
    call FILE_HISTORY_Set_Dim( "UZ", 1, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]

    dims (1) = "uzh"
    if ( PRC_IsMaster ) count(1) = UKMAX + 1
    call FILE_HISTORY_Set_Dim( "UZH", 1, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]


    ! X, Y
    dims(1) = 'lon'
    dims(2) = 'lat'
    start(1) = xs; count(1) = xc
    start(2) = ys; count(2) = yc
    call FILE_HISTORY_Set_Dim( "XY", 2, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]
    dims (3) = "height"
    start(3) = 1
    count(3) = KMAX
    call FILE_HISTORY_Set_Dim( "ZXY", 3, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]
    dims (3) = "height_xyw"
    count(3) = KMAX + 1
    call FILE_HISTORY_Set_Dim( "ZHXY", 3, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]
    dims (3) = "oz"
    count(3) = OKMAX
    call FILE_HISTORY_Set_Dim( "OXY", 3, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]
    dims (3) = "ozh"
    count(3) = OKMAX + 1
    call FILE_HISTORY_Set_Dim( "OHXY", 3, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]
    dims (3) = "lz"
    count(3) = LKMAX
    call FILE_HISTORY_Set_Dim( "LXY", 3, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]
    dims (3) = "lzh"
    count(3) = LKMAX + 1
    call FILE_HISTORY_Set_Dim( "LHXY", 3, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]
    dims (3) = "uz"
    count(3) = UKMAX
    call FILE_HISTORY_Set_Dim( "UXY", 3, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]
    dims (3) = "uzh"
    count(3) = UKMAX + 1
    call FILE_HISTORY_Set_Dim( "UHXY", 3, dims(:), nzs, zs(:), start(:), count(:) ) ! [IN]



    ! XH, Y
    dims(1) = 'lon_uy'
    dims(2) = 'lat_uy'
    if ( PRC_HAS_W ) then
       start(1) = xs+1; count(1) = xc
    else
       start(1) = xs  ; count(1) = xc+1
    endif
    call FILE_HISTORY_Set_Dim( "XHY", 2, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]
    dims(3) = "height_uyz"
    count(3) = KMAX
    call FILE_HISTORY_Set_Dim( "ZXHY", 3, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]
    dims(3) = "height_uyw"
    count(3) = KMAX + 1
    call FILE_HISTORY_Set_Dim( "ZHXHY", 3, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]


    ! X, YH
    dims(1) = 'lon_xv'
    dims(2) = 'lat_xv'
    start(1) = xs; count(1) = xc
    if ( PRC_HAS_S ) then
       start(2) = ys+1; count(2) = yc
    else
       start(2) = ys  ; count(2) = yc+1
    endif
    call FILE_HISTORY_Set_Dim( "XYH", 2, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]
    dims(3) = "height_xvz"
    count(3) = KMAX
    call FILE_HISTORY_Set_Dim( "ZXYH", 3, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]
    dims(3) = "height_xvw"
    count(3) = KMAX + 1
    call FILE_HISTORY_Set_Dim( "ZHXYH", 3, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]

    dims(1) = 'lon_uv'
    dims(2) = 'lat_uv'
    if ( PRC_HAS_W ) then
       start(1) = xs+1; count(1) = xc
    else
       start(1) = xs  ; count(1) = xc+1
    endif
    call FILE_HISTORY_Set_Dim( "XHYH", 2, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]
    dims(3) = "height_uvz"
    count(3) = KMAX
    call FILE_HISTORY_Set_Dim( "ZXHYH", 3, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]
    dims(3) = "height_uvw"
    count(3) = KMAX + 1
    call FILE_HISTORY_Set_Dim( "ZHXHYH", 3, dims(:), 1, zs(:), start(:), count(:) ) ! [IN]

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

    integer                :: ksize
    integer                :: kstart

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
    case ('LZ')
       ksize  = LKMAX
       kstart = LKS
    case ('LZH')
       ksize  = LKMAX+1
       kstart = LKS-1
    case ('OZ')
       ksize  = OKMAX
       kstart = OKS
    case ('OZH')
       ksize  = OKMAX+1
       kstart = OKS-1
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
    use scale_interpolation, only: &
       INTERP_available, &
       INTERP_vertical_xi2z, &
       INTERP_vertical_xih2zh, &
       INTERP_vertical_xi2p, &
       INTERP_vertical_xih2p
    implicit none

    real(RP),         intent(in) :: src(:,:,:)
    character(len=*), intent(in) :: dim_type
    character(len=*), intent(in) :: zcoord
    logical,          intent(in) :: fill_halo

    real(DP), intent(out) :: dst(:)

    real(RP) :: var_Z(KA,IA,JA)
    real(RP) :: var_P(FILE_HISTORY_CARTESC_PRES_nlayer,IA,JA)
    integer :: isize, jsize, ksize
    integer :: istart, jstart, kstart

    integer :: i, j, k

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
    case('L')
       ksize  = LKMAX
       kstart = LKS
    case('O')
       ksize  = OKMAX
       kstart = OKS
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
       call INTERP_vertical_xi2z( src  (:,:,:), & ! [IN]
                                  var_Z(:,:,:)  ) ! [OUT]
       call PROF_rapend  ('FILE_O_interp', 2)

       !$omp parallel do
       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = var_Z(kstart+k-1,istart+i-1,jstart+j-1)
       enddo
       enddo
       enddo

    else if( ksize == KMAX+1 .and. zcoord == "z" .and. INTERP_available ) then ! z*->z interpolation (half level)


       call PROF_rapstart('FILE_O_interp', 2)
       call INTERP_vertical_xih2zh( src  (:,:,:), & ! [IN]
                                    var_Z(:,:,:)  ) ! [OUT]
       call PROF_rapend  ('FILE_O_interp', 2)

       !$omp parallel do
       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = var_Z(kstart+k-1,istart+i-1,jstart+j-1)
       enddo
       enddo
       enddo

    elseif( ksize == KMAX .and. zcoord == "pressure" ) then ! z*->p interpolation (full level)
       ksize = FILE_HISTORY_CARTESC_PRES_nlayer

       call PROF_rapstart('FILE_O_interp', 2)
       call INTERP_vertical_xi2p( FILE_HISTORY_CARTESC_PRES_nlayer, & ! [IN]
                                  src  (:,:,:),     & ! [IN]
                                  var_P(:,:,:)      ) ! [OUT]
       call PROF_rapend  ('FILE_O_interp', 2)

       !$omp parallel do
       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = var_P(k,istart+i-1,jstart+j-1)
       enddo
       enddo
       enddo

    elseif( ksize == KMAX+1 .and. zcoord == "pressure" ) then ! z*->p interpolation (half level)
       ksize = FILE_HISTORY_CARTESC_PRES_nlayer

       call PROF_rapstart('FILE_O_interp', 2)
       call INTERP_vertical_xih2p( FILE_HISTORY_CARTESC_PRES_nlayer, & ! [IN]
                                   src  (:,:,:),     & ! [IN]
                                   var_P(:,:,:)      ) ! [OUT]
       call PROF_rapend  ('FILE_O_interp', 2)

       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          dst((k-1)*jsize*isize+(j-1)*isize+i) = var_P(k,istart+i-1,jstart+j-1)
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
    use scale_grid, only: &
       GRID_CZ,    &
       GRID_CX,    &
       GRID_CY,    &
       GRID_FZ,    &
       GRID_FX,    &
       GRID_FY,    &
       GRID_CDZ,   &
       GRID_CDX,   &
       GRID_CDY,   &
       GRID_FDZ,   &
       GRID_FDX,   &
       GRID_FDY,   &
       GRID_CBFZ,  &
       GRID_CBFX,  &
       GRID_CBFY,  &
       GRID_FBFZ,  &
       GRID_FBFX,  &
       GRID_FBFY,  &
       GRID_CXG,   &
       GRID_CYG,   &
       GRID_FXG,   &
       GRID_FYG,   &
       GRID_CDXG,  &
       GRID_CDYG,  &
       GRID_FDXG,  &
       GRID_FDYG,  &
       GRID_CBFXG, &
       GRID_CBFYG, &
       GRID_FBFXG, &
       GRID_FBFYG
    use scale_land_grid, only: &
       GRID_LCZ, &
       GRID_LFZ, &
       GRID_LCDZ
    use scale_ocean_grid, only: &
       GRID_OCZ, &
       GRID_OFZ, &
       GRID_OCDZ
    use scale_urban_grid, only: &
       GRID_UCZ, &
       GRID_UFZ, &
       GRID_UCDZ
    use scale_grid_real, only: &
       REAL_CZ,   &
       REAL_FZ,   &
       REAL_LON,  &
       REAL_LONX, &
       REAL_LONY, &
       REAL_LONXY, &
       REAL_LAT,  &
       REAL_LATX, &
       REAL_LATY, &
       REAL_LATXY
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_landuse, only: &
       LANDUSE_frac_land
    implicit none

    real(RP)         :: AXIS     (imh,jmh,0:KMAX)
    character(len=2) :: AXIS_name(3)

    integer :: k, i, j
    integer :: rankidx(2)
    integer :: start(3,4) !> 1: FF, 2: HF, 3: FH, 4: HH (x,y)
    integer :: startX, startY, startZ
    integer :: startXH, startYH
    integer :: XAG, YAG
    integer :: XAGH, YAGH
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

    ! for the shared-file I/O method, the axes are global (gsize)
    ! for one-file-per-process I/O method, the axes size is equal to the local buffer size

    call FILE_HISTORY_Set_Axis( 'z',   'Z',               'm', 'z',   GRID_CZ (KS  :KE),   gsize=KMAX   , start=startZ )
    call FILE_HISTORY_Set_Axis( 'zh',  'Z (half level)',  'm', 'zh',  GRID_FZ (KS-1:KE),   gsize=KMAX+1 , start=startZ )

    if ( FILE_HISTORY_CARTESC_PRES_nlayer > 0 ) then
       call FILE_HISTORY_Set_Axis( 'pressure', 'Pressure', 'hPa', 'pressure', FILE_HISTORY_CARTESC_PRES_val(:)/100.0_RP, &
                                    gsize=FILE_HISTORY_CARTESC_PRES_nlayer, start=startZ, down=.true. )
    endif

    call FILE_HISTORY_Set_Axis( 'lz',  'LZ',              'm', 'lz',  GRID_LCZ(LKS  :LKE), gsize=LKMAX  , start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'lzh', 'LZ (half level)', 'm', 'lzh', GRID_LFZ(LKS-1:LKE), gsize=LKMAX+1, start=startZ, down=.true. )

    call FILE_HISTORY_Set_Axis( 'oz',  'OZ',              'm', 'oz',  GRID_OCZ(OKS  :OKE), gsize=OKMAX  , start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'ozh', 'OZ (half level)', 'm', 'ozh', GRID_OFZ(OKS-1:OKE), gsize=OKMAX+1, start=startZ, down=.true. )

    call FILE_HISTORY_Set_Axis( 'uz',  'UZ',              'm', 'uz',  GRID_UCZ(UKS  :UKE), gsize=UKMAX  , start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'uzh', 'UZ (half level)', 'm', 'uzh', GRID_UFZ(UKS-1:UKE), gsize=UKMAX+1, start=startZ, down=.true. )

    call FILE_HISTORY_Set_Axis( 'x',   'X',               'm', 'x',   GRID_CX (ims :ime),  gsize=XAG , start=startX  )
    call FILE_HISTORY_Set_Axis( 'xh',  'X (half level)',  'm', 'xh',  GRID_FX (imsh:ime),  gsize=XAGH, start=startXH )

    call FILE_HISTORY_Set_Axis( 'y',   'Y',               'm', 'y',   GRID_CY (jms :jme),  gsize=YAG , start=startY  )
    call FILE_HISTORY_Set_Axis( 'yh',  'Y (half level)',  'm', 'yh',  GRID_FY (jmsh:jme),  gsize=YAGH, start=startYH )

    ! axes below always include halos when written to file regardless of PRC_PERIODIC_X/PRC_PERIODIC_Y
    call FILE_HISTORY_Set_Axis( 'CZ',   'Atmos Grid Center Position Z', 'm', 'CZ',  GRID_CZ,   gsize=KA,      start=startZ )
    call FILE_HISTORY_Set_Axis( 'FZ',   'Atmos Grid Face Position Z',   'm', 'FZ',  GRID_FZ,   gsize=KA+1,    start=startZ )
    call FILE_HISTORY_Set_Axis( 'CDZ',  'Grid Cell length Z',           'm', 'CZ',  GRID_CDZ,  gsize=KA,      start=startZ )
    call FILE_HISTORY_Set_Axis( 'FDZ',  'Grid distance Z',              'm', 'FDZ', GRID_FDZ,  gsize=KA-1,    start=startZ )
    call FILE_HISTORY_Set_Axis( 'CBFZ', 'Boundary factor Center Z',     '1', 'CZ',  GRID_CBFZ, gsize=KA,      start=startZ )
    call FILE_HISTORY_Set_Axis( 'FBFZ', 'Boundary factor Face Z',       '1', 'FZ',  GRID_FBFZ, gsize=KA+1,    start=startZ )

    call FILE_HISTORY_Set_Axis( 'LCZ',  'Land Grid Center Position Z',  'm', 'LCZ', GRID_LCZ,  gsize=LKMAX,   start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'LFZ',  'Land Grid Face Position Z',    'm', 'LFZ', GRID_LFZ,  gsize=LKMAX+1, start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'LCDZ', 'Land Grid Cell length Z',      'm', 'LCZ', GRID_LCDZ, gsize=LKMAX,   start=startZ              )

    call FILE_HISTORY_Set_Axis( 'UCZ',  'Urban Grid Center Position Z', 'm', 'UCZ', GRID_UCZ,  gsize=UKMAX,   start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'UFZ',  'Urban Grid Face Position Z',   'm', 'UFZ', GRID_UFZ,  gsize=UKMAX+1, start=startZ, down=.true. )
    call FILE_HISTORY_Set_Axis( 'UCDZ', 'Urban Grid Cell length Z',     'm', 'UCZ', GRID_UCDZ, gsize=UKMAX,   start=startZ              )

    if ( FILE_HISTORY_AGGREGATE ) then
       call FILE_HISTORY_Set_Axis( 'CX',   'Atmos Grid Center Position X', 'm', 'CX',  GRID_CXG,   gsize=IAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'CY',   'Atmos Grid Center Position Y', 'm', 'CY',  GRID_CYG,   gsize=JAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'FX',   'Atmos Grid Face Position X',   'm', 'FX',  GRID_FXG,   gsize=IAG+1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'FY',   'Atmos Grid Face Position Y',   'm', 'FY',  GRID_FYG,   gsize=JAG+1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'CDX',  'Grid Cell length X',           'm', 'CX',  GRID_CDXG,  gsize=IAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'CDY',  'Grid Cell length Y',           'm', 'CY',  GRID_CDYG,  gsize=JAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'FDX',  'Grid distance X',              'm', 'FDX', GRID_FDXG,  gsize=IAG-1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'FDY',  'Grid distance Y',              'm', 'FDY', GRID_FDYG,  gsize=JAG-1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'CBFX', 'Boundary factor Center X',     '1', 'CX',  GRID_CBFXG, gsize=IAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'CBFY', 'Boundary factor Center Y',     '1', 'CY',  GRID_CBFYG, gsize=JAG,   start=startZ )
       call FILE_HISTORY_Set_Axis( 'FBFX', 'Boundary factor Face X',       '1', 'FX',  GRID_FBFXG, gsize=IAG+1, start=startZ )
       call FILE_HISTORY_Set_Axis( 'FBFY', 'Boundary factor Face Y',       '1', 'FY',  GRID_FBFYG, gsize=JAG+1, start=startZ )
    else
       call FILE_HISTORY_Set_Axis( 'CX',   'Atmos Grid Center Position X', 'm', 'CX',  GRID_CX   )
       call FILE_HISTORY_Set_Axis( 'CY',   'Atmos Grid Center Position Y', 'm', 'CY',  GRID_CY   )
       call FILE_HISTORY_Set_Axis( 'FX',   'Atmos Grid Face Position X',   'm', 'FX',  GRID_FX   )
       call FILE_HISTORY_Set_Axis( 'FY',   'Atmos Grid Face Position Y',   'm', 'FY',  GRID_FY   )
       call FILE_HISTORY_Set_Axis( 'CDX',  'Grid Cell length X',           'm', 'CX',  GRID_CDX  )
       call FILE_HISTORY_Set_Axis( 'CDY',  'Grid Cell length Y',           'm', 'CY',  GRID_CDY  )
       call FILE_HISTORY_Set_Axis( 'FDX',  'Grid distance X',              'm', 'FDX', GRID_FDX  )
       call FILE_HISTORY_Set_Axis( 'FDY',  'Grid distance Y',              'm', 'FDY', GRID_FDY  )
       call FILE_HISTORY_Set_Axis( 'CBFX', 'Boundary factor Center X',     '1', 'CX',  GRID_CBFX )
       call FILE_HISTORY_Set_Axis( 'CBFY', 'Boundary factor Center Y',     '1', 'CY',  GRID_CBFY )
       call FILE_HISTORY_Set_Axis( 'FBFX', 'Boundary factor Face X',       '1', 'FX',  GRID_FBFX )
       call FILE_HISTORY_Set_Axis( 'FBFY', 'Boundary factor Face Y',       '1', 'FY',  GRID_FBFY )
    endif

    call FILE_HISTORY_Set_Axis('CXG',   'Grid Center Position X (global)',   'm', 'CXG',  GRID_CXG,   gsize=IAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('CYG',   'Grid Center Position Y (global)',   'm', 'CYG',  GRID_CYG,   gsize=JAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('FXG',   'Grid Face Position X (global)',     'm', 'FXG',  GRID_FXG,   gsize=IAG+1, start=startZ )
    call FILE_HISTORY_Set_Axis('FYG',   'Grid Face Position Y (global)',     'm', 'FYG',  GRID_FYG,   gsize=JAG+1, start=startZ )
    call FILE_HISTORY_Set_Axis('CDXG',  'Grid Cell length X (global)',       'm', 'CXG',  GRID_CDXG,  gsize=IAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('CDYG',  'Grid Cell length Y (global)',       'm', 'CYG',  GRID_CDYG,  gsize=JAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('FDXG',  'Grid distance X (global)',          'm', 'FDXG', GRID_FDXG,  gsize=IAG-1, start=startZ )
    call FILE_HISTORY_Set_Axis('FDYG',  'Grid distance Y (global)',          'm', 'FDYG', GRID_FDYG,  gsize=JAG-1, start=startZ )
    call FILE_HISTORY_Set_Axis('CBFXG', 'Boundary factor Center X (global)', '1', 'CXG',  GRID_CBFXG, gsize=IAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG',  GRID_CBFYG, gsize=JAG,   start=startZ )
    call FILE_HISTORY_Set_Axis('FBFXG', 'Boundary factor Face X (global)',   '1', 'FXG',  GRID_FBFXG, gsize=IAG+1, start=startZ )
    call FILE_HISTORY_Set_Axis('FBFYG', 'Boundary factor Face Y (global)',   '1', 'FYG',  GRID_FBFYG, gsize=JAG+1, start=startZ )



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
       AXIS(1:im,1:jm,k) = REAL_CZ(k+KS-1,ims:ime,jms:jme)
    enddo
    AXIS_name(1:3) = (/'x ', 'y ', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height', 'height above ground level',                        &
                                          'm', AXIS_name(1:3), AXIS(1:im,1:jm,1:KMAX), start=start(:,1) )

    do k = 0, KMAX
       AXIS(1:im,1:jm,k) = REAL_FZ(k+KS-1,ims:ime,jms:jme)
    enddo
    AXIS_name(1:3) = (/'x ', 'y ', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_xyw', 'height above ground level (half level xyw)',    &
                                          'm' , AXIS_name(1:3), AXIS(1:im,1:jm,0:KMAX), start=start(:,1) )

    do k = 1, KMAX
    do j = 1, jm
    do i = 1, min(imh,IA-imsh)
       AXIS(i,j,k) = ( REAL_CZ(k+KS-1,imsh+i-1,jms+j-1) + REAL_CZ(k+KS-1,imsh+i,jms+j-1) ) * 0.5_RP
    enddo
    enddo
    enddo
    if ( imh == IA-imsh+1 ) then
       do k = 1, KMAX
       do j = 1, jm
          AXIS(imh,j,k) = REAL_CZ(k+KS-1,imsh+imh-1,jms+j-1)
       enddo
       enddo
    endif
    AXIS_name(1:3) = (/'xh', 'y ', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_uyz', 'height above ground level (half level uyz)',    &
                                          'm', AXIS_name(1:3), AXIS(1:imh,1:jm,1:KMAX), start=start(:,2) )

    do k = 1, KMAX
    do j = 1, min(jmh,JA-jmsh)
    do i = 1, im
       AXIS(i,j,k) = ( REAL_CZ(k+KS-1,ims+i-1,jmsh+j-1) + REAL_CZ(k+KS-1,ims+i-1,jmsh+j) ) * 0.5_RP
    enddo
    enddo
    enddo
    if ( jmh == JA-jmsh+1 ) then
       do k = 1, KMAX
       do i = 1, im
          AXIS(i,jmh,k) = REAL_CZ(k+KS-1,ims+i-1,jmsh+jmh-1)
       enddo
       enddo
    endif
    AXIS_name(1:3) = (/'x ', 'yh', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_xvz', 'height above ground level (half level xvz)',    &
                                          'm', AXIS_name(1:3), AXIS(1:im,1:jmh,1:KMAX), start=start(:,3) )

    do k = 1, KMAX
    do j = 1, min(jmh,JA-jmsh)
    do i = 1, min(imh,IA-imsh)
       AXIS(i,j,k) = ( REAL_CZ(k+KS-1,imsh+i-1,jmsh+j-1) + REAL_CZ(k+KS-1,imsh+i  ,jmsh+j-1) &
                     + REAL_CZ(k+KS-1,imsh+i-1,jmsh+j  ) + REAL_CZ(k+KS-1,imsh+i  ,jmsh+j  ) ) * 0.25_RP
    enddo
    enddo
    enddo
    if ( jmh == JA-jmsh+1 ) then
       do k = 1, KMAX
       do i = 1, min(imh,IA-imsh)
          AXIS(i,jmh,k) = ( REAL_CZ(k+KS-1,imsh+i-1,jmsh+jmh-1) + REAL_CZ(k+KS-1,imsh+i,jmsh+jmh-1) ) * 0.5_RP
       enddo
       enddo
    endif
    if ( imh == IA-imsh+1 ) then
       do k = 1, KMAX
       do j = 1, min(jmh,JA-jmsh)
          AXIS(imh,j,k) = ( REAL_CZ(k+KS-1,imsh+imh-1,jmsh+j-1) + REAL_CZ(k+KS-1,imsh+imh-1,jmsh+j) ) * 0.5_RP
       enddo
       enddo
    endif
    if ( imh == IA-imsh+1 .AND. jmh == JA-jmsh+1 ) then
       do k = 1, KMAX
          AXIS(imh,jmh,k) = REAL_CZ(k+KS-1,imsh+imh-1,jmsh+jmh-1)
       enddo
    endif
    AXIS_name(1:3) = (/'xh', 'yh', 'z '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_uvz', 'height above ground level (half level uvz)',     &
                                          'm', AXIS_name(1:3), AXIS(1:imh,1:jmh,1:KMAX), start=start(:,4) )

    do k = 0, KMAX
    do j = 1, jm
    do i = 1, min(imh,IA-imsh)
       AXIS(i,j,k) = ( REAL_FZ(k+KS-1,imsh+i-1,jms+j-1) + REAL_FZ(k+KS-1,imsh+i,jms+j-1) ) * 0.5_RP
    enddo
    enddo
    enddo
    if ( imh == IA-imsh+1 ) then
       do k = 0, KMAX
       do j = 1, jm
          AXIS(imh,j,k) = REAL_FZ(k+KS-1,imsh+imh-1,jms+j-1)
       enddo
       enddo
    endif
    AXIS_name(1:3) = (/'xh', 'y ', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_uyw', 'height above ground level (half level uyw)',    &
                                          'm', AXIS_name(1:3), AXIS(1:imh,1:jm,0:KMAX), start=start(:,2) )

    do k = 0, KMAX
    do j = 1, min(jmh,JA-jmsh)
    do i = 1, im
       AXIS(i,j,k) = ( REAL_FZ(k+KS-1,ims+i-1,jmsh+j-1) + REAL_FZ(k+KS-1,ims+i-1,jmsh+j) ) * 0.5_RP
    enddo
    enddo
    enddo
    if ( jmh == JA-jmsh+1 ) then
       do k = 0, KMAX
       do i = 1, im
          AXIS(i,jmh,k) = REAL_FZ(k+KS-1,ims+i-1,jmsh+jmh-1)
       enddo
       enddo
    endif
    AXIS_name(1:3) = (/'x ', 'yh', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_xvw', 'height above ground level (half level xvw)',    &
                                          'm', AXIS_name(1:3), AXIS(1:im,1:jmh,0:KMAX), start=start(:,3) )

    do k = 0, KMAX
    do j = 1, min(jmh,JA-jmsh)
    do i = 1, min(imh,IA-imsh)
       AXIS(i,j,k) = ( REAL_FZ(k+KS-1,imsh+i-1,jmsh+j-1) + REAL_FZ(k+KS-1,imsh+i  ,jmsh+j-1) &
                     + REAL_FZ(k+KS-1,imsh+i-1,jmsh+j  ) + REAL_FZ(k+KS-1,imsh+i  ,jmsh+j  ) ) * 0.25_RP
    enddo
    enddo
    enddo
    if ( jmh == JA-jmsh+1 ) then
       do k = 0, KMAX
       do i = 1, min(imh,IA-imsh)
          AXIS(i,jmh,k) = ( REAL_FZ(k+KS-1,imsh+i-1,jmsh+jmh-1) + REAL_FZ(k+KS-1,imsh+i,jmsh+jmh-1) ) * 0.5_RP
       enddo
       enddo
    endif
    if ( imh == IA-imsh+1 ) then
       do k = 0, KMAX
       do j = 1, min(jmh,JA-jmsh)
          AXIS(imh,j,k) = ( REAL_FZ(k+KS-1,imsh+imh-1,jmsh+j-1) + REAL_FZ(k+KS-1,imsh+imh-1,jmsh+j) ) * 0.5_RP
       enddo
       enddo
    endif
    if ( imh == IA-imsh+1 .AND. jm == JA-jmsh+1 ) then
       do k = 0, KMAX
          AXIS(imh,jmh,k) = REAL_FZ(k+KS-1,imsh+imh-1,jmsh+jmh-1)
       enddo
    endif
    AXIS_name(1:3) = (/'xh', 'yh', 'zh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'height_uvw', 'height above ground level (half level uvw)',     &
                                          'm', AXIS_name(1:3), AXIS(1:imh,1:jmh,0:KMAX), start=start(:,4) )

    AXIS(1:im,1:jm,1) = REAL_LON (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lon', 'longitude',                                                 &
                                          'degrees_east', AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1) )

    AXIS(1:imh,1:jm,1) = REAL_LONX(imsh:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lon_uy', 'longitude (half level uy)',                               &
                                          'degrees_east', AXIS_name(1:2), AXIS(1:imh,1:jm,1), start=start(:,2) )

    AXIS(1:im,1:jmh,1) = REAL_LONY(ims:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'yh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lon_xv', 'longitude (half level xv)',                               &
                                          'degrees_east', AXIS_name(1:2), AXIS(1:im,1:jmh,1), start=start(:,3) )

    AXIS(1:imh,1:jmh,1) = REAL_LONXY(imsh:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'yh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lon_uv', 'longitude (half level uv)',                                &
                                          'degrees_east', AXIS_name(1:2), AXIS(1:imh,1:jmh,1), start=start(:,4) )

    AXIS(1:im,1:jm,1) = REAL_LAT (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lat', 'latitude',                                                   &
                                          'degrees_north', AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1) )

    AXIS(1:imh,1:jm,1) = REAL_LATX(imsh:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lat_uy', 'latitude (half level uy)',                                 &
                                          'degrees_north', AXIS_name(1:2), AXIS(1:imh,1:jm,1), start=start(:,2) )

    AXIS(1:im,1:jmh,1) = REAL_LATY(ims:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'yh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lat_xv', 'latitude (half level xv)',                                 &
                                          'degrees_north', AXIS_name(1:2), AXIS(1:im,1:jmh,1), start=start(:,3) )

    AXIS(1:imh,1:jmh,1) = REAL_LATXY(imsh:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'yh'/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lat_uv', 'latitude (half level uv)',                                  &
                                          'degrees_north', AXIS_name(1:2), AXIS(1:imh,1:jmh,1), start=start(:,4) )

    AXIS(1:im,1:jm,1) = TOPO_Zsfc(ims:ime,jms:jme)
    AXIS_name(1:2) = (/'x ', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'topo', 'topography',                                    &
                                          'm', AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1) )

    AXIS(1:im,1:jm,1) = LANDUSE_frac_land(ims:ime,jms:jme)
    AXIS_name(1:2) = (/'x ', 'y '/)
    call FILE_HISTORY_Set_AssociatedCoordinate( 'lsmask', 'fraction for land-sea mask',                  &
                                          '1', AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1) )

    return
  end subroutine FILE_HISTORY_CARTESC_set_axes

  !-----------------------------------------------------------------------------
  subroutine FILE_HISTORY_CARTESC_set_axes_attributes
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
    use scale_mapproj, only: &
       MPRJ_get_attributes
    implicit none

    character(len=5)  :: periodic_z, periodic_x, periodic_y
    character(len=34) :: tunits

    type(axisattinfo) :: ainfo(4) ! x, xh, y, yh
    type(mappinginfo) :: minfo
    !---------------------------------------------------------------------------

    periodic_z = "false"
    if ( PRC_PERIODIC_X ) then
       periodic_x = "true"
    else
       periodic_x = "false"
    endif
    if ( PRC_PERIODIC_Y ) then
       periodic_y = "true"
    else
       periodic_y = "false"
    endif

    if ( .NOT. FILE_HISTORY_AGGREGATE ) then
       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_rank_x", (/PRC_2Drank(PRC_myrank,1)/) ) ! [IN]
       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_rank_y", (/PRC_2Drank(PRC_myrank,2)/) ) ! [IN]

       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_num_x", (/PRC_NUM_X/) ) ! [IN]
       call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_num_y", (/PRC_NUM_Y/) ) ! [IN]
    endif

    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_periodic_z", periodic_z ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_periodic_x", periodic_x ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_prc_periodic_y", periodic_y ) ! [IN]

    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_grid_index_kmax",  (/KMAX/)  ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_grid_index_imaxg", (/IMAXG/) ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_grid_index_jmaxg", (/JMAXG/) ) ! [IN]

    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_grid_index_khalo", (/KHALO/) ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_grid_index_ihalo", (/IHALO/) ) ! [IN]
    call FILE_HISTORY_Set_Attribute( "global", "scale_cartesC_grid_index_jhalo", (/JHALO/) ) ! [IN]

    call FILE_get_CFtunits( FILE_HISTORY_CARTESCORY_STARTDATE(:), tunits )
    call FILE_HISTORY_Set_Attribute( "global", "time_units", tunits )
    call FILE_HISTORY_Set_Attribute( "global", "time_start", (/FILE_HISTORY_CARTESCORY_STARTMS/) )

    if ( PRC_PERIODIC_X ) then
       ainfo(1)%periodic = "true"
       ainfo(2)%periodic = "true"
    else
       ainfo(1)%periodic = "false"
       ainfo(2)%periodic = "false"
    endif

    if ( PRC_PERIODIC_Y ) then
       ainfo(3)%periodic = "true"
       ainfo(4)%periodic = "true"
    else
       ainfo(3)%periodic = "false"
       ainfo(4)%periodic = "false"
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
    call MPRJ_get_attributes( minfo%mapping_name,                             & ! [OUT]
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

    return
  end subroutine FILE_HISTORY_CARTESC_set_axes_attributes

end module scale_file_history_cartesC
