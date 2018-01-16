!-------------------------------------------------------------------------------
!> module HISTORY
!!
!! @par Description
!!          History output module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-05 (H.Yashiro)   [new]
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-06-11 (S.Nishizawa) [mod] use gtool_history
!!
!<
!-------------------------------------------------------------------------------
module scale_history
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_land_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: HIST_setup
  public :: HIST_switch
  public :: HIST_setpres
  public :: HIST_reg
  public :: HIST_query
  public :: HIST_put
  public :: HIST_in
  public :: HIST_get
  public :: HIST_write

  interface HIST_put
     module procedure HIST_put_0D
     module procedure HIST_put_1D
     module procedure HIST_put_2D
     module procedure HIST_put_3D
  end interface HIST_put

  interface HIST_in
     module procedure HIST_in_0D
     module procedure HIST_in_1D
     module procedure HIST_in_2D
     module procedure HIST_in_3D
  end interface HIST_in

  interface HIST_get
     module procedure HIST_get_1D
     module procedure HIST_get_2D
     module procedure HIST_get_3D
  end interface HIST_get

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: HIST_put_axes
  private :: HIST_set_axes_attributes

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter   :: I_MODEL = 0            !< model coordinate
  integer,                private, parameter   :: I_Z     = 1            !< z coordinate
  integer,                private, parameter   :: I_PRES  = 2            !< pressure coordinate

  integer,                private              :: HIST_item_limit
  integer,                private              :: HIST_item_count        !< number of the history item
  character(len=H_SHORT), private, allocatable :: HIST_item   (:)        !< name   of the history item
  integer,                private, allocatable :: HIST_variant(:)        !< number of the variants      for each history item
  integer,                private, allocatable :: HIST_zcoord (:,:)      !< vertical interpolation type for each variant of history item

  integer,                private, parameter   :: HIST_PRES_nlim   = 300 !< limit  of index size for pressure layer
  integer,                private              :: HIST_PRES_nlayer =  -1 !< Number of pressure layer
  real(RP),               private, allocatable :: HIST_PRES_val(:)       !< pressure level to output [hPa]

  logical,                private              :: enabled
  integer,                private              :: im,   jm,  km
  integer,                private              :: ims,  ime
  integer,                private              :: jms,  jme
  integer,                private              :: imh,  jmh
  integer,                private              :: imsh, jmsh

  logical, private :: HIST_BND = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine HIST_setup
    use gtool_history, only: &
       HistoryInit
    use scale_process, only: &
       PRC_MPIstop,    &
       PRC_masterrank, &
       PRC_myrank
    use scale_rm_process, only: &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y, &
       PRC_HAS_W,      &
       PRC_HAS_S
    use scale_time, only: &
       TIME_NOWDATE,     &
       TIME_NOWMS,       &
       TIME_DTSEC,       &
       TIME_STARTDAYSEC
    use scale_interpolation, only: &
       INTERP_setup_pres
    implicit none

    character(len=H_MID) :: HISTORY_H_TITLE = 'SCALE-RM HISTORY OUTPUT' !< title of the output file
    character(len=H_MID) :: HISTORY_T_SINCE

    real(RP)             :: HIST_PRES(HIST_PRES_nlim) = 0.0_RP          !< pressure level to output [hPa]

    NAMELIST / PARAM_HIST / &
       HIST_PRES_nlayer, &
       HIST_PRES,        &
       HIST_BND

    integer  :: HIST_variant_limit

    real(DP) :: start_daysec
    integer  :: ierr
    integer  :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[HISTORY] / Categ[ATMOS-RM IO] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_HIST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_HIST. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_HIST)

    ! check pressure coordinate
    if ( HIST_PRES_nlayer > 0 ) then
       if ( HIST_PRES_nlayer > 100 ) then
          write(*,*) 'xxx number of layers of pressure is larger the KMAX'
          call PRC_MPIstop
       endif

       do k = 1, HIST_PRES_nlayer
          if ( HIST_PRES(k) <= 0.0_RP ) then
             write(*,*) 'xxx Invalid value found in pressure coordinate! (k,value)=', k, HIST_PRES(k)
             call PRC_MPIstop
          elseif ( HIST_PRES(k+1) >= HIST_PRES(k) ) then
             write(*,*) 'xxx The value of pressure coordinate must be descending order! ', &
                        '(k,value[k],value[k+1])=', k, HIST_PRES(k), HIST_PRES(k+1)
             call PRC_MPIstop
          endif
       enddo
    else
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** HIST_PRES_nlayer is not set.'
       if( IO_L ) write(IO_FID_LOG,*) '*** Output with pressure coordinate is disabled'
    endif

    call PROF_rapstart('FILE_O_NetCDF', 2)

    start_daysec = TIME_STARTDAYSEC
    if ( TIME_NOWDATE(1) > 0 ) then
       write(HISTORY_T_SINCE,'(I4.4,5(A1,I2.2))')      TIME_NOWDATE(1), &
                                                  '-', TIME_NOWDATE(2), &
                                                  '-', TIME_NOWDATE(3), &
                                                  ' ', TIME_NOWDATE(4), &
                                                  ':', TIME_NOWDATE(5), &
                                                  ':', TIME_NOWDATE(6)
       start_daysec = TIME_NOWMS
    else
       HISTORY_T_SINCE = ''
    endif

    if ( HIST_BND ) then
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

    km = max( LKMAX+1, UKMAX+1, KMAX+1, HIST_PRES_nlayer )

    call HistoryInit( HIST_item_limit,                  & ! [OUT]
                      HIST_variant_limit,               & ! [OUT]
                      imh, jmh, km,                     & ! [IN]
                      PRC_masterrank,                   & ! [IN]
                      PRC_myrank,                       & ! [IN]
                      HISTORY_H_TITLE,                  & ! [IN]
                      H_SOURCE,                         & ! [IN]
                      H_INSTITUTE,                      & ! [IN]
                      time_start     = start_daysec,    & ! [IN]
                      time_interval  = TIME_DTSEC,      & ! [IN]
                      time_since     = HISTORY_T_SINCE, & ! [IN]
                      default_zcoord = 'model',         & ! [IN]
                      namelist_fid   = IO_FID_CONF      ) ! [IN]

    HIST_item_count = 0
    if ( HIST_item_limit > 0 ) then
       allocate( HIST_item   (HIST_item_limit)                    )
       allocate( HIST_variant(HIST_item_limit)                    )
       allocate( HIST_zcoord (HIST_item_limit,HIST_variant_limit) )
       HIST_item   (:)   = ''
       HIST_variant(:)   = 0
       HIST_zcoord (:,:) = 0
    endif

    if ( HIST_PRES_nlayer > 0 ) then
       allocate( HIST_PRES_val(HIST_PRES_nlayer) )

       do k = 1, HIST_PRES_nlayer
          HIST_PRES_val(k) = HIST_PRES(k) * 100.0_RP ! [hPa->Pa]
       enddo

       call INTERP_setup_pres( HIST_PRES_nlayer ) ! [IN]
    endif

    call HIST_put_axes

    enabled = .true.

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_setup

  !-----------------------------------------------------------------------------
  !> set switch
  subroutine HIST_switch( switch )
    implicit none

    logical, intent(in) :: switch
    !---------------------------------------------------------------------------

    enabled = switch

    return
  end subroutine HIST_switch

  !-----------------------------------------------------------------------------
  !> set interpolation factor for pressure coordinate
  subroutine HIST_setpres( &
       PRES,    &
       SFC_PRES )
    use scale_interpolation, only: &
       INTERP_update_pres
    implicit none

    real(RP), intent(in) :: PRES    (KA,IA,JA) ! pressure in Xi coordinate [Pa]
    real(RP), intent(in) :: SFC_PRES(   IA,JA) ! surface pressure          [Pa]
    !---------------------------------------------------------------------------

    if ( HIST_PRES_nlayer > 0 ) then
       call INTERP_update_pres( HIST_PRES_nlayer,     & ! [IN]
                                PRES         (:,:,:), & ! [IN]
                                SFC_PRES     (:,:)  , & ! [IN]
                                HIST_PRES_val(:)      ) ! [IN]
    endif

    return
  end subroutine HIST_setpres

  !-----------------------------------------------------------------------------
  !> Register/Append variable to history file
  subroutine HIST_reg( &
       itemid, &
       item,   &
       desc,   &
       unit,   &
       ndim,   &
       xdim,   &
       ydim,   &
       zdim    )
    use gtool_history, only: &
       HistoryAddVariable
    use scale_time, only: &
       TIME_NOWSTEP
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_myrank
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_HAS_W, &
       PRC_HAS_S
    use scale_mapproj, only: &
       MPRJ_get_attributes
    implicit none

    integer,          intent(out) :: itemid !< index number of the item
    character(len=*), intent(in)  :: item   !< name         of the item
    character(len=*), intent(in)  :: desc   !< description  of the item
    character(len=*), intent(in)  :: unit   !< unit         of the item
    integer,          intent(in)  :: ndim   !< dimension    of the item

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim

    character(len=H_SHORT) :: dims(3)
    logical                :: flag_half_x
    logical                :: flag_half_y
    logical                :: flag_half_z
    logical                :: atom
    character(len=H_SHORT) :: mapping_name

    integer                :: start(4), count(4)

    integer :: nvariant1, nvariant2, nvariant3
    integer :: v, id
    !---------------------------------------------------------------------------

    itemid = -1

    if( .NOT. enabled ) return

    if( HIST_item_limit == 0 ) return

    do id = 1, HIST_item_count
       if ( item == HIST_item(id) ) then ! item exists
          itemid = id
          return
       endif
    enddo

    call PROF_rapstart('FILE_O_NetCDF', 2)

    ! Try to add new item

    if ( len_trim(item) >= H_SHORT ) then
       write(*,'(1x,A,I2,A,A)') 'xxx Length of history name should be <= ', H_SHORT-1 ,' chars. STOP', trim(item)
       call PRC_MPIstop
    endif

    atom = .true.

    start(:) = 0
    count(:) = 0

    if ( ndim == 1 ) then

       ! check half/full level for vertical
       flag_half_z = .false.
       if ( present(zdim) ) then
          if( zdim == 'half' ) flag_half_z = .true.
       endif

       start(1) = 1

       if ( flag_half_z ) then
          dims (1) = "zh"
          count(1) = KMAX + 1
       else
          dims (1) = "z"
          count(1) = KMAX
       endif

       ! for shared-file parallel I/O, only rank 0 writes variables with only Z dimension
       if ( PRC_myrank > 0 ) count(1) = 0

       mapping_name = ""

    elseif ( ndim == 2 ) then

       ! check half/full level for horizontal
       flag_half_x = .false.
       if ( present(xdim) ) then
          if( xdim == 'half' ) flag_half_x = .true.
       endif

       flag_half_y = .false.
       if ( present(ydim) ) then
          if( ydim == 'half' ) flag_half_y = .true.
       endif

       if    ( flag_half_x .AND. flag_half_y ) then
          dims(1) = 'lon_uv'
          dims(2) = 'lat_uv'
       elseif( flag_half_x ) then
          dims(1) = 'lon_uy'
          dims(2) = 'lat_uy'
       elseif( flag_half_y ) then
          dims(1) = 'lon_xv'
          dims(2) = 'lat_xv'
       else
          dims(1) = 'lon'
          dims(2) = 'lat'
       endif

       call MPRJ_get_attributes( mapping_name )

    elseif ( ndim == 3 ) then

       ! check half/full level for vertical/horizontal
       flag_half_z = .false.
       if ( present(zdim) ) then
          if( zdim == 'half' ) flag_half_z = .true.
       endif

       flag_half_x = .false.
       if ( present(xdim) ) then
          if( xdim == 'half' ) flag_half_x = .true.
       endif

       flag_half_y = .false.
       if ( present(ydim) ) then
          if( ydim == 'half' ) flag_half_y = .true.
       endif

       if    ( flag_half_x .AND. flag_half_y ) then
          dims(1) = 'lon_uv'
          dims(2) = 'lat_uv'
          if ( flag_half_z ) then
             dims(3) = 'height_uvw'
          else
             dims(3) = 'height_uvz'
          endif
       elseif( flag_half_x ) then
          dims(1) = 'lon_uy'
          dims(2) = 'lat_uy'
          if ( flag_half_z ) then
             dims(3) = 'height_uyw'
          else
             dims(3) = 'height_uyz'
          endif
       elseif( flag_half_y ) then
          dims(1) = 'lon_xv'
          dims(2) = 'lat_xv'
          if ( flag_half_z ) then
             dims(3) = 'height_xvw'
          else
             dims(3) = 'height_xvz'
          endif
       else
          dims(1) = 'lon'
          dims(2) = 'lat'
          if ( flag_half_z ) then
             dims(3) = 'height_xyw'
          else
             dims(3) = 'height'
          endif
       endif

       ! start and count will be used by PnetCDF I/O
       start(3) = 1
       count(3) = KMAX

       if ( present(zdim) ) then
          if    ( zdim == 'land'      ) then
             dims (3) = 'lz'
             count(3) = LKMAX
             atom     = .false.
          elseif( zdim == 'landhalf'  ) then
             dims (3) = 'lzh'
             count(3) = LKMAX+1
             atom     = .false.
             flag_half_z = .true.
          elseif( zdim == 'urban'     ) then
             dims (3) = 'uz'
             count(3) = UKMAX
             atom     = .false.
          elseif( zdim == 'urbanhalf' ) then
             dims (3) = 'uzh'
             count(3) = UKMAX+1
             atom     = .false.
             flag_half_z = .true.
          endif
       endif

       call MPRJ_get_attributes( mapping_name )

    endif

    if ( ndim >= 2 ) then
       ! start and count will be used by PnetCDF I/O
       if ( HIST_BND ) then
          start(1) = ISGB
          start(2) = JSGB
          count(1) = IMAXB
          count(2) = JMAXB
       else
          ! for the case the shared-file contains no halos
          start(1) = PRC_2Drank(PRC_myrank,1) * IMAX + 1 ! no IHALO
          count(1) = IMAX
          if ( flag_half_x ) then
             if ( PRC_HAS_W ) then
                start(1) = start(1) + 1
             else
                count(1) = count(1) + 1
             endif
          endif

          start(2) = PRC_2Drank(PRC_myrank,2) * JMAX + 1 ! no JHALO
          count(2) = JMAX
          if ( flag_half_y ) then
             if ( PRC_HAS_S ) then
                start(2) = start(2) + 1
             else
                count(2) = count(2) + 1
             endif
          endif
       endif
    endif

    if ( atom ) then

       ! model coordinate (terrain following coordinate)
       call HistoryAddVariable( nvariant1,        & ! [OUT]
                                item,             & ! [IN]
                                dims(1:ndim),     & ! [IN]
                                desc,             & ! [IN]
                                unit,             & ! [IN]
                                mapping_name,     & ! [IN]
                                TIME_NOWSTEP,     & ! [IN]
                                zcoord = 'model', & ! [IN]
                                start  = start,   & ! [IN]
                                count  = count    ) ! [IN]

       ! absolute height coordinate
       dims(3) = 'z'
       call HistoryAddVariable( nvariant2,      & ! [OUT]
                                item,           & ! [IN]
                                dims(1:ndim),   & ! [IN]
                                desc,           & ! [IN]
                                unit,           & ! [IN]
                                mapping_name,   & ! [IN]
                                TIME_NOWSTEP,   & ! [IN]
                                zcoord = 'z',   & ! [IN]
                                start  = start, & ! [IN]
                                count  = count  ) ! [IN]

       ! pressure coordinate
       if ( HIST_PRES_nlayer > 0 ) then

          dims(3) = 'pressure'
          call HistoryAddVariable( nvariant3,           & ! [OUT]
                                   item,                & ! [IN]
                                   dims(1:ndim),        & ! [IN]
                                   desc,                & ! [IN]
                                   unit,                & ! [IN]
                                   mapping_name,        & ! [IN]
                                   TIME_NOWSTEP,        & ! [IN]
                                   zcoord = 'pressure', & ! [IN]
                                   start  = start,      & ! [IN]
                                   count  = count       ) ! [IN]

       else
          nvariant3 = 0
       endif

    else

       call HistoryAddVariable( nvariant1,     & ! [OUT]
                                item,          & ! [IN]
                                dims(1:ndim),  & ! [IN]
                                desc,          & ! [IN]
                                unit,          & ! [IN]
                                mapping_name,  & ! [IN]
                                TIME_NOWSTEP,  & ! [IN]
                                start = start, & ! [IN]
                                count = count  ) ! [IN]

       nvariant2 = 0
       nvariant3 = 0
    endif

    if ( nvariant1 + nvariant2 + nvariant3 > 0 ) then
       HIST_item_count   = HIST_item_count + 1
       itemid            = HIST_item_count
       HIST_item(itemid) = item

       do v = 1, nvariant1
          HIST_variant(itemid)                      = HIST_variant(itemid) + 1
          HIST_zcoord (itemid,HIST_variant(itemid)) = I_MODEL
       enddo

       do v = 1, nvariant2
          HIST_variant(itemid)                      = HIST_variant(itemid) + 1
          HIST_zcoord (itemid,HIST_variant(itemid)) = I_Z
       enddo

       do v = 1, nvariant3
          HIST_variant(itemid)                      = HIST_variant(itemid) + 1
          HIST_zcoord (itemid,HIST_variant(itemid)) = I_PRES
       enddo
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_reg

  !-----------------------------------------------------------------------------
  !> Check time to putting data
  subroutine HIST_query( &
       itemid, &
       answer  )
    use gtool_history, only: &
       HistoryQuery
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none

    integer,  intent(in)  :: itemid !< name of the item
    logical,  intent(out) :: answer !< is it time to store?

    integer  :: n, v, id
    !---------------------------------------------------------------------------

    answer = .false.

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    call HistoryQuery( HIST_item(itemid), & ! [IN]
                       TIME_NOWSTEP,      & ! [IN]
                       answer             ) ! [OUT]

    if ( answer ) then
       id = 0
       do n = 1, itemid-1
          id = id + HIST_variant(n)
       enddo

       do v = 1, HIST_variant(itemid)
          id = id + 1

          call HIST_checkfile( id ) ! [IN]
       enddo
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_query

  !-----------------------------------------------------------------------------
  !> Check time to switching output file
  subroutine HIST_checkfile( &
       itemid )
    use MPI, only: &
       MPI_COMM_NULL
    use gtool_history, only: &
       HistoryFileCreate
    use scale_process, only: &
       PRC_LOCAL_COMM_WORLD
    use scale_time, only: &
       TIME_NOWSTEP,      &
       TIME_gettimelabel
    implicit none

    integer,  intent(in)  :: itemid !< name of the item

    character(len=19) :: timelabel
    integer           :: comm
    logical           :: existed
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    call TIME_gettimelabel( timelabel )

    if ( IO_AGGREGATE ) then  ! user input parameter indicates to do PnetCDF I/O
       comm = PRC_LOCAL_COMM_WORLD
    else
       comm = MPI_COMM_NULL
    endif

    call HistoryFileCreate( itemid,         & ! [IN]
                            TIME_NOWSTEP,   & ! [IN]
                            timelabel,      & ! [IN]
                            comm=comm,      & ! [IN]
                            existed=existed ) ! [OUT]

    if ( .NOT. existed ) call HIST_set_axes_attributes

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_checkfile

  !-----------------------------------------------------------------------------
  !> Put 1D data to history buffer
  subroutine HIST_put_0D( &
       itemid, &
       var     )
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none

    integer,  intent(in)  :: itemid !< name of the item
    real(RP), intent(in)  :: var    !< value

    integer  :: n, v, id
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    id = 0
    do n = 1, itemid-1
       id = id + HIST_variant(n)
    enddo

    do v = 1, HIST_variant(itemid)
       id = id + 1

       call HistoryPut( id,           & ! [IN]
                        TIME_NOWSTEP, & ! [IN]
                        var           ) ! [IN]
    enddo

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_put_0D

  !-----------------------------------------------------------------------------
  !> Put 1D data to history buffer
  subroutine HIST_put_1D( &
       itemid, &
       var,    &
       zdim    )
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none

    integer,          intent(in) :: itemid !< name of the item
    real(RP),         intent(in) :: var(:) !< value

    character(len=*), intent(in), optional :: zdim

    character(len=H_SHORT) :: zd
    integer                :: ksize
    integer                :: kstart

    real(RP) :: var_trim(km)

    integer  :: n, v, id
    integer  :: k
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    zd = ''
    if( present(zdim) ) zd = zdim

    ! select dimension
    select case( zd )
      case('land')
        ksize  = LKMAX
        kstart = LKS
      case('landhalf')
        ksize  = LKMAX+1
        kstart = LKS-1
      case('urban')
        ksize  = UKMAX
        kstart = UKS
      case('urbanhalf')
        ksize  = UKMAX+1
        kstart = UKS-1
      case('half')
        ksize  = KMAX+1
        kstart = KS-1
      case default
        ksize  = KMAX
        kstart = KS
    end select

    do k = 1, ksize
       var_trim(k) = var(kstart+k-1)
    enddo

    id = 0
    do n = 1, itemid-1
       id = id + HIST_variant(n)
    enddo

    do v = 1, HIST_variant(itemid)
       id = id + 1

       call HistoryPut( id,               & ! [IN]
                        TIME_NOWSTEP,     & ! [IN]
                        var_trim(1:ksize) ) ! [IN]
    enddo

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_put_1D

  !-----------------------------------------------------------------------------
  !> Put 2D data to history buffer
  subroutine HIST_put_2D( &
       itemid, &
       var,    &
       xdim,   &
       ydim,   &
       nohalo  )
    use gtool_file, only: &
       RMISS
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none

    integer,          intent(in)  :: itemid   !< name of the item
    real(RP),         intent(in)  :: var(:,:) !< value

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    logical,          intent(in), optional :: nohalo

    character(len=H_SHORT) :: xd, yd
    integer                :: isize, jsize
    integer                :: istart, jstart

    real(RP) :: var_trim(imh*jmh)
    logical  :: nohalo_
    integer  :: s(2)

    integer  :: n, v, id
    integer  :: i, j
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    xd = ''
    yd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim

    nohalo_ = .false.
    if( present(nohalo) ) nohalo_ = nohalo

    ! select dimension
    select case( xd )
      case('half')
        isize  = imh
        istart = imsh
      case default
        isize  = im
        istart = ims
    end select

    select case( yd )
      case('half')
        jsize  = jmh
        jstart = jmsh
      case default
        jsize  = jm
        jstart = jms
    end select

    s(:) = shape(var)

    id = 0
    do n = 1, itemid-1
       id = id + HIST_variant(n)
    enddo

    do v = 1, HIST_variant(itemid)

       do j = 1, jsize
       do i = 1, isize
          var_trim((j-1)*isize+i) = var(istart+i-1,jstart+j-1)
       enddo
       enddo

       if ( nohalo_ ) then
          ! W halo
          do j = 1, jsize
          do i = 1, IS-istart
             var_trim((j-1)*isize+i) = RMISS
          enddo
          enddo
          ! E halo
          do j = 1, jsize
          do i = IE-istart+2, ime-istart+1
             var_trim((j-1)*isize+i) = RMISS
          enddo
          enddo
          ! S halo
          do j = 1, JS-jstart
          do i = 1, isize
             var_trim((j-1)*isize+i) = RMISS
          enddo
          enddo
          ! N halo
          do j = JE-jstart+2, jme-jstart+1
          do i = 1, isize
             var_trim((j-1)*isize+i) = RMISS
          enddo
          enddo
       endif

       id = id + 1

       call HistoryPut( id,                     & ! [IN]
                        TIME_NOWSTEP,           & ! [IN]
                        var_trim(1:isize*jsize) ) ! [IN]
    enddo

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_put_2D

  !-----------------------------------------------------------------------------
  !> Put 3D data to history buffer
  subroutine HIST_put_3D( &
       itemid, &
       var,    &
       xdim,   &
       ydim,   &
       zdim,   &
       nohalo  )
    use gtool_file, only: &
       RMISS
    use gtool_history, only: &
       HistoryPut
    use scale_time, only: &
       TIME_NOWSTEP
    use scale_interpolation, only: &
       INTERP_vertical_xi2z,   &
       INTERP_vertical_xi2p,   &
       INTERP_vertical_xih2zh, &
       INTERP_vertical_xih2p,  &
       INTERP_available
    implicit none

    integer,          intent(in)  :: itemid     !< name of the item
    real(RP),         intent(in)  :: var(:,:,:) !< value

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim
    logical,          intent(in), optional :: nohalo

    character(len=H_SHORT) :: xd, yd, zd
    integer                :: isize, jsize, ksize
    integer                :: istart, jstart, kstart

    real(RP) :: var_Z(KA              ,IA,JA)
    real(RP) :: var_P(HIST_PRES_nlayer,IA,JA)

    real(RP) :: var_trim(km*imh*jmh)
    logical  :: nohalo_
    integer  :: s(3)

    integer  :: n, v, id
    integer  :: i, j, k

    intrinsic shape
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    if( itemid < 0 ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    xd = ''
    yd = ''
    zd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim
    if( present(zdim) ) zd = zdim

    nohalo_ = .false.
    if( present(nohalo) ) nohalo_ = nohalo

    ! select dimension
    select case( xd )
      case('half')
        isize  = imh
        istart = imsh
      case default
        isize  = im
        istart = ims
    end select

    select case( yd )
      case('half')
        jsize  = jmh
        jstart = jmsh
      case default
        jsize  = jm
        jstart = jms
    end select

    select case( zd )
      case('land')
        ksize  = LKMAX
        kstart = LKS
      case('landhalf')
        ksize  = LKMAX+1
        kstart = LKS-1
      case('urban')
        ksize  = UKMAX
        kstart = UKS
      case('urbanhalf')
        ksize  = UKMAX+1
        kstart = UKS-1
      case('half')
        ksize  = KMAX+1
        kstart = KS-1
      case default
        ksize  = KMAX
        kstart = KS
    end select

    s(:) = shape(var)

    id = 0
    do n = 1, itemid-1
       id = id + HIST_variant(n)
    enddo

    do v = 1, HIST_variant(itemid)

       if    (       s(1)  == KA                  &
               .AND. ksize == KMAX                &
               .AND. HIST_zcoord(itemid,v) == I_Z &
               .AND. INTERP_available             ) then ! z*->z interpolation

          call PROF_rapstart('FILE_O_interp', 2)
          call INTERP_vertical_xi2z( var  (:,:,:), & ! [IN]
                                     var_Z(:,:,:)  ) ! [OUT]
          call PROF_rapend  ('FILE_O_interp', 2)

          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = var_Z(kstart+k-1,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo

       elseif(       s(1)  == KA                  &
               .AND. ksize == KMAX+1              &
               .AND. HIST_zcoord(itemid,v) == I_Z &
               .AND. INTERP_available             ) then ! z*->z interpolation

          call PROF_rapstart('FILE_O_interp', 2)
          call INTERP_vertical_xih2zh( var  (:,:,:), & ! [IN]
                                       var_Z(:,:,:)  ) ! [OUT]
          call PROF_rapend  ('FILE_O_interp', 2)

          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = var_Z(kstart+k-1,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo

       elseif(       s(1)  == KA                     &
               .AND. ksize == KMAX                   &
               .AND. HIST_zcoord(itemid,v) == I_PRES ) then ! z*->p interpolation

          ksize = HIST_PRES_nlayer

          call PROF_rapstart('FILE_O_interp', 2)
          call INTERP_vertical_xi2p( HIST_PRES_nlayer, & ! [IN]
                                     var  (:,:,:),     & ! [IN]
                                     var_P(:,:,:)      ) ! [OUT]
          call PROF_rapend  ('FILE_O_interp', 2)

          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = var_P(k,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo

       elseif(       s(1)  == KA                     &
               .AND. ksize == KMAX+1                 &
               .AND. HIST_zcoord(itemid,v) == I_PRES ) then ! z*->p interpolation

          ksize = HIST_PRES_nlayer

          call PROF_rapstart('FILE_O_interp', 2)
          call INTERP_vertical_xih2p( HIST_PRES_nlayer, & ! [IN]
                                      var  (:,:,:),     & ! [IN]
                                      var_P(:,:,:)      ) ! [OUT]
          call PROF_rapend  ('FILE_O_interp', 2)

          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = var_P(k,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo

       else ! no interpolation

          do k = 1, ksize
          do j = 1, jsize
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = var(kstart+k-1,istart+i-1,jstart+j-1)
          enddo
          enddo
          enddo

       endif

       if ( nohalo_ ) then
          ! W halo
          do k = 1, ksize
          do j = 1, jsize
          do i = 1, IS-istart
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
          enddo
          enddo
          enddo
          ! E halo
          do k = 1, ksize
          do j = 1, jsize
          do i = IE-istart+2, ime-istart+1
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
          enddo
          enddo
          enddo
          ! S halo
          do k = 1, ksize
          do j = 1, JS-jstart
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
          enddo
          enddo
          enddo
          ! N halo
          do k = 1, ksize
          do j = JE-jstart+2, jme-jstart+1
          do i = 1, isize
             var_trim((k-1)*jsize*isize+(j-1)*isize+i) = RMISS
          enddo
          enddo
          enddo
       endif

       id = id + 1

       call HistoryPut( id,                           & ! [IN]
                        TIME_NOWSTEP,                 & ! [IN]
                        var_trim(1:isize*jsize*ksize) ) ! [IN]
    enddo

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_put_3D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 0D
  subroutine HIST_in_0D( &
       var,  &
       item, &
       desc, &
       unit )
    implicit none

    real(RP),         intent(in) :: var  !< value
    character(len=*), intent(in) :: item !< name        of the item
    character(len=*), intent(in) :: desc !< description of the item
    character(len=*), intent(in) :: unit !< unit        of the item

    integer, parameter :: ndim = 0
    integer :: itemid
    logical :: do_put
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    ! Check whether the item has been already registered
    call HIST_reg  ( itemid,  & ! [OUT]
                     item,    & ! [IN]
                     desc,    & ! [IN]
                     unit,    & ! [IN]
                     ndim     ) ! [IN]

    ! Check whether it is time to input the item
    call HIST_query( itemid,  & ! [IN]
                     do_put   ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid, & ! [IN]
                      var     ) ! [IN]
    endif

    return
  end subroutine HIST_in_0D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 1D
  subroutine HIST_in_1D( &
       var,  &
       item, &
       desc, &
       unit, &
       zdim  )
    implicit none

    real(RP),         intent(in)  :: var(:) !< value
    character(len=*), intent(in)  :: item   !< name        of the item
    character(len=*), intent(in)  :: desc   !< description of the item
    character(len=*), intent(in)  :: unit   !< unit        of the item

    character(len=*), intent(in), optional :: zdim

    character(len=H_SHORT) :: zd

    integer, parameter :: ndim = 1
    integer :: itemid
    logical :: do_put
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    zd = ''
    if( present(zdim) ) zd = zdim

    ! Check whether the item has been already registered
    call HIST_reg  ( itemid,  & ! [OUT]
                     item,    & ! [IN]
                     desc,    & ! [IN]
                     unit,    & ! [IN]
                     ndim,    & ! [IN]
                     zdim=zd  ) ! [IN]

    ! Check whether it is time to input the item
    call HIST_query( itemid,  & ! [IN]
                     do_put   ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid, & ! [IN]
                      var     ) ! [IN]
    endif

    return
  end subroutine HIST_in_1D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 2D
  subroutine HIST_in_2D( &
       var,   &
       item,  &
       desc,  &
       unit,  &
       xdim,  &
       ydim,  &
       nohalo )
    implicit none

    real(RP),         intent(in)  :: var(:,:) !< value
    character(len=*), intent(in)  :: item     !< name        of the item
    character(len=*), intent(in)  :: desc     !< description of the item
    character(len=*), intent(in)  :: unit     !< unit        of the item

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    logical,          intent(in), optional :: nohalo

    character(len=H_SHORT) :: xd, yd

    integer, parameter :: ndim = 2
    integer :: itemid
    logical :: do_put
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    xd = ''
    yd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim

    ! Check whether the item has been already registered
    call HIST_reg  ( itemid,  & ! [OUT]
                     item,    & ! [IN]
                     desc,    & ! [IN]
                     unit,    & ! [IN]
                     ndim,    & ! [IN]
                     xdim=xd, & ! [IN]
                     ydim=yd  ) ! [IN]

    ! Check whether it is time to input the item
    call HIST_query( itemid,  & ! [IN]
                     do_put   ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid,       & ! [IN]
                      var,          & ! [IN]
                      nohalo=nohalo ) ! [IN]
    endif

    return
  end subroutine HIST_in_2D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of HIST_reg+HIST_put 3D
  subroutine HIST_in_3D( &
       var,   &
       item,  &
       desc,  &
       unit,  &
       xdim,  &
       ydim,  &
       zdim,  &
       nohalo )
    implicit none

    real(RP),         intent(in)  :: var(:,:,:) !< value
    character(len=*), intent(in)  :: item       !< name        of the item
    character(len=*), intent(in)  :: desc       !< description of the item
    character(len=*), intent(in)  :: unit       !< unit        of the item

    character(len=*), intent(in), optional :: xdim
    character(len=*), intent(in), optional :: ydim
    character(len=*), intent(in), optional :: zdim
    logical,          intent(in), optional :: nohalo

    character(len=H_SHORT) :: xd, yd, zd

    integer, parameter :: ndim = 3
    integer :: itemid
    logical :: do_put
    !---------------------------------------------------------------------------

    if( .NOT. enabled ) return

    xd = ''
    yd = ''
    zd = ''
    if( present(xdim) ) xd = xdim
    if( present(ydim) ) yd = ydim
    if( present(zdim) ) zd = zdim

    ! Check whether the item has been already registered
    call HIST_reg  ( itemid,  & ! [OUT]
                     item,    & ! [IN]
                     desc,    & ! [IN]
                     unit,    & ! [IN]
                     ndim,    & ! [IN]
                     xdim=xd, & ! [IN]
                     ydim=yd, & ! [IN]
                     zdim=zd  ) ! [IN]

    ! Check whether it is time to input the item
    call HIST_query( itemid,  & ! [IN]
                     do_put   ) ! [OUT]

    if ( do_put ) then
       call HIST_put( itemid,       & ! [IN]
                      var,          & ! [IN]
                      xdim=xd,      & ! [IN]
                      ydim=yd,      & ! [IN]
                      zdim=zd,      & ! [IN]
                      nohalo=nohalo ) ! [IN]
    endif

    return
  end subroutine HIST_in_3D

  !-----------------------------------------------------------------------------
  !> Get 1D data from file
  subroutine HIST_get_1D( &
       var,          &
       basename,     &
       varname,      &
       step,         &
       allow_missing )
    use gtool_history, only: &
       HistoryGet
    implicit none

    real(RP),         intent(out) :: var(:)     !< value
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    integer,          intent(in)  :: step       !< step number

    logical,          intent(in), optional :: allow_missing !< allow data is missing?

    logical :: am
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:),          & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine HIST_get_1D

  !-----------------------------------------------------------------------------
  !> Get 2D data from file
  subroutine HIST_get_2D( &
       var,          &
       basename,     &
       varname,      &
       step,         &
       allow_missing )
    use gtool_history, only: &
       HistoryGet
    implicit none

    real(RP),         intent(out) :: var(:,:)   !< value
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    integer,          intent(in)  :: step       !< step number

    logical,          intent(in), optional :: allow_missing !< allow data is missing?

    logical :: am
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:,:),        & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine HIST_get_2D

  !-----------------------------------------------------------------------------
  !> Get 3D data from file
  subroutine HIST_get_3D( &
       var,          &
       basename,     &
       varname,      &
       step,         &
       allow_missing )
    use gtool_history, only: &
       HistoryGet
    implicit none

    real(RP),         intent(out) :: var(:,:,:) !< value
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    integer,          intent(in)  :: step       !< step number

    logical,          intent(in), optional :: allow_missing !< allow data is missing?

    logical :: am
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    am = .false.
    if( present(allow_missing) ) am = allow_missing

    call HistoryGet( var(:,:,:),      & ! [OUT]
                     basename,        & ! [IN]
                     varname,         & ! [IN]
                     step,            & ! [IN]
                     allow_missing=am ) ! [IN]

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine HIST_get_3D

  !-----------------------------------------------------------------------------
  !> Flush history buffer to file
  subroutine HIST_write
    use gtool_history, only: &
       HistoryWriteAll, &
       HistoryWriteAxes
    use scale_time, only: &
       TIME_NOWSTEP
    implicit none

    logical :: axis_written_first
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    ! Note this subroutine must be called after all HIST_reg calls are completed
    ! Write registered history axes to history file
    call HistoryWriteAxes( axis_written_first )

    call HistoryWriteAll( TIME_NOWSTEP ) ![IN]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine HIST_write

  !-----------------------------------------------------------------------------
  !> Put axis coordinate to history file
  !  only register the axis and coordinate variables into internal buffers
  !  The actual write happens later when calling HIST_write
  subroutine HIST_put_axes
    use gtool_history, only: &
       HistoryPutAxis,  &
       HistoryPutAssociatedCoordinates
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

    if ( IO_AGGREGATE ) then

       startZ = 1

       if ( HIST_BND ) then
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
             XAGH  = XAGH + 1
             if( PRC_HAS_W ) startXH = startXH + 1
          endif

          if ( .NOT. PRC_PERIODIC_Y ) then
             YAGH  = YAGH + 1
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
    call HistoryPutAxis( 'z',   'Z',               'm', 'z',   GRID_CZ(KS  :KE), gsize=KMAX  , start=startZ )
    call HistoryPutAxis( 'zh',  'Z (half level)',  'm', 'zh',  GRID_FZ(KS-1:KE), gsize=KMAX+1, start=startZ )

    call HistoryPutAxis( 'lz',  'LZ',              'm', 'lz',  GRID_LCZ(LKS  :LKE), down=.true., gsize=LKMAX  , start=startZ )
    call HistoryPutAxis( 'lzh', 'LZ (half level)', 'm', 'lzh', GRID_LFZ(LKS-1:LKE), down=.true., gsize=LKMAX+1, start=startZ )

    call HistoryPutAxis( 'uz',  'UZ',              'm', 'uz',  GRID_UCZ(UKS  :UKE), down=.true., gsize=UKMAX  , start=startZ )
    call HistoryPutAxis( 'uzh', 'UZ (half level)', 'm', 'uzh', GRID_UFZ(UKS-1:UKE), down=.true., gsize=UKMAX+1, start=startZ )

    call HistoryPutAxis( 'x',   'X',               'm', 'x',   GRID_CX(ims :ime), gsize=XAG , start=startX  )
    call HistoryPutAxis( 'xh',  'X (half level)',  'm', 'xh',  GRID_FX(imsh:ime), gsize=XAGH, start=startXH )

    call HistoryPutAxis( 'y',   'Y',               'm', 'y',   GRID_CY(jms :jme), gsize=YAG , start=startY  )
    call HistoryPutAxis( 'yh',  'Y (half level)',  'm', 'yh',  GRID_FY(jmsh:jme), gsize=YAGH, start=startYH )

    if ( HIST_PRES_nlayer > 0 ) then
       call HistoryPutAxis( 'pressure', 'Pressure', 'hPa', 'pressure',         &
                            HIST_PRES_val(:)/100.0_RP, down=.true., gsize=HIST_PRES_nlayer, start=startZ )
    endif

    ! axes below always include halos when written to file regardless of PRC_PERIODIC_X/PRC_PERIODIC_Y
    call HistoryPutAxis( 'CZ',   'Atmos Grid Center Position Z', 'm', 'CZ',  GRID_CZ,   gsize=KA,   start=startZ )
    call HistoryPutAxis( 'FZ',   'Atmos Grid Face Position Z',   'm', 'FZ',  GRID_FZ,   gsize=KA+1, start=startZ )
    call HistoryPutAxis( 'CDZ',  'Grid Cell length Z',           'm', 'CZ',  GRID_CDZ,  gsize=KA,   start=startZ )
    call HistoryPutAxis( 'FDZ',  'Grid distance Z',              'm', 'FDZ', GRID_FDZ,  gsize=KA-1, start=startZ )
    call HistoryPutAxis( 'CBFZ', 'Boundary factor Center Z',     '1', 'CZ',  GRID_CBFZ, gsize=KA,   start=startZ )
    call HistoryPutAxis( 'FBFZ', 'Boundary factor Face Z',       '1', 'FZ',  GRID_FBFZ, gsize=KA+1, start=startZ )

    call HistoryPutAxis( 'LCZ',  'Land Grid Center Position Z',  'm', 'LCZ', GRID_LCZ,  down=.true., gsize=LKMAX,   start=startZ )
    call HistoryPutAxis( 'LFZ',  'Land Grid Face Position Z',    'm', 'LFZ', GRID_LFZ,  down=.true., gsize=LKMAX+1, start=startZ )
    call HistoryPutAxis( 'LCDZ', 'Land Grid Cell length Z',      'm', 'LCZ', GRID_LCDZ,              gsize=LKMAX,   start=startZ )

    call HistoryPutAxis( 'UCZ',  'Urban Grid Center Position Z', 'm', 'UCZ', GRID_UCZ,  down=.true., gsize=UKMAX,   start=startZ )
    call HistoryPutAxis( 'UFZ',  'Urban Grid Face Position Z',   'm', 'UFZ', GRID_UFZ,  down=.true., gsize=UKMAX+1, start=startZ )
    call HistoryPutAxis( 'UCDZ', 'Urban Grid Cell length Z',     'm', 'UCZ', GRID_UCDZ,              gsize=UKMAX,   start=startZ )

    if ( IO_AGGREGATE ) then
       call HistoryPutAxis( 'CX',   'Atmos Grid Center Position X', 'm', 'CX',  GRID_CXG,   gsize=IAG,   start=startZ )
       call HistoryPutAxis( 'CY',   'Atmos Grid Center Position Y', 'm', 'CY',  GRID_CYG,   gsize=JAG,   start=startZ )
       call HistoryPutAxis( 'FX',   'Atmos Grid Face Position X',   'm', 'FX',  GRID_FXG,   gsize=IAG+1, start=startZ )
       call HistoryPutAxis( 'FY',   'Atmos Grid Face Position Y',   'm', 'FY',  GRID_FYG,   gsize=JAG+1, start=startZ )
       call HistoryPutAxis( 'CDX',  'Grid Cell length X',           'm', 'CX',  GRID_CDXG,  gsize=IAG,   start=startZ )
       call HistoryPutAxis( 'CDY',  'Grid Cell length Y',           'm', 'CY',  GRID_CDYG,  gsize=JAG,   start=startZ )
       call HistoryPutAxis( 'FDX',  'Grid distance X',              'm', 'FDX', GRID_FDXG,  gsize=IAG-1, start=startZ )
       call HistoryPutAxis( 'FDY',  'Grid distance Y',              'm', 'FDY', GRID_FDYG,  gsize=JAG-1, start=startZ )
       call HistoryPutAxis( 'CBFX', 'Boundary factor Center X',     '1', 'CX',  GRID_CBFXG, gsize=IAG,   start=startZ )
       call HistoryPutAxis( 'CBFY', 'Boundary factor Center Y',     '1', 'CY',  GRID_CBFYG, gsize=JAG,   start=startZ )
       call HistoryPutAxis( 'FBFX', 'Boundary factor Face X',       '1', 'FX',  GRID_FBFXG, gsize=IAG+1, start=startZ )
       call HistoryPutAxis( 'FBFY', 'Boundary factor Face Y',       '1', 'FY',  GRID_FBFYG, gsize=JAG+1, start=startZ )
    else
       call HistoryPutAxis( 'CX',   'Atmos Grid Center Position X', 'm', 'CX',  GRID_CX   )
       call HistoryPutAxis( 'CY',   'Atmos Grid Center Position Y', 'm', 'CY',  GRID_CY   )
       call HistoryPutAxis( 'FX',   'Atmos Grid Face Position X',   'm', 'FX',  GRID_FX   )
       call HistoryPutAxis( 'FY',   'Atmos Grid Face Position Y',   'm', 'FY',  GRID_FY   )
       call HistoryPutAxis( 'CDX',  'Grid Cell length X',           'm', 'CX',  GRID_CDX  )
       call HistoryPutAxis( 'CDY',  'Grid Cell length Y',           'm', 'CY',  GRID_CDY  )
       call HistoryPutAxis( 'FDX',  'Grid distance X',              'm', 'FDX', GRID_FDX  )
       call HistoryPutAxis( 'FDY',  'Grid distance Y',              'm', 'FDY', GRID_FDY  )
       call HistoryPutAxis( 'CBFX', 'Boundary factor Center X',     '1', 'CX',  GRID_CBFX )
       call HistoryPutAxis( 'CBFY', 'Boundary factor Center Y',     '1', 'CY',  GRID_CBFY )
       call HistoryPutAxis( 'FBFX', 'Boundary factor Face X',       '1', 'FX',  GRID_FBFX )
       call HistoryPutAxis( 'FBFY', 'Boundary factor Face Y',       '1', 'FY',  GRID_FBFY )
    endif

    ! TODO: skip 8 axes below when IO_AGGREGATE is true, as all axes are now global
    call HistoryPutAxis('CXG',   'Grid Center Position X (global)',   'm', 'CXG',  GRID_CXG,   gsize=IAG,   start=startZ )
    call HistoryPutAxis('CYG',   'Grid Center Position Y (global)',   'm', 'CYG',  GRID_CYG,   gsize=JAG,   start=startZ )
    call HistoryPutAxis('FXG',   'Grid Face Position X (global)',     'm', 'FXG',  GRID_FXG,   gsize=IAG+1, start=startZ )
    call HistoryPutAxis('FYG',   'Grid Face Position Y (global)',     'm', 'FYG',  GRID_FYG,   gsize=JAG+1, start=startZ )
    call HistoryPutAxis('CBFXG', 'Boundary factor Center X (global)', '1', 'CXG',  GRID_CBFXG, gsize=IAG,   start=startZ )
    call HistoryPutAxis('CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG',  GRID_CBFYG, gsize=JAG,   start=startZ )
    call HistoryPutAxis('FBFXG', 'Boundary factor Face X (global)',   '1', 'FXG',  GRID_FBFXG, gsize=IAG+1, start=startZ )
    call HistoryPutAxis('FBFYG', 'Boundary factor Face Y (global)',   '1', 'FYG',  GRID_FBFYG, gsize=JAG+1, start=startZ )

    ! associate coordinates
    if ( IO_AGGREGATE ) then
       if ( HIST_BND ) then
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
    call HistoryPutAssociatedCoordinates( 'height', 'height above ground level',                        &
                                          'm', AXIS_name(1:3), AXIS(1:im,1:jm,1:KMAX), start=start(:,1) )

    do k = 0, KMAX
       AXIS(1:im,1:jm,k) = REAL_FZ(k+KS-1,ims:ime,jms:jme)
    enddo
    AXIS_name(1:3) = (/'x ', 'y ', 'zh'/)
    call HistoryPutAssociatedCoordinates( 'height_xyw', 'height above ground level (half level xyw)',    &
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
    call HistoryPutAssociatedCoordinates( 'height_uyz', 'height above ground level (half level uyz)',    &
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
    call HistoryPutAssociatedCoordinates( 'height_xvz', 'height above ground level (half level xvz)',    &
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
    call HistoryPutAssociatedCoordinates( 'height_uvz', 'height above ground level (half level uvz)',     &
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
    call HistoryPutAssociatedCoordinates( 'height_uyw', 'height above ground level (half level uyw)',    &
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
    call HistoryPutAssociatedCoordinates( 'height_xvw', 'height above ground level (half level xvw)',    &
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
    call HistoryPutAssociatedCoordinates( 'height_uvw', 'height above ground level (half level uvw)',     &
                                          'm', AXIS_name(1:3), AXIS(1:imh,1:jmh,0:KMAX), start=start(:,4) )

    AXIS(1:im,1:jm,1) = REAL_LON (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lon', 'longitude',                                                 &
                                          'degrees_east', AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1) )

    AXIS(1:imh,1:jm,1) = REAL_LONX(imsh:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lon_uy', 'longitude (half level uy)',                               &
                                          'degrees_east', AXIS_name(1:2), AXIS(1:imh,1:jm,1), start=start(:,2) )

    AXIS(1:im,1:jmh,1) = REAL_LONY(ims:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'yh'/)
    call HistoryPutAssociatedCoordinates( 'lon_xv', 'longitude (half level xv)',                               &
                                          'degrees_east', AXIS_name(1:2), AXIS(1:im,1:jmh,1), start=start(:,3) )

    AXIS(1:imh,1:jmh,1) = REAL_LONXY(imsh:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'yh'/)
    call HistoryPutAssociatedCoordinates( 'lon_uv', 'longitude (half level uv)',                                &
                                          'degrees_east', AXIS_name(1:2), AXIS(1:imh,1:jmh,1), start=start(:,4) )

    AXIS(1:im,1:jm,1) = REAL_LAT (ims:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lat', 'latitude',                                                   &
                                          'degrees_north', AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1) )

    AXIS(1:imh,1:jm,1) = REAL_LATX(imsh:ime,jms:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lat_uy', 'latitude (half level uy)',                                 &
                                          'degrees_north', AXIS_name(1:2), AXIS(1:imh,1:jm,1), start=start(:,2) )

    AXIS(1:im,1:jmh,1) = REAL_LATY(ims:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'x ', 'yh'/)
    call HistoryPutAssociatedCoordinates( 'lat_xv', 'latitude (half level xv)',                                 &
                                          'degrees_north', AXIS_name(1:2), AXIS(1:im,1:jmh,1), start=start(:,3) )

    AXIS(1:imh,1:jmh,1) = REAL_LATXY(imsh:ime,jmsh:jme) / D2R
    AXIS_name(1:2) = (/'xh', 'yh'/)
    call HistoryPutAssociatedCoordinates( 'lat_uv', 'latitude (half level uv)',                                  &
                                          'degrees_north', AXIS_name(1:2), AXIS(1:imh,1:jmh,1), start=start(:,4) )

    AXIS(1:im,1:jm,1) = TOPO_Zsfc(ims:ime,jms:jme)
    AXIS_name(1:2) = (/'x ', 'y '/)
    call HistoryPutAssociatedCoordinates( 'topo', 'topography',                                    &
                                          'm', AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1) )

    AXIS(1:im,1:jm,1) = LANDUSE_frac_land(ims:ime,jms:jme)
    AXIS_name(1:2) = (/'x ', 'y '/)
    call HistoryPutAssociatedCoordinates( 'lsmask', 'fraction for land-sea mask',                  &
                                          '1', AXIS_name(1:2), AXIS(1:im,1:jm,1), start=start(:,1) )

    return
  end subroutine HIST_put_axes

  !-----------------------------------------------------------------------------
  subroutine HIST_set_axes_attributes
    use gtool_history, only: &
       HistorySetGlobalAttribute, &
       HistorySetAttribute,       &
       HistorySetMapping
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
    use scale_mapproj, only: &
       MPRJ_get_attributes
    implicit none

    character(len=5) :: periodic_z, periodic_x, periodic_y

    integer :: isize, jsize
    integer :: istart, jstart
    integer :: whalo_g, ehalo_g, shalo_g, nhalo_g
    integer :: whalo_l, ehalo_l, shalo_l, nhalo_l
    integer :: isizeh, jsizeh
    integer :: istarth, jstarth
    integer :: whaloh_g, shaloh_g
    integer :: whaloh_l, shaloh_l

    character(len=5) :: logical_str

    character(len=H_SHORT) :: mapping
    real(DP) :: false_easting
    real(DP) :: false_northing
    real(DP) :: longitude_of_central_meridian
    real(DP) :: longitude_of_projection_origin
    real(DP) :: latitude_of_projection_origin
    real(DP) :: straight_vertical_longitude_from_pole
    real(DP) :: standard_parallel(2)
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

    if ( .NOT. IO_AGGREGATE ) then
       call HistorySetGlobalAttribute( "scale_rm_prc_rank_x", (/PRC_2Drank(PRC_myrank,1)/) ) ! [IN]
       call HistorySetGlobalAttribute( "scale_rm_prc_rank_y", (/PRC_2Drank(PRC_myrank,2)/) ) ! [IN]

       call HistorySetGlobalAttribute( "scale_rm_prc_num_x", (/PRC_NUM_X/) ) ! [IN]
       call HistorySetGlobalAttribute( "scale_rm_prc_num_y", (/PRC_NUM_Y/) ) ! [IN]
    endif

    call HistorySetGlobalAttribute( "scale_rm_prc_periodic_z", periodic_z ) ! [IN]
    call HistorySetGlobalAttribute( "scale_rm_prc_periodic_x", periodic_x ) ! [IN]
    call HistorySetGlobalAttribute( "scale_rm_prc_periodic_y", periodic_y ) ! [IN]

    call HistorySetGlobalAttribute( "scale_rm_grid_index_kmax",  (/KMAX/)  ) ! [IN]
    call HistorySetGlobalAttribute( "scale_rm_grid_index_imaxg", (/IMAXG/) ) ! [IN]
    call HistorySetGlobalAttribute( "scale_rm_grid_index_jmaxg", (/JMAXG/) ) ! [IN]

    call HistorySetGlobalAttribute( "scale_rm_grid_index_khalo", (/KHALO/) ) ! [IN]
    call HistorySetGlobalAttribute( "scale_rm_grid_index_ihalo", (/IHALO/) ) ! [IN]
    call HistorySetGlobalAttribute( "scale_rm_grid_index_jhalo", (/JHALO/) ) ! [IN]

    ! attributes
    if ( PRC_PERIODIC_X ) then; logical_str = "true"; else; logical_str = "false"; end if
    if( PRC_PERIODIC_X .or. .not. HIST_BND ) then
       isize = IMAX * PRC_NUM_X
       istart = IS_inG - IHALO
       whalo_g = 0
       ehalo_g = 0
       whalo_l = 0
       ehalo_l = 0
       isizeh = isize
       istarth = istart
       whaloh_g = whalo_g
       whaloh_l = whalo_l
       if ( .not. PRC_PERIODIC_X ) then
          isizeh = isize + 1
          whaloh_g = whalo_g + 1
          if ( PRC_HAS_W ) then
             istarth = istarth + 1
          else
             whaloh_l = whalo_l + 1
          end if
       end if
    else ! HIST_BND == .true.
       isize = IAG
       istart = ISGA
       whalo_g = IHALO
       ehalo_g = IHALO
       if ( IO_AGGREGATE ) then
          whalo_l = whalo_g
          ehalo_l = ehalo_g
       else
          if ( PRC_HAS_W ) then; whalo_l = 0; else; whalo_l = whalo_g; end if
          if ( PRC_HAS_E ) then; ehalo_l = 0; else; ehalo_l = ehalo_g; end if
       end if
       isizeh = isize
       istarth = istart
    end if

    call HistorySetAttribute( "x", "size_global",  (/ isize /) )
    call HistorySetAttribute( "x", "start_global", (/ istart /) )
    call HistorySetAttribute( "x", "halo_global",  (/ whalo_g, ehalo_g /) )
    call HistorySetAttribute( "x", "halo_local",   (/ whalo_l, ehalo_l /) )
    call HistorySetAttribute( "x", "periodic",     logical_str )

    call HistorySetAttribute( "xh", "size_global",  (/ isizeh /) )
    call HistorySetAttribute( "xh", "start_global", (/ istarth /) )
    call HistorySetAttribute( "xh", "halo_global",  (/ whaloh_g, ehalo_g /) )
    call HistorySetAttribute( "xh", "halo_local",   (/ whaloh_l, ehalo_l /) )
    call HistorySetAttribute( "xh", "periodic",     logical_str )


    if ( PRC_PERIODIC_Y ) then; logical_str = "true"; else; logical_str = "false"; end if
    if( PRC_PERIODIC_Y .or. .not. HIST_BND ) then
       jsize = JMAX * PRC_NUM_Y
       jstart = JS_inG - JHALO
       shalo_g = 0
       nhalo_g = 0
       shalo_l = 0
       nhalo_l = 0
       jsizeh = jsize
       jstarth = jstart
       shaloh_g = shalo_g
       shaloh_l = shalo_l
       if ( .not. PRC_PERIODIC_Y ) then
          jsizeh = jsize + 1
          shaloh_g = shalo_g + 1
          if ( PRC_HAS_S ) then
             jstarth = jstarth + 1
          else
             shaloh_l = shalo_l + 1
          end if
       end if
    else ! HIST_BND == .true.
       jsize = JAG
       jstart = JSGA
       shalo_g = JHALO
       nhalo_g = JHALO
       if ( IO_AGGREGATE ) then
          shalo_l = shalo_g
          nhalo_l = nhalo_g
       else
          if ( PRC_HAS_S ) then; shalo_l = 0; else; shalo_l = shalo_g; end if
          if ( PRC_HAS_N ) then; nhalo_l = 0; else; nhalo_l = nhalo_g; end if
       end if
       jsizeh = jsize
       jstarth = jstart
    end if

    call HistorySetAttribute( "y", "size_global",  (/ jsize /) )
    call HistorySetAttribute( "y", "start_global", (/ jstart /) )
    call HistorySetAttribute( "y", "halo_global",  (/ shalo_g, nhalo_g /) )
    call HistorySetAttribute( "y", "halo_local",   (/ shalo_l, nhalo_l /) )
    call HistorySetAttribute( "y", "periodic",     logical_str )

    call HistorySetAttribute( "yh", "size_global",  (/ jsizeh /) )
    call HistorySetAttribute( "yh", "start_global", (/ jstarth /) )
    call HistorySetAttribute( "yh", "halo_global",  (/ shaloh_g, nhalo_g /) )
    call HistorySetAttribute( "yh", "halo_local",   (/ shaloh_l, nhalo_l /) )
    call HistorySetAttribute( "yh", "periodic",     logical_str )

    call MPRJ_get_attributes( mapping,                               & ! [OUT]
                              false_easting, false_northing,         & ! [OUT]
                              longitude_of_central_meridian,         & ! [OUT]
                              longitude_of_projection_origin,        & ! [OUT]
                              latitude_of_projection_origin,         & ! [OUT]
                              straight_vertical_longitude_from_pole, & ! [OUT]
                              standard_parallel(:)                   ) ! [OUT]

    if ( mapping /= "" ) then
       call HistorySetAttribute( "x",  "standard_name", "projection_x_coordinate");
       call HistorySetAttribute( "xh", "standard_name", "projection_x_coordinate");
       call HistorySetAttribute( "y",  "standard_name", "projection_y_coordinate");
       call HistorySetAttribute( "yh", "standard_name", "projection_y_coordinate");
       call HistorySetMapping( mapping )

       if ( false_easting /= UNDEF ) then
          call HistorySetAttribute( mapping, "false_easting",  (/ false_easting /) )
       endif

       if ( false_northing /= UNDEF ) then
          call HistorySetAttribute( mapping, "false_northing", (/ false_northing /) )
       endif

       if ( longitude_of_central_meridian /= UNDEF ) then
          call HistorySetAttribute( mapping, "longitude_of_central_meridian", (/ longitude_of_central_meridian /) )
       endif

       if ( longitude_of_projection_origin /= UNDEF ) then
          call HistorySetAttribute( mapping, "longitude_of_projection_origin", (/ longitude_of_projection_origin /) )
       endif

       if ( latitude_of_projection_origin /= UNDEF ) then
          call HistorySetAttribute( mapping, "latitude_of_projection_origin", (/ latitude_of_projection_origin /) )
       endif

       if ( straight_vertical_longitude_from_pole /= UNDEF ) then
          call HistorySetAttribute( mapping, "straight_vertical_longitude_from_pole", (/ straight_vertical_longitude_from_pole /) )
       endif

       if ( standard_parallel(1) /= UNDEF ) then
          if ( standard_parallel(2) /= UNDEF ) then
             call HistorySetAttribute( mapping, "standard_parallel", standard_parallel(1:2) )
          else
             call HistorySetAttribute( mapping, "standard_parallel", standard_parallel(1:1) )
          endif
       endif
    endif

    return
  end subroutine HIST_set_axes_attributes

end module scale_history
