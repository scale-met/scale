!-------------------------------------------------------------------------------
!> module file_tiledata
!!
!! @par Description
!!          tile data reader
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_file_tiledata
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public FILE_TILEDATA_get_info
  public FILE_TILEDATA_get_data

  interface FILE_TILEDATA_get_data
     module procedure FILE_TILEDATA_get_data_real
     module procedure FILE_TILEDATA_get_data_int1
  end interface FILE_TILEDATA_get_data

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
  !> get tile information
  subroutine FILE_TILEDATA_get_info( &
       TILE_nlim,                                          &
       TILE_DLAT, TILE_DLON,                               &
       DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, &
       catalog_fname,                                      &
       GLOBAL_IA,                                          &
       TILE_nmax,                                          &
       TILE_fname, TILE_hit,                               &
       TILE_JS, TILE_JE, TILE_IS, TILE_IE,                 &
       nLATH, nLONH, jsh, jeh, ish, ieh, zonal, pole,      &
       single_fname, LATS, LATE, LONS, LONE                )
    use scale_const, only: &
       PI => CONST_PI
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: TILE_nlim
    real(RP),         intent(in)  :: TILE_DLAT, TILE_DLON
    real(RP),         intent(in)  :: DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE
    character(len=*), intent(in)  :: catalog_fname
    integer,          intent(out) :: GLOBAL_IA
    integer,          intent(out) :: TILE_nmax
    logical,          intent(out) :: TILE_hit(:)
    character(len=*), intent(out) :: TILE_fname(:)
    integer,          intent(out) :: TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:)
    integer,          intent(out) :: nLATH, nLONH
    integer,          intent(out) :: jsh, jeh, ish, ieh
    logical,          intent(out) :: zonal, pole

    character(len=*), intent(in), optional :: single_fname
    real(RP),         intent(in), optional :: LATS
    real(RP),         intent(in), optional :: LATE
    real(RP),         intent(in), optional :: LONS
    real(RP),         intent(in), optional :: LONE

    real(RP) :: TILE_LATS(TILE_nlim), TILE_LATE(TILE_nlim)
    real(RP) :: TILE_LONS(TILE_nlim), TILE_LONE(TILE_nlim)
    real(RP) :: LAT_MIN, LAT_MAX

    integer :: DOMAIN_JS, DOMAIN_JE, DOMAIN_IS, DOMAIN_IE

    call FILE_TILEDATA_get_domain_info( TILE_DLAT, TILE_DLON,                               & ! [IN]
                                        DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, & ! [IN]
                                        DOMAIN_JS,   DOMAIN_JE,   DOMAIN_IS,   DOMAIN_IE,   & ! [OUT]
                                        GLOBAL_IA                                           ) ! [OUT]


    if ( catalog_fname /= "" ) then
       LOG_INFO("FILE_TILEDATA_get_info",*) 'Input catalogue file:', trim(catalog_fname)

       call FILE_TILEDATA_read_catalog_file( TILE_nlim,                  & ! [IN]
                                             catalog_fname,              & ! [IN]
                                             TILE_DLAT, TILE_DLON,       & ! [IN]
                                             DOMAIN_IS, GLOBAL_IA,       & ! [IN]
                                             TILE_nmax,                  & ! [OUT]
                                             TILE_fname(:),              & ! [OUT]
                                             TILE_LATS(:), TILE_LATE(:), & ! [OUT]
                                             TILE_LONS(:), TILE_LONE(:)  ) ! [OUT]
       LAT_MIN = minval( TILE_LATS(1:TILE_nmax) )
       LAT_MAX = maxval( TILE_LATE(1:TILE_nmax) )

    else
       if ( .not. present(single_fname) ) then
          LOG_ERROR("FILE_TILEDATA_get_info",*) "single_fname is required if catalog_fname is empty"
          call PRC_abort
       end if
       if ( .not. present(LATS) ) then
          LOG_ERROR("FILE_TILEDATA_get_info",*) "LATS is required if catalog_fname is empty"
          call PRC_abort
       end if
       if ( .not. present(LATE) ) then
          LOG_ERROR("FILE_TILEDATA_get_info",*) "LATE is required if catalog_fname is empty"
          call PRC_abort
       end if
       if ( .not. present(LONS) ) then
          LOG_ERROR("FILE_TILEDATA_get_info",*) "LONS is required if catalog_fname is empty"
          call PRC_abort
       end if
       if ( .not. present(LONE) ) then
          LOG_ERROR("FILE_TILEDATA_get_info",*) "LONE is required if catalog_fname is empty"
          call PRC_abort
       end if

       TILE_nmax = 1
       TILE_fname(1) = single_fname
       TILE_LATS (1) = LATS
       TILE_LATE (1) = LATE
       TILE_LONS (1) = LONS
       TILE_LONE (1) = LONE

       LAT_MIN = LATS
       LAT_MAX = LATE
    end if

    zonal = ( DOMAIN_LONE - DOMAIN_LONS ) / ( 2.0_RP * PI ) > 0.9_RP

    pole =    ( DOMAIN_LATS < - PI * 0.5_RP + ( DOMAIN_LATE - DOMAIN_LATS ) * 0.1_RP ) &
         .or. ( DOMAIN_LATE >   PI * 0.5_RP - ( DOMAIN_LATE - DOMAIN_LATS ) * 0.1_RP )

    zonal = zonal .or. pole

    call FILE_TILEDATA_get_tile_info( TILE_nmax,                  & ! [IN]
                                      DOMAIN_JS, DOMAIN_JE,       & ! [IN]
                                      DOMAIN_IS, DOMAIN_IE,       & ! [IN]
                                      GLOBAL_IA,                  & ! [IN]
                                      TILE_DLAT, TILE_DLON,       & ! [IN]
                                      TILE_LATS(:), TILE_LATE(:), & ! [IN]
                                      TILE_LONS(:), TILE_LONE(:), & ! [IN]
                                      zonal,                      & ! [IN]
                                      TILE_hit(:),                & ! [OUT]
                                      TILE_JS(:), TILE_JE(:),     & ! [OUT]
                                      TILE_IS(:), TILE_IE(:),     & ! [OUT]
                                      jsh, jeh, ish, ieh,         & ! [OUT]
                                      nLATH, nLONH                ) ! [OUT]

    return
  end subroutine FILE_TILEDATA_get_info

  !-----------------------------------------------------------------------------
  !> get tile data
  subroutine FILE_TILEDATA_get_data_real( &
       nLATH, nLONH,         &
       dirname,              &
       GLOBAL_IA,            &
       TILE_nmax,            &
       TILE_DLAT, TILE_DLON, &
       TILE_fname, TILE_hit, &
       TILE_JS, TILE_JE,     &
       TILE_IS, TILE_IE,     &
       jsh, jeh, ish, ieh,   &
       data_type,            &
       DATA, LATH, LONH,     &
       min_value,            &
       yrevers               )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       PI    => CONST_PI, &
       D2R   => CONST_D2R
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: nLATH
    integer,          intent(in)  :: nLONH
    character(len=*), intent(in)  :: dirname
    integer,          intent(in)  :: GLOBAL_IA
    integer,          intent(in)  :: TILE_nmax
    real(RP),         intent(in)  :: TILE_DLAT
    real(RP),         intent(in)  :: TILE_DLON
    character(len=*), intent(in)  :: TILE_fname(:)
    logical,          intent(in)  :: TILE_hit(:)
    integer,          intent(in)  :: TILE_JS(:)
    integer,          intent(in)  :: TILE_JE(:)
    integer,          intent(in)  :: TILE_IS(:)
    integer,          intent(in)  :: TILE_IE(:)
    integer,          intent(in)  :: jsh
    integer,          intent(in)  :: jeh
    integer,          intent(in)  :: ish
    integer,          intent(in)  :: ieh
    character(len=*), intent(in)  :: data_type
    real(RP),         intent(out) :: DATA(nLONH,nLATH)
    real(RP),         intent(out) :: LATH(nLONH,nLATH)
    real(RP),         intent(out) :: LONH(nLONH,nLATH)
    real(RP),         intent(in), optional :: min_value
    logical,          intent(in), optional :: yrevers

    abstract interface
       subroutine rd( &
            jsize, isize, &
            fname,        &
            TILE_DATA,    &
            yrevers       )
         use scale_precision
         integer,          intent(in)  :: jsize
         integer,          intent(in)  :: isize
         character(len=*), intent(in)  :: fname
         real(RP),         intent(out) :: TILE_DATA(isize,jsize)
         logical, intent(in), optional :: yrevers
       end subroutine rd
    end interface

    procedure(rd), pointer :: read_data

    real(RP) :: min_value_

    character(len=H_LONG) :: fname
    real(RP), allocatable :: TILE_DATA(:,:)
    integer :: jsize, isize
    integer :: i, j, ii, jj, t

    if ( present(min_value) ) then
       min_value_ = min_value
    else
       min_value_ = - abs(UNDEF)
    end if

    select case( data_type )
    case ( "int2", "INT2" )
       read_data => FILE_TILEDATA_read_data_int2_real
    case ( "int4", "INT4" )
       read_data => FILE_TILEDATA_read_data_int4_real
    case ( "real4", "REAL4" )
       read_data => FILE_TILEDATA_read_data_real4_real
    case ( "real8", "REAL8" )
       read_data => FILE_TILEDATA_read_data_real8_real
    case default
       LOG_ERROR("FILE_TILEDATA_get_data_real",*) 'data_type is invalid: ', trim(data_type)
       call PRC_abort
    end select

    !$omp parallel do
!OCL XFILL
    do j = 1, nLATH
    do i = 1, nLONH
       DATA(i,j) = UNDEF
       LATH(i,j) = UNDEF
       LONH(i,j) = UNDEF
    end do
    end do

    do t = 1, TILE_nmax
       if ( .not. TILE_hit(t) ) cycle

       fname = trim(dirname) // '/' // trim(TILE_fname(t))

       LOG_NEWLINE
       LOG_INFO("FILE_TILEDATA_get_data_real",*) 'Input data file :', trim(fname)
       LOG_INFO_CONT(*) 'Tile   (LAT)    :', TILE_JS(t)*TILE_DLAT/D2R, (TILE_JE(t)+1)*TILE_DLAT/D2R
       LOG_INFO_CONT(*) '       (LON)    :', TILE_IS(t)*TILE_DLON/D2R, (TILE_IE(t)+1)*TILE_DLON/D2R

       isize = TILE_IE(t) - TILE_IS(t) + 1
       jsize = TILE_JE(t) - TILE_JS(t) + 1

       allocate( TILE_DATA(isize,jsize) )

       call read_data( jsize, isize,     & ! [IN]
                       fname,            & ! [IN]
                       TILE_DATA(:,:),   & ! [OUT]
                       yrevers = yrevers ) ! [IN]

       !$omp parallel do &
       !$omp private(i,j)
       do jj = 1, jsize
          j = TILE_JS(t) + jj - 1
          if ( jsh <= j .and. j <= jeh ) then
             do ii = 1, isize
                i = TILE_IS(t) + ii - 1
                if ( ish <= i .and. i <= ieh ) then
                   if ( TILE_DATA(ii,jj) < min_value_ ) then
                      DATA(i-ish+1,j-jsh+1) = UNDEF
                   else
                      DATA(i-ish+1,j-jsh+1) = TILE_DATA(ii,jj)
                   end if
                   LATH(i-ish+1,j-jsh+1) = TILE_DLAT * ( TILE_JS(t) + jj - 1 + 0.5_RP )
                   LONH(i-ish+1,j-jsh+1) = TILE_DLON * ( TILE_IS(t) + ii - 1 + 0.5_RP )
                end if
                i = i - GLOBAL_IA
                if ( ish <= i .and. i <= ieh ) then
                   if ( TILE_DATA(ii,jj) < min_value_ ) then
                      DATA(i-ish+1,j-jsh+1) = UNDEF
                   else
                      DATA(i-ish+1,j-jsh+1) = TILE_DATA(ii,jj)
                   end if
                   LATH(i-ish+1,j-jsh+1) = TILE_DLAT * ( TILE_JS(t) + jj - 1 + 0.5_RP )
                   LONH(i-ish+1,j-jsh+1) = TILE_DLON * ( TILE_IS(t) + ii - 1 + 0.5_RP ) - 2.0 * PI
                end if
             end do
          end if
       end do

       deallocate( TILE_DATA )

    enddo ! tile loop

    return
  end subroutine FILE_TILEDATA_get_data_real

  subroutine FILE_TILEDATA_get_data_int1( &
       nLATH, nLONH,         &
       dirname,              &
       GLOBAL_IA,            &
       TILE_nmax,            &
       TILE_DLAT, TILE_DLON, &
       TILE_fname, TILE_hit, &
       TILE_JS, TILE_JE,     &
       TILE_IS, TILE_IE,     &
       jsh, jeh, ish, ieh,   &
       data_type,            &
       DATA, LATH, LONH,     &
       min_value,            &
       yrevers               )
    use scale_const, only: &
       UNDEF2 => CONST_UNDEF2, &
       PI    => CONST_PI, &
       D2R   => CONST_D2R
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: nLATH
    integer,          intent(in)  :: nLONH
    character(len=*), intent(in)  :: dirname
    integer,          intent(in)  :: GLOBAL_IA
    integer,          intent(in)  :: TILE_nmax
    real(RP),         intent(in)  :: TILE_DLAT
    real(RP),         intent(in)  :: TILE_DLON
    character(len=*), intent(in)  :: TILE_fname(:)
    logical,          intent(in)  :: TILE_hit(:)
    integer,          intent(in)  :: TILE_JS(:)
    integer,          intent(in)  :: TILE_JE(:)
    integer,          intent(in)  :: TILE_IS(:)
    integer,          intent(in)  :: TILE_IE(:)
    integer,          intent(in)  :: jsh
    integer,          intent(in)  :: jeh
    integer,          intent(in)  :: ish
    integer,          intent(in)  :: ieh
    character(len=*), intent(in)  :: data_type
    integer,          intent(out) :: DATA(nLONH,nLATH)
    real(RP),         intent(out) :: LATH(nLONH,nLATH)
    real(RP),         intent(out) :: LONH(nLONH,nLATH)

    integer,          intent(in), optional :: min_value
    logical,          intent(in), optional :: yrevers

    abstract interface
       subroutine rd( &
            jsize, isize, &
            fname, &
            TILE_DATA, &
            yrevers      )
         use scale_precision
         integer,          intent(in)  :: jsize
         integer,          intent(in)  :: isize
         character(len=*), intent(in)  :: fname
         integer,          intent(out) :: TILE_DATA(isize,jsize)
         logical, intent(in), optional :: yrevers
       end subroutine rd
    end interface

    integer :: min_value_

    procedure(rd), pointer :: read_data

    character(len=H_LONG) :: fname
    integer, allocatable  :: TILE_DATA(:,:)
    integer :: jsize, isize
    integer :: i, j, ii, jj, t

    if ( present(min_value) ) then
       min_value_ = min_value
    else
       min_value_ = - abs(UNDEF2)
    end if

    select case( data_type )
    case ( "int1", "INT1" )
       read_data => FILE_TILEDATA_read_data_int1_int
    case ( "int2", "INT2" )
       read_data => FILE_TILEDATA_read_data_int2_int
    case ( "int4", "INT4" )
       read_data => FILE_TILEDATA_read_data_int4_int
    case ( "real4", "REAL4" )
       read_data => FILE_TILEDATA_read_data_real4_int
    case default
       LOG_ERROR("FILE_TILEDATA_get_data_int1",*) 'data_type is invalid: ', trim(data_type)
       call PRC_abort
    end select

    !$omp parallel do
!OCL XFILL
    do j = 1, nLATH
    do i = 1, nLONH
       DATA(i,j) = - 1
    end do
    end do

    do t = 1, TILE_nmax
       if ( .not. TILE_hit(t) ) cycle

       fname = trim(dirname) // '/' // trim(TILE_fname(t))

       LOG_NEWLINE
       LOG_INFO("FILE_TILEDATA_get_data_int1",*) 'Input data file :', trim(fname)
       LOG_INFO_CONT(*) 'Tile   (LAT)    :', TILE_JS(t)*TILE_DLAT/D2R, (TILE_JE(t)+1)*TILE_DLAT/D2R
       LOG_INFO_CONT(*) '       (LON)    :', TILE_IS(t)*TILE_DLON/D2R, (TILE_IE(t)+1)*TILE_DLON/D2R

       isize = TILE_IE(t) - TILE_IS(t) + 1
       jsize = TILE_JE(t) - TILE_JS(t) + 1

       allocate( TILE_DATA(isize,jsize) )

       call read_data( jsize, isize,     & ! [IN]
                       fname,            & ! [IN]
                       TILE_DATA(:,:), & ! [OUT]
                       yrevers = yrevers ) ! [IN]

       !$omp parallel do &
       !$omp private(i,j)
       do jj = 1, jsize
       do ii = 1, isize
          i = TILE_IS(t) + ii - 1
          j = TILE_JS(t) + jj - 1
          if ( jsh <= j .and. j <= jeh ) then
             if ( ish <= i .and. i <= ieh ) then
                if ( TILE_DATA(ii,jj) < min_value_ ) then
                   DATA(i-ish+1,j-jsh+1) = UNDEF2
                else
                   DATA(i-ish+1,j-jsh+1) = TILE_DATA(ii,jj)
                end if
                LATH  (i-ish+1,j-jsh+1) = TILE_DLAT * ( TILE_JS(t) + jj - 1 + 0.5_RP )
                LONH  (i-ish+1,j-jsh+1) = TILE_DLON * ( TILE_IS(t) + ii - 1 + 0.5_RP )
             end if
             i = i - GLOBAL_IA
             if ( ish <= i .and. i <= ieh ) then
                if ( TILE_DATA(ii,jj) < min_value_ ) then
                   DATA(i-ish+1,j-jsh+1) = UNDEF2
                else
                   DATA(i-ish+1,j-jsh+1) = TILE_DATA(ii,jj)
                end if
                LATH  (i-ish+1,j-jsh+1) = TILE_DLAT * ( TILE_JS(t) + jj - 1 + 0.5_RP )
                LONH  (i-ish+1,j-jsh+1) = TILE_DLON * ( TILE_IS(t) + ii - 1 + 0.5_RP ) - 2.0 * PI
             end if
          end if
       end do
       end do

       deallocate( TILE_DATA )

    enddo ! tile loop

    return
  end subroutine FILE_TILEDATA_get_data_int1


  ! private


  !-----------------------------------------------------------------------------
  !> read category file
  subroutine FILE_TILEDATA_read_catalog_file( &
       TILE_nlim,            &
       fname,                &
       TILE_DLAT, TILE_DLON, &
       DOMAIN_IS, GLOBAL_IA, &
       TILE_nmax,            &
       TILE_fname,           &
       TILE_LATS, TILE_LATE, &
       TILE_LONS, TILE_LONE  )
    use scale_const, only: &
         D2R => CONST_D2R
    use scale_prc, only: &
         PRC_abort
    integer,          intent(in)  :: TILE_nlim
    character(len=*), intent(in)  :: fname
    real(RP),         intent(in)  :: TILE_DLAT
    real(RP),         intent(in)  :: TILE_DLON
    integer,          intent(in)  :: DOMAIN_IS
    integer,          intent(in)  :: GLOBAL_IA
    integer,          intent(out) :: TILE_nmax
    character(len=*), intent(out) :: TILE_fname(:)
    real(RP),         intent(out) :: TILE_LATS(:)
    real(RP),         intent(out) :: TILE_LATE(:)
    real(RP),         intent(out) :: TILE_LONS(:)
    real(RP),         intent(out) :: TILE_LONE(:)

    integer :: fid, ierr
    integer :: index
    integer :: t

    fid = IO_get_available_fid()
    open( fid,                  &
          file   = fname,       &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_TILEDATA_read_catalog_file",*) 'catalogue file not found! ', trim(fname)
       call PRC_abort
    endif

    do t = 1, TILE_nlim
       read(fid,*,iostat=ierr) index, TILE_LATS(t), TILE_LATE(t), & ! South->North
                                      TILE_LONS(t), TILE_LONE(t), & ! WEST->EAST
                                      TILE_fname(t)

       if ( ierr /= 0 ) exit

       TILE_LATS(t) = TILE_LATS(t) * D2R
       TILE_LATE(t) = TILE_LATE(t) * D2R

       if ( TILE_LONS(t) > TILE_LONE(t) ) TILE_LONE(t) = TILE_LONE(t) + 360.0_RP
       TILE_LONS(t) = TILE_LONS(t) * D2R
       TILE_LONE(t) = TILE_LONE(t) * D2R
    end do

    TILE_nmax = t - 1

    close(fid)

    return
  end subroutine FILE_TILEDATA_read_catalog_file

  subroutine FILE_TILEDATA_get_domain_info( &
       TILE_DLAT, TILE_DLON,                               &
       DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, &
       DOMAIN_JS,   DOMAIN_JE,   DOMAIN_IS,   DOMAIN_IE,   &
       GLOBAL_IA                                           )
    use scale_const, only: &
       PI  => CONST_PI
    real(RP), intent(in)  :: TILE_DLAT
    real(RP), intent(in)  :: TILE_DLON
    real(RP), intent(in)  :: DOMAIN_LATS
    real(RP), intent(in)  :: DOMAIN_LATE
    real(RP), intent(in)  :: DOMAIN_LONS
    real(RP), intent(in)  :: DOMAIN_LONE
    integer,  intent(out) :: DOMAIN_JS
    integer,  intent(out) :: DOMAIN_JE
    integer,  intent(out) :: DOMAIN_IS
    integer,  intent(out) :: DOMAIN_IE
    integer,  intent(out) :: GLOBAL_IA

    DOMAIN_JS = floor  ( DOMAIN_LATS / TILE_DLAT )
    DOMAIN_JE = ceiling( DOMAIN_LATE / TILE_DLAT )
    DOMAIN_IS = floor  ( DOMAIN_LONS / TILE_DLON )
    DOMAIN_IE = ceiling( DOMAIN_LONE / TILE_DLON )

    GLOBAL_IA = nint( 2.0_RP * PI / TILE_DLON )

    return
  end subroutine FILE_TILEDATA_get_domain_info

  subroutine FILE_TILEDATA_get_tile_info( &
       TILE_nmax,            &
       DOMAIN_JS, DOMAIN_JE, &
       DOMAIN_IS, DOMAIN_IE, &
       GLOBAL_IA,            &
       TILE_DLAT, TILE_DLON, &
       TILE_LATS, TILE_LATE, &
       TILE_LONS, TILE_LONE, &
       zonal,                &
       TILE_hit,             &
       TILE_JS, TILE_JE,     &
       TILE_IS, TILE_IE,     &
       jsh, jeh, ish, ieh,   &
       nLATH, nLONH          )
    use scale_const, only: &
       PI => CONST_PI
    integer,  intent(in)  :: TILE_nmax
    integer,  intent(in)  :: DOMAIN_JS, DOMAIN_JE, DOMAIN_IS, DOMAIN_IE
    integer,  intent(in)  :: GLOBAL_IA
    real(RP), intent(in)  :: TILE_DLAT, TILE_DLON
    real(RP), intent(in)  :: TILE_LATS(:), TILE_LATE(:), TILE_LONS(:), TILE_LONE(:)
    logical,  intent(in)  :: zonal
    logical,  intent(out) :: TILE_hit(:)
    integer,  intent(out) :: TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:)
    integer,  intent(out) :: jsh, jeh, ish, ieh
    integer,  intent(out) :: nLATH, nLONH

    logical :: hit_lat, hit_lon
    integer :: nhalo
    integer :: t

    nhalo = 2

    jsh = max( DOMAIN_JS - nhalo, -floor( 0.5_RP * PI / TILE_DLAT ) )
    jeh = min( DOMAIN_JE + nhalo,  floor( 0.5_RP * PI / TILE_DLAT ) )
    ish = DOMAIN_IS - nhalo
    ieh = DOMAIN_IE + nhalo

    ! data file
    !$omp parallel do &
    !$omp private(hit_lat,hit_lon)
    do t = 1, TILE_nmax

       TILE_JS(t) = nint( TILE_LATS(t) / TILE_DLAT )
       TILE_JE(t) = nint( TILE_LATE(t) / TILE_DLAT ) - 1

       TILE_IS(t) = nint( TILE_LONS(t) / TILE_DLON )
       TILE_IE(t) = nint( TILE_LONE(t) / TILE_DLON ) - 1

       do while ( TILE_IE(t) < DOMAIN_IS )
          TILE_IS(t) = TILE_IS(t) + GLOBAL_IA
          TILE_IE(t) = TILE_IE(t) + GLOBAL_IA
       end do
       do while ( TILE_IS(t) - DOMAIN_IS >= GLOBAL_IA )
          TILE_IS(t) = TILE_IS(t) - GLOBAL_IA
          TILE_IE(t) = TILE_IE(t) - GLOBAL_IA
       end do

       if (      ( jsh <= TILE_JS(t) .AND. TILE_JS(t) <= jeh ) &
            .OR. ( jsh <= TILE_JE(t) .AND. TILE_JE(t) <= jeh ) &
            .OR. ( TILE_JS(t) <= jsh .AND. jsh <= TILE_JE(t) ) &
            .OR. ( TILE_JS(t) <= jeh .AND. jeh <= TILE_JE(t) ) ) then
          hit_lat = .true.
       else
          hit_lat = .false.
       endif

       if ( zonal ) then
          hit_lon = .true.
       else if ( ( TILE_IS(t) <= ieh             ) &
            .OR. ( ish <= TILE_IE(t) - GLOBAL_IA ) ) then
          hit_lon = .true.
       else
          hit_lon = .false.
       endif

       TILE_hit(t) = ( hit_lat .AND. hit_lon )
    end do

    if ( zonal ) then
       ish = minval(TILE_IS(1:TILE_nmax))
       ieh = maxval(TILE_IE(1:TILE_nmax))
       jsh = min( jsh, minval(TILE_JS(1:TILE_nmax)) )
       jeh = max( jeh, maxval(TILE_JE(1:TILE_nmax)) )
    end if

    nLONH = ieh - ish + 1
    nLATH = jeh - jsh + 1


    return
  end subroutine FILE_TILEDATA_get_tile_info

  subroutine FILE_TILEDATA_read_data_int2_real( &
       jsize, isize, &
       fname, &
       TILE_DATA, &
       yrevers      )
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: jsize
    integer,          intent(in)  :: isize
    character(len=*), intent(in)  :: fname
    real(RP),         intent(out) :: TILE_DATA(isize,jsize)

    logical, intent(in), optional :: yrevers

    integer(2) :: buf(isize,jsize)

    integer :: fid, ierr
    logical :: yrevers_
    integer :: i, j

    if ( present(yrevers) ) then
       yrevers_ = yrevers
    else
       yrevers_ = .false.
    end if

    fid = IO_get_available_fid()
    open( fid,                    &
          file   = fname,         &
          form   = 'unformatted', &
          access = 'direct',      &
          status = 'old',         &
          recl   = isize*jsize*2, &
          iostat = ierr            )

    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_TILEDATA_read_data_int2_real",*) 'data file not found!: ', trim(fname)
       call PRC_abort
    endif

    read(fid,rec=1) buf(:,:)
    close(fid)

    if ( yrevers_ ) then
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,jsize-j+1)
       end do
       end do
    else
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,j)
       end do
       end do
    end if

    return
  end subroutine FILE_TILEDATA_read_data_int2_real

  subroutine FILE_TILEDATA_read_data_int4_real( &
       jsize, isize, &
       fname, &
       TILE_DATA, &
       yrevers      )
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: jsize
    integer,          intent(in)  :: isize
    character(len=*), intent(in)  :: fname
    real(RP),         intent(out) :: TILE_DATA(isize,jsize)

    logical, intent(in), optional :: yrevers

    integer(4) :: buf(isize,jsize)

    integer :: fid, ierr
    logical :: yrevers_
    integer :: i, j

    if ( present(yrevers) ) then
       yrevers_ = yrevers
    else
       yrevers_ = .false.
    end if

    fid = IO_get_available_fid()
    open( fid,                    &
          file   = fname,         &
          form   = 'unformatted', &
          access = 'direct',      &
          status = 'old',         &
          recl   = isize*jsize*4, &
          iostat = ierr            )

    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_TILEDATA_read_data_int4_real",*) 'data file not found!: ', trim(fname)
       call PRC_abort
    endif

    read(fid,rec=1) buf(:,:)
    close(fid)

    if ( yrevers_ ) then
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,jsize-j+1)
       end do
       end do
    else
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,j)
       end do
       end do
    end if

    return
  end subroutine FILE_TILEDATA_read_data_int4_real

  subroutine FILE_TILEDATA_read_data_real4_real( &
       jsize, isize, &
       fname, &
       TILE_DATA, &
       yrevers      )
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: jsize
    integer,          intent(in)  :: isize
    character(len=*), intent(in)  :: fname
    real(RP),         intent(out) :: TILE_DATA(isize,jsize)

    logical, intent(in), optional :: yrevers

    real(4) :: buf(isize,jsize)

    integer :: fid, ierr
    logical :: yrevers_
    integer :: i, j

    if ( present(yrevers) ) then
       yrevers_ = yrevers
    else
       yrevers_ = .false.
    end if

    fid = IO_get_available_fid()
    open( fid,                    &
          file   = fname,         &
          form   = 'unformatted', &
          access = 'direct',      &
          status = 'old',         &
          recl   = isize*jsize*4, &
          iostat = ierr            )

    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_TILEDATA_read_data_real4_real",*) 'data file not found!: ', trim(fname)
       call PRC_abort
    endif

    read(fid,rec=1) buf(:,:)
    close(fid)

    if ( yrevers_ ) then
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,jsize-j+1)
       end do
       end do
    else
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,j)
       end do
       end do
    end if

    return
  end subroutine FILE_TILEDATA_read_data_real4_real

  subroutine FILE_TILEDATA_read_data_real8_real( &
       jsize, isize, &
       fname, &
       TILE_DATA, &
       yrevers      )
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: jsize
    integer,          intent(in)  :: isize
    character(len=*), intent(in)  :: fname
    real(RP),         intent(out) :: TILE_DATA(isize,jsize)

    logical, intent(in), optional :: yrevers

    real(8) :: buf(isize,jsize)

    integer :: fid, ierr
    logical :: yrevers_
    integer :: i, j

    if ( present(yrevers) ) then
       yrevers_ = yrevers
    else
       yrevers_ = .false.
    end if

    fid = IO_get_available_fid()
    open( fid,                    &
          file   = fname,         &
          form   = 'unformatted', &
          access = 'direct',      &
          status = 'old',         &
          recl   = isize*jsize*8, &
          iostat = ierr            )

    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_TILEDATA_read_data_real8_real",*) 'data file not found!: ', trim(fname)
       call PRC_abort
    endif

    read(fid,rec=1) buf(:,:)
    close(fid)

    if ( yrevers_ ) then
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,jsize-j+1)
       end do
       end do
    else
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,j)
       end do
       end do
    end if

    return
  end subroutine FILE_TILEDATA_read_data_real8_real

  subroutine FILE_TILEDATA_read_data_int1_int( &
       jsize, isize, &
       fname, &
       TILE_DATA, &
       yrevers      )
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: jsize
    integer,          intent(in)  :: isize
    character(len=*), intent(in)  :: fname
    integer,          intent(out) :: TILE_DATA(isize,jsize)

    logical, intent(in), optional :: yrevers

    integer(1) :: buf(isize,jsize)

    integer :: fid, ierr
    logical :: yrevers_
    integer :: i, j

    if ( present(yrevers) ) then
       yrevers_ = yrevers
    else
       yrevers_ = .false.
    end if

    fid = IO_get_available_fid()
    open( fid,                    &
          file   = fname,         &
          form   = 'unformatted', &
          access = 'direct',      &
          status = 'old',         &
          recl   = isize*jsize*1, &
          iostat = ierr            )

    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_TILEDATA_read_data_int1_int",*) 'data file not found!: ', trim(fname)
       call PRC_abort
    endif

    read(fid,rec=1) buf(:,:)
    close(fid)

    if ( yrevers_ ) then
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,jsize-j+1)
       end do
       end do
    else
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,j)
       end do
       end do
    end if

    return
  end subroutine FILE_TILEDATA_read_data_int1_int

  subroutine FILE_TILEDATA_read_data_int2_int( &
       jsize, isize, &
       fname, &
       TILE_DATA, &
       yrevers      )
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: jsize
    integer,          intent(in)  :: isize
    character(len=*), intent(in)  :: fname
    integer,          intent(out) :: TILE_DATA(isize,jsize)

    logical, intent(in), optional :: yrevers

    integer(2) :: buf(isize,jsize)

    integer :: fid, ierr
    logical :: yrevers_
    integer :: i, j

    if ( present(yrevers) ) then
       yrevers_ = yrevers
    else
       yrevers_ = .false.
    end if

    fid = IO_get_available_fid()
    open( fid,                    &
          file   = fname,         &
          form   = 'unformatted', &
          access = 'direct',      &
          status = 'old',         &
          recl   = isize*jsize*2, &
          iostat = ierr            )

    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_TILEDATA_read_data_int2_int",*) 'data file not found!: ', trim(fname)
       call PRC_abort
    endif

    read(fid,rec=1) buf(:,:)
    close(fid)

    if ( yrevers_ ) then
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,jsize-j+1)
       end do
       end do
    else
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,j)
       end do
       end do
    end if

    return
  end subroutine FILE_TILEDATA_read_data_int2_int

  subroutine FILE_TILEDATA_read_data_int4_int( &
       jsize, isize, &
       fname, &
       TILE_DATA, &
       yrevers      )
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: jsize
    integer,          intent(in)  :: isize
    character(len=*), intent(in)  :: fname
    integer,          intent(out) :: TILE_DATA(isize,jsize)

    logical, intent(in), optional :: yrevers

    integer(4) :: buf(isize,jsize)

    integer :: fid, ierr
    logical :: yrevers_
    integer :: i, j

    if ( present(yrevers) ) then
       yrevers_ = yrevers
    else
       yrevers_ = .false.
    end if

    fid = IO_get_available_fid()
    open( fid,                    &
          file   = fname,         &
          form   = 'unformatted', &
          access = 'direct',      &
          status = 'old',         &
          recl   = isize*jsize*4, &
          iostat = ierr            )

    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_TILEDATA_read_data_int4_int",*) 'data file not found!: ', trim(fname)
       call PRC_abort
    endif

    read(fid,rec=1) buf(:,:)
    close(fid)

    if ( yrevers_ ) then
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,jsize-j+1)
       end do
       end do
    else
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,j)
       end do
       end do
    end if

    return
  end subroutine FILE_TILEDATA_read_data_int4_int

  subroutine FILE_TILEDATA_read_data_real4_int( &
       jsize, isize, &
       fname, &
       TILE_DATA, &
       yrevers      )
    use scale_prc, only: &
       PRC_abort
    integer,          intent(in)  :: jsize
    integer,          intent(in)  :: isize
    character(len=*), intent(in)  :: fname
    integer,          intent(out) :: TILE_DATA(isize,jsize)

    logical, intent(in), optional :: yrevers

    real(4) :: buf(isize,jsize)

    integer :: fid, ierr
    logical :: yrevers_
    integer :: i, j

    if ( present(yrevers) ) then
       yrevers_ = yrevers
    else
       yrevers_ = .false.
    end if

    fid = IO_get_available_fid()
    open( fid,                    &
          file   = fname,         &
          form   = 'unformatted', &
          access = 'direct',      &
          status = 'old',         &
          recl   = isize*jsize*4, &
          iostat = ierr            )

    if ( ierr /= 0 ) then
       LOG_ERROR("FILE_TILEDATA_read_data_real4_int",*) 'data file not found!: ', trim(fname)
       call PRC_abort
    endif

    read(fid,rec=1) buf(:,:)
    close(fid)

    if ( yrevers_ ) then
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,jsize-j+1)
       end do
       end do
    else
       !$omp parallel do
!OCL XFILL
       do j = 1, jsize
       do i = 1, isize
          TILE_DATA(i,j) = buf(i,j)
       end do
       end do
    end if

    return
  end subroutine FILE_TILEDATA_read_data_real4_int

end module scale_file_tiledata
