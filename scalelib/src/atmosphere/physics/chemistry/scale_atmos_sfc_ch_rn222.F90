!-------------------------------------------------------------------------------
!> module atmosphere / surface / chemistry / RN222
!!
!! @par Description
!!         Surface emission component for rn222 tracer
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_sfc_ch_rn222
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
  public :: ATMOS_SFC_CH_rn222_setup
  public :: ATMOS_SFC_CH_rn222_finalize
  public :: ATMOS_SFC_CH_rn222_OCEAN_flux
  public :: ATMOS_SFC_CH_rn222_LAND_flux

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
  integer,  private, parameter :: I_ch_rn222 = 1

  character(len=H_SHORT), private :: ATMOS_SFC_CH_Rn222_emission_type = 'CONST'    ! Emission type

  ! for constant flux
  real(RP),               private :: ATMOS_SFC_CH_Rn222_const_emission_land  = 20.8E-3_RP ! Surface flux from land  [Bq/m2/s]
  real(RP),               private :: ATMOS_SFC_CH_Rn222_const_emission_ocean = 0.14E-3_RP ! Surface flux from ocean [Bq/m2/s]

  ! for flux map by Schery and Wasiolek (1998)
  character(len=H_LONG),  private :: ATMOS_SFC_CH_Rn222_SCHERY1998_dirpath = '.'

  ! for flux map by Hirao et al. (2010)
  character(len=H_LONG),  private :: ATMOS_SFC_CH_Rn222_HIRAO2010_dirpath = '.'
  integer,                private :: ATMOS_SFC_CH_Rn222_HIRAO2010_ystart  = 1979
  integer,                private :: ATMOS_SFC_CH_Rn222_HIRAO2010_yend    = 2012

  integer,                private :: ATMOS_SFC_CH_Rn222_nintrp = 5

  real(RP), private, allocatable :: emission_lat  (:,:)
  real(RP), private, allocatable :: emission_lon  (:,:)
  real(RP), private, allocatable :: emission_value(:,:,:,:)

  integer,  private, allocatable :: idx_i(:,:,:)
  integer,  private, allocatable :: idx_j(:,:,:)
  real(RP), private, allocatable :: hfact(:,:,:)

  integer,  private              :: nlon
  integer,  private              :: nlat
  integer,  private              :: nmonth
  integer,  private              :: nyear

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_SFC_CH_rn222_setup( &
       IA, JA,            &
       real_lon, real_lat )
    use scale_prc, only: &
       PRC_IsMaster, &
       PRC_abort
    use scale_const, only: &
       CONST_D2R
    use scale_comm_cartesC, only: &
       COMM_bcast
    use scale_interp, only: &
       INTERP_factor2d
    implicit none

    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    real(RP), intent(in) :: real_lon(IA,JA) ! longitude [rad]
    real(RP), intent(in) :: real_lat(IA,JA) ! latitude  [rad]

    namelist / PARAM_ATMOS_SFC_CH_RN222 / &
       ATMOS_SFC_CH_Rn222_emission_type,        &
       ATMOS_SFC_CH_Rn222_const_emission_land,  &
       ATMOS_SFC_CH_Rn222_const_emission_ocean, &
       ATMOS_SFC_CH_Rn222_SCHERY1998_dirpath,   &
       ATMOS_SFC_CH_Rn222_HIRAO2010_dirpath,    &
       ATMOS_SFC_CH_Rn222_HIRAO2010_ystart,     &
       ATMOS_SFC_CH_Rn222_HIRAO2010_yend,       &
       ATMOS_SFC_CH_Rn222_nintrp

    character(len=H_LONG) :: fname
    real(RP)              :: lon, lat, value

    integer  :: ierr, fid
    integer  :: i, j, m, y, yy
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_SFC_CH_rn222_setup",*) 'Setup'
    LOG_INFO("ATMOS_SFC_CH_rn222_setup",*) 'rn222 surface flux'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_SFC_CH_RN222,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_SFC_CH_rn222_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_SFC_CH_rn222_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_SFC_CH_RN222. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_SFC_CH_RN222)



    LOG_NEWLINE
    LOG_INFO("ATMOS_SFC_CH_rn222_setup",*) 'Type of emission of Rn222: ', trim(ATMOS_SFC_CH_RN222_emission_type)

    select case( ATMOS_SFC_CH_RN222_emission_type )
    case( 'CONST' )

       LOG_INFO_CONT('(A,ES16.6)') 'From land  [Bq/m2/s] : ', ATMOS_SFC_CH_Rn222_const_emission_land
       LOG_INFO_CONT('(A,ES16.6)') 'From ocean [Bq/m2/s] : ', ATMOS_SFC_CH_Rn222_const_emission_ocean

    case( 'SCHERY1998' )

       LOG_INFO_CONT(*) 'Flux map by Schery and Wasiolek (1998) is used'

       nlon   = 360
       nlat   = 180
       nmonth = 12
       nyear  = 1

       allocate( emission_lon  (nlon,nlat) )
       allocate( emission_lat  (nlon,nlat) )
       allocate( emission_value(nlon,nlat,nmonth,nyear) )

       if ( PRC_IsMaster ) then
          y = 1
          do m = 1, nmonth
             write(fname,'(A,A,I2.2,A)') trim(ATMOS_SFC_CH_Rn222_SCHERY1998_dirpath), "/fdh3a.", m
             LOG_INFO_CONT(*) 'Read from the ASCII file: ', trim(fname)

             fid = IO_get_available_fid()
             open( unit = fid,         &
                   file = trim(fname), &
                   status = "old",     &
                   form = "formatted"  )

                do j = 1, nlat
                do i = 1, nlon
                   lon = real(i-1,kind=RP) - 180.0_RP ! [180W-179E]
                   lat = 90.0_RP - real(j-1,kind=RP)  ! [90N-89S]

                   read(fid,*) value

                   emission_lon  (i,j)     = ( lon + 0.5_RP ) * CONST_D2R ! shift +0.5deg, [deg->rad]
                   emission_lat  (i,j)     = ( lat - 0.5_RP ) * CONST_D2R ! shift -0.5deg, [deg->rad]
                   emission_value(i,j,m,y) = value * 1.E-3_RP             ! [mBq/m2->Bq/m2]
                enddo
                enddo

             close(fid)
          enddo ! month loop
       endif

    case( 'HIRAO2010' )

       LOG_INFO_CONT(*) 'Flux map by Hirao et al. (2010) is used'
       LOG_INFO_CONT(*) 'Start year: ', ATMOS_SFC_CH_Rn222_HIRAO2010_ystart
       LOG_INFO_CONT(*) 'End   year: ', ATMOS_SFC_CH_Rn222_HIRAO2010_yend

       if (      ATMOS_SFC_CH_Rn222_HIRAO2010_ystart < 1979 &
            .OR. ATMOS_SFC_CH_Rn222_HIRAO2010_yend   > 2012 ) then
          LOG_WARN('ATMOS_SFC_CH_rn222_setup',*) 'Available period of the data is between 1979 and 2012.'
          LOG_WARN_CONT(*)                       'Please check the range of ystart and yend.'
          ATMOS_SFC_CH_Rn222_HIRAO2010_ystart = max( ATMOS_SFC_CH_Rn222_HIRAO2010_ystart, 1979 )
          ATMOS_SFC_CH_Rn222_HIRAO2010_yend   = min( ATMOS_SFC_CH_Rn222_HIRAO2010_yend,   2012 )
       endif

       nlon   = 360
       nlat   = 180
       nmonth = 12
       nyear  = ATMOS_SFC_CH_Rn222_HIRAO2010_yend - ATMOS_SFC_CH_Rn222_HIRAO2010_ystart + 1

       allocate( emission_lon  (nlon,nlat) )
       allocate( emission_lat  (nlon,nlat) )
       allocate( emission_value(nlon,nlat,nmonth,nyear) )

       if ( PRC_IsMaster ) then
          do y = 1, nyear
          do m = 1, nmonth
             yy = y+ATMOS_SFC_CH_Rn222_HIRAO2010_ystart-1
             write(fname,'(A,A,I4.4,I2.2)') trim(ATMOS_SFC_CH_Rn222_HIRAO2010_dirpath), "/flux-hra-revi", yy, m
             LOG_INFO_CONT(*) 'Read from the ASCII file: ', trim(fname)

             fid = IO_get_available_fid()
             open( unit = fid,         &
                   file = trim(fname), &
                   status = "old",     &
                   form = "formatted"  )

                do j = 1, nlat
                do i = 1, nlon
                   read(fid,*) lon, lat, value

                   emission_lon  (i,j)     = ( lon + 0.5_RP ) * CONST_D2R ! shift +0.5deg, [deg->rad]
                   emission_lat  (i,j)     = ( lat - 0.5_RP ) * CONST_D2R ! shift -0.5deg, [deg->rad]
                   emission_value(i,j,m,y) = value * 1.E-3_RP             ! [mBq/m2->Bq/m2]
                enddo
                enddo

             close(fid)
          enddo ! month loop
          enddo ! year loop
       endif

    case default
       LOG_ERROR("ATMOS_SFC_CH_rn222_setup",*) 'Not supported type of Rn222 emission! Stop.'
       call PRC_abort
    end select

    select case( ATMOS_SFC_CH_RN222_emission_type )
    case( 'SCHERY1998', 'HIRAO2010' )

       call COMM_bcast( emission_lon  (:,:),     nlon, nlat )
       call COMM_bcast( emission_lat  (:,:),     nlon, nlat )
       call COMM_bcast( emission_value(:,:,:,:), nlon, nlat, nmonth, nyear )

       allocate( idx_i(IA,JA,ATMOS_SFC_CH_Rn222_nintrp) )
       allocate( idx_j(IA,JA,ATMOS_SFC_CH_Rn222_nintrp) )
       allocate( hfact(IA,JA,ATMOS_SFC_CH_Rn222_nintrp) )

       call INTERP_factor2d( ATMOS_SFC_CH_Rn222_nintrp, & ! [IN]
                             nlon, nlat,                & ! [IN]
                             IA, JA,                    & ! [IN]
                             emission_lon(:,:),         & ! [IN]
                             emission_lat(:,:),         & ! [IN]
                             real_lon(:,:),             & ! [IN]
                             real_lat(:,:),             & ! [IN]
                             idx_i   (:,:,:),           & ! [OUT]
                             idx_j   (:,:,:),           & ! [OUT]
                             hfact   (:,:,:)            ) ! [OUT]
    end select

    return
  end subroutine ATMOS_SFC_CH_rn222_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_SFC_CH_rn222_finalize

    if ( allocated( emission_lon ) ) deallocate( emission_lon )
    if ( allocated( emission_lat ) ) deallocate( emission_lat )
    if ( allocated( emission_value ) ) deallocate( emission_value )

    if ( allocated( idx_i ) ) deallocate( idx_i )
    if ( allocated( idx_j ) ) deallocate( idx_j )
    if ( allocated( hfact ) ) deallocate( hfact )

    return
  end subroutine ATMOS_SFC_CH_rn222_finalize

  !-----------------------------------------------------------------------------
  !> Emission from the ocean surface
  subroutine ATMOS_SFC_CH_rn222_OCEAN_flux( &
       IA, IS, IE, &
       JA, JS, JE, &
       QA_CH,      &
       SFLX_QTRC   )
    implicit none

    integer,  intent(in)    :: IA, IS, IE
    integer,  intent(in)    :: JA, JS, JE
    integer,  intent(in)    :: QA_CH
    real(RP), intent(inout) :: SFLX_QTRC(IA,JA,QA_CH)

    integer :: i, j
    !---------------------------------------------------------------------------

    select case( ATMOS_SFC_CH_RN222_emission_type )
    case( 'CONST', 'SCHERY1998', 'HIRAO2010' )

       do j = JS, JE
       do i = IS, IE
          SFLX_QTRC(i,j,I_ch_rn222) = ATMOS_SFC_CH_Rn222_const_emission_ocean
       enddo
       enddo

    end select

    return
  end subroutine ATMOS_SFC_CH_rn222_OCEAN_flux

  !-----------------------------------------------------------------------------
  !> Emission from the land surface
  subroutine ATMOS_SFC_CH_rn222_LAND_flux( &
       IA, IS, IE, &
       JA, JS, JE, &
       QA_CH,      &
       NOWDATE,    &
       SFLX_QTRC   )
    use scale_prc, only: &
       PRC_abort
    use scale_interp, only: &
       INTERP_interp2d
    implicit none

    integer,  intent(in)    :: IA, IS, IE
    integer,  intent(in)    :: JA, JS, JE
    integer,  intent(in)    :: QA_CH
    integer,  intent(in)    :: NOWDATE(6)
    real(RP), intent(inout) :: SFLX_QTRC(IA,JA,QA_CH)

    integer :: i, j, m, y, yy
    !---------------------------------------------------------------------------

    select case( ATMOS_SFC_CH_RN222_emission_type )
    case( 'CONST' )

       do j = JS, JE
       do i = IS, IE
          SFLX_QTRC(i,j,I_ch_rn222) = ATMOS_SFC_CH_Rn222_const_emission_land
       enddo
       enddo

    case( 'SCHERY1998' )

       y = 1
       m = NOWDATE(2)

       call INTERP_interp2d( ATMOS_SFC_CH_Rn222_nintrp, & ! [IN]
                             nlon, nlat,                & ! [IN]
                             IA, JA,                    & ! [IN]
                             idx_i(:,:,:),              & ! [IN]
                             idx_j(:,:,:),              & ! [IN]
                             hfact(:,:,:),              & ! [IN]
                             emission_value(:,:,m,y),   & ! [IN]
                             SFLX_QTRC(:,:,I_ch_rn222)  ) ! [OUT]

    case( 'HIRAO2010' )

       yy = NOWDATE(1)
       yy = max( yy, 1979 ) ! Use flux of 1979 before 1977
       yy = min( yy, 2012 ) ! Use flux of 2012 after  2013

       y = yy-ATMOS_SFC_CH_Rn222_HIRAO2010_ystart+1

       if ( y < 1 .OR. y > nyear ) then
          LOG_ERROR("ATMOS_SFC_CH_rn222_setup",*) 'emission file does not exist for year=', yy
          call PRC_abort
       endif

       m = NOWDATE(2)

       call INTERP_interp2d( ATMOS_SFC_CH_Rn222_nintrp, & ! [IN]
                             nlon, nlat,                & ! [IN]
                             IA, JA,                    & ! [IN]
                             idx_i(:,:,:),              & ! [IN]
                             idx_j(:,:,:),              & ! [IN]
                             hfact(:,:,:),              & ! [IN]
                             emission_value(:,:,m,y),   & ! [IN]
                             SFLX_QTRC(:,:,I_ch_rn222)  ) ! [OUT]

    end select

    return
  end subroutine ATMOS_SFC_CH_rn222_LAND_flux

end module scale_atmos_sfc_ch_rn222
