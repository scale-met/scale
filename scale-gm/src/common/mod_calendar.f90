!-------------------------------------------------------------------------------
!> Module calendar
!!
!! @par Description
!!         This module contains subroutines the calendar.
!!         Orginal code is from CCSR/NIES/AGCM5.4
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_calendar
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: calendar_setup
  public :: calendar_ointvl
  public :: calendar_ss2cc
  public :: calendar_ssaft
  public :: calendar_secdy
  public :: calendar_ss2ds
  public :: calendar_ss2yd
  public :: calendar_ym2dd
  public :: calendar_ds2ss
  public :: calendar_dayyr
  public :: calendar_ss2yh
  public :: calendar_xx2ss
  public :: calendar_yh2ss
  public :: calendar_PERPR
  public :: calendar_dd2ym
  public :: calendar_daymo

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
  integer, private :: isecdy, isechr, jy, jy4, jcent, jcent4
  integer, private :: idays0, idy, ileap, id, m, idayyr, jyear, jmonth
  !
  !--- flag of automatic or not
  logical, private ::  oauto = .true.
  !<---                      yr=0-999     : 360day
  !<---                      yr=1000-1899 : 365day
  !<---                      yr=1900-     : gregorian
  !
  !--- flag of the gregorian calendar or not
  logical, private :: ogrego = .true.
  !
  !--- number of days in a month
  integer, private :: monday ( 12,2 )
  data   monday /                            &
       31,28,31,30,31,30,31,31,30,31,30,31,  &
       31,29,31,30,31,30,31,31,30,31,30,31 /
  !
  !--- flag of ideal calender (n day per month)
  logical, private :: oideal = .false.
  !------ 1 month = x days in the ideal case
  integer, private :: idaymo = 30
  !------ 1 year = x months in the ideal case
  integer, private :: imonyr = 12
  !
  !--- flag of perpetual or not
  logical, private :: operpt = .false.
  !------ perpetual date(year)
  integer, private :: iyrpp  = 0
  !------ perpetual date(month)
  integer, private :: imonpp = 3
  !------ perpetual date(day)
  integer, private :: idaypp = 21
  !
  !--- 1 minute = x sec.
  integer, private :: isecmn = 60
  !--- 1 hour = x minutes
  integer, private :: iminhr = 60
  !--- 1 day = x hours
  integer, private :: ihrday = 24
  !
  logical, private :: ooperz = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine calendar_setup
    use mod_adm, only: &
       ADM_LOG_FID, &
       ADM_CTL_FID, &
       ADM_proc_stop
    implicit none

    namelist  / NM_CALENDAR /  &
         oauto ,             &
         ogrego,             &
         oideal,             &
         operpt,             &
         idaymo,             &
         imonyr,             &
         iyrpp,              &
         imonpp,             &
         idaypp,             &
         isecmn,             &
         iminhr,             &
         ihrday

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[calendar]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=NM_CALENDAR,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** NM_CALENDAR is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist NM_CALENDAR. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist NM_CALENDAR. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=NM_CALENDAR)

    return
  end subroutine calendar_setup

  !-----------------------------------------------------------------------------
  !> fixed date
  subroutine calendar_perpt( &
       iyear ,               & !--- in : year
       imonth,               & !--- in : month
       iday                  & !--- in : day
       )
    implicit none

    integer, intent(in) :: iyear
    integer, intent(in) :: imonth
    integer, intent(in) :: iday
    !---------------------------------------------------------------------------

    operpt = .true.
    iyrpp  = iyear
    imonpp = imonth
    idaypp = iday

    return
  end subroutine calendar_perpt

  !-----------------------------------------------------------------------------
  !> calendar, refer to fixed date
  subroutine calendar_PERPR( &
       IYEAR ,               & !--- OUT
       IMONTH,               & !--- OUT
       IDAY,                 & !--- OUT
       OOPERP                & !--- OUT
       )
    implicit none

    integer, intent(out) :: IYEAR
    integer, intent(out) :: IMONTH
    integer, intent(out) :: IDAY
    logical, intent(out) :: OOPERP
    !---------------------------------------------------------------------------

    OOPERP = OPERPT
    IYEAR  = IYRPP
    IMONTH = IMONPP
    IDAY   = IDAYPP

    return
  end subroutine calendar_PERPR

  !-----------------------------------------------------------------------------
  !> fixed date(on/off)
  subroutine calendar_perpo( &
       ooperp                &
       )
    implicit none

    logical, intent(in) :: ooperp !--- in : flag of perpetual
    !---------------------------------------------------------------------------

    ooperz = ooperp

    return
  end subroutine calendar_perpo

  !-----------------------------------------------------------------------------
  !> calendar, no.of day in a month
  subroutine calendar_daymo(&
       ndaymo,              & !--- OUT : day
       iyear,               & !--- IN : year
       imonth               & !--- IN : month
       )
    implicit none

    integer, intent(out) :: ndaymo
    integer, intent(in) :: iyear
    integer, intent(in) :: imonth
    !---------------------------------------------------------------------------

    if ( oauto ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego ) then
       if ( ocleap( iyear ) ) then
          ndaymo = monday( imonth,2 )
       else
          ndaymo = monday( imonth,1 )
       endif
    elseif ( .not. oideal ) then
       ndaymo = monday( imonth,1 )
    else
       ndaymo = idaymo
    endif

    return
  end subroutine calendar_daymo

  !-----------------------------------------------------------------------------
  !> calendar, no.of day in an year
  subroutine calendar_dayyr(&
       ndayyr,              & !-- OUT : day
       iyear                & !---IN : year
       )
    implicit none

    integer, intent(out) :: ndayyr
    integer, intent(in) :: iyear
    !---------------------------------------------------------------------------

    if ( oauto ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego ) then
       if ( ocleap( iyear ) ) then
          ndayyr = 366
       else
          ndayyr = 365
       endif
    elseif ( .not. oideal ) then
       ndayyr = 365
    else
       ndayyr = idaymo*imonyr
    endif

    return
  end subroutine calendar_dayyr

  !-----------------------------------------------------------------------------
  !> calendar, no.of month in an year
  subroutine calendar_monyr( &
       nmonyr,               & !--- OUT : month
       iyear                 & !--- IN : year
       )
    implicit none

    integer, intent(out) :: nmonyr
    integer, intent(in) :: iyear
    !---------------------------------------------------------------------------

    if ( oauto  ) then
       nmonyr = 12
    else
       nmonyr = imonyr
    endif

    return
  end subroutine calendar_monyr

  !-----------------------------------------------------------------------------
  !> calendar, no.of sec. in a day
  subroutine calendar_secdy( &
       nsecdy                &
       )
    implicit none

    integer, intent(out) :: nsecdy !--- OUT : sec
    !---------------------------------------------------------------------------

    nsecdy = isecmn*iminhr*ihrday

    return
  end subroutine calendar_secdy

  !-----------------------------------------------------------------------------
  !> calendar, no of sec. in a minute
  subroutine calendar_secmi( &
       nsecmi                & !--- OUT
       )
    implicit none

    integer, intent(out) :: nsecmi
    !---------------------------------------------------------------------------

    nsecmi = isecmn

    return
  end subroutine calendar_secmi

  !-----------------------------------------------------------------------------
  !> calendar, no.of sec. in an hour
  subroutine calendar_sechr( &
       nsechr                &
       )
    implicit none

    integer, intent(out) :: nsechr
    !---------------------------------------------------------------------------

    nsechr = isecmn*iminhr

    return
  end subroutine calendar_sechr

  !-----------------------------------------------------------------------------
  !> calendar, sec. -> ddss
  subroutine calendar_ss2ds( &
       idays ,               & !--- OUT
       rsec  ,               & !--- OUT
       dsec                  & !--- IN
       )
    implicit none

    integer, intent(out) :: idays
    real(RP), intent(out) :: rsec
    real(RP), intent(in) :: dsec
    !---------------------------------------------------------------------------

    isecdy = isecmn*iminhr*ihrday
    idays  = int( dsec/real(isecdy,kind=RP) ) + 1
    rsec   = dsec - real(idays-1,kind=RP) * real(isecdy,kind=RP)

    if ( nint( rsec ) >= isecdy ) then
       idays = idays + 1
       rsec  = rsec - real(isecdy,kind=RP)
    endif

    return
  end subroutine calendar_ss2ds

  !-----------------------------------------------------------------------------
  !> calendar, ddss -> sec.
  subroutine calendar_ds2ss( &
       dsec  ,               & !--- OUT
       idays ,               & !--- IN
       rsec                  & !--- IN
       )
    implicit none

    real(RP), intent(out) :: dsec
    integer, intent(in) :: idays
    real(RP), intent(in) :: rsec
    !---------------------------------------------------------------------------

    isecdy = isecmn*iminhr*ihrday
    dsec   = real(idays-1,kind=RP) * real(isecdy,kind=RP) + real(rsec,kind=RP)

    return
  end subroutine calendar_ds2ss

  !-----------------------------------------------------------------------------
  !> calendar, sec. -> hhmmss
  subroutine calendar_rs2hm( &
       ihour ,               & !--- OUT
       imin  ,               & !--- OUT
       isec  ,               & !--- OUT
       rsec                  & !--- IN
       )
    implicit none

    integer, intent(out) :: ihour
    integer, intent(out) :: imin
    integer, intent(out) :: isec
    real(RP), intent(in) :: rsec
    !---------------------------------------------------------------------------

    isechr = isecmn*iminhr
    ihour  = int ( rsec / real(isechr,kind=RP) )
    imin   = int ( ( rsec - real(ihour*isechr,kind=RP) ) / real(isecmn,kind=RP) )
    isec   = nint( rsec - real(ihour*isechr,kind=RP) - real(imin*isecmn,kind=RP) )

    if ( isec >= isecmn ) then
       imin  = imin + 1
       isec  = isec - isecmn
    endif

    if ( imin == iminhr ) then
       ihour = ihour + 1
       imin  = imin  - iminhr
    endif

    return
  end subroutine calendar_rs2hm

  !-----------------------------------------------------------------------------
  !> calendar, hhmmss -> sec.
  subroutine calendar_hm2rs( &
       rsec  ,               & !--- OUT
       ihour ,               & !--- IN
       imin  ,               & !--- IN
       isec                  & !--- IN
       )
    implicit none

    real(RP), intent(out) :: rsec
    integer, intent(in) :: ihour
    integer, intent(in) :: imin
    integer, intent(in) :: isec
    !---------------------------------------------------------------------------

    rsec = real( ihour*isecmn*iminhr + imin*isecmn + isec, kind=RP )

    return
  end subroutine calendar_hm2rs

  !-----------------------------------------------------------------------------
  !> calendar, day -> yymmdd
  subroutine calendar_dd2ym( &
       iyear ,               & !--- OUT
       imonth,               & !--- OUT
       iday  ,               & !--- OUT
       idays                 & !--- IN
       )
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: imonth
    integer, intent(out) :: iday
    integer, intent(in) :: idays
    !---------------------------------------------------------------------------

    if ( oauto  ) then
       if ( idays >= 693961 ) then       !" 1900*365+1900/4-19+5
          ogrego = .true.
       else
          ogrego = .false.
          if ( idays >= 1000*365 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( operpt .and. ooperz ) then
       iyear  = iyrpp
       imonth = imonpp
       iday   = idaypp

       return
    endif

    if ( ogrego ) then
       jy     = int(real(idays,kind=RP)/365.24_RP)
1100   continue
       jy4    = (jy+3)/4
       jcent  = (jy+99)/100
       jcent4 = (jy+399)/400
       idays0 = jy*365+jy4-jcent+jcent4
       if ( idays <= idays0 ) then
          jy = jy -1
          if ( jy >= 0 ) goto 1100
       endif
       iyear = jy
       idy   = idays - idays0
       if ( ocleap( iyear ) ) then
          ileap  = 2
       else
          ileap  = 1
       endif
    elseif ( .not. oideal ) then
       iyear = idays/365
       idy   = idays - iyear*365
       ileap = 1
    endif

    if ( ogrego .or. .not. oideal ) then
       id = 0
       do m = 1, 12
          id = id + monday(m,ileap)
          if ( idy <= id ) then
             imonth = m
             iday   = idy-id+monday(m,ileap)
             goto 3190
          endif
       enddo
3190   continue
    else
       idayyr = idaymo*imonyr
       iyear  = ( idays-1 ) / idayyr
       imonth = ( idays-1 - iyear*idayyr )/idaymo+1
       iday   = idays - iyear*idayyr - (imonth-1)*idaymo
    endif

    return
  end subroutine calendar_dd2ym

  !-----------------------------------------------------------------------------
  !> calendar, yymmdd -> day
  subroutine calendar_ym2dd( &
       idays ,               & !--- OUT
       iyear ,               & !--- IN
       imonth,               & !--- IN
       iday                  & !--- IN
       )
    implicit none

    integer, intent(out) :: idays
    integer, intent(in) :: iyear
    integer, intent(in) :: imonth
    integer, intent(in) :: iday
    !---------------------------------------------------------------------------

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego .or. .not. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1)/12
          jmonth = mod(imonth-1,12)+1
       else
          jyear  = iyear - (-imonth)/12 - 1
          jmonth = 12-mod(-imonth,12)
       endif
    endif

    if ( ogrego ) then
       jy4    = (jyear+3)/4
       jcent  = (jyear+99)/100
       jcent4 = (jyear+399)/400
       idays0 = jyear*365+jy4-jcent+jcent4
       if ( ocleap( jyear ) ) then
          ileap = 2
       else
          ileap = 1
       endif
    elseif ( .not. oideal ) then
       idays0 = jyear*365
       ileap  = 1
    endif

    if ( ogrego .or. .not. oideal ) then
       id = 0
       do m = 1, jmonth-1
          id = id + monday(m,ileap)
       enddo
    else
       idays0 = iyear*idaymo*imonyr
       id     = (imonth-1)*idaymo
    endif

    idays = idays0 + id + iday

    return
  end subroutine calendar_ym2dd

  !-----------------------------------------------------------------------------
  !> calendar, yymmdd -> yydd
  subroutine calendar_ym2yd( &
       idaysy,               & !--- OUT
       iyear ,               & !--- IN
       imonth,               & !--- IN
       iday                  & !--- IN
       )
    implicit none

    integer, intent(out) :: idaysy
    integer, intent(in) :: iyear
    integer, intent(in) :: imonth
    integer, intent(in) :: iday
    !---------------------------------------------------------------------------

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego .or. .not. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1)/12
          jmonth = mod(imonth-1,12)+1
       else
          jyear  = iyear - (-imonth)/12 - 1
          jmonth = 12-mod(-imonth,12)
       endif
    endif

    if ( ogrego ) then
       if ( ocleap( jyear ) ) then
          ileap = 2
       else
          ileap = 1
       endif
    elseif ( .not. oideal ) then
       ileap  = 1
    endif

    if ( ogrego .or. .not. oideal ) then
       id = 0
       do m = 1, jmonth-1
          id = id + monday(m,ileap)
       enddo
    else
       id     = (imonth-1)*idaymo
    endif

    idaysy = id + iday

    return
  end subroutine calendar_ym2yd

  !-----------------------------------------------------------------------------
  !> calendar, sec. -> date
  subroutine calendar_ss2yh( &
       idate ,               & !--- OUT
       dsec                  & !--- IN
       )
    implicit none

    integer, intent(out) :: idate(6) ! yymmddhhmmss
    real(RP), intent(in) :: dsec      ! time

    integer    idays             ! serial no.of day
    real(RP)    rsec              ! no. of sec. in a day
    !---------------------------------------------------------------------------

    call calendar_ss2ds &
         ( idays , rsec  , &
         dsec            )

    call calendar_dd2ym &
         ( idate(1), idate(2), idate(3), &
         idays                        )

    call calendar_rs2hm &
         ( idate(4), idate(5), idate(6), &
         rsec                          )

    return
  end subroutine calendar_ss2yh

  !-----------------------------------------------------------------------------
  !> calendar, date -> sec.
  subroutine calendar_yh2ss( &
       dsec  ,               & !--- OUT
       idate                 & !--- IN
       )
    implicit none

    real(RP), intent(out) :: dsec    ! time
    integer, intent(in) :: idate(6) ! yymmddhhmmss

    integer    idays             ! serial no.of day
    real(RP)    rsec              ! no. of sec. in a day
    !---------------------------------------------------------------------------

    call calendar_ym2dd &
         ( idays   , &
         idate(1), idate(2), idate(3) )

    call calendar_hm2rs &
         ( rsec    , &
         idate(4), idate(5), idate(6) )

    call calendar_ds2ss &
         ( dsec  , &
         idays , rsec   )

    return
  end subroutine calendar_yh2ss

  !-----------------------------------------------------------------------------
  !> calendar, day -> yydd
  subroutine calendar_dd2yd( &
       iyear ,               & !--- OUT
       idaysy,               & !--- OUT
       idays                 & !--- IN
       )
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: idaysy
    integer, intent(in) :: idays

    integer :: imonth, iday
    !---------------------------------------------------------------------------

    call calendar_dd2ym &
         ( iyear , imonth, iday  , &
         idays                   )

    call calendar_ym2yd &
         ( idaysy, &
         iyear , imonth, iday    )

    return
  end subroutine calendar_dd2yd

  !-----------------------------------------------------------------------------
  !> calendar, sec. -> yydd
  subroutine calendar_ss2yd( &
       iyear ,               & !--- OUT
       idaysy,               & !--- OUT
       dsec                  & !--- IN
       )
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: idaysy
    real(RP), intent(in) :: dsec

    integer :: imonth, iday
    !---------------------------------------------------------------------------

    call calendar_ss2ym &
         ( iyear , imonth, iday  , &
         dsec                    )

    call calendar_ym2yd &
         ( idaysy, &
         iyear , imonth, iday    )

    return
  end subroutine calendar_ss2yd

  !-----------------------------------------------------------------------------
  !> calendar, sec. -> yymmdd
  subroutine calendar_ss2ym( &
       iyear ,               & !--- OUT
       imonth,               & !--- OUT
       iday  ,               & !--- OUT
       dsec                  & !--- IN
       )
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: imonth
    integer, intent(out) :: iday
    real(RP), intent(in) :: dsec

    integer :: idays
    real(RP) :: rsec
    !---------------------------------------------------------------------------

    call calendar_ss2ds &
         ( idays , rsec  , &
         dsec            )

    call calendar_dd2ym &
         ( iyear , imonth, iday  , &
         idays                  )

    return
  end subroutine calendar_ss2ym

  !-----------------------------------------------------------------------------
  !> calendar, hour ->sec.
  subroutine calendar_xx2ss( &
       ddsec ,               & !--- OUT
       rtdur ,               & !--- IN
       hunit,                & !--- IN
       dsec                 & !--- IN
       )
    implicit none

    real(RP), intent(out) :: ddsec
    real(RP), intent(in) :: rtdur
    character, intent(in) :: hunit *(*)
    real(RP), intent(in) :: dsec

    character(len=10) :: hunitx
    integer :: isecmi, isechr, isecdy
    integer :: iyear, imonth, iday, ndaymo, ndayyr
    !---------------------------------------------------------------------------

    hunitx = hunit

    if      ( hunitx(1:1) == 's' .or. hunitx(1:1) == 'S' ) then
       ddsec = real(rtdur,kind=RP)
    elseif ( hunitx(1:2) == 'mi' .or. hunitx(1:2) == 'MI' ) then
       call calendar_secmi( isecmi )
       ddsec = real(rtdur,kind=RP) * real(isecmi,kind=RP)
    elseif ( hunitx(1:1) == 'h' .or. hunitx(1:1) == 'H' ) then
       call calendar_sechr( isechr )
       ddsec = real(rtdur,kind=RP) * real(isechr,kind=RP)
    elseif ( hunitx(1:1) == 'd' .or. hunitx(1:1) == 'D' ) then
       call calendar_secdy( isecdy )
       ddsec = real(rtdur,kind=RP) * real(isecdy,kind=RP)
    elseif ( hunitx(1:2) == 'mo' .or. hunitx(1:2) == 'MO' ) then
       call calendar_ss2ym &
            ( iyear , imonth, iday  , &
            dsec                    )
       call calendar_daymo &
            ( ndaymo, &
            iyear , imonth  )
       call calendar_secdy( isecdy )
       ddsec = real(rtdur,kind=RP) * real(ndaymo,kind=RP) * real(isecdy,kind=RP)
    elseif ( hunitx(1:1) == 'y' .or. hunitx(1:1) == 'Y' ) then
       call calendar_ss2ym &
            ( iyear , imonth, iday  , &
            dsec                    )
       call calendar_dayyr &
            ( ndayyr, &
            iyear   )
       call calendar_secdy( isecdy )
       ddsec = real(rtdur * ndayyr * isecdy,kind=RP)
    else
       write (6,*) ' ### cxx2ss: invalid unit : ', hunit, &
            ' [sec] assumed'
       ddsec = rtdur
    endif

    return
  end subroutine calendar_xx2ss

  !-----------------------------------------------------------------------------
  !> calendar, character -> date
  subroutine calendar_cc2yh( &
       itime ,               & !--- OUT
       htime                 & !--- IN
       )
    implicit none

    integer, intent(out) :: itime ( 6 )
    character, intent(in) :: htime *(*)

    integer :: i
    !---------------------------------------------------------------------------

    read ( htime, '( i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2 )' ) (itime(i),i=1,6)

    return
  end subroutine calendar_cc2yh

  !-----------------------------------------------------------------------------
  !> calendar, date -> character
  subroutine calendar_yh2cc( &
       htime ,               & !--- OUT
       itime                 & !--- IN
       )
    implicit none

    character, intent(out) :: htime *(*)
    integer, intent(in) :: itime ( 6 )

    integer :: i
    !---------------------------------------------------------------------------

    write( htime, '( i4.4,"/",i2.2,"/",i2.2,"-",i2.2,":",i2.2,":",i2.2 )' ) (itime(i),i=1,6)

    return
  end subroutine calendar_yh2cc

  !-----------------------------------------------------------------------------
  !> calendar, sec. -> character
  subroutine calendar_ss2cc( &
       htime ,               & !--- OUT
       dsec                  & !--- IN
       )
    implicit none

    character, intent(out) :: htime *(*)
    real(RP), intent(in) :: dsec

    integer :: itime(6), i
    !---------------------------------------------------------------------------

    call calendar_ss2yh &
         ( itime , &
         dsec   )

    write ( htime, '( i4.4,"/",i2.2,"/",i2.2,"-",i2.2,":",i2.2,":",i2.2 )' ) (itime(i),i=1,6)

    return
  end subroutine calendar_ss2cc

  !-----------------------------------------------------------------------------
  !> calendar :bissextile or not
  !> @return ocleap
  function ocleap( &
       iyear       & !--- IN
       )
    implicit none

    logical :: ocleap
    integer, intent(in) :: iyear

    integer :: iy, iycen, icent
    !---------------------------------------------------------------------------

    iy     = mod(iyear,4)
    iycen  = mod(iyear,100)
    icent  = mod(iyear/100,4)

    if ( iy == 0 .and. ( iycen /= 0 .or. icent == 0 ) ) then
       ocleap = .true.
    else
       ocleap = .false.
    endif

    return
  end function ocleap

  !-----------------------------------------------------------------------------
  !> calendar, time advancing
  subroutine calendar_ssaft( &
       dseca ,               & !--- OUT
       dsec  ,               & !--- IN
       raftr ,               & !--- IN
       hunit                 & !--- IN
       )
    implicit none

    real(RP), intent(out) :: dseca
    real(RP), intent(in) :: dsec
    real(RP), intent(in) :: raftr
    character, intent(in) :: hunit *(*)

    integer :: idays, iyear, imonth, iday
    real(RP) :: rsec
    real(RP) :: ddtime
    !---------------------------------------------------------------------------

    if ( hunit(1:1) == 'y' .or. hunit(1:1) == 'Y' &
         .or. hunit(1:1) == 'mo' .or. hunit(1:2) == 'MO' ) then
       call calendar_ss2ds &
            ( idays , rsec  , &
            dsec            )
       call calendar_dd2ym &
            ( iyear , imonth, iday  , &
            idays                  )
       if ( hunit(1:1) == 'y' .or. hunit(1:1) == 'Y' ) then
          iyear  = iyear  + int(raftr)
       elseif ( hunit(1:2) == 'mo' .or. hunit(1:2) == 'MO' ) then
          imonth = imonth + int(raftr)
       endif
       call calendar_ym2dd &
            ( idays , &
            iyear , imonth, iday   )
       call calendar_ds2ss &
            ( dseca , &
            idays , rsec    )
    else
       call calendar_xx2ss &
            ( ddtime, &
            raftr , hunit , dsec  )
       dseca = dsec + ddtime
    endif

    return
  end subroutine calendar_ssaft

  !-----------------------------------------------------------------------------
  !> time step passed ?
  !> @return calendar_ointvl
  function   calendar_ointvl( &
       dtime ,                & !--- IN
       dtprev ,               & !--- IN
       dtorgn,                & !--- IN
       rintv ,                & !--- IN
       htunit                 & !--- IN
       )
    implicit none

    logical :: calendar_ointvl
    real(RP), intent(in) :: dtime
    real(RP), intent(in) :: dtprev
    real(RP), intent(in) :: dtorgn
    real(RP), intent(in) :: rintv
    character, intent(in) :: htunit *(*)

    real(RP) :: ddtime
    character(len=5) :: hunit
    integer :: iyear, imon, iday, iyearp, imonp, idayp
    integer :: iy, imo
    integer :: nmonyr, ndayyr, ndaymo
    real(RP) :: ry, rmo
    !---------------------------------------------------------------------------

    hunit = htunit

    if ( dtime == dtprev ) then
       calendar_ointvl = .true.
       return
    endif
    !
    calendar_ointvl = .false.
    call calendar_ss2ym &
         ( iyear , imon  , iday  , &
         dtime                   )
    call calendar_ss2ym &
         ( iyearp, imonp , idayp , &
         dtprev                  )
    call calendar_xx2ss &
         ( ddtime, &
         rintv , hunit , dtime )

    if ( dtime >= dtorgn ) then
       if      ( hunit(1:1) == 'y' .or. hunit(1:1) == 'Y' ) then
          call calendar_monyr( nmonyr, iyear )
          call calendar_dayyr( ndayyr, iyear )
          ry = real( iyear-iyearp, kind=RP )                           &
             + real( imon -imonp , kind=RP ) / real( nmonyr, kind=RP ) &
             + real( iday -idayp , kind=RP ) / real( ndayyr, kind=RP )
          if ( ry >= rintv ) then
             calendar_ointvl = .true.
          endif
       elseif ( hunit(1:2) == 'mo' .or. hunit(1:2) == 'MO' ) then
          imo = 0
          do iy = iyearp, iyear-1
             call calendar_monyr( nmonyr, iy )
             imo = imo + nmonyr
          enddo
          call calendar_daymo( ndaymo, iyear, imon )
          rmo = real( imon-imonp+imo, kind=RP ) &
              + real( iday-idayp    , kind=RP ) / real( ndaymo, kind=RP )
          if ( rmo >= rintv ) then
             calendar_ointvl = .true.
          endif
       elseif (      calendar_dgaus((dtime -dtorgn)/ddtime)  &
            > calendar_dgaus((dtprev-dtorgn)/ddtime) ) then
          calendar_ointvl = .true.
       endif
    endif

    return
  end function calendar_ointvl

  !-----------------------------------------------------------------------------
  !> dicard gaussian
  !> @return calendar_dgaus
  function calendar_dgaus( &
        dx                 & !--- IN
        )
    implicit none

    real(RP) :: calendar_dgaus
    real(RP), intent(in) :: dx
    !---------------------------------------------------------------------------

    calendar_dgaus = aint(dx) + aint(dx - aint(dx) + 1.0_RP) - 1.0_RP

  end function calendar_dgaus

end module mod_calendar
