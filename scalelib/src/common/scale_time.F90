!-------------------------------------------------------------------------------
!> module TIME
!!
!! @par Description
!!          general module for date/time
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_time
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
  public :: TIME_gettimelabel
  public :: TIME_time2label

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(DP), public :: TIME_DTSEC                 !< time interval of model                 [sec]

  real(DP), public :: TIME_DTSEC_ATMOS_DYN       !< time interval of dynamics              [sec]
  integer,  public :: TIME_NSTEP_ATMOS_DYN       !< small step of dynamics
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_CP    !< time interval of physics(cumulus     ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_MP    !< time interval of physics(microphysics) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_RD    !< time interval of physics(radiation   ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_SF    !< time interval of physics(surface flux) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_TB    !< time interval of physics(turbulence  ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_BL    !< time interval of physics(pbl         ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_CH    !< time interval of physics(chemistry   ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_AE    !< time interval of physics(aerosol     ) [sec]
  real(DP), public :: TIME_DTSEC_ATMOS_PHY_LT    !< time interval of physics(lightning   ) [sec]
  real(DP), public :: TIME_DTSEC_OCEAN           !< time interval of ocean step            [sec]
  real(DP), public :: TIME_DTSEC_LAND            !< time interval of land step             [sec]
  real(DP), public :: TIME_DTSEC_URBAN           !< time interval of urban step            [sec]
  real(DP), public :: TIME_DTSEC_DA              !< time interval of data-assim. step      [sec]
  real(DP), public :: TIME_DTSEC_WALLCLOCK_CHECK !< time interval of wallclock terminator  [sec]

  integer,  public :: TIME_DSTEP_ATMOS_DYN       !< step interval of dynamics
  integer,  public :: TIME_DSTEP_ATMOS_PHY_CP    !< step interval of physics(cumulus     )
  integer,  public :: TIME_DSTEP_ATMOS_PHY_MP    !< step interval of physics(microphysics)
  integer,  public :: TIME_DSTEP_ATMOS_PHY_RD    !< step interval of physics(radiation   )
  integer,  public :: TIME_DSTEP_ATMOS_PHY_SF    !< step interval of physics(surface flux)
  integer,  public :: TIME_DSTEP_ATMOS_PHY_TB    !< step interval of physics(turbulence  )
  integer,  public :: TIME_DSTEP_ATMOS_PHY_BL    !< step interval of physics(pbl         )
  integer,  public :: TIME_DSTEP_ATMOS_PHY_CH    !< step interval of physics(chemistry   )
  integer,  public :: TIME_DSTEP_ATMOS_PHY_AE    !< step interval of physics(aerosol     )
  integer,  public :: TIME_DSTEP_ATMOS_PHY_LT    !< step interval of physics(lightning   )
  integer,  public :: TIME_DSTEP_OCEAN           !< step interval of ocean step
  integer,  public :: TIME_DSTEP_LAND            !< step interval of land step
  integer,  public :: TIME_DSTEP_URBAN           !< step interval of urban step
  integer,  public :: TIME_DSTEP_DA              !< step interval of data-assim. step
  integer,  public :: TIME_DSTEP_WALLCLOCK_CHECK !< step interval of wallclock terminator

  integer,  public :: TIME_NOWDATE(6)           !< current time [YYYY MM DD HH MM SS]
  integer,  public :: TIME_NOWDAY               !< absolute day of current time [day]
  real(DP), public :: TIME_NOWSEC               !< subday part  of current time [sec]
  real(DP), public :: TIME_NOWSUBSEC            !< subsecond part of current time [sec]
  real(DP), public :: TIME_NOWDAYSEC            !< second of current time [sec]
  integer,  public :: TIME_NOWSTEP              !< current step [number]
  integer,  public :: TIME_NSTEP                !< total steps [number]

  integer,  public :: TIME_OFFSET_YEAR          !< time offset [year]
  real(DP), public :: TIME_STARTDAYSEC          !< second of start time [sec]

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
  !> generate time label
  subroutine TIME_gettimelabel( &
       timelabel )
    implicit none

    character(len=19), intent(out) :: timelabel
    !---------------------------------------------------------------------------

    ! YYYYMMDD-hhmmss.sss
    call TIME_time2label( TIME_NOWDATE(:), TIME_NOWSUBSEC, & ! [IN]
                          timelabel                        ) ! [OUT]

    return
  end subroutine TIME_gettimelabel

  !> generate time label
  subroutine TIME_time2label( &
       date, subsec, &
       timelabel )
    implicit none

    integer,  intent(in) :: date(6)
    real(DP), intent(in) :: subsec

    character(len=*), intent(out) :: timelabel
    !---------------------------------------------------------------------------

    ! YYYYMMDD-hhmmss.sss
    write(timelabel,'(I4.4,I2.2,I2.2,A1,I2.2,I2.2,I2.2,A1,I3.3)') &
         date(1:3), '-', date(4:6), '.', int(subsec*1000.0_DP)

    return
  end subroutine TIME_time2label

end module scale_time
