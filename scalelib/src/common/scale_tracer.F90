!-------------------------------------------------------------------------------
!> module TRACER
!!
!! @par Description
!!          Tracer module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-04 (S.Nishizawa)   [new]
!! @li      2016-08-02 (S.Nishizawa)   [mod] add register
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_tracer
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
  public :: TRACER_regist
  public :: TRACER_inq_id

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: QA = 0

  integer, private, parameter :: QA_MAX = 1024

  real(RP),               public :: TRACER_CP  (QA_MAX)
  real(RP),               public :: TRACER_CV  (QA_MAX)
  real(RP),               public :: TRACER_R   (QA_MAX)
  real(RP),               public :: TRACER_MASS(QA_MAX)
  logical,                public :: TRACER_ADVC(QA_MAX)
  character(len=H_SHORT), public :: TRACER_NAME(QA_MAX)
  character(len=H_MID),   public :: TRACER_DESC(QA_MAX)
  character(len=H_SHORT), public :: TRACER_UNIT(QA_MAX)

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
  !> Regist tracer
  subroutine TRACER_regist(  &
       QS,                   &
       NQ, NAME, DESC, UNIT, &
       CV, CP, R, ADVC, MASS )
    use scale_prc, only: &
      PRC_abort
    implicit none

    integer,          intent(out)          :: QS
    integer,          intent(in)           :: NQ
    character(len=*), intent(in)           :: NAME(NQ)
    character(len=*), intent(in)           :: DESC(NQ)
    character(len=*), intent(in)           :: UNIT(NQ)
    real(RP),         intent(in), optional :: CV  (NQ)
    real(RP),         intent(in), optional :: CP  (NQ)
    real(RP),         intent(in), optional :: R   (NQ)
    logical,          intent(in), optional :: ADVC(NQ) !< if .true., the tracer is advected in the dynamical process. (default is .true.)
    logical,          intent(in), optional :: MASS(NQ) !< if .true., the tracer has mass. (default is .false.)

    real(RP) :: CV_  (NQ)
    real(RP) :: CP_  (NQ)
    real(RP) :: R_   (NQ)
    logical  :: ADVC_(NQ)
    logical  :: MASS_(NQ)

    character(len=24) :: NAME_trim

    integer  :: n
    !---------------------------------------------------------------------------

    if ( QA + NQ > QA_MAX ) then
       LOG_ERROR("TRACER_regist",*) 'total number of tracer must be less or equal to ', QA_MAX
       call PRC_abort
    end if

    if ( present(CV) ) then
       CV_(:) = CV(:)
    else
       CV_(:) = 0.0_RP
    end if

    if ( present(CP) ) then
       CP_(:) = CP(:)
    else
       CP_(:) = 0.0_RP
    end if

    if ( present(R) ) then
       R_(:) = R(:)
    else
       R_(:) = 0.0_RP
    end if

    if ( present(ADVC) ) then
       ADVC_(:) = ADVC(:)
    else
       ADVC_(:) = .true.
    end if

    if ( present(MASS) ) then
       MASS_(:) = MASS(:)
    else
       MASS_(:) = .false.
    end if

    LOG_NEWLINE
    do n = 1, NQ

       NAME_trim = trim(NAME(n))

       LOG_INFO("TRACER_regist",'(1x,A,I3,A,A,A,F6.1,A,F6.1,A,L1,A,L1)') &
                                      '] Register tracer : No.', QA+n,      &
                                                    ', NAME = ', NAME_trim, &
                                                      ', CV = ', CV_  (n),  &
                                                      ', CP = ', CP_  (n),  &
                                                    ', ADVC = ', ADVC_(n),  &
                                                    ', MASS = ', MASS_(n)

       TRACER_NAME(QA+n) = NAME (n)
       TRACER_DESC(QA+n) = DESC (n)
       TRACER_UNIT(QA+n) = UNIT (n)
       TRACER_CV  (QA+n) = CV_  (n)
       TRACER_CP  (QA+n) = CP_  (n)
       TRACER_R   (QA+n) = R_   (n)
       TRACER_ADVC(QA+n) = ADVC_(n)

       if ( MASS_(n) ) then
          TRACER_MASS(QA+n) = 1.0_RP
       else
          TRACER_MASS(QA+n) = 0.0_RP
       end if
    end do

    QS = QA + 1
    QA = QA + NQ

    return
  end subroutine TRACER_regist

  !-----------------------------------------------------------------------------
  !> Inquire tracer ID
  subroutine TRACER_inq_id( &
       NAME, &
       ID    )
    implicit none
    character(len=*), intent(in)  :: NAME
    integer,          intent(out) :: ID
    integer :: iq

    ID = -1
    do iq = 1, QA
       if ( NAME == TRACER_NAME(iq) ) then
          ID = iq
          exit
       end if
    end do

    return
  end subroutine TRACER_inq_id

end module scale_tracer
