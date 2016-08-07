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
module scale_tracer
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TRACER_regist

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: QA = 0

  integer, private, parameter :: QA_MAX = 1024
  real(RP), public :: TRACER_CP(QA_MAX)
  real(RP), public :: TRACER_CV(QA_MAX)
  real(RP), public :: TRACER_R(QA_MAX)
  real(RP), public :: TRACER_MASS(QA_MAX)
  logical,  public :: TRACER_ADVC(QA_MAX)
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
    use scale_process, only: &
      PRC_MPIstop
    implicit none

    integer,  intent(out) :: QS
    integer,  intent(in)  :: NQ
    character(len=*), intent(in) :: NAME(NQ)
    character(len=*), intent(in) :: DESC(NQ)
    character(len=*), intent(in) :: UNIT(NQ)
    real(RP), intent(in), optional :: CV(NQ)
    real(RP), intent(in), optional :: CP(NQ)
    real(RP), intent(in), optional :: R(NQ)
    logical,  intent(in), optional :: ADVC(NQ) !< if .true., the tracer is advected in the dynamical process. (default is .true.)
    logical,  intent(in), optional :: MASS(NQ) !< if .true., the tracer has mass. (default is .false.)

    real(RP) :: CV_(NQ)
    real(RP) :: CP_(NQ)
    real(RP) :: R_(NQ)
    logical :: ADVC_(NQ)
    logical :: MASS_(NQ)
    integer :: n
    !---------------------------------------------------------------------------

    if ( QA + NQ > QA_MAX ) then
       write(*,*) 'xxx total number of tracer must be less or equal to ', QA_MAX
       call PRC_MPIstop
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

    do n = 1, NQ
       if( IO_L ) write(IO_FID_LOG,'(a,i4,a,a,a,f6.1,a,f6.1,a,l,a,l)') &
            'register tracer: ', &
            QA+n, ', NAME=', trim(NAME(n)), ', CV=', CV_(n), ', CP=', CP_(n), ', ADVC=', ADVC_(n), ', MASS=', MASS_(n)
       TRACER_NAME(QA+n) = NAME(n)
       TRACER_DESC(QA+n) = DESC(n)
       TRACER_UNIT(QA+n) = UNIT(n)
       TRACER_CV(QA+n) = CV_(n)
       TRACER_CP(QA+n) = CP_(n)
       TRACER_R (QA+n) = R_(n)
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

end module scale_tracer
