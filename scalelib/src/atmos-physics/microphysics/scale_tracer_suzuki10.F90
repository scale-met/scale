!-------------------------------------------------------------------------------
!> module TRACER / suzuki10
!!
!! @par Description
!!          Tracer suzuki10 module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-21 (Y.Sato)   [new] imported from inc_tracer_suzuki10.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_tracer_suzuki10
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
  public :: TRACER_suzuki10_setup

  include "inc_tracer_suzuki10.h"

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine TRACER_suzuki10_setup
    use scale_process, only: &
      PRC_MPIstop
    implicit none

    NAMELIST / PARAM_BIN / &
       nbin, &
       nccn, &
       ICEFLG, &
       kphase

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ READ BIN NUMBER'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BIN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
     if( IO_L ) write(IO_FID_LOG,*)  '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
     write(*,*) 'xxx Not appropriate names in namelist PARAM_BIN, Check!'
     call PRC_MPIstop
    end if

    if( IO_L ) write(IO_FID_LOG,nml=PARAM_BIN)

    if( ICEFLG == 0 ) then
     nspc = 1
    elseif( ICEFLG == 1 ) then
     nspc = 7
    else
     write(*,*) "ICEFLG should be 0(warm rain) or 1(mixed rain) check!!"
     call PRC_MPIstop
    endif

    !-- setup QA ...
    QA   = I_QV+nbin*nspc+nccn
    QQA  = I_QV+nbin*nspc
    QQS  = I_QV
    QQE  = I_QV+nbin*nspc
    QWS =  I_QV+1
    QWE =  I_QV+nbin
    if( ICEFLG == 0 ) then
      QIS =  0
      QIE = -1
    elseif( ICEFLG == 1 ) then
      QIS =  (I_QV+nbin+1)*ICEFLG
      QIE =  (I_QV+nbin*nspc)*ICEFLG
    endif

    allocate( AQ_NAME(QA) )
    allocate( AQ_DESC(QA) )
    allocate( AQ_UNIT(QA) )


    !-----------------------------------------------------------------------------
    !
    !++ calculate each category and aerosol
    !
    !-----------------------------------------------------------------------------
    do n = 1, QA
     write(AQ_UNIT(n),'(a)')  'kg/kg'
    enddo

    write(AQ_NAME(I_QV),'(a)') 'QV'
    do m = 1, nspc
    do n = 1, nbin
     write(AQ_NAME(I_QV+nbin*(m-1)+n),'(a,i0)') trim(namspc(m)), n
    enddo
    enddo

    do n = 1, nccn
     write(AQ_NAME(I_QV+nbin*nspc+n),'(a,i0)') trim(namspc(8)), n
    enddo

    write(AQ_DESC(I_QV),'(a)')  'Water Vapor mixing ratio'
    do m = 1, nspc
    do n = 1, nbin
     write(AQ_DESC(I_QV+nbin*(m-1)+n),'(a,i0)') trim(lnamspc(m)), n
    enddo
    enddo

    do n = 1, nccn
     write(AQ_DESC(I_QV+nbin*nspc+n),'(a,i0)') trim(lnamspc(8)), n
    enddo

    return
  end subroutine TRACER_suzuki10_setup

end module scale_tracer_suzuki10
