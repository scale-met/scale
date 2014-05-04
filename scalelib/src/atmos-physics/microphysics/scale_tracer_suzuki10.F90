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

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !------------------------------------------------------------
  !
  !++ scale-les grid parameter for inc_tracer_hbinf_h33_c20.f90
  !
  !------------------------------------------------------------
!  include "setup_bin.h"
  integer, public :: nbin = 33
  integer, public :: nspc = 7
  integer, public :: nccn = 20
  integer, public :: kphase = 0
  integer, public :: ICEFLG = 1
  integer, public :: I_QV =  1
  integer, public :: QA
  integer, public :: QQA
  integer, public :: QQS
  integer, public :: QQE
  integer, public :: QWS
  integer, public :: QWE
  integer, public :: QIS
  integer, public :: QIE
!  integer, public, parameter :: QA   = I_QV+nbin*nspc+nccn
!  integer, public, parameter :: QQA  = I_QV+nbin*nspc
!  integer, public, parameter :: QQS  = I_QV
!  integer, public, parameter :: QQE  = I_QV+nbin*nspc
!  integer, public, parameter :: QWS =  I_QV+1
!  integer, public, parameter :: QWE =  I_QV+nbin
!  integer, public, parameter :: QIS =  (I_QV+nbin+1)*ICEFLG
!  integer, public, parameter :: QIE =  (I_QV+nbin*nspc)*ICEFLG

  integer, public, parameter :: I_QC =  2
  integer, public, parameter :: I_QR =  3
  integer, public, parameter :: I_QI =  4
  integer, public, parameter :: I_QS =  5
  integer, public, parameter :: I_QG =  6
  integer, public, parameter :: I_NC =  8
  integer, public, parameter :: I_NR =  9
  integer, public, parameter :: I_NI =  10
  integer, public, parameter :: I_NS =  11
  integer, public, parameter :: I_NG =  12

  character(len=H_SHORT), public, save, allocatable :: AQ_NAME(:)
  character(len=H_MID)  , public, save, allocatable :: AQ_DESC(:)
  character(len=H_SHORT), public, save, allocatable :: AQ_UNIT(:)
!  character(len=16), public, save :: AQ_NAME(QA)
!  character(len=64), public, save :: AQ_DESC(QA)
!  character(len=16), public, save :: AQ_UNIT(QA)

  character(len=3)  :: namspc(8) =(/'Qcl','Qic','Qip','Qid','Qis','Qig','Qih','Qae'/)
  character(len=27) :: lnamspc(8) = &
                             (/'Mixing ratio of cloud   bin', &
                               'Mixing ratio of colum   bin', &
                               'Mixing ratio of plate   bin', &
                               'Mixing ratio of dendrit bin', &
                               'Mixing ratio of snow    bin', &
                               'Mixing ratio of graupel bin', &
                               'Mixing ratio of hail    bin', &
                               'Mixing ratio of aerosol bin' /)
  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_binf+AE_dummy+RD_mstrnx)
  !
  !-----------------------------------------------------------------------------

  integer, public, parameter :: MP_QA = 7 ! number of hydrometeor tracer
  integer, public, parameter :: I_mp_QC = 1
  integer, public, parameter :: I_mp_QP = 2
  integer, public, parameter :: I_mp_QCL= 3
  integer, public, parameter :: I_mp_QD = 4
  integer, public, parameter :: I_mp_QS = 5
  integer, public, parameter :: I_mp_QG = 6
  integer, public, parameter :: I_mp_QH = 7

  integer, public, save :: I_MP2ALL(MP_QA)  !--- dummy for Bin model
  data I_MP2ALL / I_mp_QC,   & ! start of Cloud bin
                  I_mp_QP,   & ! start of plate bin
                  I_mp_QCL,  & ! start of columner bin
                  I_mp_QD,   & ! start of dendrite bin
                  I_mp_QS,   & ! start of snow bin
                  I_mp_QG,   & ! start of graupel bin
                  I_mp_QH    / ! start of hail bin

  integer, public, save :: I_MP2RD(MP_QA)
  data I_MP2RD  / 1,    & ! I_mp_QC => MSTRN_nptype=1: water cloud
                  2,    & ! I_mp_QP => MSTRN_nptype=2: ice cloud (plate)
                  2,    & ! I_mp_QCL=> MSTRN_nptype=2: ice cloud (columner)
                  2,    & ! I_mp_QD => MSTRN_nptype=2: ice cloud (dendrite)
                  2,    & ! I_mp_QS => MSTRN_nptype=2: ice cloud (snow)
                  2,    & ! I_mp_QG => MSTRN_nptype=2: ice cloud (graupel)
                  2     / ! I_mp_QH => MSTRN_nptype=2: ice cloud (hail)

  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_dummy = 1

  integer, public, save :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, public, save :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust
  integer :: m, n, ierr

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
