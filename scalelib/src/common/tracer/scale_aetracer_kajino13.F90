!-------------------------------------------------------------------------------
!> module TRACER / kajino13
!!
!! @par Description
!!          Tracer kajino13 module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-03-27 (Y.Sato)   [new] imported from scale_tracer_suzuki10.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_aetracer_kajino13
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
  public :: AETRACER_kajino13_setup

  include "inc_aetracer_kajino13.h"

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine AETRACER_kajino13_setup
    use scale_process, only: &
      PRC_MPIstop
    implicit none

    integer, allocatable :: aero_idx(:,:,:,:)
    integer :: n_kap_max, n_siz_max, ncat_max
    real(RP),allocatable :: NASIZ(:), NAKAP(:)

    NAMELIST / PARAM_TRACER_KAJINO13 / &
       AE_CTG, &
       NASIZ,  &
       NAKAP 

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ READ NUMBER of TRACER for Aerosol'

    ncat_max = max( IC_MIX, IC_SEA, IC_DUS )
    allocate( NASIZ(ncat_max) )
    allocate( NAKAP(ncat_max) )

    NASIZ(:) = 64
    NAKAP(:) = 1

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TRACER_KAJINO13,iostat=ierr)

    if( ierr < 0 ) then !--- missing
     if( IO_L ) write(IO_FID_LOG,*)  '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
     write(*,*) 'xxx Not appropriate names in namelist PARAM_TRACER_KAJINO13, Check!'
     call PRC_MPIstop
    end if

    if( IO_L ) write(IO_FID_LOG,nml=PARAM_TRACER_KAJINO13)

    if( AE_CTG > ncat_max ) then
     write(*,*) 'xxx AE_CTG should be smaller than', ncat_max+1, 'stop'
     call PRC_MPIstop
    endif

    allocate( NSIZ(AE_CTG) )
    allocate( NKAP(AE_CTG) )

    NKAP(1:AE_CTG) = NAKAP(1:AE_CTG)
    NSIZ(1:AE_CTG) = NASIZ(1:AE_CTG)

    if( maxval( NKAP ) /= 1 .or. minval( NKAP ) /= 1 ) then
     write(*,*) 'xxx NKAP(:) /= 1 is not supported now, Stop!'
     call PRC_MPIstop
    end if

    QA_AE = 0
    do ia0 = 1, N_ATR
    do ic = 1, AE_CTG
    do ik = 1, NKAP(ic)
    do is0 = 1, NSIZ(ic)
      QA_AE   = QA_AE + 1
    enddo
    enddo
    enddo
    enddo
    QA_AE = QA_AE + GAS_CTG
    allocate( AQ_AE_NAME(QA_AE) )
    allocate( AQ_AE_DESC(QA_AE) )
    allocate( AQ_AE_UNIT(QA_AE) )
!    QAES = QA_MP + 1
!    QAEE = QA_MP + QA_AE

    n_siz_max = 0
    n_kap_max = 0
    do ic = 1, AE_CTG
      n_siz_max = max(n_siz_max, NSIZ(ic))
      n_kap_max = max(n_kap_max, NKAP(ic))
    enddo

    allocate( aero_idx(N_ATR,AE_CTG,n_kap_max,n_siz_max) )
    m = 0
    do ia0 = 1, N_ATR
    do ic = 1, AE_CTG
    do ik = 1, NKAP(ic)
    do is0 = 1, NSIZ(ic)
      m = m+1
      aero_idx(ia0,ic,ik,is0) = m 
    enddo
    enddo
    enddo
    enddo

    !-----------------------------------------------------------------------------
    !
    !++ calculate each category and aerosol
    !
    !-----------------------------------------------------------------------------
    do n = 1, QA_AE-GAS_CTG
     write(AQ_AE_UNIT(n),'(a)')  '#/kg'
    enddo
    ic = QA_AE-GAS_CTG+IG_H2SO4
    write(AQ_AE_UNIT(ic),'(a)')  'kg/kg'
    ic = QA_AE-GAS_CTG+IG_CGAS
    write(AQ_AE_UNIT(ic),'(a)')  'kg/kg'

!    do ic = 1, n_ctg       !aerosol category
!    do ik = 1, n_kap(ic)   !kappa bin
!    do is = 1, n_siz(ic)  
    do ia0 = 1, N_ATR
    do ic = 1, AE_CTG       !aerosol category
    do ik = 1, NKAP(ic)   !kappa bin
    do is0 = 1, NSIZ(ic) 
!     write(*,*) aero_idx(ia,ic,ik,is), ia,ic,ik,is, size(AQ_AE_NAME)
!     write(AQ_AE_NAME(aero_idx(ia,ic,ik,is)),'(4(a,i0))') 'AEROSOL(atl)', ia, 'cat', ic, 'kappa', ik, 'size', is
     write(AQ_AE_NAME(aero_idx(ia0,ic,ik,is0)),'((a,i0))') 'AEROSOL(atl)', aero_idx(ia0,ic,ik,is0)
    enddo
    enddo
    enddo
    enddo
    ic = QA_AE-GAS_CTG+IG_H2SO4
    write(AQ_AE_NAME(ic),'(a)') 'H2SO4 Gas'
    ic = QA_AE-GAS_CTG+IG_CGAS
    write(AQ_AE_NAME(ic),'(a)') 'Condensable GAS'
    

    do ia0 = 1, N_ATR
    do ic = 1, AE_CTG     !aerosol category
    do ik = 1, NKAP(ic)   !kappa bin
    do is0 = 1, NSIZ(ic) 
     write(AQ_AE_DESC(aero_idx(ia0,ic,ik,is0)),'(a,i0,a5,i0,a4,i0)') 'Mixing radio of aerosol(cat)', ic, 'kappa', ik, 'size', is0
    enddo
    enddo
    enddo
    enddo
    ic = QA_AE-GAS_CTG+IG_H2SO4
    write(AQ_AE_DESC(ic),'(a)') 'Mixing ratio of H2SO4 Gas'
    ic = QA_AE-GAS_CTG+IG_CGAS
    write(AQ_AE_DESC(ic),'(a)') 'Mixing ratio of Condensable GAS'

    deallocate(aero_idx)

    return
  end subroutine AETRACER_kajino13_setup

end module scale_aetracer_kajino13
