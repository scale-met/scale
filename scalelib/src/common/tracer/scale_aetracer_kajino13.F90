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
    character(len=H_SHORT) :: attribute, catego, aunit

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
!    do ia0 = 1, N_ATR
    do ic = 1, AE_CTG
    do ik = 1, NKAP(ic)
    do is0 = 1, NSIZ(ic)
    do ia0 = 1, N_ATR
      QA_AE   = QA_AE + 1
    enddo
    enddo
    enddo
    enddo
    QA_AE = QA_AE + GAS_CTG
    allocate( AQ_AE_NAME(QA_AE) )
    allocate( AQ_AE_DESC(QA_AE) )
    allocate( AQ_AE_UNIT(QA_AE) )

    n_siz_max = 0
    n_kap_max = 0
    do ic = 1, AE_CTG
      n_siz_max = max(n_siz_max, NSIZ(ic))
      n_kap_max = max(n_kap_max, NKAP(ic))
    enddo

    allocate( aero_idx(N_ATR,AE_CTG,n_kap_max,n_siz_max) )
    m = 0
!    do ia0 = 1, N_ATR
    do ic = 1, AE_CTG
    do ik = 1, NKAP(ic)
    do is0 = 1, NSIZ(ic)
    do ia0 = 1, N_ATR
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
    ic = QA_AE-GAS_CTG+IG_H2SO4
    write(AQ_AE_UNIT(ic),'(a)')  'kg/kg'
    ic = QA_AE-GAS_CTG+IG_CGAS
    write(AQ_AE_UNIT(ic),'(a)')  'kg/kg'

!    do ia0 = 1, N_ATR
!      if( ia0 == 1 ) then
!         write(attribute,'(a)') "Number"
!         write(aunit,'(a)') "#/kg"
!      elseif( ia0 == 2 ) then
!         write(attribute,'(a)') "Section"
!         write(aunit,'(a)') "m2/kg"
!      elseif( ia0 == 3 ) then
!         write(attribute,'(a)') "Volume"
!         write(aunit,'(a)') "m3/kg"
!      elseif( ia0 == 4 ) then
!         write(attribute,'(a)') "Mass"
!         write(aunit,'(a)') "kg/kg"
!      elseif( ia0 == 5 ) then
!         write(attribute,'(a)') "kpXmass"
!         write(aunit,'(a)') "kg/kg"
!      endif
    do ic = 1, AE_CTG       !aerosol category
    do ik = 1, NKAP(ic)   !kappa bin
    do is0 = 1, NSIZ(ic) 
    do ia0 = 1, N_ATR
      if( ia0 == 1 ) then
         write(attribute,'(a)') "Number"
         write(aunit,'(a)') "/kg"
      elseif( ia0 == 2 ) then
         write(attribute,'(a)') "Section"
         write(aunit,'(a)') "m2/kg"
      elseif( ia0 == 3 ) then
         write(attribute,'(a)') "Volume"
         write(aunit,'(a)') "m3/kg"
      elseif( ia0 == 4 ) then
         write(attribute,'(a)') "Mass"
         write(aunit,'(a)') "kg/kg"
      elseif( ia0 == 5 ) then
         write(attribute,'(a)') "kXm"
         write(aunit,'(a)') "kg/kg"
      endif
     if( ic == IC_MIX ) then
      write(catego,'(a)') "Sulf_"
     elseif( ic == IC_SEA ) then
      write(catego,'(a)') "Salt_"
     elseif( ic == IC_DUS ) then
      write(catego,'(a)') "Dust_"
     endif
     write(AQ_AE_UNIT(aero_idx(ia0,ic,ik,is0)),'(a)')  trim(aunit)
     write(AQ_AE_NAME(aero_idx(ia0,ic,ik,is0)),'(a,a,i0)') trim(catego), trim(attribute), is0
     write(AQ_AE_DESC(aero_idx(ia0,ic,ik,is0)),'(a,a,a,i0)') trim(attribute), ' mixing radio of ', trim(catego), is0
    enddo
    enddo
    enddo
    enddo
    ic = QA_AE-GAS_CTG+IG_H2SO4
    write(AQ_AE_NAME(ic),'(a)') 'H2SO4_Gas'
    ic = QA_AE-GAS_CTG+IG_CGAS
    write(AQ_AE_NAME(ic),'(a)') 'Condensable_GAS'
    
    ic = QA_AE-GAS_CTG+IG_H2SO4
    write(AQ_AE_DESC(ic),'(a)') 'Mixing ratio of H2SO4 Gas'
    ic = QA_AE-GAS_CTG+IG_CGAS
    write(AQ_AE_DESC(ic),'(a)') 'Mixing ratio of Condensable GAS'

    deallocate(aero_idx)

    return
  end subroutine AETRACER_kajino13_setup

end module scale_aetracer_kajino13
