#include "scalelib.h"
module scale_spnudge
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_CZUY => ATMOS_GRID_CARTESC_REAL_CZUY, &
       REAL_CZXV => ATMOS_GRID_CARTESC_REAL_CZXV
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  logical, public :: SPNUDGE_uv = .false.
  logical, public :: SPNUDGE_uv_divfree = .false.
  integer, public :: SPNUDGE_uv_lm = 3
  integer, public :: SPNUDGE_uv_mm = 3 
  real(RP), public :: SPNUDGE_uv_tau

  logical, public :: SPNUDGE_pt = .false.
  integer, public :: SPNUDGE_pt_lm = 3
  integer, public :: SPNUDGE_pt_mm = 3 
  real(RP), public :: SPNUDGE_pt_tau

  real(RP), allocatable, public :: SPNUDGE_u_alpha(:,:,:)
  real(RP), allocatable, public :: SPNUDGE_v_alpha(:,:,:)
  real(RP), allocatable, public :: SPNUDGE_pt_alpha(:,:,:)

  real(RP), public :: SPNUDGE_P1
  real(RP), public :: SPNUDGE_P2

  
  public :: SPNUDGE_setup
  
  contains
  
  subroutine SPNUDGE_setup(KA, KS, KE, IA, IS, IE, JA, JS, JE)
    implicit none
    
    namelist /PARAM_SPNUDGE/ &
      SPNUDGE_uv, &
      SPNUDGE_uv_divfree, &
      SPNUDGE_uv_lm, &
      SPNUDGE_uv_mm, &
      SPNUDGE_uv_tau, &
      SPNUDGE_pt, &
      SPNUDGE_pt_lm, &
      SPNUDGE_pt_mm, &
      SPNUDGE_pt_tau, &
      SPNUDGE_P1, &
      SPNUDGE_P2

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    
    real(RP) :: uv_alpha, pt_alpha
    integer :: k, i, j
    integer :: ierr
    
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SPNUDGE_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SPNUDGE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SPNUDGE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SPNUDGE_setup",*) 'Not appropriate names in namelist PARAM_SPNUDGE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_SPNUDGE)

    allocate( SPNUDGE_u_alpha(KA,IA,JA) )
    allocate( SPNUDGE_v_alpha(KA,IA,JA) )
    allocate( SPNUDGE_pt_alpha(KA,IA,JA) )
    
    if( SPNUDGE_uv_tau <= 0.0_RP ) then
        uv_alpha = 0
    else
        uv_alpha = 1.0_RP / SPNUDGE_uv_tau
    endif

    if( SPNUDGE_pt_tau <= 0.0_RP ) then
        pt_alpha = 0
    else
        pt_alpha = 1.0_RP / SPNUDGE_pt_tau
    endif

    do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             if( REAL_CZUY(k,i,j) <= SPNUDGE_P1 ) then
                SPNUDGE_u_alpha(k,i,j) = 0
             elseif( REAL_CZUY(k,i,j) <= SPNUDGE_P2 ) then
                SPNUDGE_u_alpha(k,i,j) = ( REAL_CZUY(k,i,j) - SPNUDGE_P1 ) / (SPNUDGE_P2 - SPNUDGE_P1 )*uv_alpha
             else
                SPNUDGE_u_alpha(k,i,j) = uv_alpha
             endif
             
          enddo
       enddo
     enddo

    do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             if( REAL_CZXV(k,i,j) <= SPNUDGE_P1 ) then
                SPNUDGE_v_alpha(k,i,j) = 0
             elseif( REAL_CZXV(k,i,j) <= SPNUDGE_P2 ) then
                SPNUDGE_v_alpha(k,i,j) = ( REAL_CZXV(k,i,j) - SPNUDGE_P1 ) / (SPNUDGE_P2 - SPNUDGE_P1 )*uv_alpha
             else
                SPNUDGE_v_alpha(k,i,j) = uv_alpha
             endif
             
          enddo
       enddo
    enddo     

    do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             if( REAL_CZ(k,i,j) <= SPNUDGE_P1 ) then
                SPNUDGE_pt_alpha(k,i,j) = 0
             elseif( REAL_CZ(k,i,j) <= SPNUDGE_P2 ) then
                SPNUDGE_pt_alpha(k,i,j) = ( REAL_CZ(k,i,j) - SPNUDGE_P1 ) / (SPNUDGE_P2 - SPNUDGE_P1 )*pt_alpha
             else
                SPNUDGE_pt_alpha(k,i,j) = pt_alpha
             endif
             
          enddo
       enddo
    enddo
    
  end subroutine SPNUDGE_setup
  
end module scale_spnudge
