#include "scalelib.h"
module scale_spnudge
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SPNUDGE_setup
  public :: SPNUDGE_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  logical, public :: SPNUDGE_uv         = .false.
  logical, public :: SPNUDGE_uv_divfree = .false.
  integer, public :: SPNUDGE_uv_lm      = 3
  integer, public :: SPNUDGE_uv_mm      = 3

  logical, public :: SPNUDGE_pt      = .false.
  integer, public :: SPNUDGE_pt_lm   = 3
  integer, public :: SPNUDGE_pt_mm   = 3

  logical, public :: SPNUDGE_qv      = .false.
  integer, public :: SPNUDGE_qv_lm   = 3
  integer, public :: SPNUDGE_qv_mm   = 3

  real(RP), allocatable, public :: SPNUDGE_u_alpha(:,:,:)
  real(RP), allocatable, public :: SPNUDGE_v_alpha(:,:,:)
  real(RP), allocatable, public :: SPNUDGE_pt_alpha(:,:,:)
  real(RP), allocatable, public :: SPNUDGE_qv_alpha(:,:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !

contains

  subroutine SPNUDGE_setup(KA, KS, KE, IA, IS, IE, JA, JS, JE)
    use scale_dft, only: &
       DFT_setup
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_CZUY => ATMOS_GRID_CARTESC_REAL_CZUY, &
       REAL_CZXV => ATMOS_GRID_CARTESC_REAL_CZXV
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP) :: SPNUDGE_uv_tau    = 0.0_RP
    real(RP) :: SPNUDGE_pt_tau = 0.0_RP
    real(RP) :: SPNUDGE_qv_tau = 0.0_RP

    real(RP) :: SPNUDGE_level1 = 0.0_RP  !> alpha = 0 for z < SPNUDGE_level1
    real(RP) :: SPNUDGE_level2 = 0.0_RP
    real(RP) :: SPNUDGE_level3 = 1E10_RP !> alpha = alpha for SPNUDGE_level2 <= z < SPNUDGE_level3
    real(RP) :: SPNUDGE_level4 = 1E10_RP !> alpha = 0 for z => SPNUDGE_level4

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
      SPNUDGE_qv, &
      SPNUDGE_qv_lm, &
      SPNUDGE_qv_mm, &
      SPNUDGE_qv_tau, &
      SPNUDGE_level1, &
      SPNUDGE_level2, &
      SPNUDGE_level3, &
      SPNUDGE_level4

    real(RP) :: uv_alpha, pt_alpha, qv_alpha


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

    if ( SPNUDGE_level1 > SPNUDGE_level2 ) then
       LOG_ERROR("SPNUDGE_setup",*) 'SPNUDGE_level1 must be lowere or equal to SPNUDGE_level2'
       call PRC_abort
    end if
    if ( SPNUDGE_level2 > SPNUDGE_level3 ) then
       LOG_ERROR("SPNUDGE_setup",*) 'SPNUDGE_level2 must be lowere or equal to SPNUDGE_level3'
       call PRC_abort
    end if
    if ( SPNUDGE_level3 > SPNUDGE_level4 ) then
       LOG_ERROR("SPNUDGE_setup",*) 'SPNUDGE_level3 must be lowere or equal to SPNUDGE_level4'
       call PRC_abort
    end if

    if ( SPNUDGE_uv .or. SPNUDGE_pt .or. SPNUDGE_qv ) then
       LOG_WARN("SPNUDGE_setup",*) 'Spectrul nudging is still experimental'
    end if

    allocate( SPNUDGE_u_alpha(KA,IA,JA) )
    allocate( SPNUDGE_v_alpha(KA,IA,JA) )
    allocate( SPNUDGE_pt_alpha(KA,IA,JA) )
    allocate( SPNUDGE_qv_alpha(KA,IA,JA) )

    if ( SPNUDGE_uv .and. SPNUDGE_uv_tau > 0.0_RP ) then
       uv_alpha = 1.0_RP / SPNUDGE_uv_tau

       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if ( REAL_CZUY(k,i,j) < SPNUDGE_level1 ) then
             SPNUDGE_u_alpha(k,i,j) = 0.0_RP
          elseif ( REAL_CZUY(k,i,j) < SPNUDGE_level2 ) then
             SPNUDGE_u_alpha(k,i,j) = ( REAL_CZUY(k,i,j) - SPNUDGE_level1 ) / ( SPNUDGE_level2 - SPNUDGE_level1 ) * uv_alpha
          elseif ( REAL_CZUY(k,i,j) < SPNUDGE_level3 ) then
             SPNUDGE_u_alpha(k,i,j) = uv_alpha
          elseif ( REAL_CZUY(k,i,j) < SPNUDGE_level4 ) then
             SPNUDGE_u_alpha(k,i,j) = ( REAL_CZUY(k,i,j) - SPNUDGE_level3 ) / ( SPNUDGE_level4 - SPNUDGE_level3 ) * uv_alpha
          else
             SPNUDGE_u_alpha(k,i,j) = 0.0_RP
          endif
       enddo
       enddo
       enddo

       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if ( REAL_CZXV(k,i,j) < SPNUDGE_level1 ) then
             SPNUDGE_v_alpha(k,i,j) = 0.0_RP
          elseif ( REAL_CZXV(k,i,j) < SPNUDGE_level2 ) then
             SPNUDGE_v_alpha(k,i,j) = ( REAL_CZXV(k,i,j) - SPNUDGE_level1 ) / ( SPNUDGE_level2 - SPNUDGE_level1 ) * uv_alpha
          elseif ( REAL_CZXV(k,i,j) < SPNUDGE_level3 ) then
             SPNUDGE_v_alpha(k,i,j) = uv_alpha
          elseif ( REAL_CZXV(k,i,j) < SPNUDGE_level4 ) then
             SPNUDGE_v_alpha(k,i,j) = ( REAL_CZXV(k,i,j) - SPNUDGE_level3 ) / ( SPNUDGE_level4 - SPNUDGE_level3 ) * uv_alpha
          else
             SPNUDGE_v_alpha(k,i,j) = 0.0_RP
          endif
       enddo
       enddo
       enddo

    else

       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          SPNUDGE_u_alpha(k,i,j) = 0.0_RP
          SPNUDGE_v_alpha(k,i,j) = 0.0_RP
       end do
       end do
       end do

    endif

    if ( SPNUDGE_pt .and. SPNUDGE_pt_tau > 0.0_RP ) then
        pt_alpha = 1.0_RP / SPNUDGE_pt_tau

        !$omp parallel do
        do j = JS, JE
        do i = IS, IE
        do k = KS, KE
           if ( REAL_CZ(k,i,j) < SPNUDGE_level1 ) then
              SPNUDGE_pt_alpha(k,i,j) = 0.0_RP
           elseif ( REAL_CZ(k,i,j) < SPNUDGE_level2 ) then
              SPNUDGE_pt_alpha(k,i,j) = ( REAL_CZ(k,i,j) - SPNUDGE_level1 ) / ( SPNUDGE_level2 - SPNUDGE_level1 ) * pt_alpha
           elseif ( REAL_CZ(k,i,j) < SPNUDGE_level3 ) then
              SPNUDGE_pt_alpha(k,i,j) = pt_alpha
           elseif ( REAL_CZ(k,i,j) < SPNUDGE_level4 ) then
              SPNUDGE_pt_alpha(k,i,j) = ( REAL_CZ(k,i,j) - SPNUDGE_level3 ) / ( SPNUDGE_level4 - SPNUDGE_level3 ) * pt_alpha
           else
              SPNUDGE_pt_alpha(k,i,j) = 0.0_RP
           endif
        enddo
        enddo
        enddo

     else

        !$omp parallel do
        do j = JS, JE
        do i = IS, IE
        do k = KS, KE
           SPNUDGE_pt_alpha(k,i,j) = 0.0_RP
        end do
        end do
        end do

    endif

    if ( SPNUDGE_qv .and. SPNUDGE_qv_tau > 0.0_RP ) then
        qv_alpha = 1.0_RP / SPNUDGE_qv_tau

        !$omp parallel do
        do j = JS, JE
        do i = IS, IE
        do k = KS, KE
           if ( REAL_CZ(k,i,j) < SPNUDGE_level1 ) then
              SPNUDGE_qv_alpha(k,i,j) = 0.0_RP
           elseif ( REAL_CZ(k,i,j) < SPNUDGE_level2 ) then
              SPNUDGE_qv_alpha(k,i,j) = ( REAL_CZ(k,i,j) - SPNUDGE_level1 ) / ( SPNUDGE_level2 - SPNUDGE_level1 ) * qv_alpha
           elseif ( REAL_CZ(k,i,j) < SPNUDGE_level3 ) then
              SPNUDGE_qv_alpha(k,i,j) = qv_alpha
           elseif ( REAL_CZ(k,i,j) < SPNUDGE_level4 ) then
              SPNUDGE_qv_alpha(k,i,j) = ( REAL_CZ(k,i,j) - SPNUDGE_level3 ) / ( SPNUDGE_level4 - SPNUDGE_level3 ) * qv_alpha
           else
              SPNUDGE_qv_alpha(k,i,j) = 0.0_RP
           endif
        enddo
        enddo
        enddo

     else

        !$omp parallel do
        do j = JS, JE
        do i = IS, IE
        do k = KS, KE
           SPNUDGE_qv_alpha(k,i,j) = 0.0_RP
        end do
        end do
        end do

     end if

     !$acc enter data copyin(SPNUDGE_u_alpha, SPNUDGE_v_alpha, SPNUDGE_pt_alpha, SPNUDGE_qv_alpha)

     call DFT_setup( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     max( SPNUDGE_uv_lm, SPNUDGE_pt_lm, SPNUDGE_qv_lm ), &
                     max( SPNUDGE_uv_mm, SPNUDGE_pt_mm, SPNUDGE_qv_mm ) )

     return
   end subroutine SPNUDGE_setup

   subroutine SPNUDGE_finalize
     use scale_dft, only: &
        DFT_finalize

     call DFT_finalize

    !$acc exit data delete(SPNUDGE_u_alpha, SPNUDGE_v_alpha, SPNUDGE_pt_alpha, SPNUDGE_qv_alpha)
     deallocate( SPNUDGE_u_alpha )
     deallocate( SPNUDGE_v_alpha )
     deallocate( SPNUDGE_pt_alpha )
     deallocate( SPNUDGE_qv_alpha )

     return
   end subroutine SPNUDGE_finalize

end module scale_spnudge
