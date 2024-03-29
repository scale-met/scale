!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Add the tendecy term by second diffrential order diffusion
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer

  use scale_atmos_grid_cartesC, only : &
       RCDX => ATMOS_GRID_CARTESC_RCDX, &
       RCDY => ATMOS_GRID_CARTESC_RCDY, &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
       RFDX => ATMOS_GRID_CARTESC_RFDX, &
       RFDY => ATMOS_GRID_CARTESC_RFDY, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ

  use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_finalize
  public :: USER_mkinit
  public :: USER_calc_tendency
  public :: USER_update

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
  logical,  private, save :: USER_do = .false.   !< do user step?
  real(RP), private, save :: Kdiff   = 75.E0_RP  !< Eddy diffusivity

  integer,  private              :: I_COMM_DENS = 1
  integer,  private              :: I_COMM_MOMZ = 2
  integer,  private              :: I_COMM_MOMX = 3
  integer,  private              :: I_COMM_MOMY = 4
  integer,  private              :: I_COMM_RHOT = 5

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine USER_tracer_setup

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
       Kdiff

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'
    LOG_INFO("USER_setup",*) 'User procedure in test/case/coldbubble/Straka93'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Finalization
  subroutine USER_finalize
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_finalize

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calc tendency
  subroutine USER_calc_tendency
    implicit none
    !---------------------------------------------------------------------------

    ! calculate diagnostic value and input to history buffer

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_update
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_file_history, only: &
       FILE_HISTORY_in

    use scale_time, only: &
         DTSEC => TIME_DTSEC, &
         NOWDAYSEC => TIME_NOWDAYSEC

    use scale_comm

    implicit none

    integer :: k, i, j

    !---------------------------------------------------------------------------

    if ( .not. USER_do ) return

    call FILE_HISTORY_in( (RHOT/DENS - 300.0_RP), 'PT_diff', 'PT_diff', 'K' )

    ! Consider eddy turbulent mixing with constant eddy viscosity and diffusivity

    if ( NOWDAYSEC > DTSEC ) then

       call COMM_vars8( DENS(:,:,:), I_COMM_DENS )
       call COMM_vars8( MOMZ(:,:,:), I_COMM_MOMZ )
       call COMM_vars8( MOMX(:,:,:), I_COMM_MOMX )
       call COMM_vars8( MOMY(:,:,:), I_COMM_MOMY )
       call COMM_vars8( RHOT(:,:,:), I_COMM_RHOT )

       call COMM_wait ( DENS(:,:,:), I_COMM_DENS, .false. )
       call COMM_wait ( MOMZ(:,:,:), I_COMM_MOMZ, .false. )
       call COMM_wait ( MOMX(:,:,:), I_COMM_MOMX, .false. )
       call COMM_wait ( MOMY(:,:,:), I_COMM_MOMY, .false. )
       call COMM_wait ( RHOT(:,:,:), I_COMM_RHOT, .false. )

       call append_EddyDiff_zxy( RHOT,                   & ! (inout)
            DENS        )                                  ! (in)

       call append_EddyDiff_zuy( MOMX,                   & ! (inout)
            DENS        )                                  ! (in)

       call append_EddyDiff_zxv( MOMY,                   & ! (inout)
            DENS        )                                  ! (in)

       call append_EddyDiff_wxy( MOMZ,                   & ! (inout)
            DENS        )                                  ! (in)

       call COMM_vars8( DENS(:,:,:), I_COMM_DENS )
       call COMM_vars8( MOMZ(:,:,:), I_COMM_MOMZ )
       call COMM_vars8( MOMX(:,:,:), I_COMM_MOMX )
       call COMM_vars8( MOMY(:,:,:), I_COMM_MOMY )
       call COMM_vars8( RHOT(:,:,:), I_COMM_RHOT )

       call COMM_wait ( DENS(:,:,:), I_COMM_DENS, .false. )
       call COMM_wait ( MOMZ(:,:,:), I_COMM_MOMZ, .false. )
       call COMM_wait ( MOMX(:,:,:), I_COMM_MOMX, .false. )
       call COMM_wait ( MOMY(:,:,:), I_COMM_MOMY, .false. )
       call COMM_wait ( RHOT(:,:,:), I_COMM_RHOT, .false. )

    end if

    return
  end subroutine USER_update

  !---------------------------------------------------------------------

  subroutine append_EddyDiff_zxy( RHOPHI,   & ! (inout)
       DENS                                 & ! (in)
       )
    use scale_prc

    use scale_time, only: &
       DTSEC => TIME_DTSEC

    real(RP), intent(inout) :: RHOPHI(KA,IA,JA)
    real(RP), intent(in)    :: DENS(KA,IA,JA)

    real(RP) :: DiffFlx(KA,IA,JA,3)
    real(RP) :: PHI(KA,IA,JA)
    integer :: i
    integer :: j
    integer :: k

    do j= JS-1, JE+1
    do i= IS-1, IE+1
    do k= KS-1, KE+1
       PHI(k,i,j) = RHOPHI(k,i,j) / DENS(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       DiffFlx(k,i,j,ZDIR) = 0.5_RP*(DENS(k+1,i,j) + DENS(k,i,j))*Kdiff* &
            (PHI(k+1,i,j) - PHI(k,i,j))*RFDZ(k)
    enddo
    DiffFlx(KS-1,i,j,ZDIR) = 0.0_RP
    DiffFlx(KE,i,j,ZDIR)   = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS-1, min(IE,IEH)
    do k = KS, KE
       DiffFlx(k,i,j,XDIR) = 0.5_RP*(DENS(k,i+1,j) + DENS(k,i,j))*Kdiff* &
            (PHI(k,i+1,j) - PHI(k,i,j))*RFDX(i)
    enddo
    enddo
    enddo

    do j = JS-1, min(JE,JEH)
    do i = IS, IE
    do k = KS, KE
       DiffFlx(k,i,j,YDIR) = 0.5_RP*(DENS(k,i,j+1) + DENS(k,i,j))*Kdiff* &
            (PHI(k,i,j+1) - PHI(k,i,j))*RFDY(j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOPHI(k,i,j) = RHOPHI(k,i,j) + DTSEC*( &
              (DiffFlx(k,i,j,ZDIR) - DiffFlx(k-1,i,j,ZDIR))*RCDZ(k) &
            + (DiffFlx(k,i,j,XDIR) - DiffFlx(k,i-1,j,XDIR))*RCDX(i) &
            + (DiffFlx(k,i,j,YDIR) - DiffFlx(k,i,j-1,YDIR))*RCDY(j) &
            )
    enddo
    enddo
    enddo

  end subroutine append_EddyDiff_zxy

  subroutine append_EddyDiff_zuy( RHOPHI,   & ! (inout)
       DENS                                 & ! (in)
       )

    use scale_prc
    use scale_time, only: &
       DTSEC => TIME_DTSEC

    real(RP), intent(inout) :: RHOPHI(KA,IA,JA)
    real(RP), intent(in)    :: DENS(KA,IA,JA)

    real(RP) :: DiffFlx(KA,IA,JA,3)
    real(RP) :: PHI(KA,IA,JA)
    integer :: i
    integer :: j
    integer :: k

    do j= JS-1, JE+1
    do i= IS-1, IE+1
    do k= KS-1, KE+1
       PHI(k,i,j) = RHOPHI(k,i,j) * 2.0_RP/ (DENS(k,i,j) + DENS(k,i+1,j))
    enddo
    enddo
    enddo


    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       DiffFlx(k,i,j,ZDIR) = 0.25_RP*sum(DENS(k:k+1,i:i+1,j))*Kdiff* &
            (PHI(k+1,i,j) - PHI(k,i,j))*RFDZ(k)
    enddo
    DiffFlx(KS-1,i,j,ZDIR) = 0.0_RP
    DiffFlx(KE,i,j,ZDIR)   = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE+1
    do k = KS, KE
       DiffFlx(k,i-1,j,XDIR) = DENS(k,i,j)*Kdiff* &
            (PHI(k,i,j) - PHI(k,i-1,j))*RCDX(i)
    enddo
    enddo
    enddo

    do j = JS-1, JE
    do i = IS, IE
    do k = KS, KE
       DiffFlx(k,i,j,YDIR) = 0.25_RP*sum(DENS(k,i:i+1,j:j+1))*Kdiff* &
            (PHI(k,i,j+1) - PHI(k,i,j))*RFDY(j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, min(IE, IEH)
    do k = KS, KE
       RHOPHI(k,i,j) = RHOPHI(k,i,j) + DTSEC*( &
              (DiffFlx(k,i,j,ZDIR) - DiffFlx(k-1,i,j,ZDIR))*RCDZ(k) &
            + (DiffFlx(k,i,j,XDIR) - DiffFlx(k,i-1,j,XDIR))*RFDX(i) &
            + 0.0_RP*(DiffFlx(k,i,j,YDIR) - DiffFlx(k,i,j-1,YDIR))*RCDY(j) &
            )
    enddo
    enddo
    enddo

  end subroutine append_EddyDiff_zuy

  subroutine append_EddyDiff_zxv( RHOPHI,   & ! (inout)
       DENS                                 & ! (in)
       )


    use scale_time, only: &
       DTSEC => TIME_DTSEC

    real(RP), intent(inout) :: RHOPHI(KA,IA,JA)
    real(RP), intent(in)    :: DENS(KA,IA,JA)

    real(RP) :: DiffFlx(KA,IA,JA,3)
    real(RP) :: PHI(KA,IA,JA)
    integer :: i
    integer :: j
    integer :: k

    do j= JS-1, JE+1
    do i= IS-1, IE+1
    do k= KS-1, KE+1
       PHI(k,i,j) = RHOPHI(k,i,j) * 2.0_RP/ (DENS(k,i,j) + DENS(k,i,j+1))
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       DiffFlx(k,i,j,ZDIR) = 0.25_RP*sum(DENS(k:k+1,i,j:j+1))*Kdiff* &
            (PHI(k+1,i,j) - PHI(k,i,j))*RFDZ(k)
    enddo
    DiffFlx(KS-1,i,j,ZDIR) = 0.0_RP
    DiffFlx(KE,i,j,ZDIR)   = 0.0_RP
    enddo
    enddo

    do j = JS, min(JE,JEH)
    do i = IS-1, IE
    do k = KS, KE
       DiffFlx(k,i,j,XDIR) = 0.25_RP*sum(DENS(k,i:i+1,j:j+1))*Kdiff* &
            (PHI(k,i+1,j) - PHI(k,i,j))*RFDX(i)
    enddo
    enddo
    enddo

    do j = JS, JE+1
    do i = IS, IE
    do k = KS, KE
       DiffFlx(k,i,j-1,YDIR) = DENS(k,i,j)*Kdiff* &
            (PHI(k,i,j) - PHI(k,i,j-1))*RCDY(j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOPHI(k,i,j) = RHOPHI(k,i,j) + DTSEC*( &
              (DiffFlx(k,i,j,ZDIR) - DiffFlx(k-1,i,j,ZDIR))*RCDZ(k) &
            + (DiffFlx(k,i,j,XDIR) - DiffFlx(k,i-1,j,XDIR))*RCDX(i) &
            + (DiffFlx(k,i,j,YDIR) - DiffFlx(k,i,j-1,YDIR))*RFDY(j) &
            )
    enddo
    enddo
    enddo

  end subroutine append_EddyDiff_zxv

  subroutine append_EddyDiff_wxy( RHOPHI,   & ! (inout)
       DENS                                 & ! (in)
       )


    use scale_time, only: &
       DTSEC => TIME_DTSEC

    real(RP), intent(inout) :: RHOPHI(KA,IA,JA)
    real(RP), intent(in)    :: DENS(KA,IA,JA)

    real(RP) :: DiffFlx(KA,IA,JA,3)
    real(RP) :: PHI(KA,IA,JA)
    integer :: i
    integer :: j
    integer :: k

    do j= JS-1, JE+1
    do i= IS-1, IE+1
    do k= KS-1, KE+1
       PHI(k,i,j) = RHOPHI(k,i,j) * 2.0_RP/ (DENS(k,i,j) + DENS(k+1,i,j))
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DiffFlx(k-1,i,j,ZDIR) = DENS(k,i,j)*Kdiff* &
            (PHI(k,i,j) - PHI(k-1,i,j))*RCDZ(k)
    enddo
    DiffFlx(KS-1,i,j,ZDIR) = 0.0_RP
    DiffFlx(KE,i,j,ZDIR)   = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS-1, IE
    do k = KS, KE
       DiffFlx(k,i,j,XDIR) = 0.25_RP*sum(DENS(k:k+1,i:i+1,j))*Kdiff* &
            (PHI(k,i+1,j) - PHI(k,i,j))*RFDX(i)
    enddo
    enddo
    enddo

    do j = JS-1, JE
    do i = IS, IE
    do k = KS, KE
       DiffFlx(k,i,j,YDIR) = 0.25_RP*sum(DENS(k:k+1,i,j:j+1))*Kdiff* &
            (PHI(k,i,j+1) - PHI(k,i,j))*RFDY(j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOPHI(k,i,j) = RHOPHI(k,i,j) + DTSEC*( &
              (DiffFlx(k,i,j,ZDIR) - DiffFlx(k-1,i,j,ZDIR))*RFDZ(k) &
            + (DiffFlx(k,i,j,XDIR) - DiffFlx(k,i-1,j,XDIR))*RCDX(i) &
            + (DiffFlx(k,i,j,YDIR) - DiffFlx(k,i,j-1,YDIR))*RCDY(j) &
            )
    enddo
    enddo
    enddo

  end subroutine append_EddyDiff_wxy

end module mod_user
