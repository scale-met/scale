!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Set boundary condition and add eddy diffusion for a lid-driven cavity flow test
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-09-31 (Y.Kawai)   [add comments]
!! @li      2016-08-03 (Y.Kawai)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

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
  logical,  private :: USER_do = .false. !< do user step?
  real(RP), private :: Ulid    = 1.0_RP  !< Velocity along y=+L driving interior flow
  real(RP), private :: PRES0   = 1.E4_RP !< Reference pressure
  real(RP), private :: KDiff   = 20.0_RP !< Eddy viscosity

  integer,  private :: I_COMM_DENS = 1
  integer,  private :: I_COMM_MOMZ = 2
  integer,  private :: I_COMM_MOMX = 3
  integer,  private :: I_COMM_MOMY = 4
  integer,  private :: I_COMM_RHOT = 5

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
       Ulid,    &
       PRES0,   &
       KDiff

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    call USER_step

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only : &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RCDZ => GRID_RCDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY, &
       RFDZ => GRID_RFDZ
    use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT
    use scale_rm_process, only: &
       PRC_HAS_N, &
       PRC_HAS_E, &
       PRC_HAS_S, &
       PRC_HAS_W
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_grid, only : &
       CX => GRID_CX, &
       CY => GRID_CY, &
       CZ => GRID_CZ
    use scale_time, only: &
       DTSEC => TIME_DTSEC, &
       NOWTSEC => TIME_NOWSEC
    use scale_gridtrans
    use scale_history, only: &
       HIST_in
    use scale_comm
    implicit none

    real(RP) :: one(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( .not. USER_do ) return


    ! Consider eddy turbulent mixing with constant eddy viscosity and diffusivity
    !

    if ( NOWTSEC > DTSEC ) then

       call COMM_vars8( DENS(:,:,:), I_COMM_DENS )
!!$       call COMM_vars8( MOMZ(:,:,:), I_COMM_MOMZ )
       call COMM_vars8( MOMX(:,:,:), I_COMM_MOMX )
       call COMM_vars8( MOMY(:,:,:), I_COMM_MOMY )
       call COMM_vars8( RHOT(:,:,:), I_COMM_RHOT )

       call COMM_wait ( DENS(:,:,:), I_COMM_DENS, .false. )
!!$       call COMM_wait ( MOMZ(:,:,:), I_COMM_MOMZ, .false. )
       call COMM_wait ( MOMX(:,:,:), I_COMM_MOMX, .false. )
       call COMM_wait ( MOMY(:,:,:), I_COMM_MOMY, .false. )
       call COMM_wait ( RHOT(:,:,:), I_COMM_RHOT, .false. )

       call append_EddyDiff_zxy( RHOT,                   & ! (inout)
            DENS        )                                  ! (in)

       call append_EddyDiff_zuy( MOMX,                   & ! (inout)
            DENS        )                                  ! (in)

       call append_EddyDiff_zxv( MOMY,                   & ! (inout)
            DENS        )                                  ! (in)
!!$
!!$       call append_EddyDiff_wxy( MOMZ,                   & ! (inout)
!!$            DENS        )                                  ! (in)

       one = 1.0_RP
       call append_EddyDiff_zxy( DENS,                   & ! (inout)
            one        )                                   ! (in)

       call COMM_vars8( DENS(:,:,:), I_COMM_DENS )
!!$       call COMM_vars8( MOMZ(:,:,:), I_COMM_MOMZ )
       call COMM_vars8( MOMX(:,:,:), I_COMM_MOMX )
       call COMM_vars8( MOMY(:,:,:), I_COMM_MOMY )
       call COMM_vars8( RHOT(:,:,:), I_COMM_RHOT )

       call COMM_wait ( DENS(:,:,:), I_COMM_DENS, .false. )
!!$       call COMM_wait ( MOMZ(:,:,:), I_COMM_MOMZ, .false. )
       call COMM_wait ( MOMX(:,:,:), I_COMM_MOMX, .false. )
       call COMM_wait ( MOMY(:,:,:), I_COMM_MOMY, .false. )
       call COMM_wait ( RHOT(:,:,:), I_COMM_RHOT, .false. )

    end if

    ! Apply boundary condition
    !

    if ( .NOT. PRC_HAS_N ) then
       MOMY(:,:,JE)     = 0.0_RP
       do j = 1, JHALO
          MOMX(:,:,JE+j)   = 2.0_RP * Ulid - MOMX(:,:,JE-j+1)
          MOMY(:,:,JE+j  ) = - MOMY(:,:,JE-j  )

          DENS(:,:,JE+j) = DENS(:,:,JE-j+1)
          MOMZ(:,:,JE+j) = MOMZ(:,:,JE-j+1)
          RHOT(:,:,JE+j) = RHOT(:,:,JE-j+1)
       enddo
    end if

    if ( .NOT. PRC_HAS_E ) then
       MOMX(:,IE,:)     = 0.0_RP
       do i = 1, IHALO
          MOMX(:,IE+i,:) = - MOMX(:,IE-i,:)
          MOMY(:,IE+i,:) = - MOMY(:,IE-i+1,:)

          DENS(:,IE+i,:) = DENS(:,IE-i+1,:)
          MOMZ(:,IE+i,:) = MOMZ(:,IE-i+1,:)
          RHOT(:,IE+i,:) = RHOT(:,IE-i+1,:)
       enddo
    end if

    if ( .NOT. PRC_HAS_S ) then
       MOMY(:,:,JS-1) = 0.0_RP
       do j = 1, JHALO
          MOMX(:,:,JS-j) = - MOMX(:,:,JS+j-1)

          if ( j < JHALO ) MOMY(:,:,JS-j-1) = - MOMY(:,:,JS+j-1)

          DENS(:,:,JS-j) = DENS(:,:,JS+j-1)
          MOMZ(:,:,JS-j) = MOMZ(:,:,JS+j-1)
          RHOT(:,:,JS-j) = RHOT(:,:,JS+j-1)
       enddo
    end if

    if ( .NOT. PRC_HAS_W ) then
       MOMX(:,IS-1,:)   = 0.0_RP
       do i = 1, IHALO
          if ( i < IHALO ) MOMX(:,IS-i-1,:) = - MOMX(:,IS+i-1,:)
          MOMY(:,IS-i,:)   = - MOMY(:,IS+i-1,:)

          DENS(:,IS-i,:) = DENS(:,IS+i-1,:)
          MOMZ(:,IS-i,:) = MOMZ(:,IS+i-1,:)
          RHOT(:,IS-i,:) = RHOT(:,IS+i-1,:)
       enddo
    end if

    return
  end subroutine USER_step


  !---------------------------------------------------------------------

  subroutine append_EddyDiff_zxy( RHOPHI,   & ! (inout)
       DENS                                 & ! (in)
       )
    use scale_process

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

!!$    do j = JS, JE
!!$    do i = IS, IE
!!$    do k = KS, KE-1
!!$       DiffFlx(k,i,j,ZDIR) = 0.5_RP*(DENS(k+1,i,j) + DENS(k,i,j))*Kdiff* &
!!$            (PHI(k+1,i,j) - PHI(k,i,j))*RFDZ(k)
!!$    enddo
!!$    DiffFlx(KS-1,i,j,ZDIR) = 0.0_RP
!!$    DiffFlx(KE,i,j,ZDIR)   = 0.0_RP
!!$    enddo
!!$    enddo
    DiffFlx(:,:,:,ZDIR) = 0.0_RP

    do j = JS, JE
    do i = IS-1, IE
    do k = KS, KE
       DiffFlx(k,i,j,XDIR) = 0.5_RP*(DENS(k,i+1,j) + DENS(k,i,j))*Kdiff* &
            (PHI(k,i+1,j) - PHI(k,i,j))*RFDX(i)
    enddo
    enddo
    enddo

    do j = JS-1, JE
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
            + 0d0*(DiffFlx(k,i,j,YDIR) - DiffFlx(k,i,j-1,YDIR))*RCDY(j) &
            )
    enddo
    enddo
    enddo

  end subroutine append_EddyDiff_zxy

  subroutine append_EddyDiff_zuy( RHOPHI,   & ! (inout)
       DENS                                 & ! (in)
       )

    use scale_process
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


!!$    do j = JS, JE
!!$    do i = IS, IE
!!$    do k = KS, KE-1
!!$       DiffFlx(k,i,j,ZDIR) = 0.25_RP*sum(DENS(k:k+1,i:i+1,j))*Kdiff* &
!!$            (PHI(k+1,i,j) - PHI(k,i,j))*RFDZ(k)
!!$    enddo
!!$    DiffFlx(KS-1,i,j,ZDIR) = 0.0_RP
!!$    DiffFlx(KE,i,j,ZDIR)   = 0.0_RP
!!$    enddo
!!$    enddo
    DiffFlx(:,:,:,ZDIR) = 0.0_RP

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
            + (DiffFlx(k,i,j,YDIR) - DiffFlx(k,i,j-1,YDIR))*RCDY(j) &
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

!!$    do j = JS, JE
!!$    do i = IS, IE
!!$    do k = KS, KE-1
!!$       DiffFlx(k,i,j,ZDIR) = 0.25_RP*sum(DENS(k:k+1,i,j:j+1))*Kdiff* &
!!$            (PHI(k+1,i,j) - PHI(k,i,j))*RFDZ(k)
!!$    enddo
!!$    DiffFlx(KS-1,i,j,ZDIR) = 0.0_RP
!!$    DiffFlx(KE,i,j,ZDIR)   = 0.0_RP
!!$    enddo
!!$    enddo
    DiffFlx(:,:,:,ZDIR) = 0.0_RP

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

!!$    do j = JS, JE
!!$    do i = IS, IE
!!$    do k = KS, KE
!!$       DiffFlx(k-1,i,j,ZDIR) = DENS(k,i,j)*Kdiff* &
!!$            (PHI(k,i,j) - PHI(k-1,i,j))*RCDZ(k)
!!$    enddo
!!$    DiffFlx(KS-1,i,j,ZDIR) = 0.0_RP
!!$    DiffFlx(KE,i,j,ZDIR)   = 0.0_RP
!!$    enddo
!!$    enddo
    DiffFlx(:,:,:,ZDIR) = 0.0_RP

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
