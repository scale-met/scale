!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Add the tendecy term by second diffrential order diffusion
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-06-24 (Y.Kawai)   [new]
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
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
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
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop

    use scale_const, only: &
       PI  => CONST_PI,    &
       OHM => CONST_OHM,   &
       RPlanet => CONST_RADIUS
    
    use scale_grid, only : &
       CY => GRID_CY, &
       FY => GRID_FY, &
       RFDY => GRID_RFDY
        
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
       Kdiff
    
    integer :: ierr
    integer :: i, j

    real(RP) :: f0
    real(RP) :: beta0
    real(RP) :: y0
    
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
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_USER)

    ! Set coriolis parameter in f or beta plane approximation
    

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0

    implicit none
    !---------------------------------------------------------------------------
    integer :: j

     call USER_step
    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    ! calculate diagnostic value and input to history buffer
    !call USER_step

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       GRAV  => CONST_GRAV
       
    use scale_gridtrans
    
    use scale_history, only: &
       HIST_in

    use scale_time, only: &
         DTSEC => TIME_DTSEC, &
         NOWTSEC => TIME_NOWSEC

    use scale_comm
    
    implicit none
    
    integer :: k, i, j
    
    !---------------------------------------------------------------------------

    if ( .not. USER_do ) return

    call HIST_in( (RHOT/DENS - 300.0_RP), 'PT_diff', 'PT_diff', 'K' )
    
    ! Consider eddy turbulent mixing with constant eddy viscosity and diffusivity

    if ( NOWTSEC > DTSEC ) then

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
  end subroutine USER_step

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
       DiffFlx(k,i,j,YDIR) = 0.5_RP*0.25_RP*sum(DENS(k,i:i+1,j:j+1))*Kdiff* &
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
