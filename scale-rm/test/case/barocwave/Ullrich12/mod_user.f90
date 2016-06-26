!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Set boundary condition for baroclinic wave in a channel.
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-06-26 (Y.Kawai)   [new]
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
  logical,  private, save :: USER_do = .false. !< do user step?
  real(RP), private, save :: Phi0Deg = 45.0_RP !< Central latitude for f or beta plane

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
    
    use scale_atmos_dyn, only: &
       CORIOLI
    
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
       Phi0Deg  
    
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
    
    f0 = 2.0_RP*OHM*sin( Phi0Deg*PI/180.0_RP )
    beta0 = 2.0_RP*OHM/RPlanet*cos( Phi0Deg*PI/180.0_RP )
    y0 = 0.5_RP*(FY(JE) - FY(JS-1))
    
    do j = JS-1,JE+1
    do i = IS-1,IE+1
      CORIOLI(i,j) = f0 + beta0*(CY(j) - y0)
    enddo      
    enddo

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
    use scale_rm_process, only: &
       PRC_HAS_N, &
       PRC_HAS_S
    use scale_grid, only : &
       CX => GRID_CX, &
       CY => GRID_CY, &
       CZ => GRID_CZ
    use scale_time, only: &
       DTSEC => TIME_DTSEC
    use scale_atmos_dyn, only: &
       CORIOLI
    use scale_gridtrans
    
    use scale_history, only: &
       HIST_in
    implicit none


    integer :: k, i, j
    
    !---------------------------------------------------------------------------

    if ( .not. USER_do ) return

    ! Apply boundary condition at y=+Ly and y=-Ly

    if ( .NOT. PRC_HAS_N ) then     
       do j = 1, JHALO
          MOMY(:,:,JE)   = 0d0
          MOMY(:,:,JE+j  ) = - MOMY(:,:,JE-j  )

          DENS(:,:,JE+j) = DENS(:,:,JE-j+1)
          MOMX(:,:,JE+j) = MOMX(:,:,JE-j+1)          
          MOMZ(:,:,JE+j) = MOMZ(:,:,JE-j+1)
          RHOT(:,:,JE+j) = RHOT(:,:,JE-j+1)
       enddo
    end if

    if ( .NOT. PRC_HAS_S ) then     
       do j = 1, JHALO
          MOMY(:,:,JS-1) = 0d0
          MOMY(:,:,JS-j-1) = - MOMY(:,:,JS+j-1)

          DENS(:,:,JS-j) = DENS(:,:,JS+j-1)
          MOMX(:,:,JS-j) = MOMX(:,:,JS+j-1)          
          MOMZ(:,:,JS-j) = MOMZ(:,:,JS+j-1)
          RHOT(:,:,JS-j) = RHOT(:,:,JS+j-1)
       enddo
    end if
    
    return
  end subroutine USER_step

end module mod_user
