!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Set coriolis parameter and boundary conditions for baroclinic wave in a channel.
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

  use scale_time, only: &
       DTSEC => TIME_DTSEC, &
       NOWTSEC => TIME_NOWSEC

  use scale_const, only: &
       PI  => CONST_PI,    &
       OHM => CONST_OHM,   &
       RPlanet => CONST_RADIUS, &
       Rdry => CONST_Rdry
  
  use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_pres, &
       ATMOS_REFSTATE_temp, &
       ATMOS_REFSTATE_dens, &
       ATMOS_REFSTATE_pott, &
       ATMOS_REFSTATE_qv,   &
       ATMOS_REFSTATE_write
  
  use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       PRES

  use scale_process
  
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
  logical,  private, save :: USER_do = .false. !< do user step?

  real(RP), private, allocatable :: RHOT_bc(:,:,:)
  real(RP), private, allocatable :: DENS_bc(:,:,:)
  
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
       USER_do
    
    integer :: ierr
    integer :: i, j
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
    
    !
    allocate( RHOT_bc(KA,IA,2) )
    allocate( DENS_bc(KA,IA,2) )

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0

    implicit none
    !---------------------------------------------------------------------------
    integer :: j

    ! Save some information of inital fields to set boundary conditions. 
    ! 
    RHOT_bc(:,:,1) = RHOT(:,:,JS) - 0.5_RP*(RHOT(:,:,JS+1) - RHOT(:,:,JS))
    RHOT_bc(:,:,2) = RHOT(:,:,JE-1) + 1.5_RP*(RHOT(:,:,JE) - RHOT(:,:,JE-1))

    DENS_bc(:,:,1) = DENS(:,:,JS) - 0.5_RP*(DENS(:,:,JS+1) - DENS(:,:,JS))
    DENS_bc(:,:,2) = DENS(:,:,JE-1) + 1.5_RP*(DENS(:,:,JE) - DENS(:,:,JE-1))

    call USER_step
    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume

    use scale_comm, only: &
         COMM_vars8, COMM_wait
    
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
       CZ => GRID_CZ
    use scale_time, only: &
       DTSEC => TIME_DTSEC
    use scale_gridtrans
    
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none


    integer :: k, i, j
    
    !---------------------------------------------------------------------------

    if ( .not. USER_do ) return
    
    ! Apply the boundary condition at y=+Ly and y=-Ly

    if ( .NOT. PRC_HAS_N ) then     
       MOMY(:,:,JE)   = 0.0_RP
       do j = 1, JHALO
          MOMY(:,:,JE+j  ) = - MOMY(:,:,JE-j  )
          DENS(:,:,JE+j) =  2.0_RP*DENS_bc(:,:,2) - DENS(:,:,JE-j+1)
          MOMX(:,:,JE+j) = - MOMX(:,:,JE-j+1)          
          MOMZ(:,:,JE+j) = - MOMZ(:,:,JE-j+1)
          RHOT(:,:,JE+j) = 2.0_RP*RHOT_bc(:,:,2) - RHOT(:,:,JE-j+1)
       enddo
    end if

    if ( .NOT. PRC_HAS_S ) then     
       MOMY(:,:,JS-1) = 0.0_RP
       do j = 1, JHALO
          if ( j < JHALO ) MOMY(:,:,JS-j-1) = - MOMY(:,:,JS+j-1)
          DENS(:,:,JS-j) = 2.0_RP*DENS_bc(:,:,1) - DENS(:,:,JS+j-1)
          MOMX(:,:,JS-j) = - MOMX(:,:,JS+j-1)          
          MOMZ(:,:,JS-j) = - MOMZ(:,:,JS+j-1)
          RHOT(:,:,JS-j) = 2.0_RP*RHOT_bc(:,:,1) - RHOT(:,:,JS+j-1)
       enddo
    end if
    
    return
  end subroutine USER_step

end module mod_user
