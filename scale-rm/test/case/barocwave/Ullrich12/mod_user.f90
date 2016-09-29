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
  real(RP), private, allocatable :: RHOT_bc(:,:,:)
  real(RP), private, allocatable :: DENS_bc(:,:,:)
  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    
    use scale_grid, only : &
       CY => GRID_CY,      &
       FY => GRID_FY,      &
       RFDY => GRID_RFDY,  &
       GRID_DOMAIN_CENTER_Y
    
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
    y0 = GRID_DOMAIN_CENTER_Y
    
    do j = JS-1,JE+1
    do i = IS-1,IE+1
      CORIOLI(i,j) = f0 + beta0*(CY(j) - y0)
    enddo      
    enddo

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

    ! Save some information of  boundary condition
    
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

    if ( NOWTSEC == 0.0_RP ) then
       ! Set background fields
!!$       ATMOS_REFSTATE_pres(:,:,:) = PRES(:,:,:)
!!$       ATMOS_REFSTATE_dens(:,:,:) = DENS(:,:,:)
!!$       ATMOS_REFSTATE_pott(:,:,:) = RHOT(:,:,:)/DENS(:,:,:)
!!$       ATMOS_REFSTATE_temp(:,:,:) = PRES(:,:,:)/(Rdry*DENS(:,:,:))
!!$       ATMOS_REFSTATE_qv(:,:,:)   = 0.0_RP
!!$
!!$       call COMM_vars8( ATMOS_REFSTATE_dens(:,:,:), 1 )
!!$       call COMM_vars8( ATMOS_REFSTATE_temp(:,:,:), 2 )
!!$       call COMM_vars8( ATMOS_REFSTATE_pres(:,:,:), 3 )
!!$       call COMM_vars8( ATMOS_REFSTATE_pott(:,:,:), 4 )
!!$       call COMM_vars8( ATMOS_REFSTATE_qv  (:,:,:), 5 )
!!$       call COMM_wait ( ATMOS_REFSTATE_dens(:,:,:), 1, .false. )
!!$       call COMM_wait ( ATMOS_REFSTATE_temp(:,:,:), 2, .false. )
!!$       call COMM_wait ( ATMOS_REFSTATE_pres(:,:,:), 3, .false. )
!!$       call COMM_wait ( ATMOS_REFSTATE_pott(:,:,:), 4, .false. )
!!$       call COMM_wait ( ATMOS_REFSTATE_qv  (:,:,:), 5, .false. )
!!$
!!$       call ATMOS_REFSTATE_write       
    end if
    
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

          DENS(:,:,JE+j) =  2.0_RP*DENS_bc(:,:,2) - DENS(:,:,JE-j+1)
          MOMX(:,:,JE+j) = - MOMX(:,:,JE-j+1)          
          MOMZ(:,:,JE+j) = - MOMZ(:,:,JE-j+1)
          RHOT(:,:,JE+j) = 2.0_RP*RHOT_bc(:,:,2) - RHOT(:,:,JE-j+1)
       enddo
    end if

    if ( .NOT. PRC_HAS_S ) then     
       do j = 1, JHALO
          MOMY(:,:,JS-1) = 0d0
          MOMY(:,:,JS-j-1) = - MOMY(:,:,JS+j-1)

          DENS(:,:,JS-j) = 2.0_RP*DENS_bc(:,:,1) - DENS(:,:,JS+j-1)
          MOMX(:,:,JS-j) = - MOMX(:,:,JS+j-1)          
          MOMZ(:,:,JS-j) = - MOMZ(:,:,JS+j-1)
          RHOT(:,:,JS-j) = 2.0_RP*RHOT_bc(:,:,1) - RHOT(:,:,JS+j-1)
       enddo
    end if
    
    return
  end subroutine USER_step

end module mod_user
