!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)   [new]
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
  use scale_grid
  use scale_grid_index
  use scale_tracer
  use scale_index
  use scale_time, only: &
     TIME_DTSEC
  use scale_atmos_boundary, only: &
     ATMOS_BOUNDARY_alpha, &
     ATMOS_BOUNDARY_var
  use mod_atmos_vars
  use scale_history, only: &
     HIST_in
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
  logical, private :: USER_do  = .false. !< do user step?
  logical, private :: CONST_MOMX = .false.         ! assume constant momentum-x
  integer, private :: NUM_RELAX_GRIDS = 10        ! num of relaxation grid
  real(RP), private :: UINIT  = 0.0_RP
  real(RP), private :: UEND   = 5.0_RP
  real(RP), private :: FRONT_OFFSET   = 1000.0_RP

  real(RP), private, allocatable :: DENS_INIT(:)
  real(RP), private, allocatable :: POTT_INIT(:)
  real(RP), private, allocatable :: VELX_INIT(:)
  real(RP), private, allocatable :: DENS_tend(:,:)
  real(RP), private, allocatable :: U_tend(:)

  real(RP), private :: UACCEL_per_grid
  real(RP), private :: front_posi
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_USER / &
       USER_do,             &
       NUM_RELAX_GRIDS,     &
       UINIT,               &
       UEND,                &
       FRONT_OFFSET,        &
       CONST_MOMX

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

    if( IO_L ) write(IO_FID_LOG,*) '*** This module is accelerate zonal-wind.'

    allocate( DENS_INIT(KA                  ) )
    allocate( POTT_INIT(KA                  ) )
    allocate( VELX_INIT(KA                  ) )
    allocate( DENS_tend(KA,0:NUM_RELAX_GRIDS) )
    allocate( U_tend   (   0:NUM_RELAX_GRIDS) )

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none

    integer :: i, j, k
    !---------------------------------------------------------------------------

    do k = 1, KA
       DENS_INIT(k) = DENS(k,IS,JS)
       POTT_INIT(k) = POTT(k,IS,JS)
       VELX_INIT(k) = MOMX(k,IS,JS) / ( DENS(k,IS,JS)+DENS(k,IS+1,JS) ) * 2.0_RP
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ATMOS_BOUNDARY_var(k,i,j,I_BND_DENS) = DENS_INIT(k)
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = UINIT
       ATMOS_BOUNDARY_var(k,i,j,I_BND_POTT) = POTT_INIT(k)
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 2, KA
      ATMOS_BOUNDARY_alpha(k,i,j,I_BND_DENS) = ATMOS_BOUNDARY_alpha(1,i,j,I_BND_DENS)
      ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELX) = ATMOS_BOUNDARY_alpha(1,i,j,I_BND_VELX)
      ATMOS_BOUNDARY_alpha(k,i,j,I_BND_POTT) = ATMOS_BOUNDARY_alpha(1,i,j,I_BND_POTT)
    enddo
    enddo
    enddo

    do k = 1, KA
       UACCEL_per_grid = (UEND - UINIT) / dble( NUM_RELAX_GRIDS )
    enddo

    DENS_tend(k,0) = 0.0_RP
    U_tend(0)      = 0.0_RP
    do i = 1, NUM_RELAX_GRIDS
       U_tend(i) = U_tend(i-1) + UACCEL_per_grid

       if ( CONST_MOMX ) then
          do k = 1, KA
             DENS_tend(k,i) = ((UINIT*DENS_INIT(k)) / (UINIT + U_tend(i))) - DENS_INIT(k)
          enddo
       else
          DENS_tend(:,i) = 0.0_RP
       endif
    enddo

    front_posi = -1.0_RP * FRONT_OFFSET

    do i = 1, NUM_RELAX_GRIDS
       if( IO_L ) write(IO_FID_LOG,*) '*** DENS_tend(KS)', i, DENS_tend(KS,i)
       if( IO_L ) write(IO_FID_LOG,*) '*** U_tend       ', i, U_tend(i)
    enddo

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
  !> User step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    integer :: i, j, k, ii
    integer :: front_grid
    real(RP) :: mini, diff, dist1, dist2
    real(RP) :: fact_A, fact_B

    logical  :: zintp
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       call PRC_MPIstop
    endif

    call HIST_in( ATMOS_BOUNDARY_var(:,:,:,I_BND_DENS), 'BND_DENS', 'boundary_dens', 'kg/m3' )
    call HIST_in( ATMOS_BOUNDARY_var(:,:,:,I_BND_VELX), 'BND_VELX', 'boundary_velx', 'm/s'   )

    front_posi = front_posi + UEND * TIME_DTSEC
    if( IO_L ) write(IO_FID_LOG,*) '*** front position', front_posi

    if ( front_posi > 0.0_RP .and. front_posi < GRID_FX(IA) ) then
       mini = 9.999D10
       do i = 1, IA
          diff = GRID_FX(i) - front_posi
          if ( diff > 0.0_RP .and. diff < mini ) then
             mini = diff
             front_grid = i
          endif
       enddo

       dist1 = mini
       dist2 = GRID_FDX(front_grid) - dist1
       fact_A = dist1 / GRID_FDX(front_grid)
       fact_B = dist2 / GRID_FDX(front_grid)
       if( IO_L ) write(IO_FID_LOG,*) '*** front grid', front_grid

       ii = front_grid
       do i = 1, NUM_RELAX_GRIDS
       do k = 1, KA
       do j = 1, JA
          ATMOS_BOUNDARY_var(k,ii,j,I_BND_DENS) = DENS_INIT(k) + (DENS_tend(k,i-1)*fact_A + DENS_tend(k,i)*fact_B)
          ATMOS_BOUNDARY_var(k,ii,j,I_BND_POTT) = POTT_INIT(k)
          ATMOS_BOUNDARY_var(k,ii,j,I_BND_VELX) = UINIT + (U_tend(i-1)*fact_A + U_tend(i)*fact_B)
       enddo
       enddo
       ii = ii - 1
       if ( ii < 1 ) exit
       enddo

       if ( ii > 1 ) then
          do j = 1, JA
          do i = ii, 1, -1
          do k = 1, KA
             ATMOS_BOUNDARY_var(k,i,j,I_BND_DENS) = DENS_INIT(k) + DENS_tend(k,NUM_RELAX_GRIDS)
             ATMOS_BOUNDARY_var(k,i,j,I_BND_POTT) = POTT_INIT(k)
             ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = UINIT + U_tend(NUM_RELAX_GRIDS)
          enddo
          enddo
          enddo
       endif

    elseif ( front_posi > GRID_FX(IE) ) then
          do j = 1, JA
          do i = 1, IA
          do k = 1, KA
             ATMOS_BOUNDARY_var(k,i,j,I_BND_DENS) = DENS_INIT(k) + DENS_tend(k,NUM_RELAX_GRIDS)
             ATMOS_BOUNDARY_var(k,i,j,I_BND_POTT) = POTT_INIT(k)
             ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = UEND
          enddo
          enddo
          enddo
    endif

    return
  end subroutine USER_step

end module mod_user
