!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
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

  use scale_time, only: &
     TIME_DTSEC
  use scale_atmos_boundary, only: &
     ATMOS_BOUNDARY_alpha_DENS, &
     ATMOS_BOUNDARY_alpha_VELX, &
     ATMOS_BOUNDARY_alpha_POTT, &
     ATMOS_BOUNDARY_DENS, &
     ATMOS_BOUNDARY_VELX, &
     ATMOS_BOUNDARY_POTT
  use mod_atmos_vars
  use scale_file_history, only: &
     FILE_HISTORY_in
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
       USER_do,             &
       NUM_RELAX_GRIDS,     &
       UINIT,               &
       UEND,                &
       FRONT_OFFSET,        &
       CONST_MOMX

    integer :: i, j, k
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'
    LOG_INFO("USER_setup",*) 'This module is accelerate zonal-wind.'

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

    allocate( DENS_INIT(KA                  ) )
    allocate( POTT_INIT(KA                  ) )
    allocate( VELX_INIT(KA                  ) )
    allocate( DENS_tend(KA,0:NUM_RELAX_GRIDS) )
    allocate( U_tend   (   0:NUM_RELAX_GRIDS) )



    do k = 1, KA
       DENS_INIT(k) = DENS(k,IS,JS)
       POTT_INIT(k) = POTT(k,IS,JS)
       VELX_INIT(k) = MOMX(k,IS,JS) / ( DENS(k,IS,JS)+DENS(k,IS+1,JS) ) * 2.0_RP
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ATMOS_BOUNDARY_DENS(k,i,j) = DENS_INIT(k)
       ATMOS_BOUNDARY_VELX(k,i,j) = UINIT
       ATMOS_BOUNDARY_POTT(k,i,j) = POTT_INIT(k)
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 2, KA
      ATMOS_BOUNDARY_alpha_DENS(k,i,j) = ATMOS_BOUNDARY_alpha_DENS(1,i,j)
      ATMOS_BOUNDARY_alpha_VELX(k,i,j) = ATMOS_BOUNDARY_alpha_VELX(1,i,j)
      ATMOS_BOUNDARY_alpha_POTT(k,i,j) = ATMOS_BOUNDARY_alpha_POTT(1,i,j)
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
       LOG_INFO("USER_setup",*) 'DENS_tend(KS)', i, DENS_tend(KS,i)
       LOG_INFO("USER_setup",*) 'U_tend       ', i, U_tend(i)
    enddo

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Finalization
  subroutine USER_finalize
    implicit none
    !---------------------------------------------------------------------------

    allocate( DENS_INIT )
    allocate( POTT_INIT )
    allocate( VELX_INIT )
    allocate( DENS_tend )
    allocate( U_tend )

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
  !> Calculate tendency
  subroutine USER_calc_tendency
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_update
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_FX, &
       ATMOS_GRID_CARTESC_FDX
    implicit none
    integer :: i, j, k, ii
    integer :: front_grid
    real(RP) :: mini, diff, dist1, dist2
    real(RP) :: fact_A, fact_B

    logical  :: zintp
    !---------------------------------------------------------------------------

    call FILE_HISTORY_in( ATMOS_BOUNDARY_DENS(:,:,:), 'BND_DENS', 'boundary_dens', 'kg/m3' )
    call FILE_HISTORY_in( ATMOS_BOUNDARY_VELX(:,:,:), 'BND_VELX', 'boundary_velx', 'm/s'   )

    front_posi = front_posi + UEND * TIME_DTSEC
    LOG_INFO("USER_update",*) 'front position', front_posi

    if ( front_posi > 0.0_RP .and. front_posi < ATMOS_GRID_CARTESC_FX(IA) ) then
       mini = 9.999D10
       do i = 1, IA
          diff = ATMOS_GRID_CARTESC_FX(i) - front_posi
          if ( diff > 0.0_RP .and. diff < mini ) then
             mini = diff
             front_grid = i
          endif
       enddo

       dist1 = mini
       dist2 = ATMOS_GRID_CARTESC_FDX(front_grid) - dist1
       fact_A = dist1 / ATMOS_GRID_CARTESC_FDX(front_grid)
       fact_B = dist2 / ATMOS_GRID_CARTESC_FDX(front_grid)
       LOG_INFO("USER_update",*) 'front grid', front_grid

       ii = front_grid
       do i = 1, NUM_RELAX_GRIDS
       do k = 1, KA
       do j = 1, JA
          ATMOS_BOUNDARY_DENS(k,ii,j) = DENS_INIT(k) + (DENS_tend(k,i-1)*fact_A + DENS_tend(k,i)*fact_B)
          ATMOS_BOUNDARY_POTT(k,ii,j) = POTT_INIT(k)
          ATMOS_BOUNDARY_VELX(k,ii,j) = UINIT + (U_tend(i-1)*fact_A + U_tend(i)*fact_B)
       enddo
       enddo
       ii = ii - 1
       if ( ii < 1 ) exit
       enddo

       if ( ii > 1 ) then
          do j = 1, JA
          do i = ii, 1, -1
          do k = 1, KA
             ATMOS_BOUNDARY_DENS(k,i,j) = DENS_INIT(k) + DENS_tend(k,NUM_RELAX_GRIDS)
             ATMOS_BOUNDARY_POTT(k,i,j) = POTT_INIT(k)
             ATMOS_BOUNDARY_VELX(k,i,j) = UINIT + U_tend(NUM_RELAX_GRIDS)
          enddo
          enddo
          enddo
       endif

    elseif ( front_posi > ATMOS_GRID_CARTESC_FX(IE) ) then
          do j = 1, JA
          do i = 1, IA
          do k = 1, KA
             ATMOS_BOUNDARY_DENS(k,i,j) = DENS_INIT(k) + DENS_tend(k,NUM_RELAX_GRIDS)
             ATMOS_BOUNDARY_POTT(k,i,j) = POTT_INIT(k)
             ATMOS_BOUNDARY_VELX(k,i,j) = UEND
          enddo
          enddo
          enddo
    endif

    return
  end subroutine USER_update

end module mod_user
