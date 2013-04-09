!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Boundary treatment
!!
!! @par Description
!!          Boundary treatment of model domain
!!          Additional forcing, Sponge layer, rayleigh dumping
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-07 (Y.Miyamoto) [new]
!! @li      2011-12-11 (H.Yashiro)  [mod] integrate to SCALE3
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_boundary
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_FILECHR
  use gtool_file_h, only: &
     File_HLONG, &
     File_HSHORT, &
     File_HMID,   &
     File_REAL4,  &
     File_REAL8
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_BOUNDARY_setup
  public :: ATMOS_BOUNDARY_read
  public :: ATMOS_BOUNDARY_write
  public :: ATMOS_BOUNDARY_generate
  public :: ATMOS_BOUNDARY_setalpha

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
# include "scale-les.h"
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: I_BND_VELZ = 1 ! reference velocity (z) [m/s]
  integer, public, parameter :: I_BND_VELX = 2 ! reference velocity (x) [m/s]
  integer, public, parameter :: I_BND_VELY = 3 ! reference velocity (y) [m/s]
  integer, public, parameter :: I_BND_POTT = 4 ! reference potential temperature [K]
  integer, public, parameter :: I_BND_QV   = 5 ! reference water vapor [kg/kg]

  real(RP), public, save :: ATMOS_BOUNDARY_var  (KA,IA,JA,5) !> reference container (with HALO)
  real(RP), public, save :: ATMOS_BOUNDARY_alpha(KA,IA,JA,5) ! damping coefficient [0-1]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=IO_FILECHR), private :: ATMOS_BOUNDARY_TYPE          = 'NONE'
  character(len=IO_FILECHR), private :: ATMOS_BOUNDARY_IN_BASENAME   = ''
  character(len=IO_FILECHR), private :: ATMOS_BOUNDARY_OUT_BASENAME  = ''
  character(len=FILE_HLONG), private :: ATMOS_BOUNDARY_OUT_TITLE     = 'SCALE3 BOUNDARY CONDITION'
  character(len=FILE_HLONG), private :: ATMOS_BOUNDARY_OUT_SOURCE    = 'SCALE-LES ver. '//VERSION
  character(len=FILE_HLONG), private :: ATMOS_BOUNDARY_OUT_INSTITUTE = 'AICS/RIKEN'
  logical,                   private :: ATMOS_BOUNDARY_USE_VELZ      = .false. ! read from file?
  logical,                   private :: ATMOS_BOUNDARY_USE_VELX      = .false. ! read from file?
  logical,                   private :: ATMOS_BOUNDARY_USE_VELY      = .false. ! read from file?
  logical,                   private :: ATMOS_BOUNDARY_USE_POTT      = .false. ! read from file?
  logical,                   private :: ATMOS_BOUNDARY_USE_QV        = .false. ! read from file?
  real(RP),                  private :: ATMOS_BOUNDARY_VALUE_VELZ    =  0.0_RP ! w at boundary, 0 [m/s]
  real(RP),                  private :: ATMOS_BOUNDARY_VALUE_VELX    =  5.0_RP ! u at boundary, 5 [m/s]
  real(RP),                  private :: ATMOS_BOUNDARY_VALUE_VELY    =  5.E0_RP ! v at boundary, 5 [m/s]
  real(RP),                  private :: ATMOS_BOUNDARY_VALUE_POTT    = 300.E0_RP! PT at boundary, 300 [K]
  real(RP),                  private :: ATMOS_BOUNDARY_VALUE_QV      = 1.E-3_RP ! QV at boundary, 1.d-3 [kg/kg]

  real(RP),                  private :: ATMOS_BOUNDARY_FRACZ         = 1.0_RP  ! fraction of boundary region for dumping [z]
  real(RP),                  private :: ATMOS_BOUNDARY_FRACX         = 1.0_RP  ! fraction of boundary region for dumping [x]
  real(RP),                  private :: ATMOS_BOUNDARY_FRACY         = 1.0_RP  ! fraction of boundary region for dumping [y]
  real(RP),                  private :: ATMOS_BOUNDARY_tauz          = 75.0_RP ! maximum value for damping tau (z) [s]
  real(RP),                  private :: ATMOS_BOUNDARY_taux          = 75.0_RP ! maximum value for damping tau (x) [s]
  real(RP),                  private :: ATMOS_BOUNDARY_tauy          = 75.0_RP ! maximum value for damping tau (y) [s]

  character(len=File_HSHORT), private :: REF_NAME(5)
  data REF_NAME / 'VELZ_ref','VELX_ref','VELY_ref','POTT_ref','QV_ref' /

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Initialize Boundary Treatment
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only: &
       CONST_UNDEF
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    NAMELIST / PARAM_ATMOS_BOUNDARY / &
       ATMOS_BOUNDARY_TYPE,          &
       ATMOS_BOUNDARY_IN_BASENAME,   &
       ATMOS_BOUNDARY_OUT_BASENAME,  &
       ATMOS_BOUNDARY_OUT_TITLE,     &
       ATMOS_BOUNDARY_OUT_SOURCE,    &
       ATMOS_BOUNDARY_OUT_INSTITUTE, &
       ATMOS_BOUNDARY_USE_VELZ,      &
       ATMOS_BOUNDARY_USE_VELX,      &
       ATMOS_BOUNDARY_USE_VELY,      &
       ATMOS_BOUNDARY_USE_POTT,      &
       ATMOS_BOUNDARY_USE_QV,        &
       ATMOS_BOUNDARY_VALUE_VELZ,    &
       ATMOS_BOUNDARY_VALUE_VELY,    &
       ATMOS_BOUNDARY_VALUE_VELX,    &
       ATMOS_BOUNDARY_VALUE_POTT,    &
       ATMOS_BOUNDARY_VALUE_QV,      &
       ATMOS_BOUNDARY_FRACZ,         &
       ATMOS_BOUNDARY_FRACX,         &
       ATMOS_BOUNDARY_FRACY,         &
       ATMOS_BOUNDARY_tauz,          &
       ATMOS_BOUNDARY_taux,          &
       ATMOS_BOUNDARY_tauy

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Boundary]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_BOUNDARY,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_BOUNDARY. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_BOUNDARY)

    !--- set reference field for boundary
    ATMOS_BOUNDARY_var(:,:,:,:) = CONST_UNDEF

    if     ( ATMOS_BOUNDARY_TYPE == 'NONE' ) then
       ! for backward compatibility
       if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
          write(*,*) 'xxx [obsolete] use ATMOS_BOUNDARY_TYPE'
          call ATMOS_BOUNDARY_read
       elseif( ATMOS_BOUNDARY_OUT_BASENAME /= '' ) then
          write(*,*) 'xxx [obsolete] use ATMOS_BOUNDARY_TYPE'
          call ATMOS_BOUNDARY_generate
       endif
    elseif ( ATMOS_BOUNDARY_TYPE == 'INIT' ) then
       call ATMOS_BOUNDARY_setinitval
    elseif ( ATMOS_BOUNDARY_TYPE == 'FILE' ) then
       if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
          call ATMOS_BOUNDARY_read
       else
          write(*,*) 'xxx You need specify ATMOS_BOUNDARY_IN_BASENAME'
          call PRC_MPIstop
       endif
    elseif ( ATMOS_BOUNDARY_TYPE == 'CONST' ) then
       call ATMOS_BOUNDARY_generate
    else
       write(*,*) 'xxx ATMOS_BOUNDARY_TYPE is invalid'
       call PRC_MPIstop
    endif

    if( ATMOS_BOUNDARY_OUT_BASENAME /= '' ) then
       call ATMOS_BOUNDARY_write
    endif

    call ATMOS_BOUNDARY_setalpha

    call COMM_vars8( ATMOS_BOUNDARY_var  (:,:,:,I_BND_QV), 1 )
    call COMM_vars8( ATMOS_BOUNDARY_alpha(:,:,:,i_BND_QV), 2 )
    call COMM_wait ( ATMOS_BOUNDARY_var  (:,:,:,I_BND_QV), 1 )
    call COMM_wait ( ATMOS_BOUNDARY_alpha(:,:,:,i_BND_QV), 2 )

    return
  end subroutine ATMOS_BOUNDARY_setup

  !-----------------------------------------------------------------------------
  !> Calc dumping coefficient alpha
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_setalpha
    use mod_const, only: &
       CONST_UNDEF, &
       PI => CONST_PI
    use mod_grid, only: &
       CBFZ => GRID_CBFZ, &
       CBFX => GRID_CBFX, &
       CBFY => GRID_CBFY, &
       FBFZ => GRID_FBFZ, &
       FBFX => GRID_FBFX, &
       FBFY => GRID_FBFY

    real(RP) :: coef, alpha
    real(RP) :: ee1, ee2

    integer :: i, j, k
    !---------------------------------------------------------------------------

    !--- set damping coefficient
    ATMOS_BOUNDARY_alpha(:,:,:,:) = 0.0_RP

    if ( ATMOS_BOUNDARY_tauz <= 0.0_RP ) then ! invalid tau
       coef = 0.0_RP
    else
       coef = 1.0_RP / ATMOS_BOUNDARY_tauz
    endif

    ! check invalid fraction
    if ( ATMOS_BOUNDARY_FRACZ < 0.0_RP ) ATMOS_BOUNDARY_FRACZ = 0.0_RP
    if ( ATMOS_BOUNDARY_FRACZ > 1.0_RP ) ATMOS_BOUNDARY_FRACZ = 1.0_RP
    if ( ATMOS_BOUNDARY_FRACX < 0.0_RP ) ATMOS_BOUNDARY_FRACX = 0.0_RP
    if ( ATMOS_BOUNDARY_FRACX > 1.0_RP ) ATMOS_BOUNDARY_FRACX = 1.0_RP
    if ( ATMOS_BOUNDARY_FRACY < 0.0_RP ) ATMOS_BOUNDARY_FRACY = 0.0_RP
    if ( ATMOS_BOUNDARY_FRACY > 1.0_RP ) ATMOS_BOUNDARY_FRACY = 1.0_RP

    do k = KS, KE
       ee1 = CBFZ(k)
       if ( ee1 <= 1.0_RP - ATMOS_BOUNDARY_FRACZ ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + ATMOS_BOUNDARY_FRACZ ) / ATMOS_BOUNDARY_FRACZ
       endif

       if    ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELX) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELY) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_QV  ) )
       elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELX) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELY) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_QV  ) )
       endif
    enddo

    do k = KS, KE-1
       ee2 = FBFZ(k)
       if ( ee2 <= 1.0_RP - ATMOS_BOUNDARY_FRACZ ) then
          ee2 = 0.0_RP
       else
          ee2 = ( ee2 - 1.0_RP + ATMOS_BOUNDARY_FRACZ ) / ATMOS_BOUNDARY_FRACZ
       endif

       if    ( ee2 > 0.0_RP .AND. ee2 <= 0.5_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP - cos( ee2*PI ) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELZ) )
       elseif( ee2 > 0.5_RP .AND. ee2 <= 1.0_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP + sin( (ee2-0.5_RP)*PI ) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELZ) )
       endif
    enddo

    if ( ATMOS_BOUNDARY_taux <= 0.0_RP ) then ! invalid tau
       coef = 0.0_RP
    else
       coef = 1.0_RP / ATMOS_BOUNDARY_taux
    endif

    do i = IS, IE
       ee1 = CBFX(i)
       ee2 = FBFX(i)
       if ( ee1 <= 1.0_RP - ATMOS_BOUNDARY_FRACX ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + ATMOS_BOUNDARY_FRACX ) / ATMOS_BOUNDARY_FRACX
       endif
       if ( ee2 <= 1.0_RP - ATMOS_BOUNDARY_FRACX ) then
          ee2 = 0.0_RP
       else
          ee2 = ( ee2 - 1.0_RP + ATMOS_BOUNDARY_FRACX ) / ATMOS_BOUNDARY_FRACX
       endif

       if ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELZ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELY) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_QV  ) )
       elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELZ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELY) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_QV  ) )
       endif

       if ( ee2 > 0.0_RP .AND. ee2 <= 0.5_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP - cos( ee2*PI ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELX) )
       elseif( ee2 > 0.5_RP .AND. ee2 <= 1.0_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP + sin( (ee2-0.5_RP)*PI ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELX) )
       endif
    enddo

    if ( ATMOS_BOUNDARY_tauy <= 0.0_RP ) then ! invalid tau
       coef = 0.0_RP
    else
       coef = 1.0_RP / ATMOS_BOUNDARY_tauy
    endif

    do j = JS, JE
       ee1 = CBFY(j)
       ee2 = FBFY(j)
       if ( ee1 <= 1.0_RP - ATMOS_BOUNDARY_FRACY ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + ATMOS_BOUNDARY_FRACY ) / ATMOS_BOUNDARY_FRACY
       endif
       if ( ee2 <= 1.0_RP - ATMOS_BOUNDARY_FRACY ) then
          ee2 = 0.0_RP
       else
          ee2 = ( ee2 - 1.0_RP + ATMOS_BOUNDARY_FRACY ) / ATMOS_BOUNDARY_FRACY
       endif

       if ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELZ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELX) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_QV  ) )
       elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELZ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELX) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_QV  ) )
       endif

       if ( ee2 > 0.0_RP .AND. ee2 <= 0.5_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP - cos( ee2*PI ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELY) )
       elseif( ee2 > 0.5_RP .AND. ee2 <= 1.0_RP ) then
          alpha = coef * 0.5_RP * ( 1.0_RP + sin( (ee2-0.5_RP)*PI ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELY) )
       endif
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) == CONST_UNDEF ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELX) = 0.0_RP
          endif
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_VELY) == CONST_UNDEF ) then
            ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELY) = 0.0_RP
          endif
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_POTT) == CONST_UNDEF ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_POTT) = 0.0_RP
          endif
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_QV) == CONST_UNDEF ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_QV) = 0.0_RP
          endif
       enddo

       do k = KS, KE-1
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_VELZ) == CONST_UNDEF ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELZ) = 0.0_RP
         endif
      enddo
    enddo
    enddo

    return
  end subroutine ATMOS_BOUNDARY_setalpha

  !-----------------------------------------------------------------------------
  !> Read boundary data
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_setinitval
    use mod_const, only: &
       CONST_UNDEF
    use mod_grid, only: &
       CZ_mask => GRID_CZ_mask, &
       CX_mask => GRID_CX_mask, &
       CY_mask => GRID_CX_mask
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_comm, only: &
       COMM_vars, &
       COMM_wait
    implicit none

    integer :: i, j, k, iv
    !---------------------------------------------------------------------------

    ATMOS_BOUNDARY_var(:,:,:,I_BND_VELZ) = CONST_UNDEF
    ATMOS_BOUNDARY_var(:,:,:,I_BND_VELY) = CONST_UNDEF
    ATMOS_BOUNDARY_var(:,:,:,I_BND_VELX) = CONST_UNDEF
    ATMOS_BOUNDARY_var(:,:,:,I_BND_POTT) = CONST_UNDEF
    ATMOS_BOUNDARY_var(:,:,:,I_BND_QV  ) = CONST_UNDEF

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( .not. ( CZ_mask(k) .AND. CX_mask(i) .AND. CY_mask(j) ) ) then ! Buffer Layer
          if ( ATMOS_BOUNDARY_USE_VELZ ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_VELZ) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j) + DENS(k,i,j) )
          endif
          if ( ATMOS_BOUNDARY_USE_VELX ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i+1,j) + DENS(k,i,j) )
          endif
          if ( ATMOS_BOUNDARY_USE_VELY ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_VELY) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1) + DENS(k,i,j) )
          endif
          if ( ATMOS_BOUNDARY_USE_POTT ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_POTT) = RHOT(k,i,j) / DENS(k,i,j)
          endif
          if ( ATMOS_BOUNDARY_USE_QV ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_QV  ) = QTRC(k,i,j,I_QV)
          endif
       endif
    enddo
    enddo
    enddo

    ! fill KHALO
    do iv = I_BND_VELZ, I_BND_QV
    do j  = JS, JE
    do i  = IS, IE
       ATMOS_BOUNDARY_var(   1:KS-1,i,j,iv) = ATMOS_BOUNDARY_var(KS,i,j,iv)
       ATMOS_BOUNDARY_var(KE+1:KA,  i,j,iv) = ATMOS_BOUNDARY_var(KE,i,j,iv)
    enddo
    enddo
    enddo

    ! fill IHALO & JHALO
    do iv = I_BND_VELZ, I_BND_QV
       call COMM_vars( ATMOS_BOUNDARY_var(:,:,:,iv), iv )
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_wait( ATMOS_BOUNDARY_var(:,:,:,iv), iv )
    enddo

    return
  end subroutine ATMOS_BOUNDARY_setinitval

  !-----------------------------------------------------------------------------
  !> Read boundary data
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_read
    use gtool_file, only: &
       FileRead
    use mod_comm, only: &
       COMM_vars, &
       COMM_wait
    use mod_process, only: &
         PRC_myrank
    implicit none

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname

    integer :: iv, i, j
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_IN_BASENAME

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       call FileRead( reference_atmos(:,:,:), bname, 'VELZ', 1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_VELZ) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_VELX ) then
       call FileRead( reference_atmos(:,:,:), bname, 'VELX', 1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_VELX) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_VELY ) then
       call FileRead( reference_atmos(:,:,:), bname, 'VELY', 1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_VELY) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_POTT ) then
       call FileRead( reference_atmos(:,:,:), bname, 'POTT', 1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_POTT) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_QV ) then
       call FileRead( reference_atmos(:,:,:), bname, 'QV',   1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_QV) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    ! fill IHALO & JHALO
    do iv = I_BND_VELZ, I_BND_QV
       call COMM_vars( ATMOS_BOUNDARY_var(:,:,:,iv), iv )
    enddo
    do iv = I_BND_VELZ, I_BND_QV
       call COMM_wait( ATMOS_BOUNDARY_var(:,:,:,iv), iv )
    enddo

    ! fill KHALO
    do iv = I_BND_VELZ, I_BND_QV
    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_var(   1:KS-1,i,j,iv) = ATMOS_BOUNDARY_var(KS,i,j,iv)
       ATMOS_BOUNDARY_var(KE+1:KA,  i,j,iv) = ATMOS_BOUNDARY_var(KE,i,j,iv)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_BOUNDARY_read

  !-----------------------------------------------------------------------------
  !> Write boundary data
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_write
    use mod_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use gtool_file, only: &
       FileCreate, &
       FileAddVariable, &
       FilePutAxis, &
       FileWrite, &
       FileClose
    use mod_grid, only:  &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY
    implicit none

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    logical :: fileexisted
    integer :: fid, vid
    integer :: dtype

    integer :: rankidx(2)
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_OUT_BASENAME

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)
    call FileCreate( fid, fileexisted,                          & ! (out)
         bname,                                                 & ! (in)
         ATMOS_BOUNDARY_OUT_TITLE,                              & ! (in)
         ATMOS_BOUNDARY_OUT_SOURCE,                             & ! (in)
         ATMOS_BOUNDARY_OUT_INSTITUTE,                          & ! (in)
         (/'z','x','y'/), (/KMAX,IMAX,JMAX/), (/'Z','X','Y'/),  & ! (in)
         (/'m','m','m'/), (/File_REAL4,File_REAL4,File_REAL4/), & ! (in)
         PRC_master, PRC_myrank, rankidx                        ) ! (in)

    call FilePutAxis(fid, 'z', GRID_CZ(KS:KE))
    call FilePutAxis(fid, 'x', GRID_CX(IS:IE))
    call FilePutAxis(fid, 'y', GRID_CY(JS:JE))

    if ( RP == 8 ) then
       dtype = File_REAL8
    else if ( RP == 4 ) then
       dtype = File_REAL4
    endif

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       call FileAddVariable( vid,                      & ! (out)
            fid, 'VELZ', 'Reference Velocit w', 'm/s', & ! (in)
            (/'z','x','y'/), dtype                     ) ! (in)
       reference_atmos(1:KMAX,1:IMAX,1:JMAX) = ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_VELZ)
       call FileWrite( vid, reference_atmos(:,:,:), NOWSEC, NOWSEC )
    endif

    if ( ATMOS_BOUNDARY_USE_VELX ) then
       call FileAddVariable( vid,                      & ! (out)
            fid, 'VELX', 'Reference Velocit u', 'm/s', & ! (in)
            (/'z','x','y'/), dtype                     ) ! (in)
       reference_atmos(1:KMAX,1:IMAX,1:JMAX) = ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_VELX)
       call FileWrite( vid, reference_atmos(:,:,:), NOWSEC, NOWSEC )
    endif

    if ( ATMOS_BOUNDARY_USE_VELY ) then
       call FileAddVariable( vid,                      & ! (out)
            fid, 'VELY', 'Reference Velocit v', 'm/s', & ! (in)
            (/'z','x','y'/), dtype                     ) ! (in)
       reference_atmos(1:KMAX,1:IMAX,1:JMAX) = ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_VELY)
       call FileWrite( vid, reference_atmos(:,:,:), NOWSEC, NOWSEC )
    endif

    if ( ATMOS_BOUNDARY_USE_POTT ) then
       call FileAddVariable( vid,             & ! (out)
            fid, 'POTT', 'Reference PT', 'K', & ! (in)
            (/'z','x','y'/), dtype            ) ! (in)
       reference_atmos(1:KMAX,1:IMAX,1:JMAX) = ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_POTT)
       call FileWrite( vid, reference_atmos(:,:,:), NOWSEC, NOWSEC )
    endif

    if ( ATMOS_BOUNDARY_USE_QV ) then
       call FileAddVariable( vid,                        & ! (out)
            fid, 'QV', 'Reference water vapor', 'kg/kg', & ! (in)
            (/'z','x','y'/), File_REAL8                  ) ! (in)
       reference_atmos(1:KMAX,1:IMAX,1:JMAX) = ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_QV)
       call FileWrite( vid, reference_atmos(:,:,:), NOWSEC, NOWSEC )
    endif

    call FileClose( fid )

    return
  end subroutine ATMOS_BOUNDARY_write

  !-----------------------------------------------------------------------------
  !> generate boundary data (temporal)
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_generate
    use mod_const, only: &
       CONST_UNDEF
    use mod_grid, only: &
       CZ_mask => GRID_CZ_mask, &
       CX_mask => GRID_CX_mask, &
       CY_mask => GRID_CX_mask
    use mod_comm, only: &
       COMM_vars, &
       COMM_wait
    implicit none

    integer :: i, j, k, iv
    !---------------------------------------------------------------------------

    ATMOS_BOUNDARY_var(:,:,:,I_BND_VELZ) = CONST_UNDEF
    ATMOS_BOUNDARY_var(:,:,:,I_BND_VELY) = CONST_UNDEF
    ATMOS_BOUNDARY_var(:,:,:,I_BND_VELX) = CONST_UNDEF
    ATMOS_BOUNDARY_var(:,:,:,I_BND_POTT) = CONST_UNDEF
    ATMOS_BOUNDARY_var(:,:,:,I_BND_QV  ) = CONST_UNDEF

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( .not. ( CZ_mask(k) .AND. CX_mask(i) .AND. CY_mask(j) ) ) then ! Buffer Layer
          if ( ATMOS_BOUNDARY_USE_VELZ ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_VELZ) = ATMOS_BOUNDARY_VALUE_VELZ
          endif
          if ( ATMOS_BOUNDARY_USE_VELY ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_VELY) = ATMOS_BOUNDARY_VALUE_VELY
          endif
          if ( ATMOS_BOUNDARY_USE_VELX ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = ATMOS_BOUNDARY_VALUE_VELX
          endif
          if ( ATMOS_BOUNDARY_USE_POTT ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_POTT) = ATMOS_BOUNDARY_VALUE_POTT
          endif
          if ( ATMOS_BOUNDARY_USE_QV ) then
            ATMOS_BOUNDARY_var(k,i,j,I_BND_QV  ) = ATMOS_BOUNDARY_VALUE_QV
          endif
       endif
    enddo
    enddo
    enddo

!    do j = JS-1, JE+1
!    do i = IS-1, IE+1
!       do k = KS, KE
!          if ( CZ_mask(k) .AND. CX_mask(i) ) then ! Inner Area
!             ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = CONST_UNDEF
!          else                                    ! Buffer Area
!             ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = FBFX(i) * ATMOS_BOUNDARY_VALUE_VELX &
!                                                  * ( 1.0_RP - CBFZ(k) )
!          endif
!       enddo
!    enddo
!    enddo
!    ATMOS_BOUNDARY_var(:,:,:,I_BND_VELX) = CONST_UNDEF

    ! fill KHALO
    do iv = I_BND_VELZ, I_BND_QV
    do j  = JS, JE
    do i  = IS, IE
       ATMOS_BOUNDARY_var(   1:KS-1,i,j,iv) = ATMOS_BOUNDARY_var(KS,i,j,iv)
       ATMOS_BOUNDARY_var(KE+1:KA,  i,j,iv) = ATMOS_BOUNDARY_var(KE,i,j,iv)
    enddo
    enddo
    enddo

    ! fill IHALO & JHALO
    do iv = I_BND_VELZ, I_BND_QV
       call COMM_vars( ATMOS_BOUNDARY_var(:,:,:,iv), iv )
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_wait( ATMOS_BOUNDARY_var(:,:,:,iv), iv )
    enddo

    return
  end subroutine ATMOS_BOUNDARY_generate

end module mod_atmos_boundary
