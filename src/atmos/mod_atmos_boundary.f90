!-------------------------------------------------------------------------------
!> module Atmosphere / Boundary treatment
!!
!! @par Description
!!          Boundary treatment of model domain
!!          Additional forcing, Sponge layer, rayleigh dumping
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] integrate
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_boundary
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FILECHR
  use mod_fileio_h, only: &
     FIO_HSHORT
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_BOUNDARY_setup
  public :: ATMOS_BOUNDARY
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  real(8), public, allocatable, save :: velx_ref(:,:,:) ! reference velocity (x) [m/s]
  real(8), public, allocatable, save :: vely_ref(:,:,:) ! reference velocity (y) [m/s]
  real(8), public, allocatable, save :: velz_ref(:,:,:) ! reference velocity (z) [m/s]
  real(8), public, allocatable, save :: pott_ref(:,:,:) ! reference potential temperature [K]
  real(8), public, allocatable, save :: qv_ref  (:,:,:) ! reference water vapor [kg/kg]

  real(8), public, allocatable, save :: DAMP_alphau(:,:,:) ! damping coefficient for u  [0-1]
  real(8), public, allocatable, save :: DAMP_alphav(:,:,:) ! damping coefficient for v  [0-1]
  real(8), public, allocatable, save :: DAMP_alphaw(:,:,:) ! damping coefficient for w  [0-1]
  real(8), public, allocatable, save :: DAMP_alphat(:,:,:) ! damping coefficient for pt [0-1]
  real(8), public, allocatable, save :: DAMP_alphaq(:,:,:) ! damping coefficient for qv and scalars [0-1]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  real(8),                   private, allocatable, save :: atmos_refvar(:,:,:,:)  !> reference container (with HALO)
  character(len=FIO_HSHORT), private,         parameter :: REF_NAME(5) &
                      = (/'VELX_ref','VELY_ref','VELZ_ref','POTT_ref','QV_ref'/)

  character(len=IO_FILECHR), private, save :: ATMOS_REFERENCE_IN_BASENAME = 'refvar_in'
  logical,                   private, save :: ref_velx = .false. ! read from file?
  logical,                   private, save :: ref_vely = .false. ! read from file?
  logical,                   private, save :: ref_velz = .false. ! read from file?
  logical,                   private, save :: ref_pott = .false. ! read from file?
  logical,                   private, save :: ref_qv   = .false. ! read from file?
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Initialize Boundary Treatment
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop, &
       PRC_NUM_X,   &
       PRC_NUM_Y
    use mod_const, only: &
       CONST_UNDEF8,   &
       PI => CONST_PI
    use mod_grid, only : &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       KMAX => GRID_KMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       WS   => GRID_WS,   &
       WE   => GRID_WE,   &
       GRID_DX, &
       GRID_DY, &
       GRID_DZ, &
       GRID_CX, &
       GRID_FX, &
       GRID_CY, &
       GRID_FY, &
       GRID_CZ, &
       GRID_FZ
    implicit none

    real(8) :: ATMOS_BOUNDARY_wallsponge_dx    =    0.D0 ! thickness of sponge damping layer (x) [m]
    real(8) :: ATMOS_BOUNDARY_wallsponge_dy    =    0.D0 ! thickness of sponge damping layer (y) [m]
    real(8) :: ATMOS_BOUNDARY_uppersponge_dz   = 1500.D0 ! thickness of sponge damping layer (z) [m]
    real(8) :: ATMOS_BOUNDARY_wallsponge_taux  =   75.D0 ! maximum value for damping tau (x) [s]
    real(8) :: ATMOS_BOUNDARY_wallsponge_tauy  =   75.D0 ! maximum value for damping tau (y) [s]
    real(8) :: ATMOS_BOUNDARY_uppersponge_tauz =   75.D0 ! maximum value for damping tau (z) [s]

    NAMELIST / PARAM_ATMOS_BOUNDARY / &
       ATMOS_BOUNDARY_wallsponge_dx,    &
       ATMOS_BOUNDARY_wallsponge_dy,    &
       ATMOS_BOUNDARY_uppersponge_dz,   &
       ATMOS_BOUNDARY_wallsponge_taux,  &
       ATMOS_BOUNDARY_wallsponge_tauy,  &
       ATMOS_BOUNDARY_uppersponge_tauz, &
       ATMOS_REFERENCE_IN_BASENAME

    real(8) :: ridge, sponge, coef, alpha
    real(8) :: ee1, ee2

    integer :: ierr
    integer :: i, j, k
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
    allocate( atmos_refvar(IA,JA,KA,5) ); atmos_refvar(:,:,:,:) = CONST_UNDEF8

    allocate( velx_ref(IA,JA,KA) ); velx_ref(:,:,:) = CONST_UNDEF8
    allocate( vely_ref(IA,JA,KA) ); vely_ref(:,:,:) = CONST_UNDEF8
    allocate( velz_ref(IA,JA,KA) ); velz_ref(:,:,:) = CONST_UNDEF8
    allocate( pott_ref(IA,JA,KA) ); pott_ref(:,:,:) = CONST_UNDEF8
    allocate( qv_ref  (IA,JA,KA) ); qv_ref  (:,:,:) = CONST_UNDEF8

    call ATMOS_BOUNDARY_reference_read

    !--- set damping coefficient
    allocate( DAMP_alphau(IA,JA,KA) ); DAMP_alphau(:,:,:) = 0.D0
    allocate( DAMP_alphav(IA,JA,KA) ); DAMP_alphav(:,:,:) = 0.D0
    allocate( DAMP_alphaw(IA,JA,KA) ); DAMP_alphaw(:,:,:) = 0.D0
    allocate( DAMP_alphat(IA,JA,KA) ); DAMP_alphat(:,:,:) = 0.D0
    allocate( DAMP_alphaq(IA,JA,KA) ); DAMP_alphaq(:,:,:) = 0.D0

    ridge  = GRID_DX*real(IMAX*PRC_NUM_X,kind=8)
    sponge = GRID_DX*real(IMAX*PRC_NUM_X,kind=8) - ATMOS_BOUNDARY_wallsponge_dx
    coef   = 1.D0 / ATMOS_BOUNDARY_wallsponge_taux

    do i = IS, IE+1

       if ( GRID_FX(i) <= sponge ) cycle

       ee1 = ( GRID_CX(i)-sponge ) / ( ridge-sponge )
       ee2 = ( GRID_FX(i)-sponge ) / ( ridge-sponge )

       if ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          DAMP_alphav(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphav(i,JS:JE,KS  :KE  ) )
          DAMP_alphat(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphat(i,JS:JE,KS  :KE  ) )
          DAMP_alphaq(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphaq(i,JS:JE,KS  :KE  ) )
          DAMP_alphaw(i,JS:JE,WS+1:WE-1) = max( alpha, DAMP_alphaw(i,JS:JE,WS+1:WE-1) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          DAMP_alphav(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphav(i,JS:JE,KS  :KE  ) )
          DAMP_alphat(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphat(i,JS:JE,KS  :KE  ) )
          DAMP_alphaq(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphaq(i,JS:JE,KS  :KE  ) )
          DAMP_alphaw(i,JS:JE,WS+1:WE-1) = max( alpha, DAMP_alphaw(i,JS:JE,WS+1:WE-1) )
       endif

       if ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          DAMP_alphau(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphau(i,JS:JE,KS  :KE  ) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          DAMP_alphau(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphau(i,JS:JE,KS  :KE  ) )
       endif
    enddo

    ridge  = 0
    sponge = ATMOS_BOUNDARY_wallsponge_dx
    coef   = 1.D0 / ATMOS_BOUNDARY_wallsponge_taux

    do i = IS-1, IE

       if ( GRID_FX(i) >= sponge ) cycle

       ee1 = ( GRID_CX(i)-sponge ) / ( ridge-sponge )
       ee2 = ( GRID_FX(i)-sponge ) / ( ridge-sponge )

       if ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          DAMP_alphav(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphav(i,JS:JE,KS  :KE  ) )
          DAMP_alphat(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphat(i,JS:JE,KS  :KE  ) )
          DAMP_alphaq(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphaq(i,JS:JE,KS  :KE  ) )
          DAMP_alphaw(i,JS:JE,WS+1:WE-1) = max( alpha, DAMP_alphaw(i,JS:JE,WS+1:WE-1) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          DAMP_alphav(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphav(i,JS:JE,KS  :KE  ) )
          DAMP_alphat(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphat(i,JS:JE,KS  :KE  ) )
          DAMP_alphaq(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphaq(i,JS:JE,KS  :KE  ) )
          DAMP_alphaw(i,JS:JE,WS+1:WE-1) = max( alpha, DAMP_alphaw(i,JS:JE,WS+1:WE-1) )
       endif

       if ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          DAMP_alphau(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphau(i,JS:JE,KS  :KE  ) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          DAMP_alphau(i,JS:JE,KS  :KE  ) = max( alpha, DAMP_alphau(i,JS:JE,KS  :KE  ) )
       endif
    enddo

    ridge  = GRID_DY*real(JMAX*PRC_NUM_Y,kind=8)
    sponge = GRID_DY*real(JMAX*PRC_NUM_Y,kind=8) - ATMOS_BOUNDARY_wallsponge_dy
    coef   = 1.D0 / ATMOS_BOUNDARY_wallsponge_tauy

    do j = JS, JE+1

       if ( GRID_FY(j) <= sponge ) cycle

       ee1 = ( GRID_CY(j)-sponge ) / ( ridge-sponge )
       ee2 = ( GRID_FY(j)-sponge ) / ( ridge-sponge )

       if ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          DAMP_alphau(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphau(IS:IE,j,KS  :KE  ) )
          DAMP_alphaw(IS:IE,j,WS+1:WE-1) = max( alpha, DAMP_alphaw(IS:IE,j,WS+1:WE-1) )
          DAMP_alphat(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphat(IS:IE,j,KS  :KE  ) )
          DAMP_alphaq(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphaq(IS:IE,j,KS  :KE  ) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          DAMP_alphau(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphau(IS:IE,j,KS  :KE  ) )
          DAMP_alphaw(IS:IE,j,WS+1:WE-1) = max( alpha, DAMP_alphaw(IS:IE,j,WS+1:WE-1) )
          DAMP_alphat(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphat(IS:IE,j,KS  :KE  ) )
          DAMP_alphaq(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphaq(IS:IE,j,KS  :KE  ) )
       endif

       if ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          DAMP_alphav(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphav(IS:IE,j,KS  :KE  ) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          DAMP_alphav(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphav(IS:IE,j,KS  :KE  ) )
       endif
    enddo

    ridge  = 0
    sponge = ATMOS_BOUNDARY_wallsponge_dy
    coef   = 1.D0 / ATMOS_BOUNDARY_wallsponge_tauy

    do j = JS-1, JE

       if ( GRID_FY(j) >= sponge ) cycle

       ee1 = ( GRID_CY(j)-sponge ) / ( ridge-sponge )
       ee2 = ( GRID_FY(j)-sponge ) / ( ridge-sponge )

       if ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          DAMP_alphau(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphau(IS:IE,j,KS  :KE  ) )
          DAMP_alphaw(IS:IE,j,WS+1:WE-1) = max( alpha, DAMP_alphaw(IS:IE,j,WS+1:WE-1) )
          DAMP_alphat(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphat(IS:IE,j,KS  :KE  ) )
          DAMP_alphaq(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphaq(IS:IE,j,KS  :KE  ) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          DAMP_alphau(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphau(IS:IE,j,KS  :KE  ) )
          DAMP_alphaw(IS:IE,j,WS+1:WE-1) = max( alpha, DAMP_alphaw(IS:IE,j,WS+1:WE-1) )
          DAMP_alphat(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphat(IS:IE,j,KS  :KE  ) )
          DAMP_alphaq(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphaq(IS:IE,j,KS  :KE  ) )
       endif

       if ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          DAMP_alphav(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphav(IS:IE,j,KS  :KE  ) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          DAMP_alphav(IS:IE,j,KS  :KE  ) = max( alpha, DAMP_alphav(IS:IE,j,KS  :KE  ) )
       endif
    enddo

    ridge  = GRID_DZ*real(KMAX,kind=8)
    sponge = GRID_DZ*real(KMAX,kind=8) - ATMOS_BOUNDARY_uppersponge_dz
    coef   = 1.D0 / ATMOS_BOUNDARY_uppersponge_tauz

    do k = KS, KE
       if ( GRID_FZ(k) <= sponge ) cycle

       ee1 = ( GRID_CZ(k)-sponge ) / ( ridge-sponge )

       if    ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          DAMP_alphau(IS:IE,JS:JE,k) = max( alpha, DAMP_alphau(IS:IE,JS:JE,k) )
          DAMP_alphav(IS:IE,JS:JE,k) = max( alpha, DAMP_alphav(IS:IE,JS:JE,k) )
          DAMP_alphat(IS:IE,JS:JE,k) = max( alpha, DAMP_alphat(IS:IE,JS:JE,k) )
          DAMP_alphaq(IS:IE,JS:JE,k) = max( alpha, DAMP_alphaq(IS:IE,JS:JE,k) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          DAMP_alphau(IS:IE,JS:JE,k) = max( alpha, DAMP_alphau(IS:IE,JS:JE,k) )
          DAMP_alphav(IS:IE,JS:JE,k) = max( alpha, DAMP_alphav(IS:IE,JS:JE,k) )
          DAMP_alphat(IS:IE,JS:JE,k) = max( alpha, DAMP_alphat(IS:IE,JS:JE,k) )
          DAMP_alphaq(IS:IE,JS:JE,k) = max( alpha, DAMP_alphaq(IS:IE,JS:JE,k) )
       endif
    enddo

    do k = WS+1, WE-1
       if ( GRID_FZ(k) <= sponge ) cycle

       ee2 = ( GRID_FZ(k)-sponge ) / ( ridge-sponge )

       if    ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          DAMP_alphaw(IS:IE,JS:JE,k) = max( alpha, DAMP_alphaw(IS:IE,JS:JE,k) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          DAMP_alphaw(IS:IE,JS:JE,k) = max( alpha, DAMP_alphaw(IS:IE,JS:JE,k) )
       endif
    enddo

    do k = KS,   KE
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       if ( velx_ref(i,j,k) == CONST_UNDEF8 ) then
          DAMP_alphau(i,j,k) = 0.D0
       endif
       if ( vely_ref(i,j,k) == CONST_UNDEF8 ) then
          DAMP_alphav(i,j,k) = 0.D0
       endif
       if ( pott_ref(i,j,k) == CONST_UNDEF8 ) then
          DAMP_alphat(i,j,k) = 0.D0
       endif
       if ( qv_ref  (i,j,k) == CONST_UNDEF8 ) then
          DAMP_alphaq(i,j,k) = 0.D0
       endif
    enddo
    enddo
    enddo

    do k = WS+1, WE-1
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       if ( velz_ref(i,j,k) == CONST_UNDEF8 ) then
          DAMP_alphaw(i,j,k) = 0.D0
       endif
    enddo
    enddo
    enddo

!    do k = 1, KA
!       if( IO_L ) write(IO_FID_LOG,*) 'DAMPING w(face) at k=',k,'+1/2 : ',DAMP_alphaw(50,50,k),velz_ref(50,50,k)
!    enddo
!    do I = 1, IA
!       if( IO_L ) write(IO_FID_LOG,*) 'DAMPING u       at i=',i,'     : ',DAMP_alphau(i,50,KS)
!    enddo

    return
  end subroutine ATMOS_BOUNDARY_setup

  !-----------------------------------------------------------------------------
  !> Boundary Treatment
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY( dens,   momx,   momy,   momz,   pott,  &
                                     momx_t, momy_t, momz_t, pott_t )
    use mod_grid, only : &
       IA  => GRID_IA, &
       JA  => GRID_JA, &
       KA  => GRID_KA, &
       IS  => GRID_IS, &
       IE  => GRID_IE, &
       JS  => GRID_JS, &
       JE  => GRID_JE, &
       KS  => GRID_KS, &
       KE  => GRID_KE
    implicit none

    ! prognostic value
    real(8), intent(in)    :: dens(IA,JA,KA)      ! density [kg/m3]
    real(8), intent(in)    :: momx(IA,JA,KA)      ! momentum (x) [kg/m3 * m/s]
    real(8), intent(in)    :: momy(IA,JA,KA)      ! momentum (y) [kg/m3 * m/s]
    real(8), intent(in)    :: momz(IA,JA,KA)      ! momentum (z) [kg/m3 * m/s]
    real(8), intent(in)    :: pott(IA,JA,KA)      ! potential temperature [K]

    ! prognostic tendency
    real(8), intent(inout) :: momx_t(IA,JA,KA)
    real(8), intent(inout) :: momy_t(IA,JA,KA)
    real(8), intent(inout) :: momz_t(IA,JA,KA)
    real(8), intent(inout) :: pott_t(IA,JA,KA)
    !---------------------------------------------------------------------------

    ! do nothing at here
    ! dumping is done in dynamical step

    return
  end subroutine ATMOS_BOUNDARY

  !-----------------------------------------------------------------------------
  !> Read restart of atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_reference_read
    use mod_const, only: &
       CONST_UNDEF8
    use mod_comm, only: &
       COMM_vars, &
       COMM_stats
    use mod_grid, only : &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       KMAX => GRID_KMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE
    use mod_fileio, only: &
       FIO_input
    implicit none

    real(8), allocatable :: reference_atmos(:,:,:) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=8)          :: lname

    integer :: iv, i, j
    !---------------------------------------------------------------------------

    allocate( reference_atmos(IMAX,JMAX,KMAX) ); reference_atmos(:,:,:) = CONST_UNDEF8

    bname = ATMOS_REFERENCE_IN_BASENAME
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    if ( ref_velx ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'VELX', lname, 1, KMAX, 1 )
       atmos_refvar(IS:IE,JS:JE,KS:KE,1) = reference_atmos(1:IMAX,1:JMAX,1:KMAX)
    endif

    if ( ref_vely ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'VELY', lname, 1, KMAX, 1 )
       atmos_refvar(IS:IE,JS:JE,KS:KE,2) = reference_atmos(1:IMAX,1:JMAX,1:KMAX)
    endif

    if ( ref_velz ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'VELZ', lname, 1, KMAX, 1 )
       atmos_refvar(IS:IE,JS:JE,KS:KE,3) = reference_atmos(1:IMAX,1:JMAX,1:KMAX)
    endif

    if ( ref_pott ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'POTT', lname, 1, KMAX, 1 )
       atmos_refvar(IS:IE,JS:JE,KS:KE,4) = reference_atmos(1:IMAX,1:JMAX,1:KMAX)
    endif

    if ( ref_qv ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'QV',   lname, 1, KMAX, 1 )
       atmos_refvar(IS:IE,JS:JE,KS:KE,5) = reference_atmos(1:IMAX,1:JMAX,1:KMAX)
    endif

    deallocate( reference_atmos )

    ! fill IHALO & JHALO
    call COMM_vars( atmos_refvar(:,:,:,:) )

    ! fill KHALO
    do iv = 1, 5
    do j  = 1, JA
    do i  = 1, IA
       atmos_refvar(i,j,   1:KS-1,iv) = atmos_refvar(i,j,KS,iv)
       atmos_refvar(i,j,KE+1:KA,  iv) = atmos_refvar(i,j,KE,iv)
    enddo
    enddo
    enddo
    atmos_refvar(:,:,KA-10:KA,3) = 0.D0

    call COMM_stats( atmos_refvar(:,:,:,:), REF_NAME(:) )

    velx_ref(:,:,:) = atmos_refvar(:,:,:,1)
    vely_ref(:,:,:) = atmos_refvar(:,:,:,2)
    velz_ref(:,:,:) = atmos_refvar(:,:,:,3)
    pott_ref(:,:,:) = atmos_refvar(:,:,:,4)
    qv_ref  (:,:,:) = atmos_refvar(:,:,:,5)

    return
  end subroutine ATMOS_BOUNDARY_reference_read

end module mod_atmos_boundary
