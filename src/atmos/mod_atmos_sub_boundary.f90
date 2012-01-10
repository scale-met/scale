
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
!! @li      2011-12-07 (Y.Miyamoto) [new]
!! @li      2011-12-11 (H.Yashiro)  [mod] integrate to SCALE3
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
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(8), public, allocatable, save :: atmos_refvar(:,:,:,:)  !> reference container (with HALO)

  integer, public, parameter :: I_REF_VELX = 1 ! reference velocity (x) [m/s]
  integer, public, parameter :: I_REF_VELY = 2 ! reference velocity (y) [m/s]
  integer, public, parameter :: I_REF_VELZ = 3 ! reference velocity (z) [m/s]
  integer, public, parameter :: I_REF_POTT = 4 ! reference potential temperature [K]
  integer, public, parameter :: I_REF_QV   = 5 ! reference water vapor [kg/kg]

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
  character(len=IO_FILECHR), private, save :: ATMOS_BOUNDARY_IN_BASENAME = 'refvar_in'
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
       PRC_MPIstop
    use mod_const, only: &
       CONST_UNDEF8, &
       PI => CONST_PI
    use mod_grid, only : &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KA => GRID_KA, &
       IS => GRID_IS, &
       IE => GRID_IE, &
       JS => GRID_JS, &
       JE => GRID_JE, &
       KS => GRID_KS, &
       KE => GRID_KE, &
       WS => GRID_WS, &
       WE => GRID_WE, &
       GRID_CBFX, &
       GRID_CBFY, &
       GRID_CBFZ, &
       GRID_FBFX, &
       GRID_FBFY, &
       GRID_FBFZ
    implicit none

    real(8) :: ATMOS_BOUNDARY_taux = 75.D0 ! maximum value for damping tau (x) [s]
    real(8) :: ATMOS_BOUNDARY_tauy = 75.D0 ! maximum value for damping tau (y) [s]
    real(8) :: ATMOS_BOUNDARY_tauz = 75.D0 ! maximum value for damping tau (z) [s]

    NAMELIST / PARAM_ATMOS_BOUNDARY / &
       ATMOS_BOUNDARY_taux,  &
       ATMOS_BOUNDARY_tauy,  &
       ATMOS_BOUNDARY_tauz, &
       ATMOS_BOUNDARY_IN_BASENAME

    real(8) :: coef, alpha
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
    allocate( atmos_refvar(KA,IA,JA,5) ); atmos_refvar(:,:,:,:) = CONST_UNDEF8

    call ATMOS_BOUNDARY_reference_read

    !--- set damping coefficient
    allocate( DAMP_alphau(KA,IA,JA) ); DAMP_alphau(:,:,:) = 0.D0
    allocate( DAMP_alphav(KA,IA,JA) ); DAMP_alphav(:,:,:) = 0.D0
    allocate( DAMP_alphaw(KA,IA,JA) ); DAMP_alphaw(:,:,:) = 0.D0
    allocate( DAMP_alphat(KA,IA,JA) ); DAMP_alphat(:,:,:) = 0.D0
    allocate( DAMP_alphaq(KA,IA,JA) ); DAMP_alphaq(:,:,:) = 0.D0

    coef = 1.D0 / ATMOS_BOUNDARY_taux

    do i = IS, IE
       ee1 = GRID_CBFX(i)
       ee2 = GRID_FBFX(i)

       if ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          DAMP_alphav(KS:KE,i,JS:JE) = max( alpha, DAMP_alphav(KS:KE,i,JS:JE) )
          DAMP_alphat(KS:KE,i,JS:JE) = max( alpha, DAMP_alphat(KS:KE,i,JS:JE) )
          DAMP_alphaq(KS:KE,i,JS:JE) = max( alpha, DAMP_alphaq(KS:KE,i,JS:JE) )
          DAMP_alphaw(WS:WE,i,JS:JE) = max( alpha, DAMP_alphaw(WS:WE,i,JS:JE) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          DAMP_alphav(KS:KE,i,JS:JE) = max( alpha, DAMP_alphav(KS:KE,i,JS:JE) )
          DAMP_alphat(KS:KE,i,JS:JE) = max( alpha, DAMP_alphat(KS:KE,i,JS:JE) )
          DAMP_alphaq(KS:KE,i,JS:JE) = max( alpha, DAMP_alphaq(KS:KE,i,JS:JE) )
          DAMP_alphaw(WS:WE,i,JS:JE) = max( alpha, DAMP_alphaw(WS:WE,i,JS:JE) )
       endif

       if ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          DAMP_alphau(KS:KE,i,JS:JE) = max( alpha, DAMP_alphau(KS:KE,i,JS:JE) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          DAMP_alphau(KS:KE,i,JS:JE) = max( alpha, DAMP_alphau(KS:KE,i,JS:JE) )
       endif
    enddo

    coef = 1.D0 / ATMOS_BOUNDARY_tauy

    do j = JS, JE
       ee1 = GRID_CBFY(j)
       ee2 = GRID_FBFY(j)

       if ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          DAMP_alphau(KS:KE,IS:IE,j) = max( alpha, DAMP_alphau(KS:KE,IS:IE,j) )
          DAMP_alphaw(WS:WE,IS:IE,j) = max( alpha, DAMP_alphaw(WS:WE,IS:IE,j) )
          DAMP_alphat(KS:KE,IS:IE,j) = max( alpha, DAMP_alphat(KS:KE,IS:IE,j) )
          DAMP_alphaq(KS:KE,IS:IE,j) = max( alpha, DAMP_alphaq(KS:KE,IS:IE,j) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          DAMP_alphau(KS:KE,IS:IE,j) = max( alpha, DAMP_alphau(KS:KE,IS:IE,j) )
          DAMP_alphaw(WS:WE,IS:IE,j) = max( alpha, DAMP_alphaw(WS:WE,IS:IE,j) )
          DAMP_alphat(KS:KE,IS:IE,j) = max( alpha, DAMP_alphat(KS:KE,IS:IE,j) )
          DAMP_alphaq(KS:KE,IS:IE,j) = max( alpha, DAMP_alphaq(KS:KE,IS:IE,j) )
       endif

       if ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          DAMP_alphav(KS:KE,IS:IE,j) = max( alpha, DAMP_alphav(KS:KE,IS:IE,j) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          DAMP_alphav(KS:KE,IS:IE,j) = max( alpha, DAMP_alphav(KS:KE,IS:IE,j) )
       endif
    enddo

    coef = 1.D0 / ATMOS_BOUNDARY_tauz

    do k = KS, KE
       ee1 = GRID_CBFZ(k)

       if    ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          DAMP_alphau(k,IS:IE,JS:JE) = max( alpha, DAMP_alphau(k,IS:IE,JS:JE) )
          DAMP_alphav(k,IS:IE,JS:JE) = max( alpha, DAMP_alphav(k,IS:IE,JS:JE) )
          DAMP_alphat(k,IS:IE,JS:JE) = max( alpha, DAMP_alphat(k,IS:IE,JS:JE) )
          DAMP_alphaq(k,IS:IE,JS:JE) = max( alpha, DAMP_alphaq(k,IS:IE,JS:JE) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          DAMP_alphau(k,IS:IE,JS:JE) = max( alpha, DAMP_alphau(k,IS:IE,JS:JE) )
          DAMP_alphav(k,IS:IE,JS:JE) = max( alpha, DAMP_alphav(k,IS:IE,JS:JE) )
          DAMP_alphat(k,IS:IE,JS:JE) = max( alpha, DAMP_alphat(k,IS:IE,JS:JE) )
          DAMP_alphaq(k,IS:IE,JS:JE) = max( alpha, DAMP_alphaq(k,IS:IE,JS:JE) )
       endif
    enddo

    do k = WS, WE
       ee2 = GRID_FBFZ(k)

       if    ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          DAMP_alphaw(k,IS:IE,JS:JE) = max( alpha, DAMP_alphaw(k,IS:IE,JS:JE) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          DAMP_alphaw(k,IS:IE,JS:JE) = max( alpha, DAMP_alphaw(k,IS:IE,JS:JE) )
       endif
    enddo

    do j = JS-1, JE+1
    do i = IS-1, IE+1
       do k = KS, KE
          if ( atmos_refvar(k,i,j,I_REF_VELX) == CONST_UNDEF8 ) then
             DAMP_alphau(k,i,j) = 0.D0
          endif
          if ( atmos_refvar(k,i,j,I_REF_VELY) == CONST_UNDEF8 ) then
            DAMP_alphav(k,i,j) = 0.D0
          endif
          if ( atmos_refvar(k,i,j,I_REF_POTT) == CONST_UNDEF8 ) then
             DAMP_alphat(k,i,j) = 0.D0
          endif
          if ( atmos_refvar(k,i,j,I_REF_QV) == CONST_UNDEF8 ) then
             DAMP_alphaq(k,i,j) = 0.D0
          endif
       enddo

       do k = WS, WE
          if ( atmos_refvar(k,i,j,I_REF_VELZ) == CONST_UNDEF8 ) then
             DAMP_alphaw(k,i,j) = 0.D0
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
  !> Read restart of atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_reference_read
    use mod_const, only: &
       CONST_UNDEF8
    use mod_comm, only: &
       COMM_vars, &
       COMM_wait, &
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

    real(8) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=FIO_HSHORT) :: REF_NAME(5)
    data REF_NAME / 'VELX_ref','VELY_ref','VELZ_ref','POTT_ref','QV_ref' /

    character(len=IO_FILECHR) :: bname
    character(len=8)          :: lname

    integer :: iv, i, j
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_IN_BASENAME
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    if ( ref_velx ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'VELX', lname, 1, KMAX, 1 )
       atmos_refvar(KS:KE,IS:IE,JS:JE,I_REF_VELX) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ref_vely ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'VELY', lname, 1, KMAX, 1 )
       atmos_refvar(KS:KE,IS:IE,JS:JE,I_REF_VELY) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ref_velz ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'VELZ', lname, 1, KMAX, 1 )
       atmos_refvar(KS:KE,IS:IE,JS:JE,I_REF_VELZ) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ref_pott ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'POTT', lname, 1, KMAX, 1 )
       atmos_refvar(KS:KE,IS:IE,JS:JE,I_REF_POTT) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ref_qv ) then
       call FIO_input( reference_atmos(:,:,:), bname, 'QV',   lname, 1, KMAX, 1 )
       atmos_refvar(KS:KE,IS:IE,JS:JE,I_REF_QV) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    ! fill IHALO & JHALO
    do iv = I_REF_VELX, I_REF_QV
       call COMM_vars( atmos_refvar(:,:,:,iv), iv )
    enddo

    do iv = I_REF_VELX, I_REF_QV
       call COMM_wait( iv )
    enddo

    ! fill KHALO
    do iv = I_REF_VELX, I_REF_QV
    do j  = 1, JA
    do i  = 1, IA
       atmos_refvar(   1:KS-1,i,j,iv) = atmos_refvar(KS,i,j,iv)
       atmos_refvar(KE+1:KA,  i,j,iv) = atmos_refvar(KE,i,j,iv)
    enddo
    enddo
    enddo

    call COMM_stats( atmos_refvar(:,:,:,:), REF_NAME(:) )

    return
  end subroutine ATMOS_BOUNDARY_reference_read

end module mod_atmos_boundary
