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
  include 'inc_index.h'
  include 'inc_precision.h'

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
  character(len=IO_FILECHR), private :: ATMOS_BOUNDARY_IN_BASENAME   = ''
  character(len=IO_FILECHR), private :: ATMOS_BOUNDARY_OUT_BASENAME  = ''
  character(len=FILE_HLONG), private :: ATMOS_BOUNDARY_OUT_TITLE     = 'SCALE3 BOUNDARY CONDITION'
  character(len=FILE_HLONG), private :: ATMOS_BOUNDARY_OUT_SOURCE    = 'SCALE-LES ver. 3'
  character(len=FILE_HLONG), private :: ATMOS_BOUNDARY_OUT_INSTITUTE = 'AICS/RIKEN'
  logical,                   private :: ATMOS_BOUNDARY_USE_VELZ      = .false. ! read from file?
  logical,                   private :: ATMOS_BOUNDARY_USE_VELX      = .false. ! read from file?
  logical,                   private :: ATMOS_BOUNDARY_USE_VELY      = .false. ! read from file?
  logical,                   private :: ATMOS_BOUNDARY_USE_POTT      = .false. ! read from file?
  logical,                   private :: ATMOS_BOUNDARY_USE_QV        = .false. ! read from file?
  real(RP),                  private :: ATMOS_BOUNDARY_VALUE_VELZ    =  0.D0 ! w at boundary, 0 [m/s]
  real(RP),                  private :: ATMOS_BOUNDARY_VALUE_VELX    =  5.D0 ! u at boundary, 5 [m/s]
  real(RP),                  private :: ATMOS_BOUNDARY_tauz          = 75.D0 ! maximum value for damping tau (z) [s]
  real(RP),                  private :: ATMOS_BOUNDARY_taux          = 75.D0 ! maximum value for damping tau (x) [s]
  real(RP),                  private :: ATMOS_BOUNDARY_tauy          = 75.D0 ! maximum value for damping tau (y) [s]

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
       CONST_UNDEF8
    implicit none

    NAMELIST / PARAM_ATMOS_BOUNDARY / &
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
       ATMOS_BOUNDARY_VALUE_VELX,    &
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
    ATMOS_BOUNDARY_var(:,:,:,:) = CONST_UNDEF8

    if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
       call ATMOS_BOUNDARY_read
    elseif( ATMOS_BOUNDARY_OUT_BASENAME /= '' ) then
       ATMOS_BOUNDARY_var(:,:,:,:) = 0.D0
    endif

    if ( ATMOS_BOUNDARY_OUT_BASENAME /= '' ) then
       call ATMOS_BOUNDARY_generate
       call ATMOS_BOUNDARY_write
    endif

    call ATMOS_BOUNDARY_setalpha

    return
  end subroutine ATMOS_BOUNDARY_setup

  !-----------------------------------------------------------------------------
  !> Calc dumping coefficient alpha
  !-----------------------------------------------------------------------------
  subroutine ATMOS_BOUNDARY_setalpha
    use mod_const, only: &
       CONST_UNDEF8, &
       PI => CONST_PI
    use mod_grid, only : &
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
    ATMOS_BOUNDARY_alpha(:,:,:,:) = 0.D0

    coef = 1.D0 / ATMOS_BOUNDARY_tauz

    do k = KS, KE
       ee1 = CBFZ(k)

       if    ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELX) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELY) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_QV  ) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELX) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELY) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_QV  ) )
       endif
    enddo

    do k = KS-1, KE
       ee2 = FBFZ(k)

       if    ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELZ) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(k,:,:,I_BND_VELZ) )
       endif
    enddo

    coef = 1.D0 / ATMOS_BOUNDARY_taux

    do i = IS, IE
       ee1 = CBFX(i)
       ee2 = FBFX(i)

       if ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELZ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELY) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_QV  ) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELZ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELY) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_QV  ) )
       endif

       if ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELX) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(:,i,:,I_BND_VELX) )
       endif
    enddo

    coef = 1.D0 / ATMOS_BOUNDARY_tauy

    do j = JS, JE
       ee1 = CBFY(j)
       ee2 = FBFY(j)

       if ( ee1 > 0.0D0 .AND. ee1 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee1*PI ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELZ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELX) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_QV  ) )
       elseif( ee1 > 0.5D0 .AND. ee1 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee1-0.5D0)*PI ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELZ) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELZ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELX) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELX) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_POTT) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_POTT) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_QV  ) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_QV  ) )
       endif

       if ( ee2 > 0.0D0 .AND. ee2 <= 0.5D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 - dcos( ee2*PI ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELY) )
       elseif( ee2 > 0.5D0 .AND. ee2 <= 1.0D0 ) then
          alpha = coef * 0.5D0 * ( 1.D0 + dsin( (ee2-0.5D0)*PI ) )
          ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELY) = max( alpha, ATMOS_BOUNDARY_alpha(:,:,j,I_BND_VELY) )
       endif
    enddo

    do j = JS-1, JE+1
    do i = IS-1, IE+1
       do k = KS, KE
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) == CONST_UNDEF8 ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELX) = 0.D0
          endif
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_VELY) == CONST_UNDEF8 ) then
            ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELY) = 0.D0
          endif
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_POTT) == CONST_UNDEF8 ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_POTT) = 0.D0
          endif
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_QV) == CONST_UNDEF8 ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_QV) = 0.D0
          endif
       enddo

       do k = KS-1, KE
          if ( ATMOS_BOUNDARY_var(k,i,j,I_BND_VELZ) == CONST_UNDEF8 ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELZ) = 0.D0
         endif
      enddo
    enddo
    enddo

    return
  end subroutine ATMOS_BOUNDARY_setalpha

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
       PRC_myrank
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use gtool_file, only: &
       FileCreate, &
       FileAddVariable, &
       FilePutAxis, &
       FileWrite, &
       FileClose
    use mod_grid, only :  &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY
    implicit none

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    integer :: fid, vid
    integer :: dtype
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_OUT_BASENAME

    call FileCreate( fid,                                       & ! (out)
         bname,                                                 & ! (in)
         ATMOS_BOUNDARY_OUT_TITLE,                              & ! (in)
         ATMOS_BOUNDARY_OUT_SOURCE,                             & ! (in)
         ATMOS_BOUNDARY_OUT_INSTITUTE,                          & ! (in)
         (/'z','x','y'/), (/KMAX,IMAX,JMAX/), (/'Z','X','Y'/),  & ! (in)
         (/'m','m','m'/), (/File_REAL4,File_REAL4,File_REAL4/), & ! (in)
         PRC_master, PRC_myrank                                 ) ! (in)

    call FilePutAxis(fid, 'z', GRID_CZ(KS:KE))
    call FilePutAxis(fid, 'x', GRID_CX(IS:IE))
    call FilePutAxis(fid, 'y', GRID_CY(JS:JE))

    if ( RP == 8 ) then
       dtype = File_REAL8
    else if ( RP == 4 ) then
       dtype = File_REAL4
    end if

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
       CONST_UNDEF8
    use mod_grid, only : &
       CZ_mask => GRID_CZ_mask, &
       CX_mask => GRID_CX_mask
    use mod_comm, only: &
       COMM_vars, &
       COMM_wait
    use mod_atmos_refstate, only: &
       ATMOS_REFSTATE_pott
    implicit none

    integer :: i, j, k, iv
    !---------------------------------------------------------------------------

    do k = KS, KE
       if ( CZ_mask(k) ) then ! Inner Layer
          ATMOS_BOUNDARY_var(k,:,:,I_BND_VELZ) = CONST_UNDEF8
       else                   ! Buffer Layer
          ATMOS_BOUNDARY_var(k,:,:,I_BND_VELZ) = ATMOS_BOUNDARY_VALUE_VELZ
       endif
    enddo
    ATMOS_BOUNDARY_var(:,:,:,I_BND_VELY) = CONST_UNDEF8
    ATMOS_BOUNDARY_var(:,:,:,I_BND_POTT) = CONST_UNDEF8
    ATMOS_BOUNDARY_var(:,:,:,I_BND_QV)   = CONST_UNDEF8

!    do j = JS-1, JE+1
!    do i = IS-1, IE+1
!       do k = KS, KE
!          if ( CZ_mask(k) .AND. CX_mask(i) ) then ! Inner Area
!             ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = CONST_UNDEF8
!          else                                    ! Buffer Area
!             ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = FBFX(i) * ATMOS_BOUNDARY_VALUE_VELX &
!                                                  * ( 1.D0 - CBFZ(k) )
!          endif
!       enddo
!    enddo
!    enddo
    ATMOS_BOUNDARY_var(:,:,:,I_BND_VELX) = CONST_UNDEF8

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
