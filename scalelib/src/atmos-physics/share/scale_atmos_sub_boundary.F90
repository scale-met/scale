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
!! @li      2011-12-11 (H.Yashiro)  [mod] integrate to SCALE-LES ver.3
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2014-05-18 (R.Yoshida)  [add] boudary read/update for real case
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_boundary
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use gtool_file_h, only: &
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
  public :: ATMOS_BOUNDARY_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: ATMOS_BOUNDARY_var  (:,:,:,:)         !> reference container (with HALO)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha(:,:,:,:)         !> damping coefficient [0-1]
  logical,  public              :: ATMOS_BOUNDARY_UPDATE_FLAG = .false. !> switch for real case

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_BOUNDARY_initialize
  private :: ATMOS_BOUNDARY_updatefile
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private :: ATMOS_BOUNDARY_TYPE         = 'NONE'
  character(len=H_LONG), private :: ATMOS_BOUNDARY_IN_BASENAME  = ''
  character(len=H_LONG), private :: ATMOS_BOUNDARY_OUT_BASENAME = ''
  character(len=H_MID),  private :: ATMOS_BOUNDARY_OUT_TITLE    = 'SCALE-LES BOUNDARY CONDITION' !< title of the output file
  character(len=H_MID),  private :: ATMOS_BOUNDARY_OUT_DTYPE    = 'DEFAULT'                      !< REAL4 or REAL8

  logical,               private :: ATMOS_BOUNDARY_USE_VELZ     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_VELX     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_VELY     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_POTT     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_QV       = .false. ! read from file?

  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELZ   =   0.0_RP ! velocity w      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELX   =   0.0_RP ! velocity u      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELY   =   0.0_RP ! velocity v      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_POTT   = 300.0_RP ! potential temp. at boundary, 300 [K]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_QV     =   0.0_RP ! water vapor     at boundary, 0   [kg/kg]

  real(RP),              private :: ATMOS_BOUNDARY_FRACZ        =   1.0_RP ! fraction of boundary region for dumping (z) [0-1]
  real(RP),              private :: ATMOS_BOUNDARY_FRACX        =   1.0_RP ! fraction of boundary region for dumping (x) [0-1]
  real(RP),              private :: ATMOS_BOUNDARY_FRACY        =   1.0_RP ! fraction of boundary region for dumping (y) [0-1]
  real(RP),              private :: ATMOS_BOUNDARY_tauz         =  60.0_RP ! maximum value for damping tau (z) [s]
  real(RP),              private :: ATMOS_BOUNDARY_taux         =  60.0_RP ! maximum value for damping tau (x) [s]
  real(RP),              private :: ATMOS_BOUNDARY_tauy         =  60.0_RP ! maximum value for damping tau (y) [s]

  real(RP),              private, allocatable :: ATMOS_BOUNDARY_var_ref  (:,:,:,:,:) !> reference container (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_increment(:,:,:,:)   !> damping coefficient [0-1]

  real(DP),              private :: ATMOS_BOUNDARY_UPDATE_DT    = 0.0_DP   ! inteval time of boudary data update [s]

  real(DP),              private :: last_updated      =  -999.0_DP
  integer,               private :: boundary_timestep = 0

  character(len=H_SHORT), private :: REF_NAME(5)
  data REF_NAME / 'VELZ_ref','VELX_ref','VELY_ref','POTT_ref','QV_ref' /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_BOUNDARY_setup( &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC  )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    NAMELIST / PARAM_ATMOS_BOUNDARY / &
       ATMOS_BOUNDARY_TYPE,         &
       ATMOS_BOUNDARY_IN_BASENAME,  &
       ATMOS_BOUNDARY_OUT_BASENAME, &
       ATMOS_BOUNDARY_OUT_TITLE,    &
       ATMOS_BOUNDARY_USE_VELZ,     &
       ATMOS_BOUNDARY_USE_VELX,     &
       ATMOS_BOUNDARY_USE_VELY,     &
       ATMOS_BOUNDARY_USE_POTT,     &
       ATMOS_BOUNDARY_USE_QV,       &
       ATMOS_BOUNDARY_VALUE_VELZ,   &
       ATMOS_BOUNDARY_VALUE_VELY,   &
       ATMOS_BOUNDARY_VALUE_VELX,   &
       ATMOS_BOUNDARY_VALUE_POTT,   &
       ATMOS_BOUNDARY_VALUE_QV,     &
       ATMOS_BOUNDARY_FRACZ,        &
       ATMOS_BOUNDARY_FRACX,        &
       ATMOS_BOUNDARY_FRACY,        &
       ATMOS_BOUNDARY_tauz,         &
       ATMOS_BOUNDARY_taux,         &
       ATMOS_BOUNDARY_tauy,         &
       ATMOS_BOUNDARY_UPDATE_DT

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Boundary]/Categ[ATMOS]'

    allocate( ATMOS_BOUNDARY_var  (KA,IA,JA,5) )
    allocate( ATMOS_BOUNDARY_alpha(KA,IA,JA,5) )
    ATMOS_BOUNDARY_var  (:,:,:,:) = CONST_UNDEF
    ATMOS_BOUNDARY_alpha(:,:,:,:) = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_BOUNDARY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_BOUNDARY. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML) write(IO_FID_LOG,nml=PARAM_ATMOS_BOUNDARY)

    if ( ATMOS_BOUNDARY_TYPE == 'CONST' ) then

       call ATMOS_BOUNDARY_generate

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'INIT' ) then

       call ATMOS_BOUNDARY_setinitval( DENS, & ! [IN]
                                       MOMZ, & ! [IN]
                                       MOMX, & ! [IN]
                                       MOMY, & ! [IN]
                                       RHOT, & ! [IN]
                                       QTRC  ) ! [IN]

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'FILE' ) then

       if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
          call ATMOS_BOUNDARY_read
       else
          write(*,*) 'xxx You need specify ATMOS_BOUNDARY_IN_BASENAME'
          call PRC_MPIstop
       endif

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'REAL' ) then

       allocate( ATMOS_BOUNDARY_var_ref  (KA,IA,JA,5,2) )
       allocate( ATMOS_BOUNDARY_increment(KA,IA,JA,5)   )

       if ( ATMOS_BOUNDARY_UPDATE_DT <= 0.0_DP ) then
          write(*,*) 'xxx You need specify ATMOS_BOUNDARY_UPDATE_DT as larger than 0.0'
          call PRC_MPIstop
       endif

       if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
          call ATMOS_BOUNDARY_initialize
       else
          write(*,*) 'xxx You need specify ATMOS_BOUNDARY_IN_BASENAME'
          call PRC_MPIstop
       endif

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .true.

    else
       write(*,*) 'xxx unsupported ATMOS_BOUNDARY_TYPE. Check!', trim(ATMOS_BOUNDARY_TYPE)
       call PRC_MPIstop
    endif

    if( ATMOS_BOUNDARY_OUT_BASENAME /= '' ) then
       call ATMOS_BOUNDARY_write
    endif

    return
  end subroutine ATMOS_BOUNDARY_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_BOUNDARY_var_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, iv
    !---------------------------------------------------------------------------

    do iv = I_BND_VELZ, I_BND_QV
    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_var(   1:KS-1,i,j,iv) = ATMOS_BOUNDARY_var(KS,i,j,iv)
       ATMOS_BOUNDARY_var(KE+1:KA,  i,j,iv) = ATMOS_BOUNDARY_var(KE,i,j,iv)
    enddo
    enddo
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_vars8( ATMOS_BOUNDARY_var(:,:,:,iv), iv )
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_wait ( ATMOS_BOUNDARY_var(:,:,:,iv), iv )
    enddo

    return
  end subroutine ATMOS_BOUNDARY_var_fillhalo

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_BOUNDARY_alpha_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, iv
    !---------------------------------------------------------------------------

    do iv = I_BND_VELZ, I_BND_QV
    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_alpha(   1:KS-1,i,j,iv) = ATMOS_BOUNDARY_alpha(KS,i,j,iv)
       ATMOS_BOUNDARY_alpha(KE+1:KA,  i,j,iv) = ATMOS_BOUNDARY_alpha(KE,i,j,iv)
    enddo
    enddo
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_vars8( ATMOS_BOUNDARY_alpha(:,:,:,iv), iv )
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_wait ( ATMOS_BOUNDARY_alpha(:,:,:,iv), iv )
    enddo

    return
  end subroutine ATMOS_BOUNDARY_alpha_fillhalo

  !-----------------------------------------------------------------------------
  !> Calc dumping coefficient alpha
  subroutine ATMOS_BOUNDARY_setalpha
    use scale_const, only: &
       EPS => CONST_EPS, &
       PI  => CONST_PI
    use scale_grid, only: &
       CBFZ => GRID_CBFZ, &
       CBFX => GRID_CBFX, &
       CBFY => GRID_CBFY, &
       FBFZ => GRID_FBFZ, &
       FBFX => GRID_FBFX, &
       FBFY => GRID_FBFY

    real(RP) :: coef_z, alpha_z1, alpha_z2
    real(RP) :: coef_x, alpha_x1, alpha_x2
    real(RP) :: coef_y, alpha_y1, alpha_y2
    real(RP) :: ee1, ee2

    integer :: i, j, k
    !---------------------------------------------------------------------------

    ! check invalid fraction
    ATMOS_BOUNDARY_FRACZ = max( min( ATMOS_BOUNDARY_FRACZ, 1.0_RP ), EPS )
    ATMOS_BOUNDARY_FRACX = max( min( ATMOS_BOUNDARY_FRACX, 1.0_RP ), EPS )
    ATMOS_BOUNDARY_FRACY = max( min( ATMOS_BOUNDARY_FRACY, 1.0_RP ), EPS )

    if ( ATMOS_BOUNDARY_tauz <= 0.0_RP ) then ! invalid tau
       coef_z = 0.0_RP
    else
       coef_z = 1.0_RP / ATMOS_BOUNDARY_tauz
    endif

    if ( ATMOS_BOUNDARY_taux <= 0.0_RP ) then ! invalid tau
       coef_x = 0.0_RP
    else
       coef_x = 1.0_RP / ATMOS_BOUNDARY_taux
    endif

    if ( ATMOS_BOUNDARY_tauy <= 0.0_RP ) then ! invalid tau
       coef_y = 0.0_RP
    else
       coef_y = 1.0_RP / ATMOS_BOUNDARY_tauy
    endif

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ee1 = CBFZ(k)
       if ( ee1 <= 1.0_RP - ATMOS_BOUNDARY_FRACZ ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + ATMOS_BOUNDARY_FRACZ ) / ATMOS_BOUNDARY_FRACZ
       endif

       ee2 = FBFZ(k)
       if ( ee2 <= 1.0_RP - ATMOS_BOUNDARY_FRACZ ) then
          ee2 = 0.0_RP
       else
          ee2 = ( ee2 - 1.0_RP + ATMOS_BOUNDARY_FRACZ ) / ATMOS_BOUNDARY_FRACZ
       endif

       alpha_z1 = 0.0_RP
       if    ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
          alpha_z1 = coef_z * 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
       elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
          alpha_z1 = coef_z * 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
       endif

       alpha_z2 = 0.0_RP
       if    ( ee2 > 0.0_RP .AND. ee2 <= 0.5_RP ) then
          alpha_z2 = coef_z * 0.5_RP * ( 1.0_RP - cos( ee2*PI ) )
       elseif( ee2 > 0.5_RP .AND. ee2 <= 1.0_RP ) then
          alpha_z2 = coef_z * 0.5_RP * ( 1.0_RP + sin( (ee2-0.5_RP)*PI ) )
       endif

       ee1 = CBFX(i)
       if ( ee1 <= 1.0_RP - ATMOS_BOUNDARY_FRACX ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + ATMOS_BOUNDARY_FRACX ) / ATMOS_BOUNDARY_FRACX
       endif

       ee2 = FBFX(i)
       if ( ee2 <= 1.0_RP - ATMOS_BOUNDARY_FRACX ) then
          ee2 = 0.0_RP
       else
          ee2 = ( ee2 - 1.0_RP + ATMOS_BOUNDARY_FRACX ) / ATMOS_BOUNDARY_FRACX
       endif

       alpha_x1 = 0.0_RP
       if    ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
          alpha_x1 = coef_x * 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
       elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
          alpha_x1 = coef_x * 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
       endif

       alpha_x2 = 0.0_RP
       if    ( ee2 > 0.0_RP .AND. ee2 <= 0.5_RP ) then
          alpha_x2 = coef_x * 0.5_RP * ( 1.0_RP - cos( ee2*PI ) )
       elseif( ee2 > 0.5_RP .AND. ee2 <= 1.0_RP ) then
          alpha_x2 = coef_x * 0.5_RP * ( 1.0_RP + sin( (ee2-0.5_RP)*PI ) )
       endif

       ee1 = CBFY(j)
       if ( ee1 <= 1.0_RP - ATMOS_BOUNDARY_FRACY ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + ATMOS_BOUNDARY_FRACY ) / ATMOS_BOUNDARY_FRACY
       endif

       ee2 = FBFY(j)
       if ( ee2 <= 1.0_RP - ATMOS_BOUNDARY_FRACY ) then
          ee2 = 0.0_RP
       else
          ee2 = ( ee2 - 1.0_RP + ATMOS_BOUNDARY_FRACY ) / ATMOS_BOUNDARY_FRACY
       endif

       alpha_y1 = 0.0_RP
       if    ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
          alpha_y1 = coef_y * 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
       elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
          alpha_y1 = coef_y * 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
       endif

       alpha_y2 = 0.0_RP
       if    ( ee2 > 0.0_RP .AND. ee2 <= 0.5_RP ) then
          alpha_y2 = coef_y * 0.5_RP * ( 1.0_RP - cos( ee2*PI ) )
       elseif( ee2 > 0.5_RP .AND. ee2 <= 1.0_RP ) then
          alpha_y2 = coef_y * 0.5_RP * ( 1.0_RP + sin( (ee2-0.5_RP)*PI ) )
       endif

       ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELZ) = max( alpha_z2, alpha_x1, alpha_y1 )
       ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELX) = max( alpha_z1, alpha_x2, alpha_y1 )
       ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELY) = max( alpha_z1, alpha_x1, alpha_y2 )
       ATMOS_BOUNDARY_alpha(k,i,j,I_BND_POTT) = max( alpha_z1, alpha_x1, alpha_y1 )
       ATMOS_BOUNDARY_alpha(k,i,j,I_BND_QV  ) = max( alpha_z1, alpha_x1, alpha_y1 )
    enddo
    enddo
    enddo

    if ( .NOT. ATMOS_BOUNDARY_USE_VELZ ) then
       ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELZ) = 0.0_RP
    endif
    if ( .NOT. ATMOS_BOUNDARY_USE_VELX ) then
       ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELX) = 0.0_RP
    endif
    if ( .NOT. ATMOS_BOUNDARY_USE_VELY ) then
       ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELY) = 0.0_RP
    endif
    if ( .NOT. ATMOS_BOUNDARY_USE_POTT ) then
       ATMOS_BOUNDARY_alpha(:,:,:,I_BND_POTT) = 0.0_RP
    endif
    if ( .NOT. ATMOS_BOUNDARY_USE_QV   ) then
       ATMOS_BOUNDARY_alpha(:,:,:,I_BND_QV  ) = 0.0_RP
    endif

    call ATMOS_BOUNDARY_alpha_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_setalpha

  !-----------------------------------------------------------------------------
  !> Read boundary data
  subroutine ATMOS_BOUNDARY_setinitval( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    integer :: i, j, k, iv
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELZ) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j) + DENS(k,i,j) )
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i+1,j) + DENS(k,i,j) )
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELY) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1) + DENS(k,i,j) )
       ATMOS_BOUNDARY_var(k,i,j,I_BND_POTT) = RHOT(k,i,j) / DENS(k,i,j)
       ATMOS_BOUNDARY_var(k,i,j,I_BND_QV  ) = QTRC(k,i,j,I_QV)
    enddo
    enddo
    enddo

    call ATMOS_BOUNDARY_var_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_setinitval

  !-----------------------------------------------------------------------------
  !> Read boundary data
  subroutine ATMOS_BOUNDARY_read
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    implicit none

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=H_LONG) :: bname

    integer :: iv, i, j
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_IN_BASENAME

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       call FileRead( reference_atmos(:,:,:), bname, 'VELZ', 1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_VELZ) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_VELZ', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha(KS:KE,IS:IE,JS:JE,I_BND_VELZ) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_VELX ) then
       call FileRead( reference_atmos(:,:,:), bname, 'VELX', 1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_VELX) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_VELX', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha(KS:KE,IS:IE,JS:JE,I_BND_VELX) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_VELY ) then
       call FileRead( reference_atmos(:,:,:), bname, 'VELY', 1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_VELY) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_VELY', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha(KS:KE,IS:IE,JS:JE,I_BND_VELY) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_POTT ) then
       call FileRead( reference_atmos(:,:,:), bname, 'POTT', 1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_POTT) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_POTT', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha(KS:KE,IS:IE,JS:JE,I_BND_POTT) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_QV   ) then
       call FileRead( reference_atmos(:,:,:), bname, 'QV',   1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_QV) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_QV', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha(KS:KE,IS:IE,JS:JE,I_BND_QV) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    call ATMOS_BOUNDARY_var_fillhalo
    call ATMOS_BOUNDARY_alpha_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_read

  !-----------------------------------------------------------------------------
  !> Write boundary data
  subroutine ATMOS_BOUNDARY_write
    use scale_fileio, only: &
       FILEIO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_VELZ),                  &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELZ', 'Reference Velocity w', 'm/s', 'ZXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELZ),                &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELZ', 'Alpha for w', '1', 'ZXY',               &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELX ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_VELX),                  &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELX', 'Reference Velocity u', 'm/s', 'ZXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELX),                &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELX', 'Alpha for u', '1', 'ZXY',               &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELY ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_VELY),                  &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELY', 'Reference Velocity y', 'm/s', 'ZXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELY),                &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELY', 'Alpha for v', '1', 'ZXY',               &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_POTT ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_POTT),                  &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'POTT', 'Reference PT', 'K', 'ZXY',                    &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_POTT),                &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_POTT', 'Alpha for PT', '1', 'ZXY',              &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_QV   ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_QV),                    &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'QV', 'Reference water vapor', 'kg/kg', 'ZXY',         &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_QV),                  &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_QV', 'Alpha for QV', '1', 'ZXY',                &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    return
  end subroutine ATMOS_BOUNDARY_write

  !-----------------------------------------------------------------------------
  !> generate boundary data
  subroutine ATMOS_BOUNDARY_generate
    implicit none

    integer :: i, j, k, iv
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELZ) = ATMOS_BOUNDARY_VALUE_VELZ
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELY) = ATMOS_BOUNDARY_VALUE_VELY
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = ATMOS_BOUNDARY_VALUE_VELX
       ATMOS_BOUNDARY_var(k,i,j,I_BND_POTT) = ATMOS_BOUNDARY_VALUE_POTT
       ATMOS_BOUNDARY_var(k,i,j,I_BND_QV  ) = ATMOS_BOUNDARY_VALUE_QV
    enddo
    enddo
    enddo

    call ATMOS_BOUNDARY_var_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_generate

  !-----------------------------------------------------------------------------
  !> Initialize boundary value for real case experiment
  subroutine ATMOS_BOUNDARY_initialize
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_time, only: &
       TIME_DTSEC,  &
       TIME_NOWSEC
    implicit none

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=H_LONG) :: bname

    integer  :: i, j, k, iv
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_IN_BASENAME

    ! read boundary data from input file
    boundary_timestep = 1
    call FileRead( reference_atmos(:,:,:), bname, 'VELZ', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELZ,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELX', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELX,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELY', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELY,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'POTT', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_POTT,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'QV',   boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_QV  ,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)

    boundary_timestep = 2
    call FileRead( reference_atmos(:,:,:), bname, 'VELZ', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELZ,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELX', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELX,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELY', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELY,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'POTT', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_POTT,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'QV',   boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_QV  ,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)

    do iv = I_BND_VELZ, I_BND_QV
    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_var_ref(   1:KS-1,i,j,iv,1) = ATMOS_BOUNDARY_var_ref(KS,i,j,iv,1)
       ATMOS_BOUNDARY_var_ref(KE+1:KA,  i,j,iv,1) = ATMOS_BOUNDARY_var_ref(KE,i,j,iv,1)
       ATMOS_BOUNDARY_var_ref(   1:KS-1,i,j,iv,2) = ATMOS_BOUNDARY_var_ref(KS,i,j,iv,2)
       ATMOS_BOUNDARY_var_ref(KE+1:KA,  i,j,iv,2) = ATMOS_BOUNDARY_var_ref(KE,i,j,iv,2)
    enddo
    enddo
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_vars8( ATMOS_BOUNDARY_var_ref(:,:,:,iv,1), iv )
       call COMM_vars8( ATMOS_BOUNDARY_var_ref(:,:,:,iv,2), iv )
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_wait( ATMOS_BOUNDARY_var_ref(:,:,:,iv,1), iv )
       call COMM_wait( ATMOS_BOUNDARY_var_ref(:,:,:,iv,2), iv )
    enddo

    ! set boundary data and time increment
    do iv = I_BND_VELZ, I_BND_QV
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       ATMOS_BOUNDARY_var      (k,i,j,iv) = ATMOS_BOUNDARY_var_ref(k,i,j,iv,1)
       ATMOS_BOUNDARY_increment(k,i,j,iv) = ( ATMOS_BOUNDARY_var_ref(k,i,j,iv,2) &
                                            - ATMOS_BOUNDARY_var_ref(k,i,j,iv,1) ) &
                                            / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
    enddo
    enddo
    enddo
    enddo

    last_updated = TIME_NOWSEC

    return
  end subroutine ATMOS_BOUNDARY_initialize

  !-----------------------------------------------------------------------------
  !> Update boundary value with a constant time increment
  subroutine ATMOS_BOUNDARY_update()
    use scale_time, only: &
       TIME_DTSEC, &
       TIME_NOWSEC
    implicit none

    real(DP) :: integrated_sec
    integer  :: i, j, k, iv
    !---------------------------------------------------------------------------

    integrated_sec = TIME_NOWSEC - last_updated

    if ( integrated_sec >= ATMOS_BOUNDARY_UPDATE_DT ) then
       boundary_timestep = boundary_timestep + 1

       call ATMOS_BOUNDARY_updatefile

       do iv = I_BND_VELZ, I_BND_QV
       do j  = 1, JA
       do i  = 1, IA
       do k  = 1, KA
          ATMOS_BOUNDARY_var      (k,i,j,iv) = ATMOS_BOUNDARY_var_ref(k,i,j,iv,1)
          ATMOS_BOUNDARY_increment(k,i,j,iv) = ( ATMOS_BOUNDARY_var_ref(k,i,j,iv,2) &
                                               - ATMOS_BOUNDARY_var_ref(k,i,j,iv,1) ) &
                                               / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
       enddo
       enddo
       enddo
       enddo

       last_updated = TIME_NOWSEC
    else
       do iv = I_BND_VELZ, I_BND_QV
       do j  = 1, JA
       do i  = 1, IA
       do k  = 1, KA
          ATMOS_BOUNDARY_var(k,i,j,iv) = ATMOS_BOUNDARY_var(k,i,j,iv) + ATMOS_BOUNDARY_increment(k,i,j,iv)
       enddo
       enddo
       enddo
       enddo
    endif

    return
  end subroutine ATMOS_BOUNDARY_update

  !-----------------------------------------------------------------------------
  !> Update reference boundary file
  subroutine ATMOS_BOUNDARY_updatefile()
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_time, only: &
       TIME_NOWSEC
    implicit none

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=H_LONG) :: bname

    integer :: i, j, k, iv
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_IN_BASENAME

    do iv = I_BND_VELZ, I_BND_QV
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       ATMOS_BOUNDARY_var_ref(k,i,j,iv,1) = ATMOS_BOUNDARY_var_ref(k,i,j,iv,2)
    enddo
    enddo
    enddo
    enddo

    call FileRead( reference_atmos(:,:,:), bname, 'VELZ', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELZ,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELX', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELX,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELY', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELY,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'POTT', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_POTT,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'QV',   boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_QV  ,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)

    do iv = I_BND_VELZ, I_BND_QV
    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_var_ref(   1:KS-1,i,j,iv,2) = ATMOS_BOUNDARY_var_ref(KS,i,j,iv,2)
       ATMOS_BOUNDARY_var_ref(KE+1:KA,  i,j,iv,2) = ATMOS_BOUNDARY_var_ref(KE,i,j,iv,2)
    enddo
    enddo
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_vars8( ATMOS_BOUNDARY_var_ref(:,:,:,iv,2), iv )
    enddo

    do iv = I_BND_VELZ, I_BND_QV
       call COMM_wait ( ATMOS_BOUNDARY_var_ref(:,:,:,iv,2), iv )
    enddo

    return
  end subroutine ATMOS_BOUNDARY_updatefile

end module scale_atmos_boundary
