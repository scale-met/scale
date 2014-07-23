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
  use scale_index
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
  !++ Public parameters & variables
  !
  real(RP), public :: ATMOS_BOUNDARY_SMOOTHER_FACT  =  0.2_RP ! fact for smoother to damping
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
  logical,               private :: ATMOS_BOUNDARY_USE_DENS     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_QV       = .false. ! read from file?

  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELZ   =   0.0_RP ! velocity w      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELX   =   0.0_RP ! velocity u      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELY   =   0.0_RP ! velocity v      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_POTT   = 300.0_RP ! potential temp. at boundary, 300 [K]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_QV     =   0.0_RP ! water vapor     at boundary, 0   [kg/kg]

  real(RP),              private :: ATMOS_BOUNDARY_FRACZ        =   1.0_RP ! fraction of boundary region for dumping (z) [0-1]
  real(RP),              private :: ATMOS_BOUNDARY_FRACX        =   1.0_RP ! fraction of boundary region for dumping (x) [0-1]
  real(RP),              private :: ATMOS_BOUNDARY_FRACY        =   1.0_RP ! fraction of boundary region for dumping (y) [0-1]
  real(RP),              private :: ATMOS_BOUNDARY_tauz                    ! maximum value for damping tau (z) [s]
  real(RP),              private :: ATMOS_BOUNDARY_taux                    ! maximum value for damping tau (x) [s]
  real(RP),              private :: ATMOS_BOUNDARY_tauy                    ! maximum value for damping tau (y) [s]

  real(RP),              private, allocatable :: ATMOS_BOUNDARY_var_ref  (:,:,:,:,:) !> reference container (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_increment(:,:,:,:)   !> damping coefficient [0-1]

  real(DP),              private :: ATMOS_BOUNDARY_UPDATE_DT    = 0.0_DP   ! inteval time of boudary data update [s]
  integer,               private :: ATMOS_BOUNDARY_START_DATE(6) = (/ -9999, 0, 0, 0, 0, 0 /) ! boundary initial date

  real(DP),              private :: last_updated      =  -999.0_DP
  real(DP),              private :: integrated_sec    = 0.0_DP
  integer,               private :: boundary_timestep = 0
  logical,               private :: ATMOS_BOUNDARY_LINEARZ      = .false.  ! linear or non-linear profile of relax region

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
    use scale_time, only: &
       DT => TIME_DTSEC
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
       ATMOS_BOUNDARY_USE_DENS,     &
       ATMOS_BOUNDARY_USE_QV,       &
       ATMOS_BOUNDARY_VALUE_VELZ,   &
       ATMOS_BOUNDARY_VALUE_VELY,   &
       ATMOS_BOUNDARY_VALUE_VELX,   &
       ATMOS_BOUNDARY_VALUE_POTT,   &
       ATMOS_BOUNDARY_VALUE_QV,     &
       ATMOS_BOUNDARY_SMOOTHER_FACT,&
       ATMOS_BOUNDARY_FRACZ,        &
       ATMOS_BOUNDARY_FRACX,        &
       ATMOS_BOUNDARY_FRACY,        &
       ATMOS_BOUNDARY_tauz,         &
       ATMOS_BOUNDARY_taux,         &
       ATMOS_BOUNDARY_tauy,         &
       ATMOS_BOUNDARY_UPDATE_DT,    &
       ATMOS_BOUNDARY_START_DATE, &
       ATMOS_BOUNDARY_LINEARZ

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Boundary]/Categ[ATMOS]'

    allocate( ATMOS_BOUNDARY_var  (KA,IA,JA,I_BND_SIZE) )
    allocate( ATMOS_BOUNDARY_alpha(KA,IA,JA,I_BND_SIZE) )
    ATMOS_BOUNDARY_var  (:,:,:,:) = CONST_UNDEF
    ATMOS_BOUNDARY_alpha(:,:,:,:) = 0.0_RP

    ATMOS_BOUNDARY_tauz = DT * 10.0_RP
    ATMOS_BOUNDARY_taux = DT * 10.0_RP
    ATMOS_BOUNDARY_tauy = DT * 10.0_RP

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

       allocate( ATMOS_BOUNDARY_var_ref  (KA,IA,JA,I_BND_SIZE,2) )
       allocate( ATMOS_BOUNDARY_increment(KA,IA,JA,I_BND_SIZE)   )

       if ( ATMOS_BOUNDARY_UPDATE_DT <= 0.0_DP ) then
          write(*,*) 'xxx You need specify ATMOS_BOUNDARY_UPDATE_DT as larger than 0.0'
          call PRC_MPIstop
       endif

       if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
          call ATMOS_BOUNDARY_initialize( DENS )
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

    do iv = 1, I_BND_SIZE
    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_var(   1:KS-1,i,j,iv) = ATMOS_BOUNDARY_var(KS,i,j,iv)
       ATMOS_BOUNDARY_var(KE+1:KA,  i,j,iv) = ATMOS_BOUNDARY_var(KE,i,j,iv)
    enddo
    enddo
    enddo

    do iv = 1, I_BND_SIZE
       call COMM_vars8( ATMOS_BOUNDARY_var(:,:,:,iv), iv )
    enddo

    do iv = 1, I_BND_SIZE
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

    do iv = 1, I_BND_SIZE
    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_alpha(   1:KS-1,i,j,iv) = ATMOS_BOUNDARY_alpha(KS,i,j,iv)
       ATMOS_BOUNDARY_alpha(KE+1:KA,  i,j,iv) = ATMOS_BOUNDARY_alpha(KE,i,j,iv)
    enddo
    enddo
    enddo

    do iv = 1, I_BND_SIZE
       call COMM_vars8( ATMOS_BOUNDARY_alpha(:,:,:,iv), iv )
    enddo

    do iv = 1, I_BND_SIZE
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
       if ( ATMOS_BOUNDARY_LINEARZ ) then
          alpha_z1 = coef_z * ee1
       else
          if    ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
             alpha_z1 = coef_z * 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
          elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
             alpha_z1 = coef_z * 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
          endif
       endif

       alpha_z2 = 0.0_RP
       if ( ATMOS_BOUNDARY_LINEARZ ) then
          alpha_z2 = coef_z * ee2
       else
          if    ( ee2 > 0.0_RP .AND. ee2 <= 0.5_RP ) then
             alpha_z2 = coef_z * 0.5_RP * ( 1.0_RP - cos( ee2*PI ) )
          elseif( ee2 > 0.5_RP .AND. ee2 <= 1.0_RP ) then
             alpha_z2 = coef_z * 0.5_RP * ( 1.0_RP + sin( (ee2-0.5_RP)*PI ) )
          endif
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

       alpha_x1 = coef_x * ee1
       alpha_x2 = coef_x * ee2

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

       alpha_y1 = coef_y * ee1
       alpha_y2 = coef_y * ee2


       if ( ATMOS_BOUNDARY_TYPE == 'REAL' ) then
          ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELZ) = alpha_z2
          if ( ATMOS_BOUNDARY_USE_DENS ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_DENS) = max( alpha_z1, alpha_x1, alpha_y1 )
          else
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_DENS) = max( alpha_x1, alpha_y1 )
          endif
          if ( ATMOS_BOUNDARY_USE_VELX ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELX) = max( alpha_z1, alpha_x2, alpha_y1 )
          else
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELX) = max( alpha_x2, alpha_y1 )
          endif
          if ( ATMOS_BOUNDARY_USE_VELY ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELY) = max( alpha_z1, alpha_x1, alpha_y2 )
          else
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELY) = max( alpha_x1, alpha_y2 )
          endif
          if ( ATMOS_BOUNDARY_USE_POTT ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_POTT) = max( alpha_z1, alpha_x1, alpha_y1 )
          else
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_POTT) = max( alpha_x1, alpha_y1 )
          endif
          if ( ATMOS_BOUNDARY_USE_QV   ) then
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_QV  ) = max( alpha_z1, alpha_x1, alpha_y1 )
          else
             ATMOS_BOUNDARY_alpha(k,i,j,I_BND_QV  ) = max( alpha_x1, alpha_y1 )
          endif
       else
          ATMOS_BOUNDARY_alpha(k,i,j,I_BND_DENS) = max( alpha_z1, alpha_x1, alpha_y1 )
          ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELZ) = max( alpha_z2, alpha_x1, alpha_y1 )
          ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELX) = max( alpha_z1, alpha_x2, alpha_y1 )
          ATMOS_BOUNDARY_alpha(k,i,j,I_BND_VELY) = max( alpha_z1, alpha_x1, alpha_y2 )
          ATMOS_BOUNDARY_alpha(k,i,j,I_BND_POTT) = max( alpha_z1, alpha_x1, alpha_y1 )
          ATMOS_BOUNDARY_alpha(k,i,j,I_BND_QV  ) = max( alpha_z1, alpha_x1, alpha_y1 )
       end if
    enddo
    enddo
    enddo

    if ( .NOT. ATMOS_BOUNDARY_USE_VELZ ) then
       ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELZ) = 0.0_RP
    end if
    if ( .NOT. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       if ( .NOT. ATMOS_BOUNDARY_USE_DENS ) then
          ATMOS_BOUNDARY_alpha(:,:,:,I_BND_DENS) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_VELX ) then
          ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELX) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_VELY ) then
          ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELY) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_POTT ) then
          ATMOS_BOUNDARY_alpha(:,:,:,I_BND_POTT) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_QV ) then
          ATMOS_BOUNDARY_alpha(:,:,:,I_BND_QV  ) = 0.0_RP
       end if
    end if


    call ATMOS_BOUNDARY_alpha_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_setalpha

  !-----------------------------------------------------------------------------
  !> Read boundary data
  subroutine ATMOS_BOUNDARY_setinitval( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    integer :: i, j, k
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ATMOS_BOUNDARY_var(k,i,j,I_BND_DENS) = DENS(k,i,j)
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELZ) = MOMZ(k,i,j) / ( DENS(k,i,j)+DENS(k+1,i,j) ) * 2.0_RP
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELX) = MOMX(k,i,j) / ( DENS(k,i,j)+DENS(k+1,i,j) ) * 2.0_RP
       ATMOS_BOUNDARY_var(k,i,j,I_BND_VELY) = MOMY(k,i,j) / ( DENS(k,i,j)+DENS(k+1,i,j) ) * 2.0_RP
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
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_IN_BASENAME

    if (      ATMOS_BOUNDARY_USE_DENS &
         .or. ATMOS_BOUNDARY_USE_VELZ &
         .or. ATMOS_BOUNDARY_USE_VELX &
         .or. ATMOS_BOUNDARY_USE_VELY &
         .or. ATMOS_BOUNDARY_USE_POTT &
         ) then
       call FileRead( reference_atmos(:,:,:), bname, 'DENS', 1, PRC_myrank )
       ATMOS_BOUNDARY_var(KS:KE,IS:IE,JS:JE,I_BND_DENS) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    end if
    if ( ATMOS_BOUNDARY_USE_DENS ) then
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_DENS', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha(KS:KE,IS:IE,JS:JE,I_BND_DENS) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

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

    real(RP) :: buffer(KA,IA,JA)
    !---------------------------------------------------------------------------

    if (      ATMOS_BOUNDARY_USE_DENS &
         .or. ATMOS_BOUNDARY_USE_VELZ &
         .or. ATMOS_BOUNDARY_USE_VELX &
         .or. ATMOS_BOUNDARY_USE_VELY &
         .or. ATMOS_BOUNDARY_USE_POTT &
         ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_DENS),                  &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'DENS', 'Reference Density', 'kg/m3', 'ZXY',           &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    end if
    if ( ATMOS_BOUNDARY_USE_DENS .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_DENS),                &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_DENS', 'Alpha for dens', '1', 'ZXY',            &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_VELZ), &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELZ', 'Reference Velocity w', 'm/s', 'ZXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELZ),                &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELZ', 'Alpha for w', '1', 'ZXY',               &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELX .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_VELX), &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELX', 'Reference Velocity u', 'm/s', 'ZXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELX),                &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELX', 'Alpha for u', '1', 'ZXY',               &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELY .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_VELY), &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELY', 'Reference Velocity y', 'm/s', 'ZXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_VELY),                &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELY', 'Alpha for v', '1', 'ZXY',               &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_POTT .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       call FILEIO_write( ATMOS_BOUNDARY_var(:,:,:,I_BND_POTT), &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'POTT', 'Reference PT', 'K', 'ZXY',                    &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha(:,:,:,I_BND_POTT),                &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_POTT', 'Alpha for PT', '1', 'ZXY',              &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_QV .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
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
    use scale_atmos_refstate, only: &
         ATMOS_REFSTATE_dens
    implicit none

    integer :: i, j, k
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ATMOS_BOUNDARY_var(k,i,j,I_BND_DENS) = ATMOS_REFSTATE_DENS(k,i,j)
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
  subroutine ATMOS_BOUNDARY_initialize( DENS )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_time, only: &
       TIME_NOWDATE,      &
       TIME_OFFSET_YEAR,  &
       TIME_DTSEC,        &
       TIME_NOWDAYSEC
    use scale_calendar, only: &
       CALENDAR_date2daysec,    &
       CALENDAR_combine_daysec
    implicit none
    real(RP), intent(in) :: DENS(KA,IA,JA)

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    integer :: run_time_startdate(6)
    integer :: run_time_startday
    real(RP) :: run_time_startsec
    real(RP) :: run_time_startms
    real(RP) :: run_time_nowdaysec
    integer :: boundary_time_startday
    real(RP) :: boundary_time_startsec
    real(RP) :: boundary_time_startms
    real(RP) :: boundary_time_initdaysec
    real(RP) :: boundary_diff_daysec
    real(RP) :: boundary_inc_offset
    integer  :: fillgaps_steps

    character(len=H_LONG) :: bname

    integer  :: i, j, k, iv
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_IN_BASENAME
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I5.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') '*** BOUNDARY START Date     : ', &
         ATMOS_BOUNDARY_START_DATE(1),'/',ATMOS_BOUNDARY_START_DATE(2),'/',ATMOS_BOUNDARY_START_DATE(3),' ',  &
         ATMOS_BOUNDARY_START_DATE(4),':',ATMOS_BOUNDARY_START_DATE(5),':',ATMOS_BOUNDARY_START_DATE(6)

    if ( ATMOS_BOUNDARY_START_DATE(1) == -9999 ) then
       boundary_timestep = 1
       boundary_inc_offset = 0.D0
       fillgaps_steps = 0
    else
       !--- recalculate time of the run [no offset]
       run_time_startms = 0.0_DP
       run_time_startdate(:) = TIME_NOWDATE(:)
       run_time_startdate(1) = TIME_OFFSET_YEAR
       call CALENDAR_date2daysec( run_time_startday,     & ! [OUT]
                                  run_time_startsec,     & ! [OUT]
                                  run_time_startdate(:), & ! [IN]
                                  run_time_startms       ) ! [IN]
       run_time_nowdaysec  = CALENDAR_combine_daysec( run_time_startday, run_time_startsec )

       !--- calculate time of the initial step in boundary file [no offset]
       boundary_time_startms = 0.0_DP
       call CALENDAR_date2daysec( boundary_time_startday,      & ! [OUT]
                                  boundary_time_startsec,       & ! [OUT]
                                  ATMOS_BOUNDARY_START_DATE(:), & ! [IN]
                                  boundary_time_startms         ) ! [IN]
       boundary_time_initdaysec  = CALENDAR_combine_daysec( boundary_time_startday, boundary_time_startsec )

       boundary_diff_daysec = run_time_nowdaysec - boundary_time_initdaysec
       boundary_timestep = 1 + int( boundary_diff_daysec / ATMOS_BOUNDARY_UPDATE_DT )
       boundary_inc_offset = mod( boundary_diff_daysec, ATMOS_BOUNDARY_UPDATE_DT )
       fillgaps_steps = int( boundary_inc_offset / TIME_DTSEC )
    endif

    if( IO_L ) write(IO_FID_LOG,*) '+++ BOUNDARY TIMESTEP NUMBER FOR INIT:', boundary_timestep
    if( IO_L ) write(IO_FID_LOG,*) '+++ BOUNDARY INCREMENT OFFSET:', boundary_inc_offset
    if( IO_L ) write(IO_FID_LOG,*) '+++ BOUNDARY FILLGAPS STEPS:', fillgaps_steps
    ! read boundary data from input file
    call FileRead( reference_atmos(:,:,:), bname, 'DENS', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_DENS,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELX', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELX,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELY', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELY,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'POTT', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_POTT,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'QV',   boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_QV  ,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)

    boundary_timestep = boundary_timestep + 1
    call FileRead( reference_atmos(:,:,:), bname, 'DENS', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_DENS,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELX', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELX,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELY', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELY,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'POTT', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_POTT,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'QV',   boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_QV  ,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)

    do iv = 1, I_BND_SIZE
    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_var_ref(   1:KS-1,i,j,iv,1) = ATMOS_BOUNDARY_var_ref(KS,i,j,iv,1)
       ATMOS_BOUNDARY_var_ref(KE+1:KA,  i,j,iv,1) = ATMOS_BOUNDARY_var_ref(KE,i,j,iv,1)
       ATMOS_BOUNDARY_var_ref(   1:KS-1,i,j,iv,2) = ATMOS_BOUNDARY_var_ref(KS,i,j,iv,2)
       ATMOS_BOUNDARY_var_ref(KE+1:KA,  i,j,iv,2) = ATMOS_BOUNDARY_var_ref(KE,i,j,iv,2)
    enddo
    enddo
    enddo

    do iv = 1, I_BND_SIZE
       call COMM_vars8( ATMOS_BOUNDARY_var_ref(:,:,:,iv,1), iv )
       call COMM_vars8( ATMOS_BOUNDARY_var_ref(:,:,:,iv,2), iv+I_BND_SIZE )
    enddo

    do iv = 1, I_BND_SIZE
       call COMM_wait( ATMOS_BOUNDARY_var_ref(:,:,:,iv,1), iv )
       call COMM_wait( ATMOS_BOUNDARY_var_ref(:,:,:,iv,2), iv+I_BND_SIZE )
    enddo

    ! set boundary data and time increment
    do iv = 1, I_BND_SIZE
       if ( I_BND_SIZE == I_MOMZ ) cycle
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

    ! fill in gaps of the offset
    do iv = 1, I_BND_SIZE
    if ( iv==I_BND_VELZ ) cycle
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       ATMOS_BOUNDARY_var(k,i,j,iv) = ATMOS_BOUNDARY_var(k,i,j,iv) + ATMOS_BOUNDARY_increment(k,i,j,iv) * dble(fillgaps_steps)
    enddo
    enddo
    enddo
    enddo

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       ATMOS_BOUNDARY_var(:,:,:,I_BND_VELZ) = ATMOS_BOUNDARY_VALUE_VELZ
    end if

    last_updated = TIME_NOWDAYSEC - boundary_inc_offset

    return
  end subroutine ATMOS_BOUNDARY_initialize

  !-----------------------------------------------------------------------------
  !> Update boundary value with a constant time increment
  subroutine ATMOS_BOUNDARY_update
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       epsilon => CONST_EPS
    use scale_time, only: &
       TIME_DTSEC, &
       TIME_NOWDAYSEC
    implicit none

    integer  :: i, j, k, iv
    !---------------------------------------------------------------------------

    if ( ATMOS_BOUNDARY_TYPE == 'REAL' ) then

       integrated_sec = TIME_NOWDAYSEC - last_updated

       if ( integrated_sec >= ATMOS_BOUNDARY_UPDATE_DT - epsilon ) then
          boundary_timestep = boundary_timestep + 1

          call ATMOS_BOUNDARY_updatefile

          do iv = 1, I_BND_SIZE
             if ( iv == I_BND_VELZ ) cycle
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

          last_updated = last_updated + ATMOS_BOUNDARY_UPDATE_DT
       else
          do iv = 1, I_BND_SIZE
          if ( iv==I_BND_VELZ ) cycle
          do j  = 1, JA
          do i  = 1, IA
          do k  = 1, KA
             ATMOS_BOUNDARY_var(k,i,j,iv) = ATMOS_BOUNDARY_var(k,i,j,iv) + ATMOS_BOUNDARY_increment(k,i,j,iv)
          enddo
          enddo
          enddo
          enddo
       endif

    else
       write(*,*) 'xxx [BUG] unsupported type'
       call PRC_MPIstop
    end if

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
    implicit none

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=H_LONG) :: bname

    integer :: i, j, k, iv
    !---------------------------------------------------------------------------

    if (IO_L) write(IO_FID_LOG,*)"*** Atmos Boundary: read from boundary file(timestep=", boundary_timestep, ")"

    bname = ATMOS_BOUNDARY_IN_BASENAME

    do iv = 1, I_BND_SIZE
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       ATMOS_BOUNDARY_var_ref(k,i,j,iv,1) = ATMOS_BOUNDARY_var_ref(k,i,j,iv,2)
    enddo
    enddo
    enddo
    enddo

    call FileRead( reference_atmos(:,:,:), bname, 'DENS', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_DENS,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELX', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELX,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELY', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_VELY,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'POTT', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_POTT,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'QV',   boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_var_ref(KS:KE,IS:IE,JS:JE,I_BND_QV  ,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)

    do iv = 1, I_BND_SIZE
    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_var_ref(   1:KS-1,i,j,iv,2) = ATMOS_BOUNDARY_var_ref(KS,i,j,iv,2)
       ATMOS_BOUNDARY_var_ref(KE+1:KA,  i,j,iv,2) = ATMOS_BOUNDARY_var_ref(KE,i,j,iv,2)
    enddo
    enddo
    enddo

    do iv = 1, I_BND_SIZE
       call COMM_vars8( ATMOS_BOUNDARY_var_ref(:,:,:,iv,2), iv )
    enddo

    do iv = 1, I_BND_SIZE
       call COMM_wait ( ATMOS_BOUNDARY_var_ref(:,:,:,iv,2), iv )
    enddo

    return
  end subroutine ATMOS_BOUNDARY_updatefile

end module scale_atmos_boundary
