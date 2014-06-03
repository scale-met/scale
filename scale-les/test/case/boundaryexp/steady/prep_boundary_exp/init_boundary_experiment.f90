!-------------------------------------------------------------------------------
!> Program initial boundary experiment
!!
!! @par Description
!!          This program is a tool to prepare an external initial data
!!          in order to assess correctness in initial/boundary input.
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-29 (R.Yoshida)   [new]
!!
!<
!-------------------------------------------------------------------------------
program init_boundary_experiment
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use netcdf       ! [external]
  !
  !-----------------------------------------------------------------------------
  implicit none
!  include '/usr/local/intel/netcdf-4.1.3/include/netcdf.inc'
!  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
!  private :: namelist_check
!  private :: grid_preparation
!  private :: steady_case
!  private :: steadystate
!  private :: netcdf_writeout
!  private :: put_attributes
!  private :: data_writeout
!  private :: handle_err
!  private :: vid_err
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

  integer, parameter :: SP       = 4
  integer, parameter :: DP       = 8
  integer, parameter :: CHA      = 128
  integer, parameter :: CHA_NCIO = 20
  integer, parameter :: FID_NML  = 20
  integer, parameter :: nidmax = 37

  real(DP), parameter :: CONST_PI     = 3.14159265358979_DP
  real(DP), parameter :: CONST_RADIUS = 6.37122E+6_DP       !< radius of the planet [m]
  real(DP), parameter :: CONST_OHM    = 7.2920E-5_DP        !< angular velocity of the planet [1/s]
  real(DP), parameter :: CONST_GRAV   = 9.80665_DP          !< standard acceleration of gravity [m/s2]
  real(DP), parameter :: CONST_Rdry   =  287.04_DP          !< specific gas constant (dry air)           [J/kg/K]
  real(DP), parameter :: CONST_CPdry  = 1004.64_DP          !< specific heat (dry air,constant pressure) [J/kg/K]

  real(DP), parameter :: P0           = 1000.00E+2_DP       !< standard surface pressure [Pa]
  real(DP), parameter :: U0           = 5.00_DP             !< standard zonal wind speed [m]
  real(DP), parameter :: V0           = 0.00_DP             !< standard meridional wind speed [m]
  real(DP), parameter :: W0           = 0.00_DP             !< standard vertical wind speed [m]
  real(DP), parameter :: T0           = 300.00_DP           !< standard surface temperature [K]
  real(DP), parameter :: UNDEF        = -999.d15
  
  integer :: IA    = 100
  integer :: JA    = 100
  integer :: KA    = 50
  integer :: NSTEP = 20
  integer :: SOIL  = 4

  real(SP), allocatable :: P(:,:,:)        !< pertubation pressure
  real(SP), allocatable :: PB(:,:,:)       !< base pressure
  real(SP), allocatable :: PHB(:,:,:)      !< geopotential
  real(SP), allocatable :: T(:,:,:)        !< potential temperature
  real(SP), allocatable :: W(:,:,:)        !< vertical wind
  real(SP), allocatable :: U(:,:,:)        !< zonal wind
  real(SP), allocatable :: V(:,:,:)        !< meridional wind
  real(SP), allocatable :: QVAPOR(:,:,:)   !< passive tracer
  real(SP), allocatable :: TSLB(:,:,:)     !< soil temperature
  real(SP), allocatable :: SKINTMP(:,:)    !< skin temperature
  real(SP), allocatable :: ALBEDO(:,:)     !< albedo
  real(SP), allocatable :: EMISS(:,:)      !< emissivity
  real(SP), allocatable :: DUMMY_2D(:,:)   !< dummy variable
  real(SP), allocatable :: DUMMY_3D(:,:,:) !< dummy variable
  real(SP), allocatable :: TIM_SEC(:)      !< time stamp [sec]
  real(SP), allocatable :: ZS(:)           !< depth of soil layers
  real(SP), allocatable :: DZS(:)          !< thickness of soil layers

  real(SP), allocatable :: XLAT(:,:)     !< no staggard
  real(SP), allocatable :: XLAT_U(:,:)   !< x  staggard
  real(SP), allocatable :: XLAT_V(:,:)   !< y  staggard
  real(SP), allocatable :: XLONG(:,:)    !< no staggard
  real(SP), allocatable :: XLONG_U(:,:)  !< x  staggard
  real(SP), allocatable :: XLONG_V(:,:)  !< y  staggard

  real(SP), allocatable :: CX(:)
  real(SP), allocatable :: FX(:)
  real(SP), allocatable :: CY(:)
  real(SP), allocatable :: FY(:)
  real(SP), allocatable :: CZ(:)
  real(SP), allocatable :: FZ(:)
  real(SP), allocatable :: FS(:)

  real(SP) :: lat0          = 0.00_SP        !< [degree]
  real(SP) :: lon0          = 0.00_SP        !< [degree]
  real(SP) :: dx            = 500.0_SP       !< [m]
  real(SP) :: dy            = 500.0_SP       !< [m]
  real(SP) :: height_bottom = 0.00_SP        !< [m]
  real(SP) :: height_top    = 10000.00_SP    !< [m]
  real(SP) :: dt            = 60.00_SP       !< [sec] time interval of external file

  character(len=CHA) :: testcase  = "STEADY"
  character(len=CHA) :: basename  = "boundary_exp"
  character(len=CHA) :: foutput   = "boundary_exp_00000"
  character(len=CHA) :: fnamelist = "init_boundary_exp.conf"

  character(len=CHA_NCIO) ::  nc_att_title    = "boundary experiment"
  character(len=CHA_NCIO) ::  nc_att_grid     = "lat-lon grid"
  character(len=CHA_NCIO) ::  nc_att_equation = "not available"
  character(len=CHA_NCIO) ::  nc_att_model    = "WRF-ARW wrfout"
  character(len=CHA_NCIO) ::  nc_att_descrip  = "version 4.03"

  integer :: VDIMS1DT
  integer :: VDIMS1DS(2)
  integer :: VDIMS1DZ, VDIMS1DZH
  integer :: VDIMS1DX, VDIMS1DXH
  integer :: VDIMS1DY, VDIMS1DYH
  integer :: VDIMS2D(3), VDIMS2DX(3), VDIMS2DY(3)
  integer :: VDIMS3D(4), VDIMS3DX(4), VDIMS3DY(4)
  integer :: VDIMS3DZ(4), VDIMS3DS(4)

  integer :: vid(nidmax)

  integer :: CZdimID, FZdimID, CXdimID, FXdimID
  integer :: CYdimID, FYdimID, SdimID, TdimID

  real(SP) :: cpu_t1, cpu_t2
  integer :: i ,j ,k, l
  !=============================================================================

  NAMELIST / PARAM_DIMENSION / &
     IA,             &
     JA,             &
     KA,             &
     NSTEP

  NAMELIST / PARAM_BOUNDARY_EXP / &
     lat0,           &
     lon0,           &
     dx,             &
     dy,             &
     height_bottom,  &
     height_top,     &
     testcase,       &
     basename,       &
     dt
  !-----------------------------------------------------------------------------
  call cpu_time(cpu_t1)  ! timer

  open (FID_NML, file=trim(fnamelist), status='old', delim='apostrophe')
  read (FID_NML, nml=PARAM_DIMENSION)
  read (FID_NML, nml=PARAM_BOUNDARY_EXP)
  close (FID_NML)
  call namelist_check

  allocate ( CX(IA-1) )
  allocate ( FX(IA  ) )
  allocate ( CY(JA-1) )
  allocate ( FY(JA  ) )
  allocate ( CZ(KA-1) )
  allocate ( FZ(KA  ) )
  allocate ( FS(SOIL) )
  allocate ( XLAT(IA-1,JA-1) )
  allocate ( XLAT_U(IA,JA-1) )
  allocate ( XLAT_V(IA-1,JA) )
  allocate ( XLONG(IA-1,JA-1) )
  allocate ( XLONG_U(IA,JA-1) )
  allocate ( XLONG_V(IA-1,JA) )

  call grid_preparation

  allocate ( P(IA-1,        JA-1, KA-1) )
  allocate ( PB(IA-1,       JA-1, KA-1) )
  allocate ( PHB(IA-1,      JA-1, KA  ) )
  allocate ( T(IA-1,        JA-1, KA-1) )
  allocate ( W(IA-1,        JA-1, KA  ) )
  allocate ( U(IA,          JA-1, KA-1) )
  allocate ( V(IA-1,        JA,   KA-1) )
  allocate ( QVAPOR(IA-1,   JA-1, KA-1) )
  allocate ( TSLB(IA-1,     JA-1, SOIL) )
  allocate ( SKINTMP(IA-1,  JA-1 )      )
  allocate ( ALBEDO(IA-1,   JA-1 )      )
  allocate ( EMISS(IA-1,    JA-1 )      )
  allocate ( DUMMY_2D(IA-1, JA-1 )      )
  allocate ( DUMMY_3D(IA-1, JA-1, KA-1) )
  allocate ( TIM_SEC(NSTEP)             )
  allocate ( ZS(SOIL)                   )
  allocate ( DZS(SOIL)                  )

  call time_stamp

  select case(trim(testcase))
  case('STEADY')
     write(*,*) '';
     write(*,*) '+++ STEADT STATE TEST CASE';
     call steady_case

  case('ZONALWIND')
     write(*,*) ' xxx ZONALWIND: Preparing...'; stop

  case('ADVECTION')
     write(*,*) ' xxx ADVECTION: Preparing...'; stop

  case default
     write(*,*) ' xxx Unsupported Test Case:', trim(testcase)
     stop

  endselect

  deallocate ( P )
  deallocate ( PB )
  deallocate ( PHB )
  deallocate ( T )
  deallocate ( W )
  deallocate ( U )
  deallocate ( V )
  deallocate ( QVAPOR )
  deallocate ( TSLB )
  deallocate ( SKINTMP )
  deallocate ( ALBEDO )
  deallocate ( EMISS )
  deallocate ( DUMMY_2D )
  deallocate ( DUMMY_3D )
  deallocate ( ZS )
  deallocate ( DZS )
  deallocate ( XLAT )
  deallocate ( XLAT_U )
  deallocate ( XLAT_V )
  deallocate ( XLONG )
  deallocate ( XLONG_U )
  deallocate ( XLONG_V )
  deallocate ( CX )
  deallocate ( FX )
  deallocate ( CY )
  deallocate ( FY )
  deallocate ( CZ )
  deallocate ( FZ )
  deallocate ( FS )

  call cpu_time(cpu_t2)
  write (*, '(1X,"CPU Time: ",F7.4," [sec]")') (cpu_t2 - cpu_t1)

  stop
  !=============================================================================
contains
  !-----------------------------------------------------------------------------
  !> Grid Preparation
  subroutine grid_preparation
    implicit none

    real(SP) :: dz
    real(DP) :: ddeg
    integer :: i, j, k
    !---------------------------------------------------------------------------

    FX(1) = 0.0_SP
    do i=2, IA
       FX(i) = FX(i-1) + dx
    enddo
    do i=1, IA-1
       CX(i) = (FX(i) + FX(i+1)) * 0.5_SP
    enddo

    FY(1) = 0.0_SP
    do j=2, JA
       FY(j) = FY(j-1) + dy
    enddo
    do j=1, JA-1
       CY(j) = (FY(j) + FY(j+1)) * 0.5_SP
    enddo

    dz = (height_top - height_bottom) / float(KA-1)
    FZ(1) = height_bottom
    do k=2, KA
       FZ(k) = FZ(k-1) + dz
    enddo
    do k=1, KA-1
       CZ(k) = (FZ(k) + FZ(k+1)) * 0.5_SP
    enddo

    ddeg = 360.0_DP / (2.0_DP * CONST_RADIUS * CONST_PI)
    do j=1, JA-1
    do i=1, IA-1
       XLAT(i,j)  = CY(j) * ddeg
       XLONG(i,j) = CX(i) * ddeg
    enddo
    enddo
    do j=1, JA-1
    do i=1, IA
       XLAT_U(i,j)  = CY(j) * ddeg
       XLONG_U(i,j) = FX(i) * ddeg
    enddo
    enddo
    do j=1, JA
    do i=1, IA-1
       XLAT_V(i,j)  = FY(j) * ddeg
       XLONG_V(i,j) = CX(i) * ddeg
    enddo
    enddo

    return
  end subroutine grid_preparation

  !-----------------------------------------------------------------------------
  !> Initialize STEADY CASE
  subroutine steady_case
    implicit none

    real(DP), allocatable :: pre(:)  !< pressure at steady state
    real(DP), allocatable :: tem(:)  !< potential temperature at steady state
    integer :: i, j, k, l
    !---------------------------------------------------------------------------

    allocate ( pre(KA-1) )
    allocate ( tem(KA-1) )

    call steadystate( pre,tem )

    do j=1, JA-1
    do i=1, IA-1
       PB(i,j,:) = pre(:)
       T(i,j,:) = tem(:)
    enddo
    enddo

    do k=1, KA
    do j=1, JA-1
    do i=1, IA-1
       PHB(i,j,k) = FZ(k) * CONST_GRAV
    enddo
    enddo
    enddo

    P(:,:,:)         = 0.0_SP
    U(:,:,:)         = U0
    V(:,:,:)         = V0
    W(:,:,:)         = W0
    QVAPOR(:,:,:)    = 0.0_SP

    TSLB(:,:,:)      = T0
    SKINTMP(:,:)     = T0
    ALBEDO(:,:)      = 0.2_SP
    EMISS(:,:)       = 1.0_SP - ALBEDO(:,:)

    DUMMY_2D(:,:)    = 0.0_SP
    DUMMY_3D(:,:,:)  = 0.0_SP

    !ZS(0) = 0.0_SP
       DZS(1) = 0.5_SP
    ZS(1) = 0.5_SP
       DZS(2) = 0.5_SP
    ZS(2) = 1.0_SP
       DZS(3) = 0.5_SP
    ZS(3) = 1.5_SP
       DZS(4) = 0.5_SP
    ZS(4) = 2.0_SP

    do l=1, NSTEP
       call netcdf_writeout( l )
    enddo

    deallocate( pre )
    deallocate( tem )

    return
  end subroutine steady_case

  !-----------------------------------------------------------------------------
  !> Calculation Steady State
  subroutine steadystate( &
       pre,  &
       tem   )
    implicit none

    real(DP), intent(out) :: pre(:)  !< pressure at steady state
    real(DP), intent(out) :: tem(:)  !< potential temperature at steady state

    real(DP) :: KAPPA, G, R
    real(DP) :: pre_save, f, df, dz
    real(DP), parameter :: eps = 1.0E-7_DP
    integer, parameter :: itrmax = 100
    integer :: k, itr
    !---------------------------------------------------------------------------

    KAPPA = CONST_Rdry / CONST_CPdry
    G = CONST_GRAV
    R = CONST_Rdry
    pre(1) = P0
    tem(1) = T0

    do k=2, KA-1

       ! first guess
       pre(k) = pre(k-1)
       tem(k) = 300.D0 * ( pre(k)/P0 )**KAPPA
       tem(k) = max( 200.D0, tem(k) )

       ! Newton-Lapson
       do itr = 1, itrmax
          pre_save = pre(k) ! save

          dz = FZ(k+1) - FZ(k)
          f  = log(pre(k)/pre(k-1)) / dz + G / ( R * 0.5D0 * (tem(k)+tem(k-1)) )
          df = 1.D0 / (pre(k)*dz)

          pre(k) = pre(k) - f / df
          tem(k) = 300.D0 * ( pre(k)/P0 )**KAPPA
          tem(k) = max( 200.D0, tem(k) )

          if( abs(pre_save-pre(k)) <= eps ) exit
       enddo

       if ( itr > itrmax ) then
          write(*,*) 'xxx iteration not converged!', &
                      k, pre_save-pre(k), pre(k), pre(k-1), tem(k), tem(k-1)
          stop
       endif
    enddo

    return
  end subroutine steadystate

  !-----------------------------------------------------------------------------
  !> Output Parameters as a Namelist Check
  subroutine netcdf_writeout ( &
      step   )
    use netcdf       ! [external]
    implicit none

    integer, intent(in) :: step  !< number of time step

    integer :: ncid, nid, timestep
    integer :: status, ierr
    character(5) :: num, nmax
    !---------------------------------------------------------------------------

    write(num,'(I5.5)') step-1
    write(nmax,'(I5.5)') NSTEP-1
    write(*, *) "    output netcdf file: ",num,"/",nmax
    foutput = trim(basename)//"_"//num

    ! open output file
    status = nf90_create (trim(foutput), ior(nf90_clobber,nf90_64bit_offset), ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    ! define the dimensions
    status = nf90_def_dim(ncid, 'bottom_top',       KA-1,           CZdimID)
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, 'bottom_top_stag',  KA,             FZdimID)
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, 'west_east',        IA-1,           CXdimID)
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, 'west_east_stag',   IA,             FXdimID)
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, 'south_north',      JA-1,           CYdimID)
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, 'south_north_stag', JA,             FYdimID)
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, 'soil_layers_stag', SOIL,           SdimID)
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_def_dim(ncid, 'time',             nf90_unlimited, TdimID)
    if (status .ne. nf90_noerr) call handle_err(status)

    VDIMS1DT     = TDimID
    VDIMS1DZ     = CZDimID; VDIMS1DZH    = FZDimID
    VDIMS1DX     = CXDimID; VDIMS1DXH    = FXDimID
    VDIMS1DY     = CYDimID; VDIMS1DYH    = FYDimID
    VDIMS1DS(1)  = SDimID;  VDIMS1DS(2)  = TDimID
    VDIMS2D(1)   = CXDimID; VDIMS2D(2)   = CYDimID; VDIMS2D(3)  = TDimID
    VDIMS2DX(1)  = FXDimID; VDIMS2DX(2)  = CYDimID; VDIMS2DX(3) = TDimID
    VDIMS2DY(1)  = CXDimID; VDIMS2DY(2)  = FYDimID; VDIMS2DY(3) = TDimID
    VDIMS3D(1)   = CXDimID; VDIMS3D(2)   = CYDimID; VDIMS3D(3)  = CZDimID; VDIMS3D(4)  = TDimID
    VDIMS3DX(1)  = FXDimID; VDIMS3DX(2)  = CYDimID; VDIMS3DX(3) = CZDimID; VDIMS3DX(4) = TDimID
    VDIMS3DY(1)  = CXDimID; VDIMS3DY(2)  = FYDimID; VDIMS3DY(3) = CZDimID; VDIMS3DY(4) = TDimID
    VDIMS3DZ(1)  = CXDimID; VDIMS3DZ(2)  = CYDimID; VDIMS3DZ(3) = FZDimID; VDIMS3DZ(4) = TDimID
    VDIMS3DS(1)  = CXDimID; VDIMS3DS(2)  = CYDimID; VDIMS3DS(3) = SDimID;  VDIMS3DS(4) = TDimID

    ! put title attributes
    status = nf90_put_att( ncid, nf90_global, 'title',            nc_att_title   )
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att( ncid, nf90_global, 'grid',             nc_att_grid    )
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att( ncid, nf90_global, 'equation',         nc_att_equation)
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att( ncid, nf90_global, 'test_case',        testcase       )
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att( ncid, nf90_global, 'model',            nc_att_model   )
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att( ncid, nf90_global, 'description',      nc_att_descrip )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! put variables attributes
    call put_attributes( ncid )

    ! leave define mode
    status = nf90_enddef( ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! data writeout
    timestep = 1
    call data_writeout( timestep, ncid )

    ! close output file
    status = nf90_close( ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    return
  end subroutine netcdf_writeout

  !-----------------------------------------------------------------------------
  !> Put Variables Attributes
  subroutine put_attributes ( &
      ncid   )
    use netcdf       ! [external]
    implicit none

    integer, intent(in) :: ncid  !< netcdf file id
    integer :: nid
    integer :: status, ierr
    !---------------------------------------------------------------------------

    nid = 1
    status = nf90_def_var (ncid, 'bottom_top', nf90_float, VDIMS1DZ, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'units', 'meters')
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'positive', 'up')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 2
    status = nf90_def_var (ncid, 'bottom_top_stag', nf90_float, VDIMS1DZH, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'units', 'meters')
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'positive', 'up')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 3
    status = nf90_def_var (ncid, 'west_east', nf90_float, VDIMS1DX, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'standard_name', 'west_east')
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'units', 'none')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 4
    status = nf90_def_var (ncid, 'west_east_stag', nf90_float, VDIMS1DXH, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'standard_name', 'west_east_stag')
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'units', 'none')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 5
    status = nf90_def_var (ncid, 'south_north', nf90_float, VDIMS1DY, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'standard_name', 'south_north')
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'units', 'none')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 6
    status = nf90_def_var (ncid, 'south_north_stag', nf90_float, VDIMS1DYH, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'standard_name', 'south_north_stag')
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'units', 'none')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 7
    status = nf90_def_var (ncid, 'soil_layers_stag', nf90_float, VDIMS1DS, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'standard_name', 'soil_layers_stag')
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'units', 'none')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 8
    status = nf90_def_var (ncid, 'time', nf90_float, VDIMS1DT, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'standard_name', 'time')
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'units', 'sec. since 0000-01-01 00:00:00')
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att (ncid, vid(nid), 'calender', 'noleap')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 9
    status = nf90_def_var(ncid, 'P',  nf90_float, VDIMS3D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'perturbation pressure')
    status = nf90_put_att(ncid, vid(nid), 'units', 'Pa')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 10
    status = nf90_def_var(ncid, 'PB',  nf90_float, VDIMS3D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'BASE STATE PRESSURE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'Pa')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 11
    status = nf90_def_var(ncid, 'PHB',  nf90_float, VDIMS3DZ, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'base-state geopotential')
    status = nf90_put_att(ncid, vid(nid), 'units', 'm2 s-2')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'Z')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 12
    status = nf90_def_var(ncid, 'T',  nf90_float, VDIMS3D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'perturbation potential temperature (theta-t0)')
    status = nf90_put_att(ncid, vid(nid), 'units', 'K')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 13
    status = nf90_def_var(ncid, 'U',  nf90_float, VDIMS3DX, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'x-wind component')
    status = nf90_put_att(ncid, vid(nid), 'units', 'm s-1')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'X')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 14
    status = nf90_def_var(ncid, 'V',  nf90_float, VDIMS3DY, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'y-wind component')
    status = nf90_put_att(ncid, vid(nid), 'units', 'm s-1')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'Y')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 15
    status = nf90_def_var(ncid, 'W',  nf90_float, VDIMS3DZ, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'z-wind component')
    status = nf90_put_att(ncid, vid(nid), 'units', 'm s-1')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'Z')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 16
    status = nf90_def_var(ncid, 'QVAPOR',  nf90_float, VDIMS3D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'Water vapor mixing ratio')
    status = nf90_put_att(ncid, vid(nid), 'units', 'kg kg-1')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 17
    status = nf90_def_var(ncid, 'QCLOUD',  nf90_float, VDIMS3D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'Cloud water mixing ratio')
    status = nf90_put_att(ncid, vid(nid), 'units', 'kg kg-1')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 18
    status = nf90_def_var(ncid, 'QRAIN',  nf90_float, VDIMS3D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'Rain water mixing ratio')
    status = nf90_put_att(ncid, vid(nid), 'units', 'kg kg-1')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 19
    status = nf90_def_var(ncid, 'QICE',  nf90_float, VDIMS3D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'Ice mixing ratio')
    status = nf90_put_att(ncid, vid(nid), 'units', 'kg kg-1')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 20
    status = nf90_def_var(ncid, 'QSNOW',  nf90_float, VDIMS3D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'Snow mixing ratio')
    status = nf90_put_att(ncid, vid(nid), 'units', 'kg kg-1')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 21
    status = nf90_def_var(ncid, 'QGRAUP',  nf90_float, VDIMS3D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'Graupel mixing ratio')
    status = nf90_put_att(ncid, vid(nid), 'units', 'kg kg-1')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 22
    status = nf90_def_var(ncid, 'ZS',  nf90_float, VDIMS1DS, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'Z  ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'DEPTHS OF CENTERS OF SOIL LAYERS')
    status = nf90_put_att(ncid, vid(nid), 'units', 'm')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'Z')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 23
    status = nf90_def_var(ncid, 'DZS',  nf90_float, VDIMS1DS, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'Z  ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'THICKNESSES OF SOIL LAYERS')
    status = nf90_put_att(ncid, vid(nid), 'units', 'm')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'Z')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 24
    status = nf90_def_var(ncid, 'TSLB',  nf90_float, VDIMS3DS, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XYZ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'SOIL TEMPERATURE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'K')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'Z')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 25
    status = nf90_def_var(ncid, 'SFROFF',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'SURFACE RUNOFF')
    status = nf90_put_att(ncid, vid(nid), 'units', 'mm')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 26
    status = nf90_def_var(ncid, 'SST',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'SEA SURFACE TEMPERATURE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'K')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 27
    status = nf90_def_var(ncid, 'TSK',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'SURFACE SKIN TEMPERATURE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'K')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 28
    status = nf90_def_var(ncid, 'ALBEDO',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'ALBEDO')
    status = nf90_put_att(ncid, vid(nid), 'units', '-')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 29
    status = nf90_def_var(ncid, 'EMISS',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'SURFACE EMISSIVITY')
    status = nf90_put_att(ncid, vid(nid), 'units', '-')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 30
    status = nf90_def_var(ncid, 'ZNT',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'TIME-VARYING ROUGHNESS LENGTH')
    status = nf90_put_att(ncid, vid(nid), 'units', 'm')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 31
    status = nf90_def_var(ncid, 'SNOW',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'SNOW WATER EQUIVALENT')
    status = nf90_put_att(ncid, vid(nid), 'units', 'kg m-2')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 32
    status = nf90_def_var(ncid, 'TSNAV',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'AVERAGE SNOW TEMPERATURE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'C')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 33
    status = nf90_def_var(ncid, 'XLAT',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'LATITUDE, SOUTH IS NEGATIVE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'degree_north')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 34
    status = nf90_def_var(ncid, 'XLONG',  nf90_float, VDIMS2D, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'LONGITUDE, WEST IS NEGATIVE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'degree_east')
    status = nf90_put_att(ncid, vid(nid), 'stagger', '')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 35
    status = nf90_def_var(ncid, 'XLAT_U',  nf90_float, VDIMS2DX, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'LATITUDE, SOUTH IS NEGATIVE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'degree_north')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'X')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 36
    status = nf90_def_var(ncid, 'XLONG_U',  nf90_float, VDIMS2DX, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'LONGITUDE, WEST IS NEGATIVE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'degree_east')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'X')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 37
    status = nf90_def_var(ncid, 'XLAT_V',  nf90_float, VDIMS2DY, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'LATITUDE, SOUTH IS NEGATIVE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'degree_north')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'Y')
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 38
    status = nf90_def_var(ncid, 'XLONG_V',  nf90_float, VDIMS2DY, vid(nid))
    if (status .ne. nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, vid(nid), 'FieldType', '104')
    status = nf90_put_att(ncid, vid(nid), 'MemoryOrder', 'XY ')
    status = nf90_put_att(ncid, vid(nid), 'description', 'LONGITUDE, WEST IS NEGATIVE')
    status = nf90_put_att(ncid, vid(nid), 'units', 'degree_east')
    status = nf90_put_att(ncid, vid(nid), 'stagger', 'Y')
    if (status .ne. nf90_noerr) call handle_err(status)

    return
  end subroutine put_attributes 

  !-----------------------------------------------------------------------------
  !> Data Writeout
  subroutine data_writeout ( &
      step,  &
      ncid   )
    use netcdf       ! [external]
    implicit none

    integer, intent(in) :: step  !< time step id
    integer, intent(in) :: ncid  !< netcdf file id
    integer :: nid
    integer :: status, ierr, vid_check
    !---------------------------------------------------------------------------

    nid = 1
    status = nf90_inq_varid (ncid, 'bottom_top', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), CZ, start=(/1/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 2
    status = nf90_inq_varid (ncid, 'bottom_top_stag', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), FZ, start=(/1/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 3
    status = nf90_inq_varid (ncid, 'west_east', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), CX, start=(/1/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 4
    status = nf90_inq_varid (ncid, 'west_east_stag', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), FX, start=(/1/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 5
    status = nf90_inq_varid (ncid, 'south_north', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), CY, start=(/1/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 6
    status = nf90_inq_varid (ncid, 'south_north_stag', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), FY, start=(/1/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 7
    status = nf90_inq_varid (ncid, 'soil_layers_stag', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), FS, start=(/1/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 8
    status = nf90_inq_varid (ncid, 'time', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), TIM_SEC(step), start=(/1/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 9
    status = nf90_inq_varid (ncid,'P', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), P(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 10
    status = nf90_inq_varid (ncid,'PB', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), PB(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 11
    status = nf90_inq_varid (ncid,'PHB', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), PHB(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 12
    status = nf90_inq_varid (ncid,'T', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), T(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 13
    status = nf90_inq_varid (ncid,'U', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), U(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 14
    status = nf90_inq_varid (ncid,'V', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), V(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 15
    status = nf90_inq_varid (ncid,'W', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), W(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 16
    status = nf90_inq_varid (ncid,'QVAPOR', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), QVAPOR(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 17
    status = nf90_inq_varid (ncid,'QCLOUD', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), DUMMY_3D(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 18
    status = nf90_inq_varid (ncid,'QRAIN', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), DUMMY_3D(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 19
    status = nf90_inq_varid (ncid,'QICE', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), DUMMY_3D(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 20
    status = nf90_inq_varid (ncid,'QSNOW', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), DUMMY_3D(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 21
    status = nf90_inq_varid (ncid,'QGRAUP', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), DUMMY_3D(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 22
    status = nf90_inq_varid (ncid,'ZS', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), ZS(:), start=(/1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 23
    status = nf90_inq_varid (ncid,'DZS', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), DZS(:), start=(/1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 24
    status = nf90_inq_varid (ncid,'TSLB', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), TSLB(:,:,:), start=(/1,1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 25
    status = nf90_inq_varid (ncid,'SFROFF', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), DUMMY_2D(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 26
    status = nf90_inq_varid (ncid,'SST', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), SKINTMP(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 27
    status = nf90_inq_varid (ncid,'TSK', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), SKINTMP(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 28
    status = nf90_inq_varid (ncid,'ALBEDO', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), ALBEDO(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 29
    status = nf90_inq_varid (ncid,'EMISS', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), EMISS(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 30
    status = nf90_inq_varid (ncid,'ZNT', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), DUMMY_2D(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 31
    status = nf90_inq_varid (ncid,'SNOW', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), DUMMY_2D(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 32
    status = nf90_inq_varid (ncid,'TSNAV', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), SKINTMP(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 33
    status = nf90_inq_varid (ncid,'XLAT', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), XLAT(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 34
    status = nf90_inq_varid (ncid,'XLONG', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), XLONG(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 35
    status = nf90_inq_varid (ncid,'XLAT_U', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), XLAT_U(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 36
    status = nf90_inq_varid (ncid,'XLONG_U', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), XLONG_U(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 37
    status = nf90_inq_varid (ncid,'XLAT_V', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), XLAT_V(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    nid = 38
    status = nf90_inq_varid (ncid,'XLONG_V', vid_check)
    if (status .ne. nf90_noerr) call handle_err(status)
    call vid_err ( nid,vid_check )
    status = nf90_put_var (ncid, vid(nid), XLONG_V(:,:), start=(/1,1,step/))
    if (status .ne. nf90_noerr) call handle_err(status)

    return
  end subroutine data_writeout

  !-----------------------------------------------------------------------------
  !> Output Parameters as a Namelist Check
  subroutine namelist_check
    implicit none
    !---------------------------------------------------------------------------

    write(*,*) "### SCALE-LES: A tool to prepare an external data to assess"
    write(*,*) "###            correctness in initial/boundary input."
    write(*,*) "### Output File Type: WRF-ARW history output (netcdf3)"
    write(*,*) ""
    write(*,*) ">>  namelist: PARAM_DIMENSION"
    write(*,*) "    IA = ", IA
    write(*,*) "    JA = ", JA
    write(*,*) "    KA = ", KA
    write(*,*) "    NSTEP = ", NSTEP
    write(*,*) ""
    write(*,*) ">>  namelist: PARAM_BOUNDARY_EXP"
    write(*,*) "    lat0 [deg] = ", lat0
    write(*,*) "    lon0 [deg] = ", lon0
    write(*,*) "    dx [m]  = ", dx
    write(*,*) "    dy [m]  = ", dy
    write(*,*) "    height_bottom [m] = ", height_bottom
    write(*,*) "    height_top [m] = ", height_top
    write(*,*) "    testcase = ", trim(testcase)
    write(*,*) "    basename = ", trim(basename)
    write(*,*) "    file time interval [sec] = ", dt

    return
  end subroutine namelist_check

  !-----------------------------------------------------------------------------
  !> Make Time Stamp
  subroutine time_stamp
    implicit none
    integer :: l

    TIM_SEC(1) = 0.0_SP
    do l=2, NSTEP
       TIM_SEC(l) = TIM_SEC(l-1) + dt
    enddo

    return
  end subroutine time_stamp

  !-----------------------------------------------------------------------------
  !> Printing Error Message and STOP the NETCDF sequences.
  subroutine handle_err(status)
    ! Printing Error Message and STOP sequences. 

    use netcdf       ! [external]
    implicit none
    integer :: status
    write(*,*) nf90_strerror(status)
    stop
    return
  end subroutine handle_err

  !-----------------------------------------------------------------------------
  !> Variable ID number check
  subroutine vid_err( &
     nid,    &
     check   )
    implicit none

    integer, intent(in) :: check, nid
    ! -----
    if( check .ne. vid(nid) )then
       write(*, '(A)') "||-- ERROR: not match vid between file and program"
       write(*, *) "check", check, "vid", vid(nid)
       stop
    endif

  end subroutine vid_err

end program init_boundary_experiment
