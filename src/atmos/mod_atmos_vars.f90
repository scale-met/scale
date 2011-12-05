!-------------------------------------------------------------------------------
!> module Atmospheric Variables
!!
!! @par Description
!!          Container for atmospheric variables
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-LES ver.2
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !

  use mod_stdio, only: &
     IO_SYSCHR, &
     IO_FILECHR
  use mod_fileio_h, only: &
     FIO_HSHORT, &
     FIO_HMID,   &
     FIO_REAL8

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !

  public :: ATMOS_vars_setup
  public :: ATMOS_vars_restart_read
  public :: ATMOS_vars_restart_write
  public :: ATMOS_vars_putDMP
  public :: ATMOS_vars_get
  public :: ATMOS_DMP2PVT

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer,                   public,              save :: A_VA      ! Number of Tracers + 5
  integer,                   public,              save :: A_QA      ! Number of Tracers
  character(len=FIO_HSHORT), public, allocatable, save :: A_NAME(:)
  character(len=FIO_HMID),   public, allocatable, save :: A_DESC(:)
  character(len=FIO_HSHORT), public, allocatable, save :: A_UNIT(:)

  character(len=IO_SYSCHR),  public, save :: ATMOS_TYPE_DYN    = 'NONE'
  character(len=IO_SYSCHR),  public, save :: ATMOS_TYPE_PHY_TB = 'NONE'
  character(len=IO_SYSCHR),  public, save :: ATMOS_TYPE_PHY_MP = 'NONE'
  character(len=IO_SYSCHR),  public, save :: ATMOS_TYPE_PHY_RD = 'NONE'

  logical,                   public, save :: ATMOS_sw_dyn
  logical,                   public, save :: ATMOS_sw_phy_tb
  logical,                   public, save :: ATMOS_sw_phy_mp
  logical,                   public, save :: ATMOS_sw_phy_rd

  real(8), public, allocatable, save :: atmos_var(:,:,:,:)      !> prognostics container (with HALO)
  real(8), public, allocatable, save :: atmos_diagvar(:,:,:,:)  !> diagnostics container (with HALO)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  character(len=IO_FILECHR), private, save :: ATMOS_RESTART_IN_BASENAME      = 'restart_in'
  character(len=IO_FILECHR), private, save :: ATMOS_RESTART_OUT_BASENAME     = 'restart_out'
  logical,                   private, save :: ATMOS_RESTART_IN_ALLOWMISSINGQ = .false.

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CONST_UNDEF8
    use mod_grid, only: &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KA => GRID_KA
    implicit none

    NAMELIST / PARAM_ATMOS / &
       ATMOS_TYPE_DYN,    &
       ATMOS_TYPE_PHY_TB, &
       ATMOS_TYPE_PHY_MP, &
       ATMOS_TYPE_PHY_RD

    integer                   :: ATMOS_QTRC_NMAX    = 0
    character(len=FIO_HSHORT) :: ATMOS_QTRC_VARNAME = ''
    character(len=IO_FILECHR) :: ATMOS_QTRC_VARDESC = ''
    character(len=IO_FILECHR) :: ATMOS_QTRC_VARUNIT = ''

    NAMELIST / PARAM_ATMOS_VARS / &
       ATMOS_QTRC_NMAX,            &
       ATMOS_QTRC_VARNAME,         &
       ATMOS_QTRC_VARDESC,         &
       ATMOS_QTRC_VARUNIT,         &
       ATMOS_RESTART_IN_BASENAME,  &
       ATMOS_RESTART_OUT_BASENAME, &
       ATMOS_RESTART_IN_ALLOWMISSINGQ

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Variables]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS] selected components'

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamics...'
    if ( ATMOS_TYPE_DYN == 'fent_fct' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Dynamical core   : Full-explicit No-terrain'
       if( IO_L ) write(IO_FID_LOG,*) '*** Tracer advection : FCT limitter'
       ATMOS_sw_dyn = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Dynamical core   : NONE'
       if( IO_L ) write(IO_FID_LOG,*) '*** Tracer advection : NONE'
       ATMOS_sw_dyn = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics...'
    if ( ATMOS_TYPE_PHY_TB == 'smagorinsky' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Sub-grid Turbulence : Smagorinsky'
       ATMOS_sw_phy_tb = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Sub-grid Turbulence : NONE'
       ATMOS_sw_phy_tb = .false.
    endif
    if ( ATMOS_TYPE_PHY_MP == 'NDW6' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Cloud Microphysics  : NDW6'
       ATMOS_sw_phy_mp = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Cloud Microphysics  : NONE'
       ATMOS_sw_phy_mp = .false.
    endif
    if ( ATMOS_TYPE_PHY_RD == 'mstrnX' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Radiative transfer  : mstrnX'
       ATMOS_sw_phy_rd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Radiative transfer  : NONE'
       ATMOS_sw_phy_rd = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[ATMOS VARS]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_VARS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_VARS)

    if ( ATMOS_QTRC_NMAX >= 1 ) then 

       A_VA = 5 + ATMOS_QTRC_NMAX
       allocate( A_NAME(A_VA) )
       allocate( A_DESC(A_VA) )
       allocate( A_UNIT(A_VA) )

       A_NAME( 1) = 'DENS'
       A_NAME( 2) = 'MOMX'
       A_NAME( 3) = 'MOMY'
       A_NAME( 4) = 'MOMZ'
       A_NAME( 5) = 'LWPT'

       A_DESC( 1) = 'density'
       A_DESC( 2) = 'momentum (x)'
       A_DESC( 3) = 'momentum (y)'
       A_DESC( 4) = 'momentum (z)'
       A_DESC( 5) = 'liquid water pot. temp.'

       A_UNIT( 1) = 'kg/m3'
       A_UNIT( 2) = 'kg/m2/s'
       A_UNIT( 3) = 'kg/m2/s'
       A_UNIT( 4) = 'kg/m2/s'
       A_UNIT( 5) = 'K'

       A_QA = ATMOS_QTRC_NMAX
       A_NAME(6:5+A_QA) = ATMOS_QTRC_VARNAME(1:ATMOS_QTRC_NMAX)
       A_DESC(6:5+A_QA) = ATMOS_QTRC_VARDESC(1:ATMOS_QTRC_NMAX)
       A_UNIT(6:5+A_QA) = ATMOS_QTRC_VARUNIT(1:ATMOS_QTRC_NMAX)

       if ( ATMOS_TYPE_PHY_MP == 'NDW6' ) then
          if ( ATMOS_QTRC_NMAX < 11 ) then
             write(*,*) 'xxx The number of tracers is not enough!', ATMOS_QTRC_NMAX
             call PRC_MPIstop
          endif

          A_NAME( 6) = 'QV'
          A_NAME( 7) = 'QC'
          A_NAME( 8) = 'QR'
          A_NAME( 9) = 'QI'
          A_NAME(10) = 'QS'
          A_NAME(11) = 'QG'
          A_NAME(12) = 'NC'
          A_NAME(13) = 'NR'
          A_NAME(14) = 'NI'
          A_NAME(15) = 'NS'
          A_NAME(16) = 'NG'

          A_DESC( 6) = 'Water Vapor mixing ratio'
          A_DESC( 7) = 'Cloud Water mixing ratio'
          A_DESC( 8) = 'Rain Water mixing ratio'
          A_DESC( 9) = 'Cloud Ice mixing ratio'
          A_DESC(10) = 'Snow mixing ratio'
          A_DESC(11) = 'Graupel mixing ratio'
          A_DESC(12) = 'Cloud Water Number Density'
          A_DESC(13) = 'Rain Water Number Density'
          A_DESC(14) = 'Cloud Ice Number Density'
          A_DESC(15) = 'Snow Number Density'
          A_DESC(16) = 'Graupel Number Density'

          A_UNIT( 6:11) = 'kg/kg'
          A_UNIT(12:16) = '1/m3'
       endif
    else
       write(*,*) 'xxx The number of tracers is not enough!', ATMOS_QTRC_NMAX
       A_QA = 0
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,*) &
    '***                 : VARNAME         , ', &
    'DESCRIPTION                                                     [UNIT            ]'
    do iv = 1, A_VA
       if( IO_L ) write(IO_FID_LOG,*) '*** NO.',iv,": ",A_NAME(iv),", ", A_DESC(iv),"[", A_UNIT(iv),"]"
    enddo

    allocate( atmos_var    (IA,JA,KA,A_VA) ); atmos_var(:,:,:,:)     = CONST_UNDEF8
    allocate( atmos_diagvar(IA,JA,KA,5)    ); atmos_diagvar(:,:,:,:) = CONST_UNDEF8
    
    return
  end subroutine ATMOS_vars_setup

  !-----------------------------------------------------------------------------
  !> Read restart of atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_restart_read
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
    use mod_comm, only: &
       COMM_vars, &
       COMM_stats
    use mod_fileio, only: &
       FIO_input
    implicit none

    real(8), allocatable :: restart_atmos(:,:,:) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=8)          :: lname

    integer :: iv, i, j
    !---------------------------------------------------------------------------

    allocate( restart_atmos(IMAX,JMAX,KMAX) )

    bname = ATMOS_RESTART_IN_BASENAME
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    do iv = 1, A_VA
       call FIO_input( restart_atmos(:,:,:), bname, A_NAME(iv), lname, 1, KMAX, 1 )

       atmos_var(IS:IE,JS:JE,KS:KE,iv) = restart_atmos(1:IMAX,1:JMAX,1:KMAX)
    enddo

    deallocate( restart_atmos )

    ! fill IHALO & JHALO
    call COMM_vars( atmos_var(:,:,:,:) )

    ! fill KHALO
    do iv = 1, A_VA
    do j  = 1, JA
    do i  = 1, IA
       atmos_var(i,j,   1:KS-1,iv) = atmos_var(i,j,KS,iv)
       atmos_var(i,j,KE+1:KA,  iv) = atmos_var(i,j,KE,iv)
    enddo
    enddo
    enddo

    call COMM_stats( atmos_var(:,:,:,:), A_NAME(:) )

    ! atmos_var -> atmos_diagvar
    call ATMOS_DMP2PVT

    ! fill IHALO & JHALO
    call COMM_vars( atmos_diagvar(:,:,:,:) )

    return
  end subroutine ATMOS_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Read restart of atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_restart_write
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use mod_grid, only : &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       KMAX => GRID_KMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE
    use mod_comm, only: &
       COMM_vars, &
       COMM_stats
    use mod_fileio, only: &
       FIO_output
    implicit none

    real(8), allocatable :: restart_atmos(:,:,:) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=FIO_HMID)   :: desc
    character(len=8)          :: lname

    integer :: iv
    !---------------------------------------------------------------------------

    allocate( restart_atmos(IMAX,JMAX,KMAX) )

    call COMM_stats( atmos_var(:,:,:,:), A_NAME(:) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (atmos) ***'

    write(bname,'(A,A,F15.3)') trim(ATMOS_RESTART_OUT_BASENAME), '_', NOWSEC
    desc  = 'SCALE3 PROGNOSTIC VARS.'
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    do iv = 1, A_VA
       restart_atmos(1:IMAX,1:JMAX,1:KMAX) = atmos_var(IS:IE,JS:JE,KS:KE,iv)

       call FIO_output( restart_atmos(:,:,:), bname, desc, '',       &
                        A_NAME(iv), A_DESC(iv), '', A_UNIT(iv),      &
                        FIO_REAL8, lname, 1, KMAX, 1, NOWSEC, NOWSEC )
    enddo

    deallocate( restart_atmos )

    return
  end subroutine ATMOS_vars_restart_write

  !-----------------------------------------------------------------------------
  !> Put and Communicate prognostic variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_putDMP( &
       dens, &
       momx, &
       momy, &
       momz, &
       lwpt, &
       qtrc  )
    use mod_grid, only: &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA
    use mod_comm, only: &
       COMM_vars
    implicit none

    real(8), intent(in) :: dens(IA,JA,KA)
    real(8), intent(in) :: momx(IA,JA,KA)
    real(8), intent(in) :: momy(IA,JA,KA)
    real(8), intent(in) :: momz(IA,JA,KA)
    real(8), intent(in) :: lwpt(IA,JA,KA)

    real(8), intent(in) :: qtrc(IA,JA,KA,A_QA)

    integer :: iq
    !---------------------------------------------------------------------------

    atmos_var(:,:,:,1) = dens(:,:,:)
    atmos_var(:,:,:,2) = momx(:,:,:)
    atmos_var(:,:,:,3) = momy(:,:,:)
    atmos_var(:,:,:,4) = momz(:,:,:)
    atmos_var(:,:,:,5) = lwpt(:,:,:)

    if ( A_QA > 0 ) then
       do iq = 1, A_QA
          atmos_var(:,:,:,5+iq) = qtrc(:,:,:,iq)
       enddo
    endif

    ! fill IHALO & JHALO
    call COMM_vars( atmos_var(:,:,:,:) )

    ! atmos_var -> atmos_diagvar
    call ATMOS_DMP2PVT

    ! fill IHALO & JHALO
    call COMM_vars( atmos_diagvar(:,:,:,:) )

    return
  end subroutine ATMOS_vars_putDMP

  !-----------------------------------------------------------------------------
  !> Get prognostic variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_get( &
       dens, &
       momx, &
       momy, &
       momz, &
       lwpt, &
       qtrc, &
       pres, &
       velx, &
       vely, &
       velz, &
       temp  )
    use mod_grid, only: &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KA => GRID_KA
    implicit none

    real(8), intent(out) :: dens(IA,JA,KA)
    real(8), intent(out) :: momx(IA,JA,KA)
    real(8), intent(out) :: momy(IA,JA,KA)
    real(8), intent(out) :: momz(IA,JA,KA)
    real(8), intent(out) :: lwpt(IA,JA,KA)

    real(8), intent(out) :: qtrc(IA,JA,KA,A_QA)

    real(8), intent(out) :: pres(IA,JA,KA)
    real(8), intent(out) :: velx(IA,JA,KA)
    real(8), intent(out) :: vely(IA,JA,KA)
    real(8), intent(out) :: velz(IA,JA,KA)
    real(8), intent(out) :: temp(IA,JA,KA)

    integer :: iq
    !---------------------------------------------------------------------------

    dens(:,:,:) = atmos_var(:,:,:,1)
    momx(:,:,:) = atmos_var(:,:,:,2)
    momy(:,:,:) = atmos_var(:,:,:,3)
    momz(:,:,:) = atmos_var(:,:,:,4)
    lwpt(:,:,:) = atmos_var(:,:,:,5)

    if ( A_QA > 0 ) then
       do iq = 1, A_QA
          qtrc(:,:,:,iq) = atmos_var(:,:,:,5+iq)
       enddo
    endif

    pres(:,:,:) = atmos_diagvar(:,:,:,1)
    velx(:,:,:) = atmos_diagvar(:,:,:,2)
    vely(:,:,:) = atmos_diagvar(:,:,:,3)
    velz(:,:,:) = atmos_diagvar(:,:,:,4)
    temp(:,:,:) = atmos_diagvar(:,:,:,5)

    return
  end subroutine ATMOS_vars_get

  !-----------------------------------------------------------------------------
  !> Get prognostic variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DMP2PVT
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_const, only : &
       Rair   => CONST_Rair,   &
       CPair  => CONST_CPair,  &
       RovCP  => CONST_RovCP,  &
       CPovCV => CONST_CPovCV, &
       LH0    => CONST_LH0,    &
       Pstd   => CONST_Pstd
    use mod_grid, only: &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       WS   => GRID_WS,   &
       WE   => GRID_WE
    implicit none

    real(8) :: pres(IA,JA,KA)
    real(8) :: velx(IA,JA,KA)
    real(8) :: vely(IA,JA,KA)
    real(8) :: velz(IA,JA,KA)
    real(8) :: temp(IA,JA,KA)

    real(8) :: fp, dfdp, dp
    real(8) :: eps = 1.D-10
    integer, parameter :: itmax = 20 ! max itelation cycle for pressure assumption

    integer :: i, j, k, it
    !---------------------------------------------------------------------------

    ! momentum -> velocity
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       velx(i,j,k) = 2.D0 * atmos_var(i,j,k,2) / ( atmos_var(i+1,j,  k,1)+atmos_var(i,j,k,1) )
       vely(i,j,k) = 2.D0 * atmos_var(i,j,k,3) / ( atmos_var(i,  j+1,k,1)+atmos_var(i,j,k,1) )
    enddo
    enddo
    enddo

    do k = WS+1, WE-1
    do j = JS,   JE
    do i = IS,   IE
       velz(i,j,k) = 2.D0 * atmos_var(i,j,k,4) / ( atmos_var(i,j,k+1,1)+atmos_var(i,j,k,1) )
    enddo
    enddo
    enddo
    velz(:,:,WS) = 0.D0 ! bottom boundary
    velz(:,:,WE) = 0.D0 ! top    boundary

    ! diagnose pressure, temperature
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       pres(i,j,k) = Pstd * ( atmos_var(i,j,k,1) * atmos_var(i,j,k,5) * Rair / Pstd )**CPovCV ! first guess

       do it = 1, itmax
          fp   = atmos_var(i,j,k,5) * ( pres(i,j,k)/Pstd )**RovCP          &
               + LH0 / CPair * ( atmos_var(i,j,k,7) + atmos_var(i,j,k,8) ) &
               - pres(i,j,k) / ( Rair *  atmos_var(i,j,k,1) )
          dfdp = RovCP / Pstd * atmos_var(i,j,k,5) * ( pres(i,j,k)/Pstd )**(RovCP-1) &
               - 1.D0 / ( Rair * atmos_var(i,j,k,1) )
          dp   = fp / dfdp

          pres(i,j,k) = pres(i,j,k) - dp
          if ( abs(dp) < eps ) exit
       enddo

       temp(i,j,k) = pres(i,j,k) / ( atmos_var(i,j,k,1) * Rair )
    enddo
    enddo
    enddo

    atmos_diagvar(IS:IE,JS:JE,KS:KE,1) = pres(IS:IE,JS:JE,KS:KE)
    atmos_diagvar(IS:IE,JS:JE,KS:KE,2) = velx(IS:IE,JS:JE,KS:KE)
    atmos_diagvar(IS:IE,JS:JE,KS:KE,3) = vely(IS:IE,JS:JE,KS:KE)
    atmos_diagvar(IS:IE,JS:JE,WS:WE,4) = velz(IS:IE,JS:JE,WS:WE)
    atmos_diagvar(IS:IE,JS:JE,KS:KE,5) = temp(IS:IE,JS:JE,KS:KE)

    return
  end subroutine ATMOS_DMP2PVT

end module mod_atmos_vars
