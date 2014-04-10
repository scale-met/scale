!-------------------------------------------------------------------------------
!> module LAND Variables
!!
!! @par Description
!!          Container for land variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_land_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_vars_setup
  public :: LAND_vars_fillhalo
  public :: LAND_vars_restart_read
  public :: LAND_vars_restart_write
  public :: LAND_vars_history
  public :: LAND_vars_total

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public, save :: LAND_TYPE_PHY = 'OFF' !< Land physics type
  logical,                public, save :: LAND_sw_phy           !< do land physics update?
  logical,                public, save :: LAND_sw_restart       !< output restart?

  ! prognostic variables
  real(RP), public, save, allocatable :: TG  (:,:,:) ! soil temperature [K]
  real(RP), public, save, allocatable :: STRG(:,:,:) ! water storage [kg/m2]
  real(RP), public, save, allocatable :: ROFF(:,:)   ! run-off water [kg/m2]
  real(RP), public, save, allocatable :: QVEF(:,:)   ! efficiency of evaporation [0-1]

  integer,  public, parameter :: PV_NUM = 4

  integer,  public, parameter :: I_TG   = 1
  integer,  public, parameter :: I_STRG = 2
  integer,  public, parameter :: I_ROFF = 3
  integer,  public, parameter :: I_QVEF = 4

  character(len=H_SHORT), public, save :: PV_NAME(PV_NUM) !< name  of the land variables
  character(len=H_MID),   public, save :: PV_DESC(PV_NUM) !< desc. of the land variables
  character(len=H_SHORT), public, save :: PV_UNIT(PV_NUM) !< unit  of the land variables

  data PV_NAME / 'TG',   &
                 'STRG', &
                 'ROFF', &
                 'QVEF'  /
  data PV_DESC / 'soil temperature',          &
                 'water storage',             &
                 'run-off water',             &
                 'efficiency of evaporation'  /
  data PV_UNIT / 'K',     &
                 'kg/m2', &
                 'kg/m2', &
                 '0-1'    /

  integer,  public, parameter :: LAND_PROPERTY_nmax = 8
  integer,  public, parameter :: I_STRGMAX =  1  ! maximum  water storage [kg/m2]
  integer,  public, parameter :: I_STRGCRT =  2  ! critical water storage [kg/m2]
  integer,  public, parameter :: I_TCS     =  3  ! thermal conductivity for soil [W/m/K]
  integer,  public, parameter :: I_HCS     =  4  ! heat capacity        for soil [J/K]
  integer,  public, parameter :: I_DFW     =  5  ! diffusive coefficient of soil water [m2/s]
  integer,  public, parameter :: I_Z0M     =  6  ! roughness length for momemtum [m]
  integer,  public, parameter :: I_Z0H     =  7  ! roughness length for heat     [m]
  integer,  public, parameter :: I_Z0E     =  8  ! roughness length for moisture [m]

  integer,  public, save, allocatable :: LAND_Type    (:,:)   ! land type index
  real(RP), public, save, allocatable :: LAND_PROPERTY(:,:,:) ! land surface property

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LAND_param_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private, save :: LAND_BOUNDARY_IN_BASENAME = ''                     !< basename of the boundary file

  logical,               private, save :: LAND_RESTART_OUTPUT       = .false.                !< output restart file?
  character(len=H_LONG), private, save :: LAND_RESTART_IN_BASENAME  = ''                     !< basename of the restart file
  character(len=H_LONG), private, save :: LAND_RESTART_OUT_BASENAME = ''                     !< basename of the output file
  character(len=H_MID),  private, save :: LAND_RESTART_OUT_TITLE    = 'SCALE-LES LAND VARS.' !< title    of the output file
  character(len=H_MID),  private, save :: LAND_RESTART_OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8

  logical,               private, save :: LAND_VARS_CHECKRANGE      = .false.

  integer,  private, parameter :: LAND_NUM_IDX = 2 ! # of land indices

  real(RP), private, save      :: LAND_PROPERTY_table(LAND_NUM_IDX,LAND_PROPERTY_nmax)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_LAND / &
       LAND_TYPE_PHY

    NAMELIST / PARAM_LAND_VARS /  &
       LAND_RESTART_IN_BASENAME,  &
       LAND_BOUNDARY_IN_BASENAME, &
       LAND_RESTART_OUTPUT,       &
       LAND_RESTART_OUT_BASENAME, &
       LAND_RESTART_OUT_TITLE,    &
       LAND_RESTART_OUT_DTYPE,    &
       LAND_VARS_CHECKRANGE

    integer :: ierr
    integer :: ip
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[LAND VARS]/Categ[LAND]'

    allocate( TG  (LKS:LKE,IA,JA) )
    allocate( STRG(LKS:LKE,IA,JA) )
    allocate( ROFF(IA,JA)     )
    allocate( QVEF(IA,JA)     )

    allocate( LAND_Type    (IA,JA) )
    allocate( LAND_PROPERTY(IA,JA,LAND_PROPERTY_nmax) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_LAND)

    if( IO_L ) write(IO_FID_LOG,*) '*** [LAND] selected components'

    if ( LAND_TYPE_PHY /= 'OFF' .AND. LAND_TYPE_PHY /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Land physics : ON'
       LAND_sw_phy = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Land physics : OFF'
       LAND_sw_phy = .false.
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_VARS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_LAND_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [LAND] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, PV_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(PV_NAME(ip)),'|', PV_DESC(ip),'[', PV_UNIT(ip),']'
    enddo

    ! restart switch
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( LAND_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Land restart output : YES'
       LAND_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Land restart output : NO'
       LAND_sw_restart = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    call LAND_param_read

    return
  end subroutine LAND_vars_setup

  !-----------------------------------------------------------------------------
  !> fill HALO region of land variables
  subroutine LAND_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO
    do k = LKS, LKE
      call COMM_vars8( TG  (k,:,:), 1 )
      call COMM_vars8( STRG(k,:,:), 2 )

      call COMM_wait ( TG  (k,:,:), 1 )
      call COMM_wait ( STRG(k,:,:), 2 )
    end do

    ! 2D variable
    call COMM_vars8( ROFF(:,:), 3 )
    call COMM_vars8( QVEF(:,:), 4 )

    call COMM_wait ( ROFF(:,:), 3 )
    call COMM_wait ( QVEF(:,:), 4 )
    return
  end subroutine LAND_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read land restart
  subroutine LAND_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_const, only: &
       CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, v
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (land) ***'

    call PROF_rapstart('FILE I NetCDF')

    if ( LAND_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( TG(:,:,:),                                       & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'TG',   'Land', step=1 ) ! [IN]
       call FILEIO_read( STRG(:,:,:),                                     & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'STRG', 'Land', step=1 ) ! [IN]
       call FILEIO_read( ROFF(:,:),                                       & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'ROFF', 'XY',   step=1 ) ! [IN]
       call FILEIO_read( QVEF(:,:),                                       & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'QVEF', 'XY',   step=1 ) ! [IN]

       call LAND_vars_fillhalo

       call LAND_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for land is not specified.'

       TG  (:,:,:) = CONST_UNDEF
       STRG(:,:,:) = CONST_UNDEF
       ROFF(:,:)   = CONST_UNDEF
       QVEF(:,:)   = CONST_UNDEF
    endif

    LAND_PROPERTY(:,:,:) = CONST_UNDEF

    ! tentative
    ! this will be merged into the landuse module
    if ( LAND_BOUNDARY_IN_BASENAME /= '' ) then

       LAND_Type(:,:) = 1
!       call FILEIO_read( LAND_Type(:,:),                                      & ! [OUT]
!                         LAND_BOUNDARY_IN_BASENAME, 'LAND_Type', 'XY', step=1 ) ! [IN]

       do v = 1,  LAND_PROPERTY_nmax
       do j = JS, JE
       do i = IS, IE
          LAND_PROPERTY(i,j,v) = LAND_PROPERTY_table( LAND_Type(i,j),v )
       enddo
       enddo
       enddo

       do v = 1,  LAND_PROPERTY_nmax
          call COMM_vars8( LAND_PROPERTY(:,:,v), v )
       enddo
       do v = 1,  LAND_PROPERTY_nmax
          call COMM_wait ( LAND_PROPERTY(:,:,v), v )
       enddo
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** boundary file for land is not specified.'
    endif

    call PROF_rapend  ('FILE I NetCDF')

    return
  end subroutine LAND_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write land restart
  subroutine LAND_vars_restart_write
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=H_LONG) :: basename

    integer :: n
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF')

    if ( LAND_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (land) ***'

       basename = ''
       write(basename(1:15), '(F15.3)') NOWSEC
       do n = 1, 15
          if ( basename(n:n) == ' ' ) basename(n:n) = '0'
       enddo
       basename = trim(LAND_RESTART_OUT_BASENAME) // '_' // trim(basename)

       call FILEIO_write( TG(:,:,:),   basename,                                     LAND_RESTART_OUT_TITLE,        & ! [IN]
                          PV_NAME(I_TG),   PV_DESC(I_TG),   PV_UNIT(I_TG),   'Land', LAND_RESTART_OUT_DTYPE, .true. ) ! [IN]
       call FILEIO_write( STRG(:,:,:), basename,                                     LAND_RESTART_OUT_TITLE,        & ! [IN]
                          PV_NAME(I_STRG), PV_DESC(I_STRG), PV_UNIT(I_STRG), 'Land', LAND_RESTART_OUT_DTYPE, .true. ) ! [IN]
       call FILEIO_write( ROFF(:,:),   basename,                                     LAND_RESTART_OUT_TITLE,        & ! [IN]
                          PV_NAME(I_ROFF), PV_DESC(I_ROFF), PV_UNIT(I_ROFF), 'XY',   LAND_RESTART_OUT_DTYPE, .true. ) ! [IN]
       call FILEIO_write( QVEF(:,:),   basename,                                     LAND_RESTART_OUT_TITLE,        & ! [IN]
                          PV_NAME(I_QVEF), PV_DESC(I_QVEF), PV_UNIT(I_QVEF), 'XY',   LAND_RESTART_OUT_DTYPE, .true. ) ! [IN]

    endif

    call PROF_rapend  ('FILE O NetCDF')

    call LAND_vars_total

    return
  end subroutine LAND_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for land variables
  subroutine LAND_vars_history
    use scale_time, only: &
       TIME_DTSEC_LAND
    use scale_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( LAND_VARS_CHECKRANGE ) then
       call VALCHECK( TG  (:,:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_TG)  , __FILE__, __LINE__ )
       call VALCHECK( STRG(:,:,:), 0.0_RP, 1000.0_RP, PV_NAME(I_STRG), __FILE__, __LINE__ )
       call VALCHECK( ROFF(:,:),   0.0_RP, 1000.0_RP, PV_NAME(I_ROFF), __FILE__, __LINE__ )
       call VALCHECK( QVEF(:,:),   0.0_RP,    2.0_RP, PV_NAME(I_QVEF), __FILE__, __LINE__ )
    endif

    call HIST_in( TG  (:,:,:), 'TG',   PV_DESC(I_TG),   PV_UNIT(I_TG),   TIME_DTSEC_LAND, zdim ='land' )
    call HIST_in( STRG(:,:,:), 'STRG', PV_DESC(I_STRG), PV_UNIT(I_STRG), TIME_DTSEC_LAND, zdim ='land' )
    call HIST_in( ROFF(:,:),   'ROFF', PV_DESC(I_ROFF), PV_UNIT(I_ROFF), TIME_DTSEC_LAND )
    call HIST_in( QVEF(:,:),   'QVEF', PV_DESC(I_QVEF), PV_UNIT(I_QVEF), TIME_DTSEC_LAND )

    return
  end subroutine LAND_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_vars_total
    use scale_stats, only: &
       STAT_checktotal, &
       STAT_total
    implicit none

    !real(RP) :: total
    !---------------------------------------------------------------------------

    if ( STAT_checktotal ) then

!       call STAT_total( total, TG(:,:,:),   PV_NAME(I_TG)   )
!       call STAT_total( total, STRG(:,:,:), PV_NAME(I_STRG) )
!       call STAT_total( total, ROFF(:,:),   PV_NAME(I_ROFF) )
!       call STAT_total( total, QVEF(:,:),   PV_NAME(I_QVEF) )

    endif

    return
  end subroutine LAND_vars_total

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_param_read
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer                :: index
    character(len=H_SHORT) :: description
    real(RP)               :: STRGMAX
    real(RP)               :: STRGCRT
    real(RP)               :: TCS
    real(RP)               :: HCS
    real(RP)               :: DFW
    real(RP)               :: Z0M
    real(RP)               :: Z0H
    real(RP)               :: Z0E

    NAMELIST / PARAM_LAND_DATA / &
       index,       &
       description, &
       STRGMAX,     &
       STRGCRT,     &
       TCS,         &
       HCS,         &
       DFW,         &
       Z0M,         &
       Z0H,         &
       Z0E

    integer :: n
    integer :: ierr
    !---------------------------------------------------------------------------

    LAND_PROPERTY_table(:,:) = CONST_UNDEF

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [LAND ] vegetation parameters'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,11(1x,A))') &
               '***        ',  'description ', &
                               'Max Stg.', &
                               'CRT Stg.', &
                               'T condu.', &
                               'H capac.', &
                               'DFC Wat.', &
                               '   Z0(m)', &
                               '   Z0(h)', &
                               '   Z0(e)'

    !--- read namelist
    rewind(IO_FID_CONF)
    do n = 1, LAND_NUM_IDX
       ! undefined roughness length
       Z0H = -1.0_RP
       Z0E = -1.0_RP

       read(IO_FID_CONF,nml=PARAM_LAND_DATA,iostat=ierr)

       if ( ierr < 0 ) then !--- no more data
          exit
       elseif( ierr > 0 ) then !--- fatal error
          write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND. Check!'
          call PRC_MPIstop
       endif

       if( Z0H < 0.0_RP ) then
         Z0H = Z0M / 7.4_RP ! defined by Garratt and Francey (1978)
       endif
       if( Z0E < 0.0_RP ) then
         Z0E = Z0M / 7.4_RP ! defined by Garratt and Francey (1978)
       endif

       LAND_PROPERTY_table(index,I_STRGMAX) = STRGMAX
       LAND_PROPERTY_table(index,I_STRGCRT) = STRGCRT
       LAND_PROPERTY_table(index,I_TCS    ) = TCS
       LAND_PROPERTY_table(index,I_HCS    ) = HCS
       LAND_PROPERTY_table(index,I_DFW    ) = DFW
       LAND_PROPERTY_table(index,I_Z0M    ) = Z0M
       LAND_PROPERTY_table(index,I_Z0H    ) = Z0H
       LAND_PROPERTY_table(index,I_Z0E    ) = Z0E

       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,1x,A12,5(1x,F8.2),1x,F8.0,4(1x,F8.2))') &
                  '*** IDX=', index, trim(description), &
                                     STRGMAX, &
                                     STRGCRT, &
                                     TCS,     &
                                     HCS,     &
                                     DFW,     &
                                     Z0M,     &
                                     Z0H,     &
                                     Z0E
    enddo

    return
  end subroutine LAND_param_read

end module mod_land_vars
