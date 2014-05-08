!-------------------------------------------------------------------------------
!> module LAND Variables
!!
!! @par Description
!!          Container for land variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
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
  logical,                public, save :: LAND_do          = .true. ! main switch for the model

  character(len=H_SHORT), public, save :: LAND_TYPE        = 'OFF' !< Land physics type

  logical,                public, save :: LAND_sw_phy
  logical,                public, save :: LAND_sw_restart

  ! prognostic variables
  real(RP), public, save, allocatable :: TG  (:,:,:) ! soil temperature [K]
  real(RP), public, save, allocatable :: STRG(:,:,:) ! water storage [kg/m2]
  real(RP), public, save, allocatable :: ROFF(:,:)   ! run-off water [kg/m2]
  real(RP), public, save, allocatable :: QVEF(:,:)   ! efficiency of evaporation [0-1]

  integer,  public, parameter :: LAND_PROPERTY_nmax = 8
  integer,  public, parameter :: I_STRGMAX          = 1 ! maximum  water storage [kg/m2]
  integer,  public, parameter :: I_STRGCRT          = 2 ! critical water storage [kg/m2]
  integer,  public, parameter :: I_TCS              = 3 ! thermal conductivity for soil [W/m/K]
  integer,  public, parameter :: I_HCS              = 4 ! heat capacity        for soil [J/K]
  integer,  public, parameter :: I_DFW              = 5 ! diffusive coefficient of soil water [m2/s]
  integer,  public, parameter :: I_Z0M              = 6 ! roughness length for momemtum [m]
  integer,  public, parameter :: I_Z0H              = 7 ! roughness length for heat     [m]
  integer,  public, parameter :: I_Z0E              = 8 ! roughness length for moisture [m]

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
  logical,               private, save :: LAND_RESTART_OUTPUT       = .false.        !< output restart file?
  character(len=H_LONG), private, save :: LAND_RESTART_IN_BASENAME  = ''             !< basename of the restart file
  character(len=H_LONG), private, save :: LAND_RESTART_OUT_BASENAME = ''             !< basename of the output file
  character(len=H_MID),  private, save :: LAND_RESTART_OUT_TITLE    = 'LAND restart' !< title    of the output file
  character(len=H_MID),  private, save :: LAND_RESTART_OUT_DTYPE    = 'DEFAULT'      !< REAL4 or REAL8

  logical,               private, save :: LAND_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX   = 4     !< number of the variables
  integer,                private, parameter :: I_TG   = 1
  integer,                private, parameter :: I_STRG = 2
  integer,                private, parameter :: I_ROFF = 3
  integer,                private, parameter :: I_QVEF = 4

  character(len=H_SHORT), private, save      :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private, save      :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private, save      :: VAR_UNIT(VMAX) !< unit  of the variables
 !< unit  of the land variables

  data VAR_NAME / 'TG',   &
                  'STRG', &
                  'ROFF', &
                  'QVEF'  /
  data VAR_DESC / 'soil temperature',          &
                  'water storage',             &
                  'run-off water',             &
                  'efficiency of evaporation'  /
  data VAR_UNIT / 'K',     &
                  'kg/m2', &
                  'kg/m2', &
                  '0-1'    /

  integer,  private, parameter :: LAND_NUM_IDX = 2 ! # of land indices

  real(RP), private, save      :: LAND_PROPERTY_table(LAND_NUM_IDX,LAND_PROPERTY_nmax)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_landuse, only: &
       LANDUSE_index_vegetation
    implicit none

    NAMELIST / PARAM_LAND / &
       LAND_do,  &
       LAND_TYPE

    NAMELIST / PARAM_LAND_VARS /  &
       LAND_RESTART_IN_BASENAME,  &
       LAND_RESTART_OUTPUT,       &
       LAND_RESTART_OUT_BASENAME, &
       LAND_RESTART_OUT_TITLE,    &
       LAND_RESTART_OUT_DTYPE,    &
       LAND_VARS_CHECKRANGE

    integer :: ierr
    integer :: i, j, v, ip
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[LAND VARS]/Categ[LAND]'
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

    if ( LAND_TYPE /= 'OFF' .AND. LAND_TYPE /= 'NONE' ) then
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
    do ip = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(VAR_NAME(ip)),'|', VAR_DESC(ip),'[', VAR_UNIT(ip),']'
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

    allocate( TG  (LKS:LKE,IA,JA) )
    allocate( STRG(LKS:LKE,IA,JA) )
    allocate( ROFF(        IA,JA) )
    allocate( QVEF(        IA,JA) )



    call LAND_param_read

    allocate( LAND_PROPERTY(IA,JA,LAND_PROPERTY_nmax) )
    LAND_PROPERTY(:,:,:) = CONST_UNDEF
    LANDUSE_index_vegetation(:,:,:) = 1 ! tentative!

    do v = 1, LAND_PROPERTY_nmax
    do j = JS, JE
    do i = IS, IE
       LAND_PROPERTY(i,j,v) = LAND_PROPERTY_table( LANDUSE_index_vegetation(i,j,1),v )
    enddo
    enddo
    enddo

    do v = 1, LAND_PROPERTY_nmax
       call COMM_vars8( LAND_PROPERTY(:,:,v), v )
    enddo
    do v = 1, LAND_PROPERTY_nmax
       call COMM_wait ( LAND_PROPERTY(:,:,v), v )
    enddo

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
    enddo

    call COMM_vars8( ROFF(:,:), 1 )
    call COMM_vars8( QVEF(:,:), 2 )

    call COMM_wait ( ROFF(:,:), 1 )
    call COMM_wait ( QVEF(:,:), 2 )

    return
  end subroutine LAND_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read land restart
  subroutine LAND_vars_restart_read
    use scale_const, only: &
       CONST_UNDEF
    use scale_fileio, only: &
       FILEIO_read
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (land) ***'

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

    return
  end subroutine LAND_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write land restart
  subroutine LAND_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( LAND_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(LAND_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (land) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** filename: ', trim(basename)

       call LAND_vars_total

       call FILEIO_write( TG(:,:,:),   basename,                                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TG),   VAR_DESC(I_TG),   VAR_UNIT(I_TG),   'Land', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( STRG(:,:,:), basename,                                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_STRG), VAR_DESC(I_STRG), VAR_UNIT(I_STRG), 'Land', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ROFF(:,:),   basename,                                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_ROFF), VAR_DESC(I_ROFF), VAR_UNIT(I_ROFF), 'XY',   LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( QVEF(:,:),   basename,                                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_QVEF), VAR_DESC(I_QVEF), VAR_UNIT(I_QVEF), 'XY',   LAND_RESTART_OUT_DTYPE  ) ! [IN]

    endif

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
       call VALCHECK( TG  (:,:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TG)  , __FILE__, __LINE__ )
       call VALCHECK( STRG(:,:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_STRG), __FILE__, __LINE__ )
       call VALCHECK( ROFF(:,:),   0.0_RP, 1000.0_RP, VAR_NAME(I_ROFF), __FILE__, __LINE__ )
       call VALCHECK( QVEF(:,:),   0.0_RP,    2.0_RP, VAR_NAME(I_QVEF), __FILE__, __LINE__ )
    endif

    call HIST_in( TG  (:,:,:), 'TG',   VAR_DESC(I_TG),   VAR_UNIT(I_TG),   TIME_DTSEC_LAND, zdim ='land' )
    call HIST_in( STRG(:,:,:), 'STRG', VAR_DESC(I_STRG), VAR_UNIT(I_STRG), TIME_DTSEC_LAND, zdim ='land' )
    call HIST_in( ROFF(:,:),   'ROFF', VAR_DESC(I_ROFF), VAR_UNIT(I_ROFF), TIME_DTSEC_LAND )
    call HIST_in( QVEF(:,:),   'QVEF', VAR_DESC(I_QVEF), VAR_UNIT(I_QVEF), TIME_DTSEC_LAND )

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

!       call STAT_total( total, TG  (:,:,:), VAR_NAME(I_TG)   )
!       call STAT_total( total, STRG(:,:,:), VAR_NAME(I_STRG) )
!       call STAT_total( total, ROFF(:,:),   VAR_NAME(I_ROFF) )
!       call STAT_total( total, QVEF(:,:),   VAR_NAME(I_QVEF) )

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
               '***        ',  ' description', &
                               ' Max Stg.', &
                               ' CRT Stg.', &
                               ' T condu.', &
                               ' H capac.', &
                               ' DFC Wat.', &
                               '    Z0(m)', &
                               '    Z0(h)', &
                               '    Z0(e)'

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

       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,1x,A12,8(1x,F9.2))') &
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
