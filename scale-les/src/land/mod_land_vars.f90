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
  public :: LAND_vars_external_in

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,  public :: LAND_sw_restart

  ! prognostic variables
  real(RP), public, allocatable :: TG  (:,:,:) ! soil temperature [K]
  real(RP), public, allocatable :: STRG(:,:,:) ! soil water [m3/m3]

  real(RP), public, allocatable :: LAND_PROPERTY(:,:,:) ! land surface property

  integer,  public, parameter   :: LAND_PROPERTY_nmax = 8
  integer,  public, parameter   :: I_STRGMAX          = 1 ! maximum  soil water [m3/m3]
  integer,  public, parameter   :: I_STRGCRT          = 2 ! critical soil water [m3/m3]
  integer,  public, parameter   :: I_TCS              = 3 ! thermal conductivity for soil [W/m/K]
  integer,  public, parameter   :: I_HCS              = 4 ! heat capacity        for soil [J/K]
  integer,  public, parameter   :: I_DFW              = 5 ! diffusive coefficient of soil water [m2/s]
  integer,  public, parameter   :: I_Z0M              = 6 ! roughness length for momemtum [m]
  integer,  public, parameter   :: I_Z0H              = 7 ! roughness length for heat     [m]
  integer,  public, parameter   :: I_Z0E              = 8 ! roughness length for moisture [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LAND_param_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: LAND_RESTART_OUTPUT       = .false.        !< output restart file?
  character(len=H_LONG),  private :: LAND_RESTART_IN_BASENAME  = ''             !< basename of the restart file
  character(len=H_LONG),  private :: LAND_RESTART_OUT_BASENAME = ''             !< basename of the output file
  character(len=H_MID),   private :: LAND_RESTART_OUT_TITLE    = 'LAND restart' !< title    of the output file
  character(len=H_MID),   private :: LAND_RESTART_OUT_DTYPE    = 'DEFAULT'      !< REAL4 or REAL8

  logical,                private :: LAND_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX   = 2     !< number of the variables
  integer,                private, parameter :: I_TG   = 1
  integer,                private, parameter :: I_STRG = 2

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'TG',   &
                  'STRG'  /
  data VAR_DESC / 'soil temperature', &
                  'soil water'        /
  data VAR_UNIT / 'K',     &
                  'm3/m3'  /

  integer,  private, parameter :: LAND_NUM_IDX = 2 ! # of land indices

  real(RP), private            :: LAND_PROPERTY_table(LAND_NUM_IDX,LAND_PROPERTY_nmax)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_landuse, only: &
       LANDUSE_index_PFT
    implicit none

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
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[LAND] / Origin[SCALE-LES]'

    allocate( TG  (LKS:LKE,IA,JA) )
    allocate( STRG(LKS:LKE,IA,JA) )

    TG  (:,:,:) = UNDEF
    STRG(:,:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** List of prognostic variables (LAND) ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A16,A,A32,3(A))') &
               '***       |','         VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A16,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(VAR_NAME(ip)),'|', VAR_DESC(ip),'[', VAR_UNIT(ip),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( LAND_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : ', trim(LAND_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       LAND_RESTART_OUTPUT             &
         .AND. LAND_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : ', trim(LAND_RESTART_OUT_BASENAME)
       LAND_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       LAND_RESTART_OUTPUT = .false.
       LAND_sw_restart = .false.
    endif



    call LAND_param_read

    allocate( LAND_PROPERTY(IA,JA,LAND_PROPERTY_nmax) )

    ! tentative, mosaic is off
    do v = 1, LAND_PROPERTY_nmax
    do j = JS, JE
    do i = IS, IE
       LAND_PROPERTY(i,j,v) = LAND_PROPERTY_table( LANDUSE_index_PFT(i,j,1),v )
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
  !> HALO Communication
  subroutine LAND_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: k

    real(RP) :: tmp1(IA,JA,KA)
    real(RP) :: tmp2(IA,JA,KA)
    !---------------------------------------------------------------------------

    do k = LKS, LKE
      tmp1(:,:,k) = TG  (k,:,:)
      tmp2(:,:,k) = STRG(k,:,:)
    end do

    do k = LKS, LKE
      call COMM_vars8( tmp1(:,:,k), k     )
      call COMM_vars8( tmp2(:,:,k), k+LKE )
    end do

    do k = LKS, LKE
      call COMM_wait ( tmp1(:,:,k), k     )
      call COMM_wait ( tmp2(:,:,k), k+LKE )
    end do

    do k = LKS, LKE
      TG  (k,:,:) = tmp1(:,:,k)
      STRG(k,:,:) = tmp2(:,:,k)
    end do

    return
  end subroutine LAND_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read land restart
  subroutine LAND_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (LAND) ***'

    if ( LAND_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(LAND_RESTART_IN_BASENAME)

       call FILEIO_read( TG(:,:,:),                                       & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'TG',   'Land', step=1 ) ! [IN]
       call FILEIO_read( STRG(:,:,:),                                     & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'STRG', 'Land', step=1 ) ! [IN]

       call LAND_vars_fillhalo

       call LAND_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for land is not specified.'
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
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (LAND) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call LAND_vars_total

       call FILEIO_write( TG(:,:,:),   basename,                                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TG),   VAR_DESC(I_TG),   VAR_UNIT(I_TG),   'Land', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( STRG(:,:,:), basename,                                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_STRG), VAR_DESC(I_STRG), VAR_UNIT(I_STRG), 'Land', LAND_RESTART_OUT_DTYPE  ) ! [IN]

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
    endif

    call HIST_in( TG  (:,:,:), 'TG',   VAR_DESC(I_TG),   VAR_UNIT(I_TG),   TIME_DTSEC_LAND, zdim='land' )
    call HIST_in( STRG(:,:,:), 'STRG', VAR_DESC(I_STRG), VAR_UNIT(I_STRG), TIME_DTSEC_LAND, zdim='land' )

    return
  end subroutine LAND_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_vars_total
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total
    integer  :: k
    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       do k = LKS, LKE
          call STAT_total( total, TG  (k,:,:), VAR_NAME(I_TG)   )
          call STAT_total( total, STRG(k,:,:), VAR_NAME(I_STRG) )
       enddo

    endif

    return
  end subroutine LAND_vars_total

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine LAND_vars_external_in( &
      tg_in,     & ! (in)
      strg_in    ) ! (in)
    implicit none

    real(RP), intent(in) :: tg_in(:,:,:)
    real(RP), intent(in) :: strg_in(:,:,:)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input (land) ***'

    TG(:,:,:)   = tg_in(:,:,:)
    STRG(:,:,:) = strg_in(:,:,:)

    call LAND_vars_fillhalo

    call LAND_vars_total

    return
  end subroutine LAND_vars_external_in

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
