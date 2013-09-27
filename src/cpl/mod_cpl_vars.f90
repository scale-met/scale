!-------------------------------------------------------------------------------
!> module COUPLER Variables
!!
!! @par Description
!!          Container for coupler variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_cpl_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use gtool_file_h, only: &
     File_HSHORT, &
     File_HMID,   &
     File_HLONG
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR,  &
     IO_FILECHR
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
!  include 'inc_land.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_vars_setup
  public :: CPL_vars_fillhalo
  public :: CPL_vars_restart_read
  public :: CPL_vars_restart_write
  public :: CPL_vars_history
  public :: CPL_vars_total

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: LST  (IA,JA) ! Land Surface Temperature [K]
  real(RP), public, save :: SST  (IA,JA) ! Sea Surface Temperature [K]
  real(RP), public, save :: SkinT(IA,JA) ! Ground Skin Temperature [K]
  real(RP), public, save :: SkinW(IA,JA) ! Ground Skin Water       [kg/m2]
  real(RP), public, save :: SnowQ(IA,JA) ! Ground Snow amount      [kg/m2]
  real(RP), public, save :: SnowT(IA,JA) ! Ground Snow Temperature [K]

  integer,                    public, save :: I_LST   = 1
  integer,                    public, save :: I_SST   = 2
  integer,                    public, save :: I_SkinT = 3
  integer,                    public, save :: I_SkinW = 4
  integer,                    public, save :: I_SnowQ = 5
  integer,                    public, save :: I_SnowT = 6
  character(len=File_HSHORT), public, save :: LP_NAME(6) !< name  of the coupler variables
  character(len=File_HMID),   public, save :: LP_DESC(6) !< desc. of the coupler variables
  character(len=File_HSHORT), public, save :: LP_UNIT(6) !< unit  of the coupler variables

  data LP_NAME / 'LST',   &
                 'SST',   &
                 'SkinT', &
                 'SkinW', &
                 'SnowQ', &
                 'SnowT'  /
  data LP_DESC / 'land surface temp.', &
                 'sea surface temp.',  &
                 'ground skin temp.',  &
                 'ground skin water',  &
                 'ground snow amount', &
                 'ground snow temp.'   /
  data LP_UNIT / 'K',     &
                 'K',     &
                 'K',     &
                 'kg/m2', &
                 'kg/m2', &
                 'K'      /

  character(len=IO_SYSCHR),  public, save :: CPL_TYPE_AtmLnd = 'OFF' !< atmos-land coupler type
  logical,                   public, save :: CPL_sw_AtmLnd           !< do atmos-land coupler calculation?
  logical,                   public, save :: CPL_sw_AtmLnd_restart   !< output atmos-land coupler restart?

  character(len=IO_FILECHR), public, save :: CPL_RESTART_IN_BASENAME = '' !< basename of the input file

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                   private, save :: CPL_RESTART_OUTPUT       = .false.             !< output restart file?
  character(len=IO_FILECHR), private, save :: CPL_RESTART_OUT_BASENAME = 'restart_out'       !< basename of the output file
  character(len=IO_SYSCHR),  private, save :: CPL_RESTART_OUT_TITLE    = 'SCALE3 CPL VARS.'  !< title    of the output file
  character(len=IO_SYSCHR),  private, save :: CPL_RESTART_OUT_DTYPE    = 'DEFAULT'           !< REAL4 or REAL8

  logical,                   private, save :: CPL_VARS_CHECKRANGE      = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_vars_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_CPL / &
       CPL_TYPE_AtmLnd

    integer :: ierr
    integer :: ip
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CPL VARS]/Categ[CPL]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CPL)

    if( IO_L ) write(IO_FID_LOG,*) '*** [CPL] selected components'

    if ( CPL_TYPE_AtmLnd /= 'OFF' .AND. CPL_TYPE_AtmLnd /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Land Coupler : ON'
       CPL_sw_AtmLnd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Land Coupler : OFF'
       CPL_sw_AtmLnd = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CPL VARS]/Categ[CPL]'

    return
  end subroutine CPL_vars_setup

  !-----------------------------------------------------------------------------
  !> fill HALO region of coupler variables
  subroutine CPL_vars_fillhalo
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO
    call COMM_vars8( LST  (:,:), 1 )
    call COMM_vars8( SST  (:,:), 2 )
    call COMM_vars8( SkinT(:,:), 3 )
    call COMM_vars8( SkinW(:,:), 4 )
    call COMM_vars8( SnowQ(:,:), 5 )
    call COMM_vars8( SnowT(:,:), 6 )

    call COMM_wait ( LST  (:,:), 1 )
    call COMM_wait ( SST  (:,:), 2 )
    call COMM_wait ( SkinT(:,:), 3 )
    call COMM_wait ( SkinW(:,:), 4 )
    call COMM_wait ( SnowQ(:,:), 5 )
    call COMM_wait ( SnowT(:,:), 6 )

    return
  end subroutine CPL_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read coupler restart
  subroutine CPL_vars_restart_read
    use mod_fileio, only: &
       FILEIO_read
    use mod_const, only: &
       CONST_UNDEF
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (coupler) ***'

    call TIME_rapstart('FILE I NetCDF')

    if ( CPL_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( LST(:,:),                                      & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'LST',   'XY', step=1 ) ![IN]
       call FILEIO_read( SST(:,:),                                      & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SST',   'XY', step=1 ) ![IN]
       call FILEIO_read( SkinT(:,:),                                    & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SkinT', 'XY', step=1 ) ![IN]
       call FILEIO_read( SkinW(:,:),                                    & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SkinW', 'XY', step=1 ) ![IN]
       call FILEIO_read( SnowQ(:,:),                                    & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SnowQ', 'XY', step=1 ) ![IN]
       call FILEIO_read( SnowT(:,:),                                    & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SnowT', 'XY', step=1 ) ![IN]

       call CPL_vars_fillhalo

       call CPL_vars_total

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for coupler is not specified.'

!       LST  (:,:) = CONST_UNDEF
       LST  (:,:) = 300.0_RP
       SST  (:,:) = CONST_UNDEF
       SkinT(:,:) = CONST_UNDEF
       SkinW(:,:) = CONST_UNDEF
       SnowQ(:,:) = CONST_UNDEF
       SnowT(:,:) = CONST_UNDEF

    endif

    call TIME_rapend  ('FILE I NetCDF')

    return
  end subroutine CPL_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write coupler restart
  subroutine CPL_vars_restart_write
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use mod_fileio, only: &
       FILEIO_write
    implicit none

    character(len=IO_FILECHR) :: bname

    integer :: n
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

    if ( CPL_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (coupler) ***'

       bname = ''
       write(bname(1:15), '(F15.3)') NOWSEC
       do n = 1, 15
          if ( bname(n:n) == ' ' ) bname(n:n) = '0'
       enddo
       write(bname,'(A,A,A)') trim(CPL_RESTART_OUT_BASENAME), '_', trim(bname)

       call FILEIO_write( LST(:,:),   bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_LST),   LP_DESC(I_LST),   LP_UNIT(I_LST),   'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SST(:,:),   bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SST),   LP_DESC(I_SST),   LP_UNIT(I_SST),   'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SkinT(:,:), bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SkinT), LP_DESC(I_SkinT), LP_UNIT(I_SkinT), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SkinW(:,:), bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SkinW), LP_DESC(I_SkinW), LP_UNIT(I_SkinW), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SnowQ(:,:), bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SnowQ), LP_DESC(I_SnowQ), LP_UNIT(I_SnowQ), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SnowT(:,:), bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SnowT), LP_DESC(I_SnowT), LP_UNIT(I_SnowT), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    call TIME_rapend  ('FILE O NetCDF')

    call CPL_vars_total

    return
  end subroutine CPL_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for coupler variables
  subroutine CPL_vars_history
    use mod_misc, only: &
       MISC_valcheck
    use mod_time, only: &
       TIME_DTSEC_CPL_AtmLnd
    use mod_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( CPL_VARS_CHECKRANGE ) then
       call MISC_valcheck( LST  (:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_LST)  )
       call MISC_valcheck( SST  (:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SST)  )
       call MISC_valcheck( SkinT(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SkinT))
       call MISC_valcheck( SkinW(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SkinW))
       call MISC_valcheck( SnowQ(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SnowQ))
       call MISC_valcheck( SnowT(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SnowT))
    endif

    call HIST_in( LST  (:,:),   'L_LST',   LP_DESC(I_LST),   LP_UNIT(I_LST),   TIME_DTSEC_CPL_AtmLnd )
    call HIST_in( SST  (:,:),   'L_SST',   LP_DESC(I_SST),   LP_UNIT(I_SST),   TIME_DTSEC_CPL_AtmLnd ) !!! DTSEC will be modified
    call HIST_in( SkinT(:,:),   'L_SkinT', LP_DESC(I_SkinT), LP_UNIT(I_SkinT), TIME_DTSEC_CPL_AtmLnd ) !!! DTSEC will be modified
    call HIST_in( SkinW(:,:),   'L_SkinW', LP_DESC(I_SkinW), LP_UNIT(I_SkinW), TIME_DTSEC_CPL_AtmLnd ) !!! DTSEC will be modified
    call HIST_in( SnowQ(:,:),   'L_SnowQ', LP_DESC(I_SnowQ), LP_UNIT(I_SnowQ), TIME_DTSEC_CPL_AtmLnd ) !!! DTSEC will be modified
    call HIST_in( SnowT(:,:),   'L_SnowT', LP_DESC(I_SnowT), LP_UNIT(I_SnowT), TIME_DTSEC_CPL_AtmLnd ) !!! DTSEC will be modified

    return
  end subroutine CPL_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for coupler
  subroutine CPL_vars_total
    use mod_comm, only: &
       COMM_total_doreport, &
       COMM_total
    implicit none

    real(RP) :: total ! dummy
    !---------------------------------------------------------------------------

!    if ( COMM_total_doreport ) then
!
!       call COMM_total( total, LST(:,:),     LP_NAME(I_LST)   )
!       call COMM_total( total, SST(:,:),     LP_NAME(I_SST)   )
!       call COMM_total( total, SkinT(:,:),   LP_NAME(I_SkinT) )
!       call COMM_total( total, SkinW(:,:),   LP_NAME(I_SkinW) )
!       call COMM_total( total, SnowQ(:,:),   LP_NAME(I_SnowQ) )
!       call COMM_total( total, SnowT(:,:),   LP_NAME(I_SnowT) )
!
!    endif

    return
  end subroutine CPL_vars_total

end module mod_CPL_vars
