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
  include 'inc_land.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_vars_setup

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

  integer,                    public, save :: I_SkinT = 1
  integer,                    public, save :: I_SkinW = 1
  integer,                    public, save :: I_SnowQ = 1
  integer,                    public, save :: I_SnowT = 1
  character(len=File_HSHORT), public, save :: LP_NAME(4) !< name  of the coupler variables
  character(len=File_HMID),   public, save :: LP_DESC(4) !< desc. of the coupler variables
  character(len=File_HSHORT), public, save :: LP_UNIT(4) !< unit  of the coupler variables

  data LP_NAME / 'SkinT', &
                 'SkinW', &
                 'SnowQ', &
                 'SnowT'  /
  data LP_DESC / 'ground skin temp.',  &
                 'ground skin water',  &
                 'ground snow amount', &
                 'ground snow temp.'   /
  data LP_UNIT / 'K',     &
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

    LST  (:,:) = 300.0_RP
    SST  (:,:) = 300.0_RP

    SkinT(:,:) = 300.0_RP
    SkinW(:,:) = 0.0_RP
    SnowQ(:,:) = 0.0_RP
    SnowT(:,:) = 0.0_RP

    return
  end subroutine CPL_vars_setup

end module mod_CPL_vars
