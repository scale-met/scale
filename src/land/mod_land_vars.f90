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
  character(len=IO_SYSCHR),  public, save :: LAND_TYPE_PHY = 'OFF' !< Land physics type
  logical,                   public, save :: LAND_sw_phy           !< do land physics update?
  logical,                   public, save :: LAND_sw_restart       !< output restart?

  ! land-atmosphere flux
  real(RP), public, save :: SFLX_GH   (IA,JA) ! ground heat flux (upward positive) [W/m2]
  real(RP), public, save :: SFLX_PREC (IA,JA) ! precipitation flux [kg/m2/s]
  real(RP), public, save :: SFLX_QVLnd(IA,JA) ! moisture flux [kg/m2/s]

  ! prognostic variables
  real(RP), public, save :: TG   (IA,JA)      ! soil temperature [K]
  real(RP), public, save :: QvEfc(IA,JA)      ! efficiency of evaporation [0-1]
  real(RP), public, save :: ROFF (IA,JA)      ! run-off water [kg/m2]
  real(RP), public, save :: STRG (IA,JA)      ! water storage [kg/m2]

  integer,  public, save :: I_TG    = 1
  integer,  public, save :: I_QvEfc = 2
  integer,  public, save :: I_ROFF  = 3
  integer,  public, save :: I_STRG  = 4

  character(len=File_HSHORT), public, save :: LP_NAME(4) !< name  of the land variables
  character(len=File_HMID),   public, save :: LP_DESC(4) !< desc. of the land variables
  character(len=File_HSHORT), public, save :: LP_UNIT(4) !< unit  of the land variables

  data LP_NAME / 'TG',    &
                 'QvEfc', &
                 'ROFF',  &
                 'STRG'   /
  data LP_DESC / 'soil temperature',          &
                 'efficiency of evaporation', &
                 'run-off water',             &
                 'water storage'              /
  data LP_UNIT / 'K',     &
                 '0-1',   &
                 'kg/m2', &
                 'kg/m2'  /

  integer,  public, parameter :: LAND_PROPERTY_nmax = 10
  integer,  public, parameter :: I_STRGMAX =  1             ! maximum  water storage [kg/m2]
  integer,  public, parameter :: I_STRGCRT =  2             ! critical water storage [kg/m2]
  integer,  public, parameter :: I_EMIT    =  3             ! surface emissivity in long-wave  radiation [0-1]
  integer,  public, parameter :: I_ALB     =  4             ! surface albedo     in short-wave radiation [0-1]
  integer,  public, parameter :: I_TCS     =  5             ! thermal conductivity for soil [W/m/K]
  integer,  public, parameter :: I_HCS     =  6             ! heat capacity        for soil [J/K]
  integer,  public, parameter :: I_DZg     =  7             ! soil depth [m]
  integer,  public, parameter :: I_Z0M     =  8             ! roughness length for momemtum [m]
  integer,  public, parameter :: I_Z0H     =  9             ! roughness length for heat     [m]
  integer,  public, parameter :: I_Z0E     = 10             ! roughness length for moisture [m]

  integer,  public,      save :: LAND_Type    (IA,JA)                    ! land type index
  real(RP), public,      save :: LAND_PROPERTY(IA,JA,LAND_PROPERTY_nmax) ! land surface property

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LAND_param_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=IO_FILECHR), private, save :: LAND_BOUNDARY_IN_BASENAME = '' !< basename of the boundary file

  logical,                   private, save :: LAND_RESTART_OUTPUT       = .false.             !< output restart file?
  character(len=IO_FILECHR), private, save :: LAND_RESTART_IN_BASENAME  = 'restart_in'        !< basename of the restart file
  character(len=IO_FILECHR), private, save :: LAND_RESTART_OUT_BASENAME = 'restart_out'       !< basename of the output file
  character(len=IO_SYSCHR),  private, save :: LAND_RESTART_OUT_TITLE    = 'SCALE3 LAND VARS.' !< title    of the output file
  character(len=IO_SYSCHR),  private, save :: LAND_RESTART_OUT_DTYPE    = 'DEFAULT'           !< REAL4 or REAL8

  logical,                   private, save :: LAND_VARS_CHECKRANGE      = .false.

  real(RP), private, save :: LAND_PROPERTY_table(LAND_NUM_IDX,LAND_PROPERTY_nmax)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_vars_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[LAND VARS]/Categ[LAND]'

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
    if( IO_L ) write(IO_FID_LOG,*) '*** [LAND ] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, 7
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(LP_NAME(ip)),'|', LP_DESC(ip),'[', LP_UNIT(ip),']'
    enddo

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
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO
    call COMM_vars8( TG   (:,:), 1 )
    call COMM_vars8( QvEfc(:,:), 2 )
    call COMM_vars8( ROFF (:,:), 3 )
    call COMM_vars8( STRG (:,:), 4 )

    call COMM_wait ( TG   (:,:), 1 )
    call COMM_wait ( QvEfc(:,:), 2 )
    call COMM_wait ( ROFF (:,:), 3 )
    call COMM_wait ( STRG (:,:), 4 )

    return
  end subroutine LAND_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read land restart
  subroutine LAND_vars_restart_read
    use mod_fileio, only: &
       FILEIO_read
    use mod_const, only: &
       CONST_UNDEF
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, v
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (land) ***'

    call TIME_rapstart('FILE I NetCDF')

    if ( LAND_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( TG(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'TG',    'XY', step=1 ) ! [IN]
       call FILEIO_read( QvEfc(:,:),                                     & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'QvEfc', 'XY', step=1 ) ! [IN]
       call FILEIO_read( ROFF(:,:),                                      & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'ROFF',  'XY', step=1 ) ! [IN]
       call FILEIO_read( STRG(:,:),                                      & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'STRG',  'XY', step=1 ) ! [IN]

       call LAND_vars_fillhalo

       call LAND_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for land is not specified.'
       TG   (:,:) = 300.0_RP
       QvEfc(:,:) =   1.0_RP
       ROFF (:,:) =   0.0_RP
       STRG (:,:) = 100.0_RP
    endif

    LAND_PROPERTY(:,:,:) = CONST_UNDEF

    if ( LAND_BOUNDARY_IN_BASENAME /= '' ) then

       LAND_Type(:,:)       = 2
       LAND_Type(:,JS:JS+1) = 1
       LAND_Type(:,JE-1:JE) = 1
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

    call TIME_rapend  ('FILE I NetCDF')

    return
  end subroutine LAND_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write land restart
  subroutine LAND_vars_restart_write
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use mod_fileio, only: &
       FILEIO_write
    implicit none

    character(len=IO_FILECHR) :: basename

    integer :: n
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

    if ( LAND_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (land) ***'

       basename = ''
       write(basename(1:15), '(F15.3)') NOWSEC
       do n = 1, 15
          if ( basename(n:n) == ' ' ) basename(n:n) = '0'
       enddo
       basename = trim(LAND_RESTART_OUT_BASENAME) // '_' // trim(basename)

       call FILEIO_write( TG(:,:),    basename,                                       LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(I_TG),    LP_DESC(I_TG),    LP_UNIT(I_TG),    'XY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( QvEfc(:,:), basename,                                       LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(I_QvEfc), LP_DESC(I_QvEfc), LP_UNIT(I_QvEfc), 'XY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ROFF(:,:),  basename,                                       LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(I_ROFF),  LP_DESC(I_ROFF),  LP_UNIT(I_ROFF),  'XY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( STRG(:,:),  basename,                                       LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(I_STRG),  LP_DESC(I_STRG),  LP_UNIT(I_STRG),  'XY', LAND_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    call TIME_rapend  ('FILE O NetCDF')

    call LAND_vars_total

    return
  end subroutine LAND_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for land variables
  subroutine LAND_vars_history
    use mod_misc, only: &
       MISC_valcheck
    use mod_time, only: &
       TIME_DTSEC_LAND
    use mod_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( LAND_VARS_CHECKRANGE ) then
       call MISC_valcheck( TG   (:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_TG)    )
       call MISC_valcheck( QvEfc(:,:), 0.0_RP,    2.0_RP, LP_NAME(I_QvEfc) )
       call MISC_valcheck( ROFF (:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_ROFF)  )
       call MISC_valcheck( STRG (:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_STRG)  )
    endif

    call HIST_in( TG   (:,:),   'L_TG',    LP_DESC(I_TG),    LP_UNIT(I_TG),    TIME_DTSEC_LAND )
    call HIST_in( QvEfc(:,:),   'L_QvEfc', LP_DESC(I_QvEfc), LP_UNIT(I_QvEfc), TIME_DTSEC_LAND )
    call HIST_in( ROFF (:,:),   'L_ROFF',  LP_DESC(I_ROFF),  LP_UNIT(I_ROFF),  TIME_DTSEC_LAND )
    call HIST_in( STRG  (:,:),  'L_STRG',  LP_DESC(I_STRG),  LP_UNIT(I_STRG),  TIME_DTSEC_LAND )

    return
  end subroutine LAND_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_vars_total
    use mod_comm, only: &
       COMM_total_doreport, &
       COMM_total
    implicit none

    real(RP) :: total ! dummy
    !---------------------------------------------------------------------------

    if ( COMM_total_doreport ) then

    endif

    return
  end subroutine LAND_vars_total

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_param_read
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only: &
       CONST_UNDEF
    implicit none

    integer                  :: index
    character(len=IO_SYSCHR) :: description
    real(RP)                 :: STRGMAX
    real(RP)                 :: STRGCRT
    real(RP)                 :: EMIT
    real(RP)                 :: ALB
    real(RP)                 :: TCS
    real(RP)                 :: HCS
    real(RP)                 :: DZg
    real(RP)                 :: Z0M
    real(RP)                 :: Z0H
    real(RP)                 :: Z0E

    NAMELIST / PARAM_LAND_DATA / &
       index,       &
       description, &
       STRGMAX,     &
       STRGCRT,     &
       EMIT,        &
       ALB,         &
       TCS,         &
       HCS,         &
       DZg,         &
       Z0M,         &
       Z0H,         &
       Z0E

    integer :: n, v
    integer :: ierr
    !---------------------------------------------------------------------------

    LAND_PROPERTY_table(:,:) = CONST_UNDEF

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [LAND ] vegetation parameters'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,11(1x,A))') &
               '***        ',  'description ', &
                               'Max Stg.', &
                               'CRT Stg.', &
                               'Emissiv.', &
                               '  Albedo', &
                               'T condu.', &
                               'H capac.', &
                               '   Depth', &
                               '   Z0(m)', &
                               '   Z0(h)', &
                               '   Z0(e)'

    !--- read namelist
    rewind(IO_FID_CONF)
    do n = 1, LAND_NUM_IDX
       read(IO_FID_CONF,nml=PARAM_LAND_DATA,iostat=ierr)

       if ( ierr < 0 ) then !--- no more data
          exit
       elseif( ierr > 0 ) then !--- fatal error
          write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND. Check!'
          call PRC_MPIstop
       endif

       LAND_PROPERTY_table(index,I_STRGMAX) = STRGMAX
       LAND_PROPERTY_table(index,I_STRGCRT) = STRGCRT
       LAND_PROPERTY_table(index,I_EMIT   ) = EMIT
       LAND_PROPERTY_table(index,I_ALB    ) = ALB
       LAND_PROPERTY_table(index,I_TCS    ) = TCS
       LAND_PROPERTY_table(index,I_HCS    ) = HCS
       LAND_PROPERTY_table(index,I_DZg    ) = DZg
       LAND_PROPERTY_table(index,I_Z0M    ) = Z0M
       LAND_PROPERTY_table(index,I_Z0H    ) = Z0H
       LAND_PROPERTY_table(index,I_Z0E    ) = Z0E

       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,1x,A12,5(1x,F8.2),1x,F8.0,4(1x,F8.2))') &
                  '*** IDX=', index, trim(description), &
                                     STRGMAX, &
                                     STRGCRT, &
                                     EMIT,    &
                                     ALB,     &
                                     TCS,     &
                                     HCS,     &
                                     DZg,     &
                                     Z0M,     &
                                     Z0H,     &
                                     Z0E
    enddo

    return
  end subroutine LAND_param_read

end module mod_land_vars
