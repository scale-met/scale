!-------------------------------------------------------------------------------
!> module MONITOR
!!
!! @par Description
!!          Monitor output module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_monitor
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: MONITOR_setup
  public :: MONITOR_set_dim
  public :: MONITOR_reg
  public :: MONITOR_put
  public :: MONITOR_in
  public :: MONITOR_write
  public :: MONITOR_finalize

  interface MONITOR_in
     module procedure MONITOR_in_2D
     module procedure MONITOR_in_3D
  end interface MONITOR_in

  interface MONITOR_put
     module procedure MONITOR_put_2D
     module procedure MONITOR_put_3D
  end interface MONITOR_put

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: MONITOR_writeheader

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer                :: MONITOR_FID = -1                    !< fileID for monitor output file

  character(len=H_LONG)  :: MONITOR_OUT_BASENAME  = 'monitor' !< filename of monitor output
  logical                :: MONITOR_USEDEVIATION  = .true.    !< use deviation from first step?
  integer                :: MONITOR_STEP_INTERVAL = 1         !< step interval
  logical                :: MONITOR_GLOBAL_SUM    = .true.    !< global or local sum

  real(DP)               :: MONITOR_dt

  integer, parameter     :: MONITOR_req_max = 1000      !< number limit for item request
  integer                :: MONITOR_nreqs   = 0         !< number of requested item
  character(len=H_SHORT) :: MONITOR_reqs(MONITOR_req_max) !< name of requested monitor item

  type item
     character(len=H_SHORT) :: name     !< name
     character(len=H_MID)   :: desc     !< description
     character(len=H_SHORT) :: unit     !< unit
     real(DP)               :: var      !< value
     real(DP)               :: var0     !< value at first time
     logical                :: first    !< first time?
     logical                :: tendency !< integrate value?
     integer                :: dimid    !< dimension type
  end type item
  integer                 :: MONITOR_nitems = 0 !< number of item to output
  type(item), allocatable :: MONITOR_items(:)

  type dim_type
     character(len=H_SHORT) :: name
     integer                :: KA, KS, KE
     integer                :: IA, IS, IE
     integer                :: JA, JS, JE
     integer                :: dim_size
     real(RP), allocatable  :: area(:,:)
     real(RP)               :: total_area
     real(RP), allocatable  :: volume(:,:,:)
     real(RP)               :: total_volume
  end type dim_type
  integer, parameter :: MONITOR_dim_max = 30
  integer            :: MONITOR_ndims = 0
  type(dim_type)     :: MONITOR_dims(MONITOR_dim_max)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MONITOR_setup( dt )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(DP), intent(in) :: dt

    namelist / PARAM_MONITOR / &
       MONITOR_OUT_BASENAME, &
       MONITOR_USEDEVIATION, &
       MONITOR_GLOBAL_SUM,   &
       MONITOR_STEP_INTERVAL

    character(len=H_SHORT) :: NAME  !> name of monitor item

    namelist / MONITOR_ITEM / &
       NAME

    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MONITOR_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MONITOR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO('MONITOR_setup',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR('MONITOR_setup',*) 'Not appropriate names in namelist PARAM_MONITOR. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MONITOR)

    ! listup monitor request
    rewind(IO_FID_CONF)
    do n = 1, MONITOR_req_max
       read(IO_FID_CONF,nml=MONITOR_ITEM,iostat=ierr)
       if( ierr /= 0 ) exit
    enddo
    MONITOR_nreqs = n - 1

    if    ( MONITOR_nreqs > MONITOR_req_max ) then
       LOG_ERROR('MONITOR_setup',*) 'request of monitor file is exceed! n >', MONITOR_req_max
       call PRC_abort
    elseif( MONITOR_nreqs == 0 ) then
       LOG_INFO('MONITOR_setup',*) 'No monitor file specified.'
       return
    else
       LOG_INFO('MONITOR_setup',*) 'Number of requested monitor item : ', MONITOR_nreqs
       LOG_INFO('MONITOR_setup',*) 'Monitor output interval   [step] : ', MONITOR_STEP_INTERVAL
       LOG_INFO('MONITOR_setup',*) 'Use deviation from first step?   : ', MONITOR_USEDEVIATION
    endif

    allocate( MONITOR_items(MONITOR_nreqs) )

    rewind(IO_FID_CONF)
    do n = 1, MONITOR_nreqs
       ! set default
       NAME  = 'unknown'

       read(IO_FID_CONF,nml=MONITOR_ITEM,iostat=ierr)
       if( ierr /= 0 ) exit

       if ( IO_FID_NML /= IO_FID_LOG ) then
          LOG_NML(MONITOR_ITEM)
       end if

       MONITOR_reqs(n) = NAME
    enddo


    MONITOR_dt = dt

    return
  end subroutine MONITOR_setup

  !-----------------------------------------------------------------------------
  !> Set area and volume
  subroutine MONITOR_set_dim( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       dim_type, dim_size, &
       area, total_area, &
       volume, total_volume )
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    character(len=*), intent(in) :: dim_type
    integer,          intent(in) :: dim_size
    real(RP),         intent(in), optional :: area(IA,JA)
    real(RP),         intent(in), optional :: total_area
    real(RP),         intent(in), optional :: volume(KA,IA,JA)
    real(RP),         intent(in), optional :: total_volume

    integer :: n

    MONITOR_ndims = MONITOR_ndims + 1
    n = MONITOR_ndims

    MONITOR_dims(n)%name = dim_type
    MONITOR_dims(n)%dim_size = dim_size

    MONITOR_dims(n)%KA = KA
    MONITOR_dims(n)%KS = KS
    MONITOR_dims(n)%KE = KE
    MONITOR_dims(n)%IA = IA
    MONITOR_dims(n)%IS = IS
    MONITOR_dims(n)%IE = IE
    MONITOR_dims(n)%JA = JA
    MONITOR_dims(n)%JS = JS
    MONITOR_dims(n)%JE = JE

    if ( dim_size >= 2 ) then
       allocate( MONITOR_dims(n)%area(IA,JA) )
       MONITOR_dims(n)%area(:,:) = area(:,:)
       MONITOR_dims(n)%total_area = total_area
    end if

    if ( dim_size >= 3 ) then
       allocate( MONITOR_dims(n)%volume(KA,IA,JA) )
       MONITOR_dims(n)%volume(:,:,:) = volume(:,:,:)
       MONITOR_dims(n)%total_volume = total_volume
    end if

    return
  end subroutine MONITOR_set_dim

  !-----------------------------------------------------------------------------
  !> Search existing item, or matching check between requested and registered item
  subroutine MONITOR_reg( &
       name, desc, unit, &
       itemid,           &
       ndims, dim_type,  &
       is_tendency       )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in)  :: name !< name         of the item
    character(len=*), intent(in)  :: desc !< description  of the item
    character(len=*), intent(in)  :: unit !< unit         of the item

    integer, intent(out) :: itemid !< index number of the item

    integer,          intent(in), optional :: ndims       !< # of dimension
    character(len=*), intent(in), optional :: dim_type    !< dimension type
    logical,          intent(in), optional :: is_tendency !< need to integrate value?

    integer :: n, reqid, dimid
    !---------------------------------------------------------------------------

    !--- search existing item
    do itemid = 1, MONITOR_nitems
       if ( name == MONITOR_items(itemid)%name ) return ! match existing item
    enddo

    do reqid = 1, MONITOR_nreqs
       if ( name == MONITOR_reqs(reqid) ) then

          call PROF_rapstart('Monit', 2)

          MONITOR_nitems = MONITOR_nitems + 1
          itemid = MONITOR_nitems

          ! new file registration
          MONITOR_items(itemid)%name = name
          MONITOR_items(itemid)%desc = desc
          MONITOR_items(itemid)%unit = unit

          dimid = -1
          if ( present(dim_type) ) then
             do n = 1, MONITOR_ndims
                if ( MONITOR_dims(n)%name == dim_type ) then
                   dimid = n
                   exit
                end if
             end do
             if ( dimid < 0 ) then
                LOG_ERROR('MONITOR_reg',*) 'dim_type (', trim(dim_type), ') must be registerd by MONITOR_set_dim'
                call PRC_abort
             end if
          else if ( present(ndims) ) then
             do n = 1, MONITOR_ndims
                if ( MONITOR_dims(n)%dim_size == ndims ) then
                   dimid = n
                   exit
                end if
             end do
             if ( dimid == -1 ) then
                LOG_ERROR('MONITOR_reg','(a,i1,a)') 'dim_type of ', ndims, 'D must be registerd with MONITOR_set_dim'
                call PRC_abort
             end if
          else
             ! ndims = 3 is assumed as default
             do n = 1, MONITOR_ndims
                if ( MONITOR_dims(n)%dim_size == 3 ) then
                   dimid = n
                   exit
                end if
             end do
             if ( dimid == -1 ) then
                LOG_ERROR('MONITOR_reg',*) 'dim_type or ndims must be specified'
                call PRC_abort
             end if
          end if

          MONITOR_items(itemid)%dimid = dimid

          MONITOR_items(itemid)%var   = 0.0_DP
          MONITOR_items(itemid)%var0  = 0.0_DP
          MONITOR_items(itemid)%first = .true.
          if ( present(is_tendency) ) then
             MONITOR_items(itemid)%tendency = is_tendency
          else
             MONITOR_items(itemid)%tendency = .false.
          end if

          LOG_NEWLINE
          LOG_INFO('MONOTOR_reg','(A,I3)') ' Item registration No.= ', itemid
          LOG_INFO_CONT(*) 'Name            : ', trim(MONITOR_items(itemid)%name)
          LOG_INFO_CONT(*) 'Description     : ', trim(MONITOR_items(itemid)%desc)
          LOG_INFO_CONT(*) 'Unit            : ', trim(MONITOR_items(itemid)%unit)
          LOG_INFO_CONT(*) 'Dimension type  : ', trim(MONITOR_dims(MONITOR_items(itemid)%dimid)%name)
          LOG_INFO_CONT(*) 'Integ. with dt? : ', MONITOR_items(itemid)%tendency

          call PROF_rapend('Monit', 2)

          return
       end if
    end do

    itemid = -1 ! not found

    return
  end subroutine MONITOR_reg

  !-----------------------------------------------------------------------------
  !> Put total value to the monitor buffer
  subroutine MONITOR_put_2D( &
      itemid, var )
    use scale_statistics, only: &
       STATISTICS_total
    implicit none
    integer,  intent(in) :: itemid     !< index number of the item
    real(RP), intent(in) :: var(:,:)   !< value

    integer :: dimid
    real(DP) :: total
    !---------------------------------------------------------------------------

    if( itemid <= 0 ) return

    call PROF_rapstart('Monit', 2)

    dimid = MONITOR_items(itemid)%dimid

    call STATISTICS_total( MONITOR_dims(dimid)%IA, MONITOR_dims(dimid)%IS, MONITOR_dims(dimid)%IE, &
                           MONITOR_dims(dimid)%JA, MONITOR_dims(dimid)%JS, MONITOR_dims(dimid)%JE, &
                           var(:,:), MONITOR_items(itemid)%name,                          & ! (in)
                           MONITOR_dims(dimid)%area(:,:), MONITOR_dims(dimid)%total_area, & ! (in)
                           log_suppress = .true., global = MONITOR_GLOBAL_SUM,            & ! (in)
                           sum = total                                                    ) ! (out)

    if ( MONITOR_items(itemid)%tendency ) then
       if ( MONITOR_items(itemid)%first ) then
          MONITOR_items(itemid)%var = 0.0_RP
          MONITOR_items(itemid)%first = .false.
       else
          MONITOR_items(itemid)%var = MONITOR_items(itemid)%var + total * MONITOR_dt ! integrate by last put
       endif
    else
       if ( MONITOR_USEDEVIATION ) then
          if ( MONITOR_items(itemid)%first ) then
             MONITOR_items(itemid)%var  = 0.0_RP
             MONITOR_items(itemid)%var0 = total
             MONITOR_items(itemid)%first = .false.
          else
             MONITOR_items(itemid)%var = total - MONITOR_items(itemid)%var0 ! overwrite by last put
          endif
       else
          MONITOR_items(itemid)%var = total ! overwrite by last put
       endif
    endif

    call PROF_rapend('Monit', 2)

    return
  end subroutine MONITOR_put_2D

  !-----------------------------------------------------------------------------
  !> Put total value to the monitor buffer
  subroutine MONITOR_put_3D( &
      itemid, var )
    use scale_statistics, only: &
       STATISTICS_total
    implicit none

    integer,  intent(in) :: itemid     !< index number of the item
    real(RP), intent(in) :: var(:,:,:) !< value

    integer :: dimid

    real(DP) :: total
    !---------------------------------------------------------------------------

    if( itemid <= 0 ) return

    call PROF_rapstart('Monit', 2)

    dimid = MONITOR_items(itemid)%dimid


    call STATISTICS_total( MONITOR_dims(dimid)%KA, MONITOR_dims(dimid)%KS, MONITOR_dims(dimid)%KE, &
                           MONITOR_dims(dimid)%IA, MONITOR_dims(dimid)%IS, MONITOR_dims(dimid)%IE, &
                           MONITOR_dims(dimid)%JA, MONITOR_dims(dimid)%JS, MONITOR_dims(dimid)%JE, &
                           var(:,:,:), MONITOR_items(itemid)%name,                              & ! (in)
                           MONITOR_dims(dimid)%volume(:,:,:), MONITOR_dims(dimid)%total_volume, & ! (in)
                           log_suppress = .true., global = MONITOR_GLOBAL_SUM,                  & ! (in)
                           sum = total                                                          ) ! (out)

    if ( MONITOR_items(itemid)%tendency ) then
       if ( MONITOR_items(itemid)%first ) then
          MONITOR_items(itemid)%var   = total * MONITOR_dt ! first put
          MONITOR_items(itemid)%first = .false.
       else
          MONITOR_items(itemid)%var = MONITOR_items(itemid)%var + total * MONITOR_dt ! integrate by last put
       endif
    else
       if ( MONITOR_USEDEVIATION ) then
          if ( MONITOR_items(itemid)%first ) then
             MONITOR_items(itemid)%var   = 0.0_RP
             MONITOR_items(itemid)%var0  = total
             MONITOR_items(itemid)%first = .false.
          else
             MONITOR_items(itemid)%var = total - MONITOR_items(itemid)%var0 ! overwrite by last put
          endif
       else
          MONITOR_items(itemid)%var = total ! overwrite by last put
       endif
    endif

    call PROF_rapend('Monit', 2)

    return
  end subroutine MONITOR_put_3D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of MONITOR_reg+MONITOR_put
  subroutine MONITOR_in_2D( &
      var,              &
      name, desc, unit, &
      ndims, dim_type,  &
      is_tendency       )
    implicit none

    real(RP),         intent(in) :: var(:,:)   !< value
    character(len=*), intent(in) :: name       !< name
    character(len=*), intent(in) :: desc       !< description
    character(len=*), intent(in) :: unit       !< unit

    integer,          intent(in), optional :: ndims       !< # of dimension
    character(len=*), intent(in), optional :: dim_type    !< dimension type
    logical,          intent(in), optional :: is_tendency !< need to integrate values?

    integer :: itemid
    !---------------------------------------------------------------------------

    call MONITOR_reg( name, desc, unit,             & ! (in)
                    itemid,                         & ! (out)
                    ndims=ndims, dim_type=dim_type, & ! (in)
                    is_tendency=is_tendency         ) ! (in)
    call MONITOR_put( itemid, var(:,:) )

    return
  end subroutine MONITOR_in_2D

  !-----------------------------------------------------------------------------
  !> Wrapper routine of MONITOR_reg+MONITOR_put
  subroutine MONITOR_in_3D( &
      var,              &
      name, desc, unit, &
      ndims, dim_type,  &
      is_tendency       )
    implicit none

    real(RP),         intent(in) :: var(:,:,:) !< value
    character(len=*), intent(in) :: name       !< name
    character(len=*), intent(in) :: desc       !< description
    character(len=*), intent(in) :: unit       !< unit

    integer,          intent(in), optional :: ndims       !< # of dimension
    character(len=*), intent(in), optional :: dim_type    !< dimension type
    logical,          intent(in), optional :: is_tendency !< need to integrate values?

    integer :: itemid
    !---------------------------------------------------------------------------

    call MONITOR_reg( name, desc, unit,             & ! (in)
                    itemid,                         & ! (out)
                    ndims=ndims, dim_type=dim_type, & ! (in)
                    is_tendency=is_tendency         ) ! (in)
    call MONITOR_put( itemid, var(:,:,:) )

    return
  end subroutine MONITOR_in_3D

  !-----------------------------------------------------------------------------
  !> Flush monitor buffer to formatted file
  subroutine MONITOR_write( memo, nowstep )
    implicit none
    character(len=*), intent(in) :: memo !< note
    integer         , intent(in) :: nowstep

    logical, save :: firsttime = .true.

    integer :: n
    !---------------------------------------------------------------------------

    if( MONITOR_nitems == 0 ) return

    call PROF_rapstart('Monit', 2)

    if (firsttime) then
       firsttime = .false.
       call MONITOR_writeheader
    endif

    if ( MONITOR_FID > 0 ) then

       if ( mod(NOWSTEP-1,MONITOR_STEP_INTERVAL) == 0 ) then
          LOG_PROGRESS(*) 'output monitor'

          write(MONITOR_FID,'(A,i7,A,A4,A)',advance='no') 'STEP=',NOWSTEP,' (',memo,')'
          do n = 1, MONITOR_nitems
             write(MONITOR_FID,'(A,ES15.8)',advance='no') ' ', MONITOR_items(n)%var
          enddo
          write(MONITOR_FID,*)
       endif

    endif

    call PROF_rapend('Monit', 2)

    return
  end subroutine MONITOR_write

  !-----------------------------------------------------------------------------
  !> Open file and write header at the first time
  subroutine MONITOR_writeheader
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank, &
       PRC_IsMaster
    implicit none

    character(len=H_LONG) :: fname !< name of monitor file for each process

    logical :: MONITOR_L
    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO('MONITOR_writeheader',*) 'Output item list '
    LOG_INFO_CONT(*) 'Number of monitor item :', MONITOR_nreqs
    LOG_INFO_CONT('(1x,2A)') 'NAME                   :description                                    ', &
                             ':UNIT           :dimension_type'
    LOG_INFO_CONT('(1x,2A)') '=======================================================================', &
                             '==============================='
    do n = 1, MONITOR_nitems
       LOG_INFO_CONT('(1x,A24,A48,A16,A16)') MONITOR_items(n)%name, MONITOR_items(n)%desc, MONITOR_items(n)%unit, MONITOR_dims(MONITOR_items(n)%dimid)%name
    enddo
    LOG_INFO_CONT('(1x,2A)') '=======================================================================', &
                             '==============================='

    if ( PRC_IsMaster ) then ! master node
       MONITOR_L = .true.
    else
       MONITOR_L = IO_LOG_ALLNODE .and. ( .not. MONITOR_GLOBAL_SUM )
    endif

    if ( MONITOR_L ) then

       !--- Open logfile
       MONITOR_FID = IO_get_available_fid()
       if ( MONITOR_GLOBAL_SUM ) then
          call IO_get_fname(fname, MONITOR_OUT_BASENAME, rank=-1)
       else
          call IO_get_fname(fname, MONITOR_OUT_BASENAME, rank=PRC_myrank)
       end if
       open( unit   = MONITOR_FID,  &
             file   = trim(fname),  &
             form   = 'formatted',  &
             iostat = ierr          )
       if ( ierr /= 0 ) then
          LOG_ERROR('MONITOR_writeheader',*) 'File open error! :', trim(fname)
          call PRC_abort
       endif

       LOG_NEWLINE
       LOG_INFO('MONITOR_writeheader',*) 'Open ASCII file for monitor, name : ', trim(fname)

       write(MONITOR_FID,'(A)',advance='no') '                   '
       do n = 1, MONITOR_nitems
          write(MONITOR_FID,'(A16)',advance='no') MONITOR_items(n)%name
       enddo
       write(MONITOR_FID,*)

    endif

    return
  end subroutine MONITOR_writeheader

  !-----------------------------------------------------------------------------
  !> Close file
  subroutine MONITOR_finalize
    use scale_prc, only: &
       PRC_myrank
    implicit none

    character(len=H_LONG) :: fname !< name of monitor file for each process

    integer :: n
    !---------------------------------------------------------------------------

    if ( MONITOR_FID > 0 ) then
       LOG_NEWLINE
       LOG_INFO('MONITOR_finalize',*) 'Close monitor file'

       close(MONITOR_FID)
    endif

    do n = 1, MONITOR_ndims
       if ( MONITOR_dims(n)%dim_size >= 2 ) deallocate( MONITOR_dims(n)%area )
       if ( MONITOR_dims(n)%dim_size >= 3 ) deallocate( MONITOR_dims(n)%volume )
    end do
    MONITOR_ndims = 0

    if ( allocated(MONITOR_items) ) deallocate( MONITOR_items )
    MONITOR_nitems = 0

    return
  end subroutine MONITOR_finalize

end module scale_monitor
