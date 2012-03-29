!-------------------------------------------------------------------------------
!> module Monitor
!!
!! @par Description
!!          Monitor output module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-03-22 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_monitor
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_FILECHR
  use mod_fileio_h, only: &
     FIO_HSHORT, &
     FIO_HMID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: MONITOR_setup
  public :: MONITOR_reg
  public :: MONITOR_put
  public :: MONITOR_in

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_index.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  character(len=IO_FILECHR), private,      save :: MONITOR_OUT_BASENAME = 'monitor'
  integer,                   private,      save :: MONITOR_DTYPE

  integer,                   private, parameter :: MONITOR_req_limit = 1000 !> number limit for monitor item request
  character(len=FIO_HSHORT), private,      save :: MONITOR_req_item   (MONITOR_req_limit)
  real(8),                   private,      save :: MONITOR_req_tintsec(MONITOR_req_limit)
  logical,                   private,      save :: MONITOR_req_tavg   (MONITOR_req_limit)

  integer,                   private,              save :: MONITOR_req_nmax = 0 !> number of requested item
  character(len=FIO_HSHORT), private, allocatable, save :: MONITOR_item   (:)
  character(len=FIO_HMID),   private, allocatable, save :: MONITOR_desc   (:)
  character(len=FIO_HSHORT), private, allocatable, save :: MONITOR_unit   (:)
  character(len=FIO_HSHORT), private, allocatable, save :: MONITOR_ktype  (:)
  integer,                   private, allocatable, save :: MONITOR_kmax   (:)
  real(8),                   private, allocatable, save :: MONITOR_tintsec(:)
  logical,                   private, allocatable, save :: MONITOR_tavg   (:)

  real(8),                   private, allocatable, save :: MONITOR_varsum (:,:,:,:)
  integer,                   private, allocatable, save :: MONITOR_step   (:)
  real(8),                   private, allocatable, save :: MONITOR_tstrsec(:)
  real(8),                   private, allocatable, save :: MONITOR_tsumsec(:)

  integer,                   private,              save :: MONITOR_id_count = 1 !> number of registered item

  real(8), private, parameter :: eps = 1.D-10 !> epsilon for timesec

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  subroutine MONITOR_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_MONITOR / &
       MONITOR_OUT_BASENAME

    character(len=FIO_HSHORT) :: ITEM  !> name of monitor item

    NAMELIST / MONITITEM / &
       ITEM

    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[MONITOR]/Categ[IO]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MONITOR,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MONITOR. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MONITOR)

    ! listup monitor request
    rewind(IO_FID_CONF)
    do n = 1, MONITOR_req_limit
       read(IO_FID_CONF,nml=MONITITEM,iostat=ierr)
       if( ierr /= 0 ) exit
    enddo
    MONITOR_req_nmax = n - 1

    if    ( MONITOR_req_nmax > MONITOR_req_limit ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** request of monitor file is exceed! n >', MONITOR_req_limit
    elseif( MONITOR_req_nmax == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** No monitor file specified.'
       return
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Number of requested monitor item             : ', MONITOR_req_nmax
    endif

    allocate( MONITOR_item   (MONITOR_req_nmax) )
    allocate( MONITOR_desc   (MONITOR_req_nmax) )
    allocate( MONITOR_unit   (MONITOR_req_nmax) )
    allocate( MONITOR_ktype  (MONITOR_req_nmax) )
    allocate( MONITOR_kmax   (MONITOR_req_nmax) )
    allocate( MONITOR_step   (MONITOR_req_nmax) )

    rewind(IO_FID_CONF)
    do n = 1, MONITOR_req_nmax
       ! set default
       ITEM  = 'unknown'

       read(IO_FID_CONF,nml=MONITITEM,iostat=ierr)
       if( ierr /= 0 ) exit

       MONITOR_req_item(n) = ITEM
    enddo

    return
  end subroutine MONITOR_setup

  !-----------------------------------------------------------------------------
  subroutine MONITOR_reg( &
      itemid, &
      item,   &
      desc,   &
      unit,   &
      ktype   )
    implicit none

    integer,          intent(out) :: itemid
    character(len=*), intent( in) :: item
    character(len=*), intent( in) :: desc
    character(len=*), intent( in) :: unit
    character(len=*), intent( in) :: ktype

    character(len=8) :: lname

    integer :: n, nmax, reqid
    !---------------------------------------------------------------------------

    !--- search existing item
    itemid = -1
    nmax = min( MONITOR_id_count, MONITOR_req_nmax )
    do n = 1, nmax
       if ( trim(item) == trim(MONITOR_item(n)) ) then ! match existing item
          itemid = n
          return
       endif
    enddo

    if ( itemid < 0 ) then ! request-register matching check
       do n = 1, MONITOR_req_nmax
          if ( trim(item) == MONITOR_req_item(n) ) then
             itemid = MONITOR_id_count
             reqid  = n
             MONITOR_id_count = MONITOR_id_count + 1

             ! new file registration
             MONITOR_item(itemid) = trim(item)
             MONITOR_desc(itemid) = trim(desc)
             MONITOR_unit(itemid) = trim(unit)
             if    ( trim(ktype) == '2D' ) then
                MONITOR_ktype(itemid) = 'ZSFC'
                MONITOR_kmax (itemid) = 1
             elseif( trim(ktype) == '3D' ) then
                write(lname,'(A,I4.4)') 'ZDEF', KMAX
                MONITOR_ktype(itemid) = lname
                MONITOR_kmax (itemid) = KMAX
             endif
             MONITOR_step        (itemid) = 1

             if( IO_L ) write(IO_FID_LOG,*) '*** [MONIT] Item registration No.= ', itemid
             if( IO_L ) write(IO_FID_LOG,*) '] Name           : ', trim(MONITOR_item (itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] Description    : ', trim(MONITOR_desc (itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] Unit           : ', trim(MONITOR_unit (itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] Vert. type     : ', trim(MONITOR_ktype(itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] # of layer     : ', MONITOR_kmax(itemid)
          endif
       enddo
    endif

    return
  end subroutine MONITOR_reg

  !-----------------------------------------------------------------------------
  subroutine MONITOR_put( &
      itemid, &
      var     )
    use mod_geometrics, only: &
       area => GEOMETRICS_area, &
       vol  => GEOMETRICS_vol
    implicit none

    integer, intent(in) :: itemid
    real(8), intent(in) :: var(:,:,:)

    real(8) :: total

    integer :: k, i, j
    !---------------------------------------------------------------------------

    total = 0.D0
    if ( MONITOR_kmax(itemid) == KMAX ) then ! 3D
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          total = total + var(k,i,j) * vol(k,i,j)
       enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
          total = total + var(1,i,j) * area(1,i,j)
       enddo
       enddo
    endif

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,i4,A,f13.3,A,A,A)') &
               ']',trim(MONITOR_item(itemid)),', step=',MONITOR_step(itemid),':',total,'[',trim(MONITOR_unit(itemid)),']'

    MONITOR_step(itemid) = MONITOR_step(itemid) + 1

    return
  end subroutine MONITOR_put

  !-----------------------------------------------------------------------------
  subroutine MONITOR_in( &
      var,   &
      item,  &
      desc,  &
      unit,  &
      ktype  )
    implicit none

    real(8),          intent(in) :: var(:,:,:)
    character(len=*), intent(in) :: item
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: unit
    character(len=*), intent(in) :: ktype

    integer :: itemid
    !---------------------------------------------------------------------------

    call MONITOR_reg( itemid, item, desc, unit, ktype )
    
    if ( itemid > 0 ) then
       call MONITOR_put( itemid, var(:,:,:) )
    endif

    return
  end subroutine MONITOR_in

end module mod_monitor
!-------------------------------------------------------------------------------
