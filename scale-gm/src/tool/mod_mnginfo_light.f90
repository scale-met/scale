!-------------------------------------------------------------------------------
!
!+  Module manageinfo light
!
!-------------------------------------------------------------------------------
module mod_mnginfo_light
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !      mnginfo reader (extraced from mod_adm)
  !
  !++ Current Corresponding Author: H.Yashiro
  ! 
  !++ History: 
  !      Version   Date      Comment 
  !      -----------------------------------------------------------------------
  !      0.90      11-09-01  H.Yashiro : [NEW]
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_misc, only : &
    MISC_get_available_fid
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ public param & variable
  !
  public :: MNG_mnginfo_input
  public :: MNG_mnginfo_noinput
  !-----------------------------------------------------------------------------
  !
  !++ public param & variable
  !
  integer, public              :: MNG_PALL

  integer, allocatable, public :: MNG_prc_rnum(:)
  integer, allocatable, public :: MNG_prc_tab (:,:)
  integer, allocatable, public :: MNG_rgn2prc (:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> read mnginfo (light ver.)
  !-----------------------------------------------------------------------------
  subroutine MNG_mnginfo_input( rlevel,fname )
    implicit none

    integer,          intent(in) :: rlevel
    character(len=*), intent(in) :: fname

    integer, parameter :: PRC_RGN_NMAX = 2560

    integer :: num_of_rgn
    namelist /rgn_info/ &
         num_of_rgn             !--- number of region

    integer :: num_of_proc
    namelist /proc_info/ &
         num_of_proc             !--- number of run-processes

    integer :: peid
    integer :: num_of_mng
    integer :: mng_rgnid(PRC_RGN_NMAX)
    namelist /rgn_mng_info/ &
         peid,              & !--- process ID
         num_of_mng,        & !--- number of regions be managed
         mng_rgnid            !--- managed region ID

    integer :: fid, ierr
    integer :: m, n
    integer :: lall
    !---------------------------------------------------------------------------

    lall = 10 * (4**rlevel)

    mng_rgnid(:)=-1

    fid = MISC_get_available_fid()
    open(fid,file=trim(fname),status='old',form='formatted',iostat=ierr)
       if (ierr /= 0) then
          write(*,*) "cannot read mnginfo file :",trim(fname)
          stop
       endif

       read(fid,nml=rgn_info) 
       if ( num_of_rgn /= lall ) then ! [add] H.Yashiro 20120621
          write(*,*) "Inconsintent between rlevel and mnginfo."
          write(*,*) "rlevel,  num_of_mng=", rlevel, lall
          write(*,*) "mnginfo, num_of_rgn=", trim(fname), num_of_rgn
          stop
       endif

       read(fid,nml=proc_info) 
       MNG_PALL = num_of_proc

       allocate( MNG_prc_rnum(num_of_proc) )
       allocate( MNG_prc_tab (num_of_rgn,num_of_proc) )
       allocate( MNG_rgn2prc (num_of_rgn) )
       MNG_prc_tab(:,:) = -1

       do m = 1, num_of_proc
          read(fid,nml=rgn_mng_info) 

          MNG_prc_rnum(m)                = num_of_mng
          MNG_prc_tab(1:num_of_mng,peid) = mng_rgnid(1:num_of_mng)

          do n = 1, num_of_mng
             MNG_rgn2prc(mng_rgnid(n)) = peid
          enddo
       enddo

    close(fid)

    return
  end subroutine MNG_mnginfo_input

  !-----------------------------------------------------------------------------
  !> read mnginfo (light ver.)
  !-----------------------------------------------------------------------------
  subroutine MNG_mnginfo_noinput(rlevel)
    implicit none

    integer, intent(in) :: rlevel

    integer :: n
    integer :: num_of_mng
    !---------------------------------------------------------------------------

    num_of_mng = 10 * (4**rlevel)
    MNG_PALL   = 1

    allocate( MNG_prc_rnum(1) )
    allocate( MNG_prc_tab (num_of_mng,1) )
    allocate( MNG_rgn2prc (num_of_mng) )
    MNG_prc_tab(:,:) = -1

    MNG_prc_rnum(1) = num_of_mng

    do n = 1, num_of_mng
       MNG_prc_tab(n,1) = n
    enddo

    MNG_rgn2prc(:) = 1

    return
  end subroutine MNG_mnginfo_noinput

end module mod_mnginfo_light
!-----------------------------------------------------------------------------

