!-------------------------------------------------------------------------------
!
!+  communication module
!
!-------------------------------------------------------------------------------
module mod_comm
  !-----------------------------------------------------------------------------
  !
  !++ description:
  !       this module is for the communication based on mpi library.
  !
  !++ Current Corresponding Author : K.Goto, H.Tomita
  !
  !++ History:
  !      Version    Date      Comment
  !      -----------------------------------------------------------------------
  !      0.00       04-02-17  Imported from igdc-4.33
  !                 06-09-??  K.Goto bug fix
  !                 06-10-08  S.Iga  add namelist (&COMMPARAM  max_varmax)
  !                 07-11-07  T.Mitsui add varmax check option(opt_check_varmax)
  !                 09-03-10  H.Tomita : Transplanting COMM_data_transfer2 from
  !                                      mod[mod_varcomm].
  !                 09-03-10  H.Tomita : rename COMM_data_transfer2 to COMM_var.
  !                 09-09-17  S.Iga : Add debug option and barrier option
  !                 10-06-07  S.Iga: new grid is implemented
  !                              (only the attribute of max_comm_xxx and
  !                              max_comm is changed from parameter to variable)
  !                 11-01-24  C.Kodama: Reduce memory usage in large rlevel.
  !                              (provided by Terai-san @ RIKEN)
  !                              Modified line: (20101207 teraim)
  !                 11-04-26  C.Kodama: default value of opt_check_varmax is changed to .true.
  !                                     and modify its behavior to abort when cmax exceeds max_maxvar*ADM_kall.
  !                 11-05-06  Y.Yamada: Merge tuning code with original code
  !                              (provided by Yamamoto-san @ NEC)
  !                              Modified line: !=org=
  !                 11-07-21  T.Ohno: A public variable 'comm_pl' is added.
  !                           If 'comm_pl' is false, pole data is not used in
  !                           COMM_data_transfer and COMM_var.
  !                 11-11-30  S.Iga (commit) : Modification around COMM_var,
  !                              suggested and modified by T.Inoue on 11-10-24
  !                 11-12-14  T.Seiki : allocatable variables are not permitted in type structure.
  !                              allocatable => pointer  (only @ SR16000 and ES)
  !                 12-03-26  T.Seiki : bug-fix if opt_comm_dbg=.true.
  !                 12-06-27  T.Ohno : bug fix for simulations at which 'comm_pl' is false
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ public procedure
  !
  public ::  COMM_setup
  public ::  COMM_data_transfer
  public ::  COMM_var
  public ::  COMM_Stat_sum
  public ::  COMM_Stat_sum_eachlayer
  public ::  COMM_Stat_avg
  public ::  COMM_Stat_max
  public ::  COMM_Stat_min

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: comm_pl = .true. ! T.Ohno 110721

  integer, public :: COMM_datatype

  ! for send
  integer, public, allocatable :: nsmax(:,:)
  integer, public, allocatable :: sendinfo(:,:,:,:)
  integer, public, allocatable :: sendlist(:,:,:,:)
  integer, public, allocatable :: nsmax_pl(:,:)
  integer, public, allocatable :: sendinfo_pl(:,:,:,:)
  integer, public, allocatable :: sendlist_pl(:,:,:,:)

  ! for recv
  integer, public, allocatable :: nrmax(:,:)
  integer, public, allocatable :: recvinfo(:,:,:,:)
  integer, public, allocatable :: recvlist(:,:,:,:)
  integer, public, allocatable :: nrmax_pl(:,:)
  integer, public, allocatable :: recvinfo_pl(:,:,:,:)
  integer, public, allocatable :: recvlist_pl(:,:,:,:)

  ! for copy
  integer, public, allocatable :: ncmax_r2r(:)
  integer, public, allocatable :: copyinfo_r2r(:,:,:)
  integer, public, allocatable :: recvlist_r2r(:,:,:)
  integer, public, allocatable :: sendlist_r2r(:,:,:)

  integer, public, allocatable :: ncmax_r2p(:)
  integer, public, allocatable :: copyinfo_r2p(:,:,:)
  integer, public, allocatable :: recvlist_r2p(:,:,:)
  integer, public, allocatable :: sendlist_r2p(:,:,:)

  integer, public, allocatable :: ncmax_p2r(:)
  integer, public, allocatable :: copyinfo_p2r(:,:,:)
  integer, public, allocatable :: recvlist_p2r(:,:,:)
  integer, public, allocatable :: sendlist_p2r(:,:,:)

  integer, public, allocatable :: ncmax_sgp(:)
  integer, public, allocatable :: copyinfo_sgp(:,:,:)
  integer, public, allocatable :: recvlist_sgp(:,:,:)
  integer, public, allocatable :: sendlist_sgp(:,:,:)

#ifdef _ACCCUDA
  real(RP), public, allocatable, pinned :: sendbuf(:,:)
  real(RP), public, allocatable, pinned :: recvbuf(:,:)
#else
  real(RP), public, allocatable :: sendbuf(:,:)
  real(RP), public, allocatable :: recvbuf(:,:)
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: ptr_prcid       = 1
  integer, private, parameter :: ptr_lrgnid      = 2
  !
  integer, private, parameter :: elemsize_comm   = 3
  integer, private, parameter :: SIZE_COMM       = 1
  integer, private, parameter :: LRGNID_COMM     = 2
  integer, private, parameter :: BASE_COMM       = 3
  !
  integer, private, parameter :: elemsize_copy   = 3
  integer, private, parameter :: SIZE_COPY       = 1
  integer, private, parameter :: LRGNID_COPY     = 2
  integer, private, parameter :: SRC_LRGNID_COPY = 3

  integer, private, parameter :: max_comm_r2r    = 9

  integer, private :: max_comm_r2p
  integer, private :: max_comm_p2r
  integer, private :: max_comm
  integer, private :: max_varmax = 100
  logical, private :: opt_check_varmax = .true.
  logical, private :: opt_comm_dbg = .false.
  logical, private :: opt_comm_barrier = .false.
  real(RP), private :: dbg_sendbuf_init
  real(RP), private :: dbg_recvbuf_init
  integer,private,allocatable :: dbg_areq_save(:,:)

  integer, private :: rank_me
  integer, private :: max_comm_prc
  integer, private :: maxdatasize_s
  integer, private :: maxdatasize_r

  integer, private :: maxn
  integer, private :: maxm
  integer, private :: maxl

  integer, private :: maxn_pl
  integer, private :: maxm_pl
  integer, private :: maxl_pl

  integer, private :: maxn_r2r
  integer, private :: maxm_r2r
  integer, private :: maxl_r2r

  integer, private :: maxn_r2p
  integer, private :: maxm_r2p
  integer, private :: maxl_r2p

  integer, private :: maxn_p2r
  integer, private :: maxm_p2r
  integer, private :: maxl_p2r

  integer, private :: maxn_sgp
  integer, private :: maxm_sgp
  integer, private :: maxl_sgp

  integer, private, allocatable :: prc_tab_rev(:,:)

  integer, private, allocatable :: clist(:)

  integer, private, allocatable :: temp_sendorder(:,:)
  integer, private, allocatable :: temp_recvorder(:,:)
  integer, private, allocatable :: temp_dest_rgn(:,:,:)
  integer, private, allocatable :: temp_src_rgn(:,:,:)
  integer, private, allocatable :: temp_dest_rgn_pl(:,:,:)
  integer, private, allocatable :: temp_src_rgn_pl(:,:,:)
  integer, private, allocatable :: tsb(:)

  integer, private, allocatable :: ssize(:,:)
  integer, private, allocatable :: sendtag(:,:)
  integer, private, allocatable :: somax(:)
  integer, private, allocatable :: destrank(:,:)
  integer, private, allocatable :: rsize(:,:)
  integer, private, allocatable :: recvtag(:,:)
  integer, private, allocatable :: romax(:)
  integer, private, allocatable :: sourcerank(:,:)

  integer, private, allocatable :: n_nspl(:,:)

  integer, private, allocatable :: n_hemisphere_copy(:,:,:)
  integer, private, allocatable :: s_hemisphere_copy(:,:,:)

  integer, private, allocatable :: tempbuf2D(:,:)
  integer, private, allocatable :: tempbuf3D(:,:,:)
  integer, private, allocatable :: tempbuf4D(:,:,:,:)
  integer, private, allocatable :: tempbuf3D2(:,:,:)

  integer, private, allocatable :: rsize_r2r(:,:,:)
  integer, private, allocatable :: ssize_r2r(:,:,:)
  integer, private, allocatable :: sourceid_r2r(:,:,:)
  integer, private, allocatable :: destid_r2r(:,:,:)
  integer, private, allocatable :: msend_r2r(:,:,:)
  integer, private, allocatable :: maxcommrecv_r2r(:,:)
  integer, private, allocatable :: maxcommsend_r2r(:,:)
  integer, private, allocatable :: rlist_r2r(:,:,:,:)
  integer, private, allocatable :: qlist_r2r(:,:,:,:)
  integer, private, allocatable :: slist_r2r(:,:,:,:)

  integer, private              :: max_datasize_r2r
  real(RP), private, allocatable :: recvbuf_r2r(:,:,:)
  real(RP), private, allocatable :: sendbuf_r2r(:,:,:)

  integer, private, allocatable :: rsize_r2p(:,:,:)
  integer, private, allocatable :: ssize_r2p(:,:,:)
  integer, private, allocatable :: source_prc_r2p(:,:,:)
  integer, private, allocatable :: source_rgn_r2p(:,:,:)
  integer, private, allocatable :: dest_prc_r2p(:,:,:)
  integer, private, allocatable :: maxcommrecv_r2p(:,:)
  integer, private, allocatable :: maxcommsend_r2p(:,:)
  integer, private, allocatable :: recvtag_r2p(:,:,:)
  integer, private, allocatable :: sendtag_r2p(:,:,:)
  integer, private, allocatable :: rlist_r2p(:,:,:,:)
  integer, private, allocatable :: qlist_r2p(:,:,:,:)
  integer, private, allocatable :: slist_r2p(:,:,:,:)

  integer,private               ::  max_datasize_r2p
  real(RP), private, allocatable :: recvbuf_r2p(:,:,:)
  real(RP), private, allocatable :: sendbuf_r2p(:,:,:)

  integer, private, allocatable :: recvtag_p2r(:,:)
  integer, private, allocatable :: sendtag_p2r(:,:)
  real(RP), private, allocatable :: sendbuf_p2r(:,:)
  real(RP), private, allocatable :: recvbuf_p2r(:,:)

  integer, private, allocatable :: dest_rank_all(:,:,:)
  integer, private, allocatable :: src_rank_all(:,:,:)

  integer, private, allocatable :: imin(:), imax(:)
  integer, private, allocatable :: jmin(:), jmax(:)
  integer, private, allocatable :: gmin(:), gmax(:),gall(:)

  integer, private, allocatable :: nmin_nspl(:), nmax_nspl(:)
  integer, private, allocatable :: pmin_nspl(:), pmax_nspl(:)
  integer, private, allocatable :: lmin_nspl(:), lmax_nspl(:)
  integer, private, allocatable :: gmin_nspl(:), gmax_nspl(:)
  integer, private, allocatable :: gall_nspl(:)

  integer, private, allocatable :: pl_index(:,:,:,:)

  integer, private :: halomax
  integer, private :: kmax

  type type_tempsb
    integer        ::  num
    integer, pointer ::  col(:)
    integer, pointer ::  val(:)
  end type

  type(type_tempsb),allocatable :: tempsb(:)

  integer, parameter ::  max_size = 10 ! This is not optimal value.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !(20101207) added by teraim
  subroutine init_tempsb
    use mod_adm, only : ADM_rgn_nmax
    implicit none
    integer :: i
    !
    allocate(tempsb(ADM_rgn_nmax+2))
    tempsb(:)%num=0
    !
    do i=1, ADM_rgn_nmax+2
      allocate(tempsb(i)%col(max_size))
      allocate(tempsb(i)%val(max_size))
      tempsb(i)%col(:)=-1
      tempsb(i)%val(:)=-1
    enddo
  end subroutine
  !(20101207) added by teraim
  subroutine finalize_tempsb
    use mod_adm, only : ADM_rgn_nmax
    implicit none
    integer :: i
    !
    do i=1, ADM_rgn_nmax+2
      deallocate(tempsb(i)%col)
      deallocate(tempsb(i)%val)
    enddo
    deallocate(tempsb)
  end subroutine
  !(20101207) added by teraim
  subroutine add_tempsb(icol, irow, ival)
    implicit none
    integer,intent(in) :: icol,irow,ival
    !
    if(ival > 0) then
      if(tempsb(irow)%num < max_size) then
        tempsb(irow)%num = tempsb(irow)%num + 1
        tempsb(irow)%col(tempsb(irow)%num)=icol
        tempsb(irow)%val(tempsb(irow)%num)=ival
      else
        write(*,*)"range of list is over."
        stop
      endif
    endif
  end subroutine

  !(20101207) added by teraim
  subroutine get_tempsb(icol, irow, ret)
    implicit none
    integer,intent(in) :: icol,irow
    integer,intent(out) :: ret
    integer :: i
    !
    ret = 0
    do i=1, max_size
      if(tempsb(irow)%col(i) == icol) then
        ret = tempsb(irow)%val(i)
        exit
      endif
    enddo
  end subroutine

  !-----------------------------------------------------------------------------
  subroutine COMM_setup( &
       max_hallo_num,    & !--- IN : number of hallo regions
       debug             ) !--- IN : debug flag
    use mod_adm, only :    &
       ADM_w,           &
       ADM_e,           &
       ADM_n,           &
       ADM_s,           &
       ADM_sw,          &
       ADM_nw,          &
       ADM_ne,          &
       ADM_se,          &
       ADM_rid,         &
       ADM_dir,         &
       ADM_vlink_nmax,  &
       ADM_rgn_nmax_pl, &
       ADM_npl,         &
       ADM_spl,         &
       ADM_gslf_pl,     &
       ADM_prc_all,     &
       ADM_prc_rnum,    &
       ADM_prc_tab,     &
       ADM_prc_me,      &
       ADM_rgn_nmax,    &
       ADM_rgn_etab,    &
       ADM_rgn_vnum,    &
       ADM_rgn_vtab,    &
       ADM_rgn_vtab_pl, &
       ADM_gmin,        &
       ADM_gmax,        &
       ADM_gall_1d,     &
       ADM_lall,        &
       ADM_kall,        &
       ADM_prc_nspl,    &
       ADM_COMM_world,  &
       ADM_rgn2prc,     &
       ADM_CTL_FID,     &
       ADM_LOG_FID,     &
       ADM_proc_stop
    implicit none

    integer,intent(in),optional ::  max_hallo_num
    logical,intent(in),optional ::  debug

    integer ::  i,j,l,n,m,p,q
    integer ::  rgnid

    integer :: lr,mr
    integer :: ls,ms
    integer :: nr,nc,ns
    integer :: rs,cs,ss
    integer :: nd,ld,pl,halo
    integer :: in,jn
    integer :: srgnid,rrgnid
    integer :: ck
    integer :: srank,drank
    integer :: ro,so
    !
    integer :: suf,g_1d
    suf(i,j,g_1d)=(g_1d)*((j)-1)+(i)
    !
    integer :: rgnid1,rgnid2,ret !(20101207) added by teraim
    !
    ! Iga(061008) ==>
    namelist / COMMPARAM /   &
         max_varmax,         & ! max number of communication variables
         opt_check_varmax,   & ! check option of varmax [Add] T.Mitsui 07/11/07
         opt_comm_dbg,       & ! debug option of comm_data_transfer [Add] S.Iga 0909XX
         opt_comm_barrier      ! debug option of comm_data_transfer [Add] S.Iga 0909XX

    integer ::  ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[comm]/Category[common share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=COMMPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** COMMPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist COMMPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist COMMPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=COMMPARAM)

    if ( RP == DP ) then
       COMM_datatype = MPI_DOUBLE_PRECISION
    elseif( RP == SP ) then
       COMM_datatype = MPI_REAL
    else
       write(*,*) 'xxx precision is not supportd'
       call ADM_proc_stop
    endif

    if (present(max_hallo_num)) then
       halomax=max_hallo_num
    else
       halomax=1
    endif

    max_comm_r2p=ADM_vlink_nmax*2!S.Iga100607
    max_comm_p2r=ADM_vlink_nmax*2!S.Iga100607
    max_comm=max_comm_r2r+max_comm_r2p+max_comm_p2r!S.Iga100607

    kmax=ADM_kall

    allocate(prc_tab_rev(ptr_prcid:ptr_lrgnid,ADM_rgn_nmax))

    do p=1,ADM_prc_all
       do n=1,ADM_prc_rnum(p)
          prc_tab_rev(ptr_prcid,ADM_prc_tab(n,p))=p
          prc_tab_rev(ptr_lrgnid,ADM_prc_tab(n,p))=n
       enddo
    enddo

    if(ADM_prc_nspl(ADM_npl) < 0 .and. ADM_prc_nspl(ADM_spl) <0 ) comm_pl = .false. ! T.Ohno 110721

    !
    allocate(imin(halomax))
    allocate(imax(halomax))
    allocate(jmin(halomax))
    allocate(jmax(halomax))
    allocate(gmin(halomax))
    allocate(gmax(halomax))
    allocate(gall(halomax))
    !
    !(20101207) changed by teraim
    allocate(rsize_r2r(max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(ssize_r2r(max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(sourceid_r2r(max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(destid_r2r(max_comm_r2r,halomax,ADM_rgn_nmax))
    !allocate(mrecv_r2r(ADM_rgn_nmax,halomax,ADM_rgn_nmax))
    !allocate(msend_r2r(ADM_rgn_nmax,halomax,ADM_rgn_nmax))
    rgnid1=ADM_prc_tab(1,ADM_prc_me)
    rgnid2=ADM_prc_tab(ADM_prc_rnum(ADM_prc_me),ADM_prc_me)
    allocate(msend_r2r(ADM_rgn_nmax,halomax,rgnid1:rgnid2))
    allocate(maxcommrecv_r2r(halomax,ADM_rgn_nmax))
    allocate(maxcommsend_r2r(halomax,ADM_rgn_nmax))
    !allocate(recvtag_r2r(ADM_rgn_nmax,halomax,ADM_rgn_nmax))
    !allocate(sendtag_r2r(ADM_rgn_nmax,halomax,ADM_rgn_nmax))
    !
    imin(halomax)=(ADM_gmin-1)+halomax
    imax(halomax)=(ADM_gmax-1)+halomax
    jmin(halomax)=(ADM_gmin-1)+halomax
    jmax(halomax)=(ADM_gmax-1)+halomax
    gmin(halomax)=(ADM_gmin-1)+halomax
    gmax(halomax)=(ADM_gmax-1)+halomax
    gall(halomax)=(ADM_gall_1d-2)+2*halomax
    !
    max_datasize_r2r=(gmax(halomax)-gmin(halomax)+1)*halomax
    allocate(rlist_r2r(max_datasize_r2r,max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(qlist_r2r(max_datasize_r2r,max_comm_r2r,halomax,ADM_rgn_nmax))
    allocate(slist_r2r(max_datasize_r2r,max_comm_r2r,halomax,ADM_rgn_nmax))
    !
    allocate(recvbuf_r2r(max_datasize_r2r*kmax*max_varmax  &
         ,ADM_prc_rnum(ADM_prc_me),max_comm_r2r))
    allocate(sendbuf_r2r(max_datasize_r2r*kmax*max_varmax  &
         ,ADM_prc_rnum(ADM_prc_me),max_comm_r2r))

    allocate(n_hemisphere_copy(ADM_w:ADM_s,halomax,ADM_rgn_nmax))
    allocate(s_hemisphere_copy(ADM_w:ADM_s,halomax,ADM_rgn_nmax))

    allocate(tempbuf2D(halomax,ADM_prc_rnum(ADM_prc_me)))
    allocate(tempbuf3D(max_comm_r2r,halomax,ADM_prc_rnum(ADM_prc_me)))
    allocate(tempbuf4D(max_datasize_r2r,max_comm_r2r,halomax,ADM_prc_rnum(ADM_prc_me)))
    allocate(tempbuf3D2(ADM_w:ADM_s,halomax,ADM_prc_rnum(ADM_prc_me)))

    !
    rsize_r2r(:,:,:)=0
    ssize_r2r(:,:,:)=0
    sourceid_r2r(:,:,:)=-1
    destid_r2r(:,:,:)=-1
    !mrecv_r2r(:,:,:)=-1
    msend_r2r(:,:,:)=-1
    maxcommrecv_r2r(:,:)=max_comm_r2r
    maxcommsend_r2r(:,:)=max_comm_r2r
    !(20101207) removed by teraim
    !recvtag_r2r(:,:,:)=-1
    !sendtag_r2r(:,:,:)=-1
    !
    rlist_r2r(:,:,:,:)=-1
    qlist_r2r(:,:,:,:)=-1
    slist_r2r(:,:,:,:)=-1
    !
    n_hemisphere_copy(:,:,:)=0
    s_hemisphere_copy(:,:,:)=0
    !
    allocate(rsize_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(ssize_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(source_prc_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(source_rgn_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(dest_prc_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(maxcommrecv_r2p(ADM_npl:ADM_spl,halomax))
    allocate(maxcommsend_r2p(ADM_npl:ADM_spl,halomax))
    allocate(recvtag_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(sendtag_r2p(max_comm_r2p,ADM_npl:ADM_spl,halomax))
    !
    max_datasize_r2p=halomax*(halomax+1)/2
    !
    allocate(rlist_r2p(max_datasize_r2p,max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(qlist_r2p(max_datasize_r2p,max_comm_r2p,ADM_npl:ADM_spl,halomax))
    allocate(slist_r2p(max_datasize_r2p,max_comm_r2p,ADM_npl:ADM_spl,halomax))
    !
    allocate(recvbuf_r2p(max_datasize_r2p*kmax*max_varmax &
         ,max_comm_r2p,ADM_npl:ADM_spl))
    allocate(sendbuf_r2p(max_datasize_r2p*kmax*max_varmax &
         ,max_comm_r2p,ADM_npl:ADM_spl))
    !
    rsize_r2p(:,:,:)=0
    ssize_r2p(:,:,:)=0
    source_prc_r2p(:,:,:)=-1
    source_rgn_r2p(:,:,:)=-1
    dest_prc_r2p(:,:,:)=-1
    maxcommrecv_r2p(:,:)=max_comm_r2p
    maxcommsend_r2p(:,:)=max_comm_r2p
    recvtag_r2p(:,:,:)=-1
    sendtag_r2p(:,:,:)=-1
    !
    rlist_r2p(:,:,:,:)=-1
    qlist_r2p(:,:,:,:)=-1
    slist_r2p(:,:,:,:)=-1
    !
!!!!!!!!!!!!!!!
    allocate(nmin_nspl(1:halomax))
    allocate(nmax_nspl(1:halomax))
    allocate(pmin_nspl(1:halomax))
    allocate(pmax_nspl(1:halomax))
    allocate(lmin_nspl(1:halomax))
    allocate(lmax_nspl(1:halomax))
    allocate(gmin_nspl(1:halomax))
    allocate(gmax_nspl(1:halomax))
    allocate(gall_nspl(1:halomax))
    nmin_nspl(halomax)=1
    nmax_nspl(halomax)=halomax+1
    pmin_nspl(halomax)=1
    pmax_nspl(halomax)=ADM_vlink_nmax
    lmin_nspl(halomax)=1
    lmax_nspl(halomax)=halomax
    gmin_nspl(halomax)=2
    gmax_nspl(halomax)=1+5*halomax*(halomax+1)/2
    gall_nspl(halomax)=1+5*halomax*(halomax+1)/2
    allocate(pl_index(nmin_nspl(halomax):nmax_nspl(halomax) &
         ,pmin_nspl(halomax):pmax_nspl(halomax) &
         ,lmin_nspl(halomax):lmax_nspl(halomax),halomax))
    !
    do halo=1,halomax
       !
       imin(halo)=(ADM_gmin-1)+halo
       imax(halo)=(ADM_gmax-1)+halo
       jmin(halo)=(ADM_gmin-1)+halo
       jmax(halo)=(ADM_gmax-1)+halo
       gmin(halo)=(ADM_gmin-1)+halo
       gmax(halo)=(ADM_gmax-1)+halo
       gall(halo)=(ADM_gall_1d-2)+2*halo
       !
       nmin_nspl(halo)=1
       nmax_nspl(halo)=halo+1
       pmin_nspl(halo)=1
       pmax_nspl(halo)=ADM_vlink_nmax
       lmin_nspl(halo)=1
       lmax_nspl(halo)=halo
       gmin_nspl(halo)=2
       gmax_nspl(halo)=1+5*halo*(halo+1)/2
       gall_nspl(halo)=1+5*halo*(halo+1)/2
       pl_index(:,:,:,halo)=-1
       do l=lmin_nspl(halo),lmax_nspl(halo)
          do p=pmin_nspl(halo),pmax_nspl(halo)
             do n=nmin_nspl(halo),l+1
                pl_index(n,p,l,halo)=n+(p-1)*l+(1+5*(l-1)*l/2)
             enddo
          enddo
          pl_index(l+1,pmax_nspl(halo),l,halo)=nmin_nspl(halo) &
               +(pmin_nspl(halo)-1)*l &
               +(1+5*(l-1)*l/2)
       enddo
       !
       ! --- r2p ----
       if(comm_pl)then
        do p=pmin_nspl(halo),pmax_nspl(halo)
           rsize_r2p(p,ADM_npl,halo)=halo*(halo+1)/2
           source_prc_r2p(p,ADM_npl,halo)=prc_tab_rev(ptr_prcid &
                ,ADM_rgn_vtab_pl(ADM_rid,ADM_npl,p))
           source_rgn_r2p(p,ADM_npl,halo)=prc_tab_rev(ptr_lrgnid &
                ,ADM_rgn_vtab_pl(ADM_rid,ADM_npl,p))
        enddo
        do p=pmin_nspl(halo),pmax_nspl(halo)
           q=0
           do ld=lmin_nspl(halo),lmax_nspl(halo)
              do nd=nmin_nspl(halo),ld
                 q=q+1
                 in=-nd+ld+nmin_nspl(halo)-lmin_nspl(halo)+imin(halo)
                 jn=-nd+nmin_nspl(halo)+(jmax(halo)-jmin(halo))+jmin(halo)
                 rlist_r2p(q,mod(p,ADM_vlink_nmax)+1,ADM_npl,halo) &
                      =pl_index(nd+1,p,ld,halo)
                 qlist_r2p(q,mod(p,ADM_vlink_nmax)+1,ADM_npl,halo)=suf(in,jn,gall(halo))
              enddo
           enddo
        enddo
        do p=pmin_nspl(halo),pmax_nspl(halo)
           rsize_r2p(p,ADM_spl,halo)=halo*(halo+1)/2
           source_prc_r2p(p,ADM_spl,halo)=prc_tab_rev(ptr_prcid &
                ,ADM_rgn_vtab_pl(ADM_rid,ADM_spl,p))
           source_rgn_r2p(p,ADM_spl,halo)=prc_tab_rev(ptr_lrgnid &
                ,ADM_rgn_vtab_pl(ADM_rid,ADM_spl,p))
        enddo
        do p=pmin_nspl(halo),pmax_nspl(halo)
           q=0
           do ld=lmin_nspl(halo),lmax_nspl(halo)
              do nd=nmin_nspl(halo),ld
                 q=q+1
                 in=nd-ld-nmin_nspl(halo) &
                      +lmin_nspl(halo)+(imax(halo)-imin(halo))+imin(halo)
                 jn=nd-nmin_nspl(halo)+jmin(halo)
                 rlist_r2p(q,p,ADM_spl,halo)=pl_index(nd,p,ld,halo)
                 qlist_r2p(q,p,ADM_spl,halo)=suf(in,jn,gall(halo))
              enddo
           enddo
        enddo
        maxcommrecv_r2p(ADM_npl,halo)=(pmax_nspl(halo)-pmin_nspl(halo)+1)
        maxcommrecv_r2p(ADM_spl,halo)=(pmax_nspl(halo)-pmin_nspl(halo)+1)
        !
        do pl=ADM_npl,ADM_spl
           do p=1,maxcommrecv_r2p(pl,halo)
              if (ADM_prc_me==source_prc_r2p(p,pl,halo)) then
                 dest_prc_r2p(p,pl,halo)=ADM_prc_nspl(pl)
                 ssize_r2p(p,pl,halo)=rsize_r2p(p,pl,halo)
                 do q=1,ssize_r2p(p,pl,halo)
                    slist_r2p(q,p,pl,halo)=qlist_r2p(q,p,pl,halo)
                 enddo
              endif
              sendtag_r2p(p,pl,halo)=pl+(ADM_spl-ADM_npl+1)*(p-1) &
                   +ADM_rgn_nmax**2+ADM_vlink_nmax*2
              recvtag_r2p(p,pl,halo)=sendtag_r2p(p,pl,halo)

!              write(*,*) 'sendtag_r2p',ADM_prc_me,p,pl,halo,sendtag_r2p(p,pl,halo)
           enddo
           maxcommsend_r2p(pl,halo)=(pmax_nspl(halo)-pmin_nspl(halo)+1)
        enddo
       endif
       !
       ! --- r2r ----
       do l=1,ADM_prc_rnum(ADM_prc_me)
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          m=0
          if (ADM_rgn_etab(ADM_dir,ADM_sw,rgnid)==ADM_ne) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_sw,rgnid)
                ! mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                n=0
                do j=jmin(halo)-halo,jmin(halo)-1
                   do i=imin(halo),imax(halo)
                      n=n+1
                      in=i
                      jn=j+jmax(halo)+1-jmin(halo)
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          endif
          if ((ADM_rgn_vnum(ADM_w,rgnid)==3)) then
             if (ADM_rgn_etab(ADM_dir,ADM_sw,rgnid)==ADM_se) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_sw,rgnid)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207)removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)+1,imax(halo)+(j-(jmin(halo)-1))
                         n=n+1
                         in=j+jmax(halo)-2*jmin(halo)+imin(halo)+1
                         jn=-i+j+imin(halo)+jmax(halo)-jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                      do i=imin(halo)-halo,imin(halo)-1+(j-(jmin(halo)-1))
                         n=n+1
                         in=i+imax(halo)+1-imin(halo)
                         jn=j+jmax(halo)+1-jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
          endif
          !
          if (ADM_rgn_etab(ADM_dir,ADM_nw,rgnid)==ADM_ne) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_nw,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                n=0
                do i=imin(halo)-halo,imin(halo)-1
                   do j=jmin(halo)+(i-(imin(halo)-1)),jmax(halo)+(i-(imin(halo)-1))
                      n=n+1
                      in=i-j+imax(halo)-imin(halo)+jmin(halo)+1
                      jn=i+imax(halo)-2*imin(halo)+jmin(halo)+1
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          elseif (ADM_rgn_etab(ADM_dir,ADM_nw,rgnid)==ADM_se) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_nw,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                n=0
                do j=jmin(halo),jmax(halo)
                   do i=imin(halo)-halo,imin(halo)-1
                      n=n+1
                      in=i+imax(halo)-imin(halo)+1
                      jn=j
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          endif
          !
!!!!!
          if ((ADM_rgn_vnum(ADM_n,rgnid)==5)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,2)==ADM_n) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo+1)/2-1
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imin(halo)-halo,imin(halo)-1
                      do j=jmax(halo)+1+(i-(imin(halo)-1)) &
                           ,min(jmax(halo)+1,jmax(halo)+(halo-1)+(i-(imin(halo)-1)))
                         n=n+1
                         in=-j+jmax(halo)+imin(halo)+1
                         jn=i-j+imax(halo)-2*imin(halo)+jmax(halo)+jmin(halo)+2
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             if (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,3)==ADM_n) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,3)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmax(halo)+2,jmax(halo)+halo
                      do i=imin(halo)+1,imin(halo)+1+(j-(jmax(halo)+2))
                         n=n+1
                         in=-i+j-2*jmax(halo)+jmin(halo)+imax(halo)+imin(halo)-1
                         jn=-i+imax(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
!!!!!
          !
          if (ADM_rgn_etab(ADM_dir,ADM_ne,rgnid)==ADM_nw) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_ne,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                n=0
                do j=jmax(halo)+1,jmax(halo)+halo
                   do i=imin(halo)+1+(j-(jmax(halo)+1)),imax(halo)+1+(j-(jmax(halo)+1))
                      n=n+1
                      in=j-jmax(halo)+imin(halo)-1
                      jn=-i+j+imax(halo)-jmax(halo)+jmin(halo)
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          elseif (ADM_rgn_etab(ADM_dir,ADM_ne,rgnid)==ADM_sw) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_ne,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207)removed by teraim
                n=0
                do j=jmax(halo)+1,jmax(halo)+halo
                   do i=imin(halo),jmax(halo)
                      n=n+1
                      in=i
                      jn=j-jmax(halo)-1+jmin(halo)
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          endif
          !
          if (ADM_rgn_etab(ADM_dir,ADM_se,rgnid)==ADM_nw) then
             if (halo>=1) then
                m=m+1
                rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_se,rgnid)
                !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207)removed by teraim
                n=0
                do j=jmin(halo),jmax(halo)
                   do i=imax(halo)+1,imax(halo)+halo
                      n=n+1
                      in=i-imax(halo)+imin(halo)-1
                      jn=j
                      rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                      qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                   enddo
                enddo
             endif
          endif
          if ((ADM_rgn_vnum(ADM_e,rgnid)==3)) then
             if (ADM_rgn_etab(ADM_dir,ADM_se,rgnid)==ADM_sw) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_se,rgnid)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+1,imax(halo)+halo
                      do j=jmin(halo)+1+(i-(imax(halo)+1)),jmax(halo)
                         n=n+1
                         in=i-j-imax(halo)+jmax(halo)+imin(halo)
                         jn=i-imax(halo)-1+jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                      do j=jmax(halo)+1+(i-(imax(halo)+1)),jmax(halo)+halo
                         n=n+1
                         in=i-imax(halo)-1+imin(halo)
                         jn=j-jmax(halo)-1+jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
          endif
          !
!!!!!
          if ((ADM_rgn_vnum(ADM_s,rgnid)==5)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,2)==ADM_s) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+2,imax(halo)+halo
                      do j=jmin(halo)+1,jmin(halo)+1+(i-(imax(halo)+2))
                         n=n+1
                         in=-j+jmax(halo)+imin(halo)+1
                         jn=i-j-2*imax(halo)+imin(halo)+jmax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,3)==ADM_s) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo+1)/2-1
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,3)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imax(halo)+1+(j-(jmin(halo)-1)) &
                           ,min(imax(halo)+1,imax(halo)+(halo-1)+(j-(jmin(halo)-1)))
                         n=n+1
                         in=-i+j+imax(halo)+jmax(halo)+imin(halo)-2*jmin(halo)+2
                         jn=-i+imax(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
          !!
          !
          !!
          if ((ADM_rgn_vnum(ADM_w,rgnid)==4)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_w,rgnid,2)==ADM_n) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo+1)/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_w,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)-1+(j-(jmin(halo)-1)),imin(halo)-1
                         n=n+1
                         in=i-j+imax(halo)-imin(halo)-jmax(halo)+2*jmin(halo)
                         jn=i+imax(halo)-2*imin(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_w,rgnid,2)==ADM_e) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_w,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)-halo,imin(halo)-1
                         n=n+1
                         in=i+imax(halo)-imin(halo)+1
                         jn=j+jmax(halo)-jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_w,rgnid,2)==ADM_s) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo+1)/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_w,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)-halo,imin(halo)-halo+(j-(jmin(halo)-halo))
                         !do i=imin(halo)-halo,imin(halo)-1
                         !  do j=i+jmin(halo)-imin(halo),jmin(halo)-1
                         n=n+1
                         in=j+jmax(halo)+1-2*jmin(halo)+imin(halo)
                         jn=-i+j-imax(halo)+2*imin(halo)+jmax(halo)-jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             if (ADM_rgn_etab(ADM_dir,ADM_sw,rgnid)==ADM_se) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_sw,rgnid)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imin(halo)+(j-(jmin(halo)-1)),imax(halo)+(j-(jmin(halo)-1))
                         n=n+1
                         in=j+jmax(halo)-2*jmin(halo)+imin(halo)+1
                         jn=-i+j+imin(halo)+jmax(halo)-jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
          !!
          !
          !!
          if ((ADM_rgn_vnum(ADM_n,rgnid)==4)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,2)==ADM_e) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo-1)
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imin(halo)-halo,imin(halo)-1
                      do j=jmax(halo)+1+(i-(imin(halo)-1)) &
                           ,jmax(halo)+(halo-1)+(i-(imin(halo)-1))
                         n=n+1
                         in=i-j+imax(halo)-imin(halo)+jmax(halo)+2
                         jn=i+imax(halo)-2*imin(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,2)==ADM_s) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imin(halo)-(halo-1),imin(halo)-1
                      do j=jmax(halo)+1,jmax(halo)+1+(i-(imin(halo)-(halo-1)))
                         n=n+1
                         in=i+imax(halo)-imin(halo)+1
                         jn=j-jmax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_n,rgnid,2)==ADM_w) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_n,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmax(halo)+1,jmax(halo)+halo
                      do i=imin(halo)-(halo-1)+(j-(jmax(halo)+1)) &
                           ,imin(halo)+(j-(jmax(halo)+1))
                         n=n+1
                         in=j-jmax(halo)+imin(halo)-1
                         jn=-i+j+imin(halo)-jmax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
          !!
          !
          !!
          if ((ADM_rgn_vnum(ADM_e,rgnid)==4)) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_e,rgnid,2)==ADM_n) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_e,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+2,imax(halo)+halo
                      do j=jmax(halo)+1,jmax(halo)+1+(i-(imax(halo)+2))
                         n=n+1
                         in=j-jmax(halo)-1+imin(halo)
                         jn=-i+j+2*imax(halo)-imin(halo)-jmax(halo)+jmin(halo)+1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_e,rgnid,2)==ADM_w) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_e,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmax(halo)+1,jmax(halo)+halo
                      do i=imax(halo)+1,imax(halo)+halo
                         n=n+1
                         in=i-imax(halo)+imin(halo)-1
                         jn=j-jmax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_e,rgnid,2)==ADM_s) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_e,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmax(halo)+2,jmax(halo)+halo
                      do i=imax(halo)+1,imax(halo)+1+(j-(jmax(halo)+2))
                         n=n+1
                         in=i-j-imax(halo)+2*jmax(halo)+imin(halo)-jmin(halo)+1
                         jn=i-imax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             if (ADM_rgn_etab(ADM_dir,ADM_se,rgnid)==ADM_sw) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(gmax(halo)-gmin(halo)+1)*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_etab(ADM_rid,ADM_se,rgnid)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+1,imax(halo)+halo
                      do j=jmin(halo)+1+(i-(imax(halo)+1)),jmax(halo)+1+(i-(imax(halo)+1))
                         n=n+1
                         in=i-j-imax(halo)+jmax(halo)+imin(halo)
                         jn=i-imax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
                !
             endif
          endif
          !!
          !
          !!
          if (ADM_rgn_vnum(ADM_s,rgnid)==4) then
             !
             if (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,2)==ADM_n) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=(halo-1)*halo/2
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do i=imax(halo)+1,imax(halo)+(halo-1)
                      do j=jmin(halo)-(halo-1)+(i-(imax(halo)+1)),jmin(halo)-1
                         n=n+1
                         in=i-imax(halo)-1+imin(halo)
                         jn=j+jmax(halo)+1-jmin(halo)
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,2)==ADM_e) then
                if (halo>=2) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*(halo-1)
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207) removed by teraim
                   n=0
                   do j=jmin(halo)-halo,jmin(halo)-1
                      do i=imax(halo)+1+(j-(jmin(halo)-1)) &
                           ,imax(halo)+(halo-1)+(j-(jmin(halo)-1))
                         n=n+1
                         in=j+jmax(halo)-2*jmin(halo)+imin(halo)+1
                         jn=-i+j+imax(halo)+jmax(halo)-jmin(halo)+2
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             elseif (ADM_rgn_vtab(ADM_dir,ADM_s,rgnid,2)==ADM_w) then
                if (halo>=1) then
                   m=m+1
                   rsize_r2r(m,halo,rgnid)=halo*halo
                   sourceid_r2r(m,halo,rgnid)=ADM_rgn_vtab(ADM_rid,ADM_s,rgnid,2)
                   !mrecv_r2r(sourceid_r2r(m,halo,rgnid),halo,rgnid)=m !(20101207)removed by teraim
                   n=0
                   do i=imax(halo)+1,imax(halo)+halo
                      do j=jmin(halo)-(halo-1)+(i-(imax(halo)+1)) &
                           ,jmin(halo)+(i-(imax(halo)+1))
                         n=n+1
                         in=i-j-imax(halo)+jmin(halo)+imin(halo)-1
                         jn=i-imax(halo)+jmin(halo)-1
                         rlist_r2r(n,m,halo,rgnid)=suf(i,j,gall(halo))
                         qlist_r2r(n,m,halo,rgnid)=suf(in,jn,gall(halo))
                      enddo
                   enddo
                endif
             endif
             !
          endif
          !
          maxcommrecv_r2r(halo,rgnid)=m
          !
          if ((ADM_rgn_vnum(ADM_w,rgnid)==3)) then
             if ((ADM_rgn_etab(ADM_dir,ADM_nw,rgnid)==ADM_ne)) then
                n_hemisphere_copy(ADM_w,halo,rgnid)=1
             elseif ((ADM_rgn_etab(ADM_dir,ADM_nw,rgnid)==ADM_se)) then
                s_hemisphere_copy(ADM_w,halo,rgnid)=1
             endif
          endif
          if ((ADM_rgn_vnum(ADM_n,rgnid)==5)) then
             n_hemisphere_copy(ADM_n,halo,rgnid)=1
          endif
          if ((ADM_rgn_vnum(ADM_s,rgnid)==3)) then
             n_hemisphere_copy(ADM_s,halo,rgnid)=1
          endif
          if ((ADM_rgn_vnum(ADM_e,rgnid)==3)) then
             if ((ADM_rgn_etab(ADM_dir,ADM_ne,rgnid)==ADM_nw)) then
                n_hemisphere_copy(ADM_e,halo,rgnid)=1
             elseif ((ADM_rgn_etab(ADM_dir,ADM_ne,rgnid)==ADM_sw)) then
                s_hemisphere_copy(ADM_e,halo,rgnid)=1
             endif
          endif
          if ((ADM_rgn_vnum(ADM_s,rgnid)==5)) then
             s_hemisphere_copy(ADM_s,halo,rgnid)=1
          endif
          if ((ADM_rgn_vnum(ADM_n,rgnid)==3)) then
             s_hemisphere_copy(ADM_n,halo,rgnid)=1
          endif
          !
       enddo !loop l
       !
       !(20101207) removed by teraim
       !do rrgnid=1,ADM_rgn_nmax
       !   do srgnid=1,ADM_rgn_nmax
       !      sendtag_r2r(rrgnid,halo,srgnid)=rrgnid+ADM_rgn_nmax*(srgnid-1)
       !      recvtag_r2r(srgnid,halo,rrgnid)=sendtag_r2r(rrgnid,halo,srgnid)
!      !       write(*,*) 'sendtag_r2r',ADM_prc_me,rrgnid,srgnid,sendtag_r2r(rrgnid,halo,srgnid)
       !   enddo
       !enddo
       !
    enddo !loop halo



    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf3D(:,:,l) = rsize_r2r(:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf3D,                                     &
                        max_comm_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                   &
                        rsize_r2r,                                     &
                        max_comm_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                   &
                        ADM_COMM_world,                                &
                        ierr                                           )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf3D(:,:,l) = sourceid_r2r(:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf3D,                                     &
                        max_comm_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                   &
                        sourceid_r2r,                                  &
                        max_comm_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                   &
                        ADM_COMM_world,                                &
                        ierr                                           )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf2D(:,l) = maxcommrecv_r2r(:,rgnid)
    enddo

    call MPI_Allgather( tempbuf2D,                        &
                        halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                      &
                        maxcommrecv_r2r,                  &
                        halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                      &
                        ADM_COMM_world,                   &
                        ierr                              )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf4D(:,:,:,l) = rlist_r2r(:,:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf4D,                                                      &
                        max_comm_r2r*max_datasize_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                                    &
                        rlist_r2r,                                                      &
                        max_comm_r2r*max_datasize_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                                    &
                        ADM_COMM_world,                                                 &
                        ierr                                                            )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf4D(:,:,:,l) = qlist_r2r(:,:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf4D,                                                      &
                        max_comm_r2r*max_datasize_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                                    &
                        qlist_r2r,                                                      &
                        max_comm_r2r*max_datasize_r2r*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                                                    &
                        ADM_COMM_world,                                                 &
                        ierr                                                            )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf3D2(:,:,l) = n_hemisphere_copy(:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf3D2,                         &
                        4*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                        &
                        n_hemisphere_copy,                  &
                        4*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                        &
                        ADM_COMM_world,                     &
                        ierr                                )

    do l = 1, ADM_prc_rnum(ADM_prc_me)
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       tempbuf3D2(:,:,l) = s_hemisphere_copy(:,:,rgnid)
    enddo

    call MPI_Allgather( tempbuf3D2,                         &
                        4*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                        &
                        s_hemisphere_copy,                  &
                        4*halomax*ADM_prc_rnum(ADM_prc_me), &
                        MPI_INTEGER,                        &
                        ADM_COMM_world,                     &
                        ierr                                )

    call MPI_Barrier(ADM_COMM_world,ierr)

    do halo=1,halomax
       do ls=1,ADM_prc_rnum(ADM_prc_me)
          srgnid=ADM_prc_tab(ls,ADM_prc_me)
          ms=0
          do lr=1,ADM_rgn_nmax
             rrgnid=lr
             do mr=1,maxcommrecv_r2r(halo,rrgnid)
                if (srgnid==sourceid_r2r(mr,halo,rrgnid)) then
                   ms=ms+1
                   !
                   !(20101207)added by teraim
                   if(ADM_rgn2prc(srgnid)==ADM_prc_me) then
                     msend_r2r(rrgnid,halo,srgnid)=ms
                   else
                     write(*,*)"This process is abort because irregular access in msend_r2r."
                     exit
                   endif
                   !
                   destid_r2r(ms,halo,srgnid)=rrgnid
                   ssize_r2r(ms,halo,srgnid)=rsize_r2r(mr,halo,rrgnid)
                   do n=1,rsize_r2r(mr,halo,rrgnid)
                      slist_r2r(n,ms,halo,srgnid)=qlist_r2r(n,mr,halo,rrgnid)
                   enddo
                endif
             enddo
          enddo
          maxcommsend_r2r(halo,srgnid)=ms
       enddo
    enddo !loop halo
    !
    call MPI_Barrier(ADM_COMM_world,ierr)
    do l=1,ADM_rgn_nmax
       call MPI_bcast(                  &
            destid_r2r(1,1,l),          &
            max_comm_r2r*halomax,       &
            MPI_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_COMM_world,             &
            ierr)
    enddo
    do l=1,ADM_rgn_nmax
       call MPI_bcast(                  &
            ssize_r2r(1,1,l),           &
            max_comm_r2r*halomax,       &
            MPI_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_COMM_world,             &
            ierr)
    enddo
    !(20101207)removed by teraim
    !do l=1,ADM_rgn_nmax
    !   call MPI_bcast(                  &
    !        msend_r2r(1,1,l),            &
    !        ADM_rgn_nmax*halomax,  &
    !        MPI_integer,                &
    !        prc_tab_rev(ptr_prcid,l)-1, &
    !        ADM_COMM_world,             &
    !        ierr)
    !enddo
    do l=1,ADM_rgn_nmax
       call MPI_bcast(                             &
            slist_r2r(1,1,1,l),                    &
            max_comm_r2r*max_datasize_r2r*halomax, &
            MPI_integer,                           &
            prc_tab_rev(ptr_prcid,l)-1,            &
            ADM_COMM_world,                        &
            ierr)
    enddo
    do l=1,ADM_rgn_nmax
       call MPI_bcast(                  &
            maxcommsend_r2r(1,l),       &
            1*halomax,                  &
            MPI_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_COMM_world,             &
            ierr)
    enddo
    !
    allocate(sendbuf_p2r(kmax*max_varmax*2, &
         ADM_rgn_nmax_pl))
    allocate(recvbuf_p2r(kmax*max_varmax*2, &
         ADM_rgn_nmax_pl))
    allocate(recvtag_p2r(max_comm_p2r,ADM_npl:ADM_spl))
    allocate(sendtag_p2r(max_comm_p2r,ADM_npl:ADM_spl))
    do pl=ADM_npl,ADM_spl
       do p=1,ADM_vlink_nmax
          recvtag_p2r(p,pl)=ADM_rgn_nmax*ADM_rgn_nmax+p+ADM_vlink_nmax*(pl-1)
          sendtag_p2r(p,pl)=ADM_rgn_nmax*ADM_rgn_nmax+p+ADM_vlink_nmax*(pl-1)

!          write(*,*) 'sendtag_p2r',ADM_prc_me,p,pl,halo,sendtag_p2r(p,pl)

       enddo
    enddo
    !
    allocate(clist(max_varmax))
    !
    call MPI_Barrier(ADM_COMM_world,ierr)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  re-setup comm_table !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rank_me=ADM_prc_me-1
    max_comm_prc=min(ADM_prc_all,max_comm_r2r*ADM_lall+2*max_comm_r2p)
    !
    allocate(n_nspl(ADM_npl:ADM_spl,halomax))
    do halo=1,halomax
       n_nspl(ADM_npl,halo)=suf(imin(halo)+0,jmax(halo)+1,gall(halo))
       n_nspl(ADM_spl,halo)=suf(imax(halo)+1,jmin(halo)+0,gall(halo))
    enddo
    !
    allocate(temp_sendorder(0:ADM_prc_all-1,halomax))
    allocate(temp_recvorder(0:ADM_prc_all-1,halomax))
    !
    !--------------------------------------------------
    allocate(romax(halomax))
    allocate(somax(halomax))
    allocate(sourcerank(max_comm_prc,halomax))
    allocate(destrank(max_comm_prc,halomax))
    allocate(rsize(max_comm_prc,halomax))
    allocate(ssize(max_comm_prc,halomax))
    romax(:)=0
    somax(:)=0
    sourcerank(:,:)=-1
    destrank(:,:)=-1
    rsize(:,:)=0
    ssize(:,:)=0
    !--------------------------------------------------
    !
    maxn=((gmax(halomax)-gmin(halomax)+1)+2)*halomax
    maxm=max_comm_r2r+1
    maxl=ADM_lall+2
    !----
    maxn_pl=halomax*(halomax+1)/2
    maxm_pl=ADM_vlink_nmax
    maxl_pl=(ADM_spl-ADM_npl+1)
    !----
    maxn_r2r=(gmax(halomax)-gmin(halomax)+1)*halomax
    maxm_r2r=max_comm_r2r
    maxl_r2r=ADM_lall
    !----
    maxn_r2p=halomax*(halomax+1)/2
    maxm_r2p=ADM_vlink_nmax
    maxl_r2p=(ADM_spl-ADM_npl+1)
    !----
    maxn_p2r=1
    maxm_p2r=ADM_vlink_nmax
    maxl_p2r=(ADM_spl-ADM_npl+1)
    !----
    maxn_sgp=halomax
    maxm_sgp=4
    maxl_sgp=12
    !
    !--------------------------------------------------
    !  for send
    !--------------------------------------------------
    allocate(nsmax(max_comm_prc,halomax))
    allocate(sendinfo(elemsize_comm,maxm*maxl,max_comm_prc,halomax))
    allocate(sendlist(maxn,maxm*maxl,max_comm_prc,halomax))
    nsmax(:,:)=0
    sendinfo(:,:,:,:)=0
    sendlist(:,:,:,:)=0
    allocate(nsmax_pl(max_comm_prc,halomax))
    allocate(sendinfo_pl(elemsize_comm,maxm_pl*maxl_pl,max_comm_prc,halomax))
    allocate(sendlist_pl(maxn_pl,maxm_pl*maxl_pl,max_comm_prc,halomax))
    nsmax_pl(:,:)=0
    sendinfo_pl(:,:,:,:)=0
    sendlist_pl(:,:,:,:)=0
    !--------------------------------------------------
    !
    !--------------------------------------------------
    !  for copy
    !--------------------------------------------------
    allocate(ncmax_r2r(halomax))
    allocate(copyinfo_r2r(elemsize_copy,maxm_r2r*maxl_r2r,halomax))
    allocate(recvlist_r2r(maxn_r2r,maxm_r2r*maxl_r2r,halomax))
    allocate(sendlist_r2r(maxn_r2r,maxm_r2r*maxl_r2r,halomax))
    ncmax_r2r(:)=0
    copyinfo_r2r(:,:,:)=0
    recvlist_r2r(:,:,:)=0
    sendlist_r2r(:,:,:)=0
    !--------------------------------------------------
    allocate(ncmax_r2p(halomax))
    allocate(copyinfo_r2p(elemsize_copy,maxm_r2p*maxl_r2p,halomax))
    allocate(recvlist_r2p(maxn_r2p,maxm_r2p*maxl_r2p,halomax))
    allocate(sendlist_r2p(maxn_r2p,maxm_r2p*maxl_r2p,halomax))
    ncmax_r2p(:)=0
    copyinfo_r2p(:,:,:)=0
    recvlist_r2p(:,:,:)=0
    sendlist_r2p(:,:,:)=0
    !--------------------------------------------------
    allocate(ncmax_p2r(halomax))
    allocate(copyinfo_p2r(elemsize_copy,maxm_p2r*maxl_p2r,halomax))
    allocate(recvlist_p2r(maxn_p2r,maxm_p2r*maxl_p2r,halomax))
    allocate(sendlist_p2r(maxn_p2r,maxm_p2r*maxl_p2r,halomax))
    ncmax_p2r(:)=0
    copyinfo_p2r(:,:,:)=0
    recvlist_p2r(:,:,:)=0
    sendlist_p2r(:,:,:)=0
    !--------------------------------------------------
    !
    !--------------------------------------------------
    !  for recv
    !--------------------------------------------------
    allocate(nrmax(max_comm_prc,halomax))
    allocate(recvinfo(elemsize_comm,maxm*maxl,max_comm_prc,halomax))
    allocate(recvlist(maxn,maxm*maxl,max_comm_prc,halomax))
    nrmax(:,:)=0
    recvinfo(:,:,:,:)=0
    recvlist(:,:,:,:)=0
    allocate(nrmax_pl(max_comm_prc,halomax))
    allocate(recvinfo_pl(elemsize_comm,maxm_pl*maxl_pl,max_comm_prc,halomax))
    allocate(recvlist_pl(maxn_pl,maxm_pl*maxl_pl,max_comm_prc,halomax))
    nrmax_pl(:,:)=0
    recvinfo_pl(:,:,:,:)=0
    recvlist_pl(:,:,:,:)=0
    !--------------------------------------------------
    allocate(temp_dest_rgn(maxm*maxl,max_comm_prc,halomax))
    allocate(temp_src_rgn(maxm*maxl,max_comm_prc,halomax))
    allocate(temp_dest_rgn_pl(maxm_pl*maxl_pl,max_comm_prc,halomax))
    allocate(temp_src_rgn_pl(maxm_pl*maxl_pl,max_comm_prc,halomax))
    temp_dest_rgn(:,:,:)=0
    temp_dest_rgn_pl(:,:,:)=0
    temp_src_rgn(:,:,:)=0
    temp_src_rgn_pl(:,:,:)=0
    !--------------------------------------------------
    !
    do halo=1,halomax
       do l=1,ADM_lall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          do m=1,maxcommrecv_r2r(halo,rgnid)
             srank=prc_tab_rev(ptr_prcid,sourceid_r2r(m,halo,rgnid))-1
             if (srank/=rank_me) then
                ck=0
                loop_ro1:do ro=1,romax(halo)
                   if (srank==sourcerank(ro,halo)) exit loop_ro1
                   ck=ck+1
                enddo loop_ro1
                if (ck==romax(halo)) then
                   romax(halo)=romax(halo)+1
                   ro=romax(halo)
                   sourcerank(ro,halo)=srank
                   temp_recvorder(srank,halo)=ro
                endif
                ro=temp_recvorder(srank,halo)
                nrmax(ro,halo)=nrmax(ro,halo)+1
                nr=nrmax(ro,halo)
                recvinfo(SIZE_COMM,nr,ro,halo)=rsize_r2r(m,halo,rgnid)
                recvinfo(LRGNID_COMM,nr,ro,halo)=l
                temp_src_rgn(nr,ro,halo)=sourceid_r2r(m,halo,rgnid)
                rs=recvinfo(SIZE_COMM,nr,ro,halo)
                rsize(ro,halo)=rsize(ro,halo)+rs
                do n=1,rs
                   recvlist(n,nr,ro,halo)=rlist_r2r(n,m,halo,rgnid)
                enddo
             else
                ncmax_r2r(halo)=ncmax_r2r(halo)+1
                nc=ncmax_r2r(halo)
                copyinfo_r2r(SIZE_COPY,nc,halo)=rsize_r2r(m,halo,rgnid)
                copyinfo_r2r(LRGNID_COPY,nc,halo)=l
                copyinfo_r2r(SRC_LRGNID_COPY,nc,halo) &
                     =prc_tab_rev(ptr_lrgnid,sourceid_r2r(m,halo,rgnid))
                cs=copyinfo_r2r(SIZE_COPY,nc,halo)
                srgnid=sourceid_r2r(m,halo,rgnid)
                do n=1,cs
                   recvlist_r2r(n,nc,halo)=rlist_r2r(n,m,halo,rgnid)
                   !
                   !(20101207)added by teraim
                   if(ADM_rgn2prc(srgnid)==ADM_prc_me) then
                     sendlist_r2r(n,nc,halo)=slist_r2r(n,msend_r2r(rgnid,halo,srgnid),halo,srgnid)
                   else
                     write(*,*)"This process is abort because irregular access is msend_r2r."
                     exit
                   endif
                   !
                enddo
             endif
          enddo !loop m
          !enddo !loop l
          !!
          !do l=1,ADM_lall
          !  rgnid=ADM_prc_tab(l,ADM_prc_me)
          do m=1,maxcommsend_r2r(halo,rgnid)
             drank=prc_tab_rev(ptr_prcid,destid_r2r(m,halo,rgnid))-1
             if (drank/=rank_me) then
                ck=0
                loop_so1:do so=1,somax(halo)
                   if (drank==destrank(so,halo)) exit loop_so1
                   ck=ck+1
                enddo loop_so1
                if (ck==somax(halo)) then
                   somax(halo)=somax(halo)+1
                   so=somax(halo)
                   destrank(so,halo)=drank
                   temp_sendorder(drank,halo)=so
                endif
                so=temp_sendorder(drank,halo)
                nsmax(so,halo)=nsmax(so,halo)+1
                ns=nsmax(so,halo)
                sendinfo(SIZE_COMM,ns,so,halo)=ssize_r2r(m,halo,rgnid)
                sendinfo(LRGNID_COMM,ns,so,halo)=l
                temp_dest_rgn(ns,so,halo)=destid_r2r(m,halo,rgnid)
                ss=sendinfo(SIZE_COMM,ns,so,halo)
                ssize(so,halo)=ssize(so,halo)+ss
                do n=1,ss
                   sendlist(n,ns,so,halo)=slist_r2r(n,m,halo,rgnid)
                enddo
             endif
          enddo !loop m
       enddo !loop l
       !enddo !loop halo
       !do halo=1,halomax
       if(comm_pl) call re_setup_pl_comm_info ! T.Ohno 110721
    enddo !loop halo
    deallocate(temp_sendorder)
    deallocate(temp_recvorder)
    !
    !allocate(temp_sb(ADM_rgn_nmax+2,halomax,ADM_rgn_nmax+2)) !(20101207) removed by teraim
    allocate(tsb(somax(halomax)))
    !temp_sb(:,:,:)=0 !(20101207) removed by teraim

    call init_tempsb !(20101207) added by teraim

    do halo=1,halomax
       tsb(:)=0
       do so=1,somax(halo)
          do ns=1,nsmax(so,halo)
             ss=sendinfo(SIZE_COMM,ns,so,halo)
             srgnid=ADM_prc_tab(sendinfo(LRGNID_COMM,ns,so,halo),ADM_prc_me)
             rrgnid=temp_dest_rgn(ns,so,halo)
             sendinfo(BASE_COMM,ns,so,halo)=tsb(so)
             !temp_sb(rrgnid,halo,srgnid)=tsb(so) !(20101207)removed by teraim
             call add_tempsb(rrgnid, srgnid, tsb(so)) !(20101207)added by teraim
             tsb(so)=tsb(so)+ss
          enddo
          do ns=1,nsmax_pl(so,halo)
             ss=sendinfo_pl(SIZE_COMM,ns,so,halo)
             pl=sendinfo_pl(LRGNID_COMM,ns,so,halo)
             srgnid=ADM_rgn_nmax+pl
             rrgnid=temp_dest_rgn_pl(ns,so,halo)
             sendinfo_pl(BASE_COMM,ns,so,halo)=tsb(so)
             !temp_sb(rrgnid,halo,srgnid)=tsb(so) !(20101207)removed by teraim
             call add_tempsb(rrgnid, srgnid, tsb(so)) !(20101207)added by teraim
             tsb(so)=tsb(so)+ss
          enddo
       enddo
    enddo
    deallocate(tsb)
    !
    !(20101207)removed by teraim
    !call MPI_Barrier(ADM_COMM_world,ierr)
    !do l=1,ADM_rgn_nmax
    !   call MPI_bcast(                  &
    !        temp_sb(1,1,l),        &
    !        (ADM_rgn_nmax+2)*halomax,      &
    !        MPI_integer,                &
    !        prc_tab_rev(ptr_prcid,l)-1, &
    !        ADM_COMM_world,             &
    !        ierr)
    !enddo
    !do pl=ADM_npl,ADM_spl
    !   call MPI_bcast(                  &
    !        temp_sb(1,1,ADM_rgn_nmax+pl),       &
    !        (ADM_rgn_nmax+2)*halomax,      &
    !        MPI_integer,                &
    !        ADM_prc_nspl(pl)-1,         &
    !        ADM_COMM_world,             &
    !        ierr)
    !enddo
    !call MPI_barrier(ADM_COMM_world,ierr)
    !
    !(20101207)added by teraim
    call MPI_Barrier(ADM_COMM_world,ierr)
    do l=1,ADM_rgn_nmax
       call MPI_bcast(                  &
            tempsb(l)%col,              &
            max_size,                   &
            MPI_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_COMM_world,             &
            ierr                        )
       call MPI_bcast(                  &
            tempsb(l)%val,              &
            max_size,                   &
            MPI_integer,                &
            prc_tab_rev(ptr_prcid,l)-1, &
            ADM_COMM_world,             &
            ierr                        )
    enddo
    if(comm_pl) then ! T.Ohno 110721
      do pl=ADM_npl,ADM_spl
         call MPI_bcast(                    &
              tempsb(ADM_rgn_nmax+pl)%col,  &
              max_size,                     &
              MPI_integer,                  &
              ADM_prc_nspl(pl)-1,           &
              ADM_COMM_world,               &
              ierr                          )
         call MPI_bcast(                    &
              tempsb(ADM_rgn_nmax+pl)%val,  &
              max_size,                     &
              MPI_integer,                  &
              ADM_prc_nspl(pl)-1,           &
              ADM_COMM_world,               &
              ierr                          )
      enddo
    endif ! T.Ohno 110721
    call MPI_Barrier(ADM_COMM_world,ierr)
    !
    !call show_tempsb !(20101209) added by teraim
    !
    do halo=1,halomax
       do ro=1,romax(halo)
          do nr=1,nrmax(ro,halo)
             rrgnid=ADM_prc_tab(recvinfo(LRGNID_COMM,nr,ro,halo),ADM_prc_me)
             srgnid=temp_src_rgn(nr,ro,halo)
             !recvinfo(BASE_COMM,nr,ro,halo)=temp_sb(rrgnid,halo,srgnid) !(20101207)removed by teraim
             !(20101207) added by teraim
             call get_tempsb(rrgnid,srgnid,ret)
             recvinfo(BASE_COMM,nr,ro,halo)=ret
          enddo
          do nr=1,nrmax_pl(ro,halo)
             pl=recvinfo_pl(LRGNID_COMM,nr,ro,halo)
             rrgnid=pl+ADM_rgn_nmax
             srgnid=temp_src_rgn_pl(nr,ro,halo)
             !recvinfo_pl(BASE_COMM,nr,ro,halo)=temp_sb(rrgnid,halo,srgnid) !(20101207)removed by teraim
             !(20101207) added by teraim
             call get_tempsb(rrgnid,srgnid,ret)
             recvinfo_pl(BASE_COMM,nr,ro,halo)=ret
          enddo
       enddo !loop ro
    enddo !loop halo
    deallocate(temp_dest_rgn)
    deallocate(temp_dest_rgn_pl)
    deallocate(temp_src_rgn)
    deallocate(temp_src_rgn_pl)
    !deallocate(temp_sb) !(20101207)removed by teraim
    call finalize_tempsb !(20101207)added by teraim
    !
    allocate(recvtag(romax(halomax),halomax))
    allocate(sendtag(somax(halomax),halomax))
    recvtag(:,:)=-1
    sendtag(:,:)=-1
    do halo=1,halomax
       do ro=1,romax(halo)
          recvtag(ro,halo)=rank_me
       enddo
       do so=1,somax(halo)
          sendtag(so,halo)=destrank(so,halo)
       enddo
    enddo
!    maxdatasize=(max_comm_r2r*(gmax(halomax)-gmin(halomax)+1)+2*max_comm_r2p*(halomax+1)/2)*halomax*kmax*max_varmax
!    maxdatasize=(maxn_r2r*maxm_r2r*maxl_r2r+maxn_r2p*maxm_r2p*maxl_r2p+maxn_p2r*maxm_p2r*maxl_p2r)*kmax*max_varmax
    maxdatasize_s=0
    do so=1,somax(halomax)
      maxdatasize_s=maxdatasize_s+ssize(so,halomax)*kmax*max_varmax
    enddo
    maxdatasize_r=0
    do ro=1,romax(halomax)
      maxdatasize_r=maxdatasize_r+rsize(ro,halomax)*kmax*max_varmax
    enddo
    allocate(recvbuf(maxdatasize_r,romax(halomax)))
    allocate(sendbuf(maxdatasize_s,somax(halomax)))
    recvbuf(:,:)=0.0_RP
    sendbuf(:,:)=0.0_RP

!!    allocate(comm_dbg_recvbuf(maxdatasize_r,romax(halomax),2)) !iga
!!    allocate(comm_dbg_sendbuf(maxdatasize_s,somax(halomax),2)) !iga
!!    comm_dbg_recvbuf=CNST_UNDEF !iga
!!    comm_dbg_sendbuf=CNST_UNDEF !iga

    !
    allocate(ncmax_sgp(halomax))
    allocate(copyinfo_sgp(elemsize_copy,maxm_sgp*maxl_sgp,halomax))
    allocate(recvlist_sgp(maxn_sgp,maxm_sgp*maxl_sgp,halomax))
    allocate(sendlist_sgp(maxn_sgp,maxm_sgp*maxl_sgp,halomax))
    ncmax_sgp(:)=0
    copyinfo_sgp(:,:,:)=0
    recvlist_sgp(:,:,:)=0
    sendlist_sgp(:,:,:)=0
    do halo=1,halomax
       ncmax_sgp(halo)=0
       do l=1,ADM_lall
          rgnid=ADM_prc_tab(l,ADM_prc_me)
          if (n_hemisphere_copy(ADM_w,halo,rgnid)==1) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmin(halo)-n,gmin(halo)-n,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmin(halo),gmin(halo)-n,gall(halo))
             enddo
          endif
          if ((n_hemisphere_copy(ADM_n,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmin(halo),gmax(halo)+n+1,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmin(halo)-n,gmax(halo)+1,gall(halo))
             enddo
          endif
          if ((n_hemisphere_copy(ADM_e,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmax(halo)+1,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmax(halo)+n+1,gall(halo))
             enddo
          endif
          if ((n_hemisphere_copy(ADM_s,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmax(halo)+1,gmin(halo)-n,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmin(halo),gall(halo))
             enddo
          endif
          if (s_hemisphere_copy(ADM_w,halo,rgnid)==1) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmin(halo),gmin(halo)-n,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmin(halo)-n,gmin(halo)-n,gall(halo))
             enddo
          endif
          if ((s_hemisphere_copy(ADM_n,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmin(halo)-n,gmax(halo)+1,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmin(halo),gmax(halo)+n+1,gall(halo))
             enddo
          endif
          if ((s_hemisphere_copy(ADM_e,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmax(halo)+1,gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmax(halo)+n+1,gall(halo))
             enddo
          endif
          if ((s_hemisphere_copy(ADM_s,halo,rgnid)==1)) then
             ncmax_sgp(halo)=ncmax_sgp(halo)+1
             nc=ncmax_sgp(halo)
             copyinfo_sgp(SIZE_COPY,nc,halo)=halo-1
             copyinfo_sgp(LRGNID_COPY,nc,halo)=l
             copyinfo_sgp(SRC_LRGNID_COPY,nc,halo)=l
             cs=copyinfo_sgp(SIZE_COPY,nc,halo)
             do n=1,cs
                recvlist_sgp(n,nc,halo)=suf(gmax(halo)+n+1,gmin(halo),gall(halo))
                sendlist_sgp(n,nc,halo)=suf(gmax(halo)+1,gmin(halo)-n,gall(halo))
             enddo
          endif
       enddo !loop l
    enddo !loop halo
    !
    !-- for output_info  ---
    allocate( src_rank_all(max_comm_prc,halomax,ADM_prc_all))
    allocate(dest_rank_all(max_comm_prc,halomax,ADM_prc_all))
    src_rank_all(:,:,:)=-1
    dest_rank_all(:,:,:)=-1
    src_rank_all(:,:,ADM_prc_me)=sourcerank(:,:)
    dest_rank_all(:,:,ADM_prc_me)=destrank(:,:)
    call MPI_Barrier(ADM_COMM_world,ierr)
    do l=1,ADM_prc_all
       call MPI_bcast(                  &
            src_rank_all(1,1,l),        &
            max_comm_prc*halomax,       &
            MPI_integer,                &
            l-1,                        &
            ADM_COMM_world,             &
            ierr)
       call MPI_bcast(                  &
            dest_rank_all(1,1,l),       &
            max_comm_prc*halomax,       &
            MPI_integer,                &
            l-1,                        &
            ADM_COMM_world,             &
            ierr)
    enddo
    !
    call MPI_Barrier(ADM_COMM_WORLD,ierr)
    !
    !--- output for debug
    if(present(debug)) then
       if(debug) call output_info
    endif
    !
    ! <== iga for dbg 090917
    if (opt_comm_dbg) then
!       dbg_sendbuf_init = -1d66  * (ADM_prc_me+1000)
       dbg_sendbuf_init = -888E+30_RP
       dbg_recvbuf_init = -777E+30_RP
       allocate(dbg_areq_save(2*(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4),4))
       dbg_areq_save(:,:) = -999 ! [Add] 12/03/26 T.Seiki
    endif
    ! iga for dbg 090916 ==>
    contains
    !
    subroutine re_setup_pl_comm_info
       do pl=ADM_npl,ADM_spl
          if (ADM_prc_me==ADM_prc_nspl(pl)) then
             do p=1,maxcommrecv_r2p(pl,halo)
                srank=source_prc_r2p(p,pl,halo)-1
                if (srank/=rank_me) then
                   ck=0
                   loop_ro2:do ro=1,romax(halo)
                      if (srank==sourcerank(ro,halo)) exit loop_ro2
                      ck=ck+1
                   enddo loop_ro2
                   if (ck==romax(halo)) then
                      romax(halo)=romax(halo)+1
                      ro=romax(halo)
                      sourcerank(ro,halo)=srank
                      temp_recvorder(srank,halo)=ro
                   endif
                   ro=temp_recvorder(srank,halo)
                   nrmax_pl(ro,halo)=nrmax_pl(ro,halo)+1
                   nr=nrmax_pl(ro,halo)
                   recvinfo_pl(SIZE_COMM,nr,ro,halo)=rsize_r2p(p,pl,halo)
                   recvinfo_pl(LRGNID_COMM,nr,ro,halo)=pl
                   temp_src_rgn_pl(nr,ro,halo)=ADM_prc_tab(source_rgn_r2p(p,pl,halo),srank+1)
                   rs=recvinfo_pl(SIZE_COMM,nr,ro,halo)
                   rsize(ro,halo)=rsize(ro,halo)+rs
                   do n=1,rs
                      recvlist_pl(n,nr,ro,halo)=rlist_r2p(n,p,pl,halo)
                   enddo
                else
                   ncmax_r2p(halo)=ncmax_r2p(halo)+1
                   nc=ncmax_r2p(halo)
                   copyinfo_r2p(SIZE_COPY,nc,halo)=rsize_r2p(p,pl,halo)
                   copyinfo_r2p(LRGNID_COPY,nc,halo)=pl
                   copyinfo_r2p(SRC_LRGNID_COPY,nc,halo)=source_rgn_r2p(p,pl,halo)
                   cs=copyinfo_r2p(SIZE_COPY,nc,halo)
                   do n=1,cs
                      recvlist_r2p(n,nc,halo)=rlist_r2p(n,p,pl,halo)
                      sendlist_r2p(n,nc,halo)=slist_r2p(n,p,pl,halo)
                   enddo
                endif
             enddo !loop p
             !
             do p=1,ADM_vlink_nmax
                rgnid=ADM_rgn_vtab_pl(ADM_rid,pl,p)
                drank=prc_tab_rev(ptr_prcid,rgnid)-1
                if (drank/=rank_me) then
                   ck=0
                   loop_so2:do so=1,somax(halo)
                      if (drank==destrank(so,halo)) exit loop_so2
                      ck=ck+1
                   enddo loop_so2
                   if (ck==somax(halo)) then
                      somax(halo)=somax(halo)+1
                      so=somax(halo)
                      destrank(so,halo)=drank
                      temp_sendorder(drank,halo)=so
                   endif
                   so=temp_sendorder(drank,halo)
                   nsmax_pl(so,halo)=nsmax_pl(so,halo)+1
                   ns=nsmax_pl(so,halo)
                   sendinfo_pl(SIZE_COMM,ns,so,halo)=1
                   sendinfo_pl(LRGNID_COMM,ns,so,halo)=pl
                   temp_dest_rgn_pl(ns,so,halo)=rgnid
                   ss=sendinfo_pl(SIZE_COMM,ns,so,halo)
                   ssize(so,halo)=ssize(so,halo)+ss
                   do n=1,ss
                      sendlist_pl(n,ns,so,halo)=ADM_gslf_pl
                   enddo
                endif
             enddo !loop p
          endif
          !
          do p=1,ADM_vlink_nmax
             rgnid=ADM_rgn_vtab_pl(ADM_rid,pl,p)
             drank=prc_tab_rev(ptr_prcid,rgnid)-1
             if (rank_me==drank) then
                srank=ADM_prc_nspl(pl)-1
                if (srank/=rank_me) then
                   ck=0
                   loop_ro3:do ro=1,romax(halo)
                      if (srank==sourcerank(ro,halo)) exit loop_ro3
                      ck=ck+1
                   enddo loop_ro3
                   if (ck==romax(halo)) then
                      romax(halo)=romax(halo)+1
                      ro=romax(halo)
                      sourcerank(ro,halo)=srank
                      temp_recvorder(srank,halo)=ro
                   endif
                   ro=temp_recvorder(srank,halo)
                   nrmax(ro,halo)=nrmax(ro,halo)+1
                   nr=nrmax(ro,halo)
                   recvinfo(SIZE_COMM,nr,ro,halo)=1
                   recvinfo(LRGNID_COMM,nr,ro,halo)=prc_tab_rev(ptr_lrgnid,rgnid)
                   temp_src_rgn(nr,ro,halo)=ADM_rgn_nmax+pl
                   rs=recvinfo(SIZE_COMM,nr,ro,halo)
                   rsize(ro,halo)=rsize(ro,halo)+rs
                   do n=1,rs
                      recvlist(n,nr,ro,halo)=n_nspl(pl,halo)
                   enddo
                else
                   ncmax_p2r(halo)=ncmax_p2r(halo)+1
                   nc=ncmax_p2r(halo)
                   copyinfo_p2r(SIZE_COPY,nc,halo)=1
                   copyinfo_p2r(LRGNID_COPY,nc,halo)=prc_tab_rev(ptr_lrgnid,rgnid)
                   copyinfo_p2r(SRC_LRGNID_COPY,nc,halo)=pl
                   cs=copyinfo_p2r(SIZE_COPY,nc,halo)
                   do n=1,cs
                      recvlist_p2r(n,nc,halo)=n_nspl(pl,halo)
                      sendlist_p2r(n,nc,halo)=ADM_gslf_pl
                   enddo
                endif
             endif
          enddo !loop p
          !
          do p=1,maxcommsend_r2p(pl,halo)
             srank=source_prc_r2p(p,pl,halo)-1
             if (rank_me==srank) then
                rgnid=ADM_rgn_vtab_pl(ADM_rid,pl,p)
                drank=ADM_prc_nspl(pl)-1
                if (drank/=rank_me) then
                   ck=0
                   loop_so3:do so=1,somax(halo)
                      if (drank==destrank(so,halo)) exit loop_so3
                      ck=ck+1
                   enddo loop_so3
                   if (ck==somax(halo)) then
                      somax(halo)=somax(halo)+1
                      so=somax(halo)
                      destrank(so,halo)=drank
                      temp_sendorder(drank,halo)=so
                   endif
                   so=temp_sendorder(drank,halo)
                   nsmax(so,halo)=nsmax(so,halo)+1
                   ns=nsmax(so,halo)
                   sendinfo(SIZE_COMM,ns,so,halo)=ssize_r2p(p,pl,halo)
                   sendinfo(LRGNID_COMM,ns,so,halo)=prc_tab_rev(ptr_lrgnid,rgnid)
                   temp_dest_rgn(ns,so,halo)=ADM_rgn_nmax+pl
                   ss=sendinfo(SIZE_COMM,ns,so,halo)
                   ssize(so,halo)=ssize(so,halo)+ss
                   do n=1,ss
                      sendlist(n,ns,so,halo)=slist_r2p(n,p,pl,halo)
                   enddo
                endif
             endif
          enddo !loop p
       enddo !loop pl
    end subroutine re_setup_pl_comm_info

  end subroutine COMM_setup
  !-----------------------------------------------------------------------------
  subroutine output_info
    use mod_adm, only :    &
         ADM_log_fid,    &
         ADM_prc_all
    implicit none

    integer :: halo
    !
    integer ::  varmax
    integer ::  cmax
    !
    !integer ::  srgnid,rrgnid
    !
    integer ::  ns,nr
    integer ::  so,sl,sb,ss
    integer ::  ro,rl,rb,rs
    integer ::  l,n
    !
    write(ADM_log_fid,*)
    write(ADM_log_fid,*) &
         'msg : sub[output_info]/mod[comm]'
    write(ADM_log_fid,*) &
         'version : comm.f90.test5.2.1_wtime'
    write(ADM_log_fid,*) &
         '---------------------------------------&
         &       commnication table  start       &
         &---------------------------------------'
    !
    varmax=1
    cmax=kmax*varmax
    do halo=1,halomax
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &       halo region =',halo,'           &
            &---------------------------------------'
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &                count                  &
            &---------------------------------------'
       write(ADM_log_fid,*) &
            'romax =',romax(halo) &
            ,'somax =',somax(halo)
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &                send                   &
            &---------------------------------------'
       do so=1,somax(halo)
          write(ADM_log_fid,*) &
               'so =',so   &
               ,'mrank =',rank_me   &
               ,'drank =',destrank(so,halo)
       enddo
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &                recv                   &
            &---------------------------------------'
       do ro=1,romax(halo)
          write(ADM_log_fid,*) &
               'ro =',ro   &
               ,'mrank =',rank_me   &
               ,'srank =',sourcerank(ro,halo)
       enddo
       write(ADM_log_fid,*) &
            '---------------------------------------&
            &                table                   &
            &---------------------------------------'
       do l=1,ADM_prc_all
          do n=1,max_comm_prc
             if (dest_rank_all(n,halo,l)==-1) cycle
             write(ADM_log_fid,*) &
                  'n =',n   &
                  ,'rank =',l-1   &
                  ,'dest_rank =',dest_rank_all(n,halo,l) &
                  ,' src_rank =', src_rank_all(n,halo,l)
          enddo
       enddo
       do so=1,somax(halo)
          do ns=1,nsmax(so,halo)
             ss=sendinfo(SIZE_COMM,ns,so,halo)
             sl=sendinfo(LRGNID_COMM,ns,so,halo)
             sb=sendinfo(BASE_COMM,ns,so,halo)*cmax
             write(ADM_log_fid,*) &
                  'so =',so   &
                  ,'rank =',rank_me   &
                  ,' dest_rank =', destrank(so,halo) &
                  ,' ss =', ss &
                  ,' sl =', sl &
                  ,' sb =', sb
          enddo
          do ns=1,nsmax_pl(so,halo)
             ss=sendinfo_pl(SIZE_COMM,ns,so,halo)
             sl=sendinfo_pl(LRGNID_COMM,ns,so,halo)
             sb=sendinfo_pl(BASE_COMM,ns,so,halo)*cmax
             write(ADM_log_fid,*) &
                  'so =',so   &
                  ,'rank =',rank_me   &
                  ,' dest_rank =', destrank(so,halo) &
                  ,' ss =', ss &
                  ,' sl =', sl &
                  ,' sb =', sb
          enddo
       enddo !loop so
       do ro=1,romax(halo)
          do nr=1,nrmax(ro,halo)
             rs=recvinfo(SIZE_COMM,nr,ro,halo)
             rl=recvinfo(LRGNID_COMM,nr,ro,halo)
             rb=recvinfo(BASE_COMM,nr,ro,halo)*cmax
             write(ADM_log_fid,*) &
                  'ro =',ro   &
                  ,'rank =',rank_me   &
                  ,' src_rank =', sourcerank(ro,halo) &
                  ,' rs =', rs &
                  ,' rl =', rl &
                  ,' rb =', rb
          enddo
          do nr=1,nrmax_pl(ro,halo)
             rs=recvinfo_pl(SIZE_COMM,nr,ro,halo)
             rl=recvinfo_pl(LRGNID_COMM,nr,ro,halo)
             rb=recvinfo_pl(BASE_COMM,nr,ro,halo)*cmax
             write(ADM_log_fid,*) &
                  'ro =',ro   &
                  ,'rank =',rank_me   &
                  ,' src_rank =', sourcerank(ro,halo) &
                  ,' rs =', rs &
                  ,' rl =', rl &
                  ,' rb =', rb
          enddo
       enddo !loop ro
    enddo !loop halo
    !
    write(ADM_log_fid,*) &
         '---------------------------------------&
         &       commnication table  end         &
         &---------------------------------------'
    !
    !call ADM_proc_stop
    return
    !
  end subroutine output_info

  !-----------------------------------------------------------------------------
  subroutine COMM_data_transfer(&
       var,   &
       var_pl )
    use mod_adm, only: &
       ADM_COMM_world, &
       ADM_proc_stop,  &
       ADM_vlink_nmax, &
       ADM_lall,       &
       ADM_kall
    implicit none

    real(RP), intent(inout) ::  var   (:,:,:,:)
    real(RP), intent(inout) ::  var_pl(:,:,:,:)

    integer ::  shp(4)
    integer ::  cmax, kmax, varmax

    integer ::  acount
    integer ::  areq(2*(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4))
    integer ::  stat(MPI_status_size,2*(ADM_lall*max_comm_r2r+ADM_vlink_nmax*4))
    integer ::  ierr

    integer ::  k, m, n
    integer ::  nr, nc, ns
    integer ::  sl, so, sb, ss
    integer ::  rl, ro, rb, rs
    integer ::  cl, scl, cs

    integer ::  max_ssize
    integer ::  max_rsize
    integer ::  max_nsmax,max_nsmax_pl
    integer ::  max_nrmax,max_nrmax_pl
    integer ::  max_ss, max_ss_pl
    integer ::  max_rs, max_rs_pl
    integer ::  max_cs_r2r, max_cs_r2p, max_cs_p2r, max_cs_sgp
    integer ::  cur_somax, cur_romax
    integer ::  cur_ncmax_r2r, cur_ncmax_r2p, cur_ncmax_p2r, cur_ncmax_sgp
    integer ::  cur_ssize, cur_rsize
    integer ::  cur_nsmax, cur_nsmax_pl
    integer ::  cur_nrmax, cur_nrmax_pl
    !---------------------------------------------------------------------------

    if ( opt_comm_barrier ) then
       call DEBUG_rapstart('COMM_barrier')
       call MPI_Barrier( ADM_COMM_world, ierr )
       call DEBUG_rapend  ('COMM_barrier')
    endif

    !$acc wait

    call DEBUG_rapstart('COMM_data_transfer')

    shp    = shape(var)
    kmax   = shp(2)
    varmax = shp(4)

    cmax = varmax * kmax

    if( opt_check_varmax ) then
       if ( cmax > max_varmax * ADM_kall ) then
          write(ADM_LOG_FID,*)  'error: cmax >  max_varmax * ADM_kall, stop!'
          write(ADM_LOG_FID,*)  'cmax=', cmax, 'max_varmax*ADM_kall=', max_varmax*ADM_kall
          call ADM_proc_stop
       endif
    endif

    max_rsize    = 1
    max_nrmax    = 0
    max_nrmax_pl = 0
    max_rs       = 0
    max_rs_pl    = 0
    do ro = 1, romax(1)
       max_rsize    = max( max_rsize,    rsize   (ro,1) )
       max_nrmax    = max( max_nrmax,    nrmax   (ro,1) )
       max_nrmax_pl = max( max_nrmax_pl, nrmax_pl(ro,1) )

       do nr = 1, nrmax(ro,1)
          max_rs = max( max_rs, recvinfo(SIZE_COMM,nr,ro,1) )
       enddo

       do nr = 1, nrmax_pl(ro,1)
          max_rs_pl = max( max_rs_pl, recvinfo_pl(SIZE_COMM,nr,ro,1) )
       enddo
    enddo
    max_rsize = max_rsize * cmax

    max_ssize    = 1
    max_nsmax    = 0
    max_nsmax_pl = 0
    max_ss       = 0
    max_ss_pl    = 0
    do so = 1, somax(1)
       max_ssize    = max( max_ssize,    ssize   (so,1) )
       max_nsmax    = max( max_nsmax,    nsmax   (so,1) )
       max_nsmax_pl = max( max_nsmax_pl, nsmax_pl(so,1) )

       do ns = 1, nsmax(so,1)
          max_ss = max( max_ss, sendinfo(SIZE_COMM,ns,so,1) )
       enddo

       do ns = 1, nsmax_pl(so,1)
          max_ss_pl = max( max_ss_pl, sendinfo_pl(SIZE_COMM,ns,so,1) )
       enddo
    enddo
    max_ssize = max_ssize * cmax

    max_cs_r2r = 0
    max_cs_r2p = 0
    max_cs_p2r = 0
    max_cs_sgp = 0
    do nc = 1, ncmax_r2r(1)
       max_cs_r2r = max( max_cs_r2r, copyinfo_r2r(SIZE_COPY,nc,1) )
    enddo
    do nc = 1, ncmax_r2p(1)
       max_cs_r2p = max( max_cs_r2p, copyinfo_r2p(SIZE_COPY,nc,1) )
    enddo
    do nc = 1, ncmax_p2r(1)
       max_cs_p2r = max( max_cs_p2r, copyinfo_p2r(SIZE_COPY,nc,1) )
    enddo
    do nc = 1, ncmax_sgp(1)
       max_cs_sgp = max( max_cs_sgp, copyinfo_sgp(SIZE_COPY,nc,1) )
    enddo

    cur_somax     = somax(1)
    cur_romax     = romax(1)
    cur_ncmax_r2r = ncmax_r2r(1)
    cur_ncmax_r2p = ncmax_r2p(1)
    cur_ncmax_p2r = ncmax_p2r(1)
    cur_ncmax_sgp = ncmax_sgp(1)

    !$acc data present(var) pcopy(var_pl) pcopyin(clist) async(0)

    !-----------------------------------------
    ! call MPI_IRECV
    !-----------------------------------------
    do ro = 1, romax(1)
       call MPI_IRECV( recvbuf(1,ro),    &
                       rsize(ro,1)*cmax, &
                       COMM_datatype,    &
                       sourcerank(ro,1), &
                       recvtag(ro,1),    &
                       ADM_COMM_WORLD,   &
                       areq(ro),         &
                       ierr              )
    enddo

    !$acc data &
    !$acc& pcopyin(sendlist,sendlist_pl) &
    !$acc& pcopyin(sendinfo,sendinfo_pl) &
    !$acc& pcopyin(nsmax,nsmax_pl) async(0)

    !$acc wait

    do so = 1, cur_somax
       cur_ssize    = ssize   (so,1) * cmax
       cur_nsmax    = nsmax   (so,1)
       cur_nsmax_pl = nsmax_pl(so,1)

       !$acc kernels copyout(sendbuf(1:cur_ssize,so)) async(so)

       !-----------------------------------------
       ! var -> sendbuf
       !-----------------------------------------
       !$acc loop independent gang
       do m = 1, varmax
       !$acc loop independent gang vector(8)
       do k = 1, kmax
       !$acc loop independent gang
       do ns = 1, cur_nsmax
       !$acc loop independent vector(32)
       do n = 1, max_ss
          if ( ns <= nsmax(so,1) ) then
             ss = sendinfo(SIZE_COMM  ,ns,so,1)
             sl = sendinfo(LRGNID_COMM,ns,so,1)
             sb = sendinfo(BASE_COMM  ,ns,so,1) * cmax
             if ( n <= ss ) then
                sendbuf(n+(k-1)*ss+(m-1)*ss*kmax+sb,so) = var(sendlist(n,ns,so,1),k,sl,m)
             endif
          endif
       enddo
       enddo
       enddo
       enddo

       !-----------------------------------------
       !  var_pl -> sendbuf
       !-----------------------------------------
       !$acc loop independent gang
       do m = 1, varmax
       !$acc loop independent gang vector(8)
       do k = 1, kmax
       !$acc loop independent gang
       do ns = 1, cur_nsmax_pl
       !$acc loop independent vector(32)
       do n = 1, max_ss_pl
          if ( ns <= nsmax_pl(so,1) ) then
             ss = sendinfo_pl(SIZE_COMM  ,ns,so,1)
             sl = sendinfo_pl(LRGNID_COMM,ns,so,1)
             sb = sendinfo_pl(BASE_COMM  ,ns,so,1) * cmax
             if ( n <= ss ) then
                sendbuf(n+(k-1)*ss+(m-1)*ss*kmax+sb,so) = var_pl(sendlist_pl(n,ns,so,1),k,sl,m)
             endif
          endif
       enddo
       enddo
       enddo
       enddo

       !$acc end kernels
    enddo

    !$acc end data

    !-----------------------------------------
    ! call MPI_ISEND
    !-----------------------------------------
    do so = 1, somax(1)
       !$acc wait(so)
       call MPI_ISEND( sendbuf(1,so),     &
                       ssize(so,1)*cmax,  &
                       COMM_datatype,     &
                       destrank(so,1),    &
                       sendtag(so,1),     &
                       ADM_COMM_WORLD,    &
                       areq(so+romax(1)), &
                       ierr               )
    enddo

    !$acc wait

    !---------------------------------------------------
    !  var -> var (region to region copy in same rank)
    !---------------------------------------------------
    !$acc kernels present(var) &
    !$acc& pcopyin(recvlist_r2r,sendlist_r2r,copyinfo_r2r) async(0)
    !$acc loop independent gang
    do m = 1, varmax
    !$acc loop independent gang
    do nc = 1, cur_ncmax_r2r
    !$acc loop independent gang vector(8)
    do k = 1, kmax
    !$acc loop independent vector(32)
    do n = 1, max_cs_r2r
       cs  = copyinfo_r2r(SIZE_COPY      ,nc,1)
       cl  = copyinfo_r2r(LRGNID_COPY    ,nc,1)
       scl = copyinfo_r2r(SRC_LRGNID_COPY,nc,1)
       if ( n <= cs ) then
          var(recvlist_r2r(n,nc,1),k,cl ,m) = var(sendlist_r2r(n,nc,1),k,scl,m)
       endif
    enddo
    enddo
    enddo
    enddo
    !$acc end kernels

    !------------------------------------------
    !  var -> var_pl ( data copy in same rank)
    !------------------------------------------
    !$acc kernels present(var,var_pl) &
    !$acc& pcopyin(recvlist_r2p,sendlist_r2p,copyinfo_r2p) async(0)
    !$acc loop independent gang
    do m = 1, varmax
    !$acc loop independent gang
    do nc = 1, cur_ncmax_r2p
    !$acc loop independent gang vector(8)
    do k = 1, kmax
    !$acc loop independent vector(32)
    do n = 1, max_cs_r2p
       cs  = copyinfo_r2p(SIZE_COPY      ,nc,1)
       cl  = copyinfo_r2p(LRGNID_COPY    ,nc,1)
       scl = copyinfo_r2p(SRC_LRGNID_COPY,nc,1)
       if ( n <= cs ) then
          var_pl(recvlist_r2p(n,nc,1),k,cl,m) = var(sendlist_r2p(n,nc,1),k,scl,m)
       endif
    enddo
    enddo
    enddo
    enddo
    !$acc end kernels

    !-----------------------------------------
    !  var_pl -> var (data copy in same rank)
    !-----------------------------------------
    !$acc kernels present(var,var_pl) &
    !$acc& pcopyin(recvlist_p2r,sendlist_p2r,copyinfo_p2r) async(0)
    !$acc loop independent gang
    do m = 1, varmax
    !$acc loop independent gang
    do nc = 1, cur_ncmax_p2r
    !$acc loop independent gang vector(8)
    do k = 1, kmax
    !$acc loop independent vector(32)
    do n = 1, max_cs_p2r
       cs  = copyinfo_p2r(SIZE_COPY      ,nc,1)
       cl  = copyinfo_p2r(LRGNID_COPY    ,nc,1)
       scl = copyinfo_p2r(SRC_LRGNID_COPY,nc,1)
       if ( n <= cs ) then
          var(recvlist_p2r(n,nc,1),k,cl,m) = var_pl(sendlist_p2r(n,nc,1),k,scl,m)
       endif
    enddo
    enddo
    enddo
    enddo
    !$acc end kernels

    acount = romax(1) + somax(1)

    call MPI_WAITALL(acount,areq,stat,ierr)

    if ( opt_comm_barrier ) then
       call MPI_Barrier(ADM_COMM_world,ierr)
    endif

    !$acc data pcopyin(recvlist,recvlist_pl) &
    !$acc& pcopyin(recvinfo,recvinfo_pl) &
    !$acc& pcopyin(nrmax,nrmax_pl) async(0)

    !$acc wait

    do ro = 1, cur_romax
       cur_rsize    = rsize(ro,1)*cmax
       cur_nrmax    = nrmax(ro,1)
       cur_nrmax_pl = nrmax_pl(ro,1)

       !$acc kernels copyin(recvbuf(1:cur_rsize,ro)) async(ro)

       !-----------------------------------------
       !  recvbuf -> var ( recieve in region )
       !-----------------------------------------
       !$acc loop independent gang
       do m = 1, varmax
       !$acc loop independent gang vector(8)
       do k = 1, kmax
       !$acc loop independent gang
       do nr = 1, cur_nrmax
       !$acc loop independent vector(32)
       do n = 1, max_rs
          if ( nr <= nrmax(ro,1) ) then
             rs = recvinfo(SIZE_COMM  ,nr,ro,1)
             rl = recvinfo(LRGNID_COMM,nr,ro,1)
             rb = recvinfo(BASE_COMM  ,nr,ro,1) * cmax
             if ( n <= rs ) then
                var(recvlist(n,nr,ro,1),k,rl,m) = recvbuf(n+(k-1)*rs+(m-1)*rs*kmax+rb,ro)
             endif
          endif
       enddo
       enddo
       enddo
       enddo

       !-----------------------------------------
       !  recvbuf -> var_pl ( recieve in pole )
       !-----------------------------------------
       !$acc loop independent gang
       do m = 1, varmax
       !$acc loop independent gang vector(8)
       do k = 1, kmax
       !$acc loop independent gang
       do nr = 1, cur_nrmax_pl
       !$acc loop independent vector(32)
       do n = 1, max_rs_pl
          if ( nr <= nrmax_pl(ro,1) ) then
             rs = recvinfo_pl(SIZE_COMM  ,nr,ro,1)
             rl = recvinfo_pl(LRGNID_COMM,nr,ro,1)
             rb = recvinfo_pl(BASE_COMM  ,nr,ro,1) * cmax
             if ( n <= rs ) then
                var_pl(recvlist_pl(n,nr,ro,1),k,rl,m) = recvbuf(n+(k-1)*rs+(m-1)*rs*kmax+rb,ro)
             endif
          endif
       enddo
       enddo
       enddo
       enddo

       !$acc end kernels
    enddo

    !$acc wait
    !$acc end data

    !-----------------------------------------
    !  copy data around singular point
    !-----------------------------------------
    !$acc kernels present(var) &
    !$acc& pcopyin(recvlist_sgp,sendlist_sgp,copyinfo_sgp) async(0)
    !$acc loop independent gang
    do m = 1, varmax
    !$acc loop independent gang
    do nc = 1, cur_ncmax_sgp
    !$acc loop independent gang vector(8)
    do k = 1, kmax
    !$acc loop independent vector(32)
    do n = 1, max_cs_sgp
       cs  = copyinfo_sgp(SIZE_COPY      ,nc,1)
       cl  = copyinfo_sgp(LRGNID_COPY    ,nc,1)
       scl = copyinfo_sgp(SRC_LRGNID_COPY,nc,1)
       if ( n <= cs ) then
          var(recvlist_sgp(n,nc,1),k,cl ,m) = var(sendlist_sgp(n,nc,1),k,scl,m)
       endif
    enddo
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    !$acc wait

    call DEBUG_rapend('COMM_data_transfer')

    return
  end subroutine COMM_data_transfer

  !-----------------------------------------------------------------------------
  subroutine COMM_var( &
       var,    &
       var_pl, &
       knum,   &
       nnum    )
    use mod_adm, only: &
       ADM_COMM_WORLD,     &
       ADM_prc_tab,        &
       ADM_rgn2prc,        &
       ADM_prc_me,         &
       ADM_NPL,            &
       ADM_SPL,            &
       ADM_prc_npl,        &
       ADM_prc_spl,        &
       ADM_rgnid_npl_mng,  &
       ADM_rgnid_spl_mng,  &
       ADM_gall,           &
       ADM_gall_pl,        &
       ADM_lall,           &
       ADM_lall_pl,        &
       ADM_gall_1d,        &
       ADM_gmin,           &
       ADM_gmax,           &
       ADM_GSLF_PL
    implicit none

    integer, intent(in)  ::  knum
    integer, intent(in)  ::  nnum
    real(RP), intent(inout) ::  var   (ADM_gall,   knum,ADM_lall,   nnum)
    real(RP), intent(inout) ::  var_pl(ADM_gall_pl,knum,ADM_lall_pl,nnum)

    real(RP) ::  v_npl_send(knum,nnum)
    real(RP) ::  v_spl_send(knum,nnum)
    real(RP) ::  v_npl_recv(knum,nnum)
    real(RP) ::  v_spl_recv(knum,nnum)

    integer ::  ireq(4)
    integer ::  istat(MPI_STATUS_SIZE)

    integer ::  ierr
    integer ::  k, l, n, rgnid

    integer ::  i,j,suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    if ( opt_comm_barrier ) then
       call DEBUG_rapstart('COMM_barrier')
       call MPI_Barrier( ADM_COMM_world, ierr )
       call DEBUG_rapend  ('COMM_barrier')
    endif

    call DEBUG_rapstart('COMM_var')

    !$acc data present(var)

    if( comm_pl ) then ! T.Ohno 110721

       !--- recv pole value
       !--- north pole
       if ( ADM_prc_me == ADM_prc_npl ) then
          call MPI_IRECV( v_npl_recv,                       &
                          knum * nnum,                      &
                          COMM_datatype,                    &
                          ADM_rgn2prc(ADM_rgnid_npl_mng)-1, &
                          ADM_NPL,                          &
                          ADM_COMM_WORLD,                   &
                          ireq(3),                          &
                          ierr                              )
       endif

       !--- south pole
       if ( ADM_prc_me == ADM_prc_spl ) then
          call MPI_IRECV( v_spl_recv,                       &
                          knum * nnum,                      &
                          COMM_datatype,                    &
                          ADM_rgn2prc(ADM_rgnid_spl_mng)-1, &
                          ADM_SPL,                          &
                          ADM_COMM_WORLD,                   &
                          ireq(4),                          &
                          ierr                              )
       endif

       !--- send pole value
       do l = 1, ADM_lall
          rgnid = ADM_prc_tab(l,ADM_prc_me)

          !--- north pole
          if ( rgnid == ADM_rgnid_npl_mng ) then
             !$acc kernels copyout(v_npl_send) pcopyin(var) async(0)
             do n = 1, nnum
             do k = 1, knum
                v_npl_send(k,n) = var(suf(ADM_gmin,ADM_gmax+1),k,l,n)
             enddo
             enddo
             !$acc end kernels

             !$acc wait

             call MPI_ISEND( v_npl_send,     &
                             knum * nnum,    &
                             COMM_datatype,  &
                             ADM_prc_npl-1,  &
                             ADM_NPL,        &
                             ADM_COMM_WORLD, &
                             ireq(1),        &
                             ierr            )
          endif

          !--- south pole
          if ( rgnid == ADM_rgnid_spl_mng ) then
             !$acc kernels copyout(v_spl_send) pcopyin(var) async(0)
             do n = 1, nnum
             do k = 1, knum
                v_spl_send(k,n) = var(suf(ADM_gmax+1,ADM_gmin),k,l,n)
             enddo
             enddo
             !$acc end kernels

             !$acc wait

             call MPI_ISEND( v_spl_send,     &
                             knum * nnum,    &
                             COMM_datatype,  &
                             ADM_prc_spl-1,  &
                             ADM_SPL,        &
                             ADM_COMM_WORLD, &
                             ireq(2),        &
                             ierr            )
          endif

       enddo

       do l = 1, ADM_lall
          rgnid = ADM_prc_tab(l,ADM_prc_me)

          !--- north pole
          if( rgnid == ADM_rgnid_npl_mng ) call MPI_WAIT(ireq(1),istat,ierr)

          !--- south pole
          if( rgnid == ADM_rgnid_spl_mng ) call MPI_WAIT(ireq(2),istat,ierr)
       enddo

       if ( ADM_prc_me == ADM_prc_npl ) then
          call MPI_WAIT(ireq(3),istat,ierr)

          do n = 1, nnum
          do k = 1, knum
             var_pl(ADM_GSLF_PL,k,ADM_NPL,n) = v_npl_recv(k,n)
          enddo
          enddo
       endif

       if ( ADM_prc_me == ADM_prc_spl ) then
          call MPI_WAIT(ireq(4),istat,ierr)

          do n = 1, nnum
          do k = 1, knum
             var_pl(ADM_GSLF_PL,k,ADM_SPL,n) = v_spl_recv(k,n)
          enddo
          enddo
       endif

    endif

    !$acc end data

    call COMM_data_transfer(var,var_pl)

    var(suf(ADM_gmax+1,ADM_gmin-1),:,:,:) = var(suf(ADM_gmax+1,ADM_gmin),:,:,:)
    var(suf(ADM_gmin-1,ADM_gmax+1),:,:,:) = var(suf(ADM_gmin,ADM_gmax+1),:,:,:)

    call DEBUG_rapend('COMM_var')

    return
  end subroutine COMM_var

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_sum( localsum, globalsum )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all
    implicit none

    real(RP), intent(in) ::  localsum
    real(RP), intent(out) ::  globalsum

    real(RP) ::  sendbuf(1)
    real(RP) ::  recvbuf(ADM_prc_all)

    integer ::  ierr
    !---------------------------------------------------------------------------

    if ( COMM_pl ) then
       sendbuf(1) = localsum

       call MPI_Allgather( sendbuf,        &
                           1,              &
                           COMM_datatype,  &
                           recvbuf,        &
                           1,              &
                           COMM_datatype,  &
                           ADM_COMM_WORLD, &
                           ierr            )

       globalsum = sum( recvbuf(:) )
    else
       globalsum = localsum
    endif

    return
  end subroutine COMM_Stat_sum

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_sum_eachlayer( kall, localsum, globalsum )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all
    implicit none

    integer, intent(in) ::  kall
    real(RP), intent(in) ::  localsum (kall)
    real(RP), intent(out) ::  globalsum(kall)

    real(RP) ::  sendbuf(kall)
    integer ::  displs (ADM_prc_all)
    integer ::  counts (ADM_prc_all)
    real(RP) ::  recvbuf(kall,ADM_prc_all)

    integer ::  ierr
    integer ::  k, p
    !---------------------------------------------------------------------------

    do p = 1, ADM_prc_all
       displs(p) = (p-1) * kall
       counts(p) = kall
    enddo

    if ( COMM_pl ) then
       sendbuf(:) = localsum(:)

       call MPI_Allgatherv( sendbuf,        &
                            kall,           &
                            COMM_datatype,  &
                            recvbuf,        &
                            counts,         &
                            displs,         &
                            COMM_datatype,  &
                            ADM_COMM_WORLD, &
                            ierr            )

       do k = 1, kall
          globalsum(k) = sum( recvbuf(k,:) )
       enddo
    else
       do k = 1, kall
          globalsum(k) = localsum(k)
       enddo
    endif

    return
  end subroutine COMM_Stat_sum_eachlayer

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_avg( localavg, globalavg )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all
    implicit none

    real(RP), intent(in) ::  localavg
    real(RP), intent(out) ::  globalavg

    real(RP) ::  sendbuf(1)
    real(RP) ::  recvbuf(ADM_prc_all)

    integer ::  ierr
    !---------------------------------------------------------------------------

    if ( COMM_pl ) then
       sendbuf(1) = localavg

       call MPI_Allgather( sendbuf,        &
                           1,              &
                           COMM_datatype,  &
                           recvbuf,        &
                           1,              &
                           COMM_datatype,  &
                           ADM_COMM_WORLD, &
                           ierr            )

       globalavg = sum( recvbuf(:) ) / real(ADM_prc_all,kind=RP)
    else
       globalavg = localavg
    endif

    !write(ADM_LOG_FID,*) 'COMM_Stat_avg', sendbuf(1), recvbuf(:)

    return
  end subroutine COMM_Stat_avg

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_max( localmax, globalmax )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all
    implicit none

    real(RP), intent(in) ::  localmax
    real(RP), intent(out) ::  globalmax

    real(RP) ::  sendbuf(1)
    real(RP) ::  recvbuf(ADM_prc_all)

    integer ::  ierr
    !---------------------------------------------------------------------------

    sendbuf(1) = localmax

    call MPI_Allgather( sendbuf,        &
                        1,              &
                        COMM_datatype,  &
                        recvbuf,        &
                        1,              &
                        COMM_datatype,  &
                        ADM_COMM_WORLD, &
                        ierr            )

    globalmax = maxval( recvbuf(:) )

    !write(ADM_LOG_FID,*) 'COMM_Stat_max', sendbuf(1), recvbuf(:)

    return
  end subroutine COMM_Stat_max

  !-----------------------------------------------------------------------------
  subroutine COMM_Stat_min( localmin, globalmin )
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all
    implicit none

    real(RP), intent(in) ::  localmin
    real(RP), intent(out) ::  globalmin

    real(RP) ::  sendbuf(1)
    real(RP) ::  recvbuf(ADM_prc_all)

    integer ::  ierr
    !---------------------------------------------------------------------------

    sendbuf(1) = localmin

    call MPI_Allgather( sendbuf,        &
                        1,              &
                        COMM_datatype,  &
                        recvbuf,        &
                        1,              &
                        COMM_datatype,  &
                        ADM_COMM_WORLD, &
                        ierr            )

    globalmin = minval( recvbuf(:) )

    !write(ADM_LOG_FID,*) 'COMM_Stat_min', sendbuf(1), recvbuf(:)

    return
  end subroutine COMM_Stat_min

end module mod_comm
!-------------------------------------------------------------------------------
