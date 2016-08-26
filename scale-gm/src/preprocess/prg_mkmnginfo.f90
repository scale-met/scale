!-------------------------------------------------------------------------------
!> Program mkmnginfo
!!
!! @par Description
!!          Making information file for paralell computation
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
program mkmnginfo
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof
  use mod_adm, only: &
     ADM_RID, &
     ADM_DIR, &
     ADM_SW,  &
     ADM_NW,  &
     ADM_NE,  &
     ADM_SE
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  integer, parameter :: fid=10
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  integer                :: rlevel                   ! region division level
  integer                :: prc_num                  ! process number
  character(len=H_SHORT) :: HGRID_SYSTEM = 'ICO'     ! grid system( default ico)
                                         ! 'LCP'
                                         ! 'MLCP'
                                         ! 'MLCP-OLD'
                                         ! 'PERIODIC-1DMD'  : T.Ohno 110721
                                         ! '1DMD-ON-SPHERE' : M.Hara 110721
  character(len=H_SHORT) :: MAPPING_TYPE = ''        ! mapping method [add] C.Kodama 2011/12/14
                                         ! ''        : standard
                                         ! 'K-TERAI' : TERAI Mapping for K-Computer
  integer                :: XTMS_K       = 6         ! (not used for icosahedral)
  integer                :: XTMS_MLCP_S  = 1         ! only for MLCP
  character(len=H_LONG)  :: output_fname = ''        ! output region-management filename

  namelist / mkmnginfo_cnf / &
       rlevel,       &
       prc_num,      &
       hgrid_system, &
       MAPPING_TYPE, &
       XTMS_K,       &
       XTMS_MLCP_S,  &
       output_fname

  integer :: ierr
  !=============================================================================

  call MPI_Init(ierr)

  open( unit   = fid,             &
        file   = 'mkmnginfo.cnf', &
        status = 'old',           &
        form   = 'formatted'      )

     read(fid,nml=mkmnginfo_cnf)

  close(fid)

  if    ( HGRID_SYSTEM == 'LCP') then!

     call generate_mngtab_lcp(rlevel,prc_num,output_fname)

  elseif( HGRID_SYSTEM == 'MLCP') then

     call generate_mngtab_mlcp(rlevel,prc_num,output_fname)

  elseif( HGRID_SYSTEM == 'MLCP-OLD') then

     call generate_mngtab_mlcp_old(rlevel,prc_num,output_fname)

  elseif( HGRID_SYSTEM == 'PERIODIC-1DMD') then

     call generate_mngtab_periodic_1dmd(rlevel,prc_num,output_fname)

  elseif( HGRID_SYSTEM == '1DMD-ON-SPHERE') then

     call generate_mngtab_1dmd_on_sphere(rlevel,prc_num,output_fname)

  else ! default is icosahedral

     call generate_mngtab(rlevel,prc_num,output_fname)

  endif

  call MPI_Finalize(ierr)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine generate_mngtab( rl, nmax_prc, fname )
    use mod_adm, only: &
       nmax_mng => ADM_l_limit
    implicit none

    integer,          intent(in) :: rl
    integer,          intent(in) :: nmax_prc
    character(len=*), intent(in) :: fname

    integer :: i,j,d
    integer :: i_nb,j_nb,d_nb,edgid_nb
    integer :: l,l_nb
    integer :: k,m,p
    integer :: rgnlen
    integer :: tmp, tmp_4r, tmp_m
    integer, parameter :: nmax_dmd = 10

    integer :: all_rgn

    integer, allocatable :: rgn_tab(:,:,:)
    integer, allocatable :: mngrgn (:)
    integer, allocatable :: prc_tab(:,:)

    integer, parameter   :: fid = 20

    integer :: num_of_rgn
    namelist / rgn_info / num_of_rgn

    integer :: rgnid
    integer :: &
         sw(ADM_RID:ADM_DIR),&
         nw(ADM_RID:ADM_DIR),&
         ne(ADM_RID:ADM_DIR),&
         se(ADM_RID:ADM_DIR)
    namelist / rgn_link_info / rgnid, sw, nw, ne, se

    integer :: num_of_proc
    namelist /proc_info/ num_of_proc

    integer :: peid
    integer :: num_of_mng
    integer :: mng_rgnid(nmax_mng)
    namelist /rgn_mng_info/ peid, num_of_mng,mng_rgnid

    integer :: dmd_data(ADM_SW:ADM_SE,nmax_dmd)
    !---------------------------------------------------------------------------

    dmd_data(ADM_SW:ADM_SE, 1)=(/ 6, 5, 2,10/)
    dmd_data(ADM_SW:ADM_SE, 2)=(/10, 1, 3, 9/)
    dmd_data(ADM_SW:ADM_SE, 3)=(/ 9, 2, 4, 8/)
    dmd_data(ADM_SW:ADM_SE, 4)=(/ 8, 3, 5, 7/)
    dmd_data(ADM_SW:ADM_SE, 5)=(/ 7, 4, 1, 6/)
    dmd_data(ADM_SW:ADM_SE, 6)=(/ 7, 5, 1,10/)
    dmd_data(ADM_SW:ADM_SE, 7)=(/ 8, 4, 5, 6/)
    dmd_data(ADM_SW:ADM_SE, 8)=(/ 9, 3, 4, 7/)
    dmd_data(ADM_SW:ADM_SE, 9)=(/10, 2, 3, 8/)
    dmd_data(ADM_SW:ADM_SE,10)=(/ 6, 1, 2, 9/)

    rgnlen = 2**rl
    all_rgn = nmax_dmd*rgnlen*rgnlen

    allocate(rgn_tab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,all_rgn))

    do d=1,nmax_dmd
       do i=1,rgnlen
          do j=1,rgnlen
             !
             l=(rgnlen*rgnlen)*(d-1)+rgnlen*(j-1)+i
             !
             do k=ADM_SW,ADM_SE
                select case(k)
                case(ADM_SW)
                   if(j==1) then
                      if(d<=5) then
                         i_nb=i
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_NE
                      else
                         i_nb=rgnlen
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_SE
                      endif
                   else
                      i_nb=i
                      j_nb=j-1
                      d_nb=d
                      edgid_nb=ADM_NE
                   endif
                case(ADM_NW)
                   if(i==1) then
                      if(d<=5) then
                         i_nb=rgnlen+1-j
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_NE
                      else
                         i_nb=rgnlen
                         j_nb=j
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_SE
                      endif
                   else
                      i_nb=i-1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_SE
                   endif
                case(ADM_NE)
                   if(j==rgnlen) then
                      if(d<=5) then
                         i_nb=1
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_NW
                      else
                         i_nb=i
                         j_nb=1
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_SW
                      endif
                   else
                      i_nb=i
                      j_nb=j+1
                      d_nb=d
                      edgid_nb=ADM_SW
                   endif
                case(ADM_SE)
                   if(i==rgnlen) then
                      if(d<=5) then
                         i_nb=1
                         j_nb=j
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_NW
                      else
                         i_nb=rgnlen+1-j
                         j_nb=1
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_SW
                      endif
                   else
                      i_nb=i+1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_NW
                   endif
                end select
                !
                l_nb=(rgnlen*rgnlen)*(d_nb-1)+rgnlen*(j_nb-1)+i_nb
                rgn_tab(ADM_RID,k,l)=l_nb
                rgn_tab(ADM_DIR,k,l)=edgid_nb
                !
             enddo
          enddo
       enddo
    enddo

    allocate( mngrgn (nmax_prc)          )
    allocate( prc_tab(nmax_mng,nmax_prc) )

    if ( mod(all_rgn,nmax_prc) /= 0 ) then
       write(*,*) 'Invalid number of process!'
       stop
    else
       mngrgn(:) = all_rgn / nmax_prc
    endif

    do m = 1, nmax_prc
       prc_tab(1:nmax_mng,m) = -1

       do p = 1, mngrgn(m)

          tmp_m = (m-1) * mngrgn(m) + p

          if ( MAPPING_TYPE == 'K-TERAI' ) then
             ! mngrgn(m) = 1 if one region per one process
             ! m: peid (-> rank)
             ! p: 1 or 2 if one or two regions per one process
             ! prc_tab(p,m): region number
             !
             ! mapping between m of Nreg-1prc and m of 1reg-1prc (N=1,2)

             if( p >= 3 ) then ! = m (1reg-1prc),  = 2m-1 or 2m (2reg-1prc)
                write(*,*) 'More than two regions is not allowed when MAPING_TYPE = K-TERAI'
                stop
             endif

             tmp_4r = 4**rlevel
             if ( mod(tmp_m-1,2*tmp_4r) < tmp_4r ) then
                tmp =      (       tmp_m-1-mod(tmp_m-1,tmp_4r)) / (2*tmp_4r) + 1
             else
                tmp = 10 - (tmp_4r+tmp_m-1-mod(tmp_m-1,tmp_4r)) / (2*tmp_4r) + 1
             endif

             prc_tab(p,m) = (tmp-1)*tmp_4r + mod(tmp_m-1,tmp_4r) + 1

             write(*,*) 'peid=', m, 'regid=', prc_tab(p,m)

          else

             prc_tab(p,m) = tmp_m ! default

          endif

       enddo
    enddo

    open( unit = fid,         &
          file = trim(fname), &
          form = 'formatted'  )

       num_of_rgn = all_rgn
       write(fid,nml=rgn_info)

       do l = 1, all_rgn
          rgnid = l
          sw    = rgn_tab(:,ADM_SW,l)
          nw    = rgn_tab(:,ADM_NW,l)
          ne    = rgn_tab(:,ADM_NE,l)
          se    = rgn_tab(:,ADM_SE,l)

          write(fid,nml=rgn_link_info)
       enddo

       num_of_proc = nmax_prc
       write(fid,nml=proc_info)

       do m = 1, nmax_prc
          peid       = m
          num_of_mng = mngrgn(m)
          mng_rgnid  = prc_tab(:,m)

          write(fid,nml=rgn_mng_info)
       enddo

    close(fid)

    return
  end subroutine generate_mngtab

!---------------- for LCP
  subroutine generate_mngtab_lcp( rl, nmax_prc, fname )
    use mod_adm, only: &
         nmax_mng => ADM_l_limit
    implicit none

    integer, intent(in) :: rl
    integer, intent(in) :: nmax_prc
    character(len=*), intent(in) :: fname
    !
    integer :: i,j,d
    integer :: i_nb,j_nb,d_nb,edgid_nb
    integer :: l,l_nb
    integer :: k,m,p
    integer :: rgnlen
    integer :: nmax_dmd=-1
    !
    integer :: all_rgn
    !
    integer, allocatable :: rgn_tab(:,:,:)
    integer, allocatable :: mngrgn(:)
    integer, allocatable :: prc_tab(:,:)
    !
    integer, parameter :: fid=20
    !
    integer :: num_of_rgn
    namelist / rgn_info / num_of_rgn
    !
    integer :: rgnid
    integer :: &
         sw(ADM_RID:ADM_DIR),&
         nw(ADM_RID:ADM_DIR),&
         ne(ADM_RID:ADM_DIR),&
         se(ADM_RID:ADM_DIR)
    namelist / rgn_link_info / rgnid, sw, nw, ne, se
    !
    integer :: num_of_proc
    namelist /proc_info/ num_of_proc
    !
    integer :: peid
    integer :: num_of_mng
    integer :: mng_rgnid(nmax_mng)
    namelist /rgn_mng_info/ peid, num_of_mng,mng_rgnid

    integer,allocatable :: dmd_data(:,:)
    !---------------------------------------------------------------------------

    nmax_dmd = XTMS_K * 4
    write(*,*) nmax_dmd, XTMS_K
    allocate(dmd_data(ADM_SW:ADM_SE,nmax_dmd))

    do d=1,XTMS_K
       dmd_data(ADM_SW:ADM_SE, d)=(/ &
            d+XTMS_K, &
            mod(d-2+XTMS_K,XTMS_K)+1  , &
            mod(d,XTMS_K)+1  , &
            d+XTMS_K*2  /)
    enddo
    do d=XTMS_K+1,XTMS_K*2
       dmd_data(ADM_SW:ADM_SE, d)=(/ &
            mod(d-2,XTMS_K)+XTMS_K*3+1   , &
            mod(d-2,XTMS_K)+XTMS_K*2+1  , &
            d-XTMS_K,&
            d+XTMS_K  /)
    enddo
    do d=XTMS_K*2+1,XTMS_K*3
       dmd_data(ADM_SW:ADM_SE, d)=(/ &
            d-XTMS_K  , &
            d-XTMS_K*2  , &
            mod(d,XTMS_K)+XTMS_K*1+1   , &
            d+XTMS_K  /)
    enddo
    do d=XTMS_K*3+1,XTMS_K*4
       dmd_data(ADM_SW:ADM_SE, d)=(/ &
            mod(d-2,XTMS_K)+XTMS_K*3+1, &
            d-XTMS_K, &
            mod(d,XTMS_K)+XTMS_K*1+1, &
            mod(d,XTMS_K)+XTMS_K*3+1   /)
    enddo

    rgnlen=2**rl
    all_rgn=nmax_dmd*rgnlen*rgnlen
    !
    allocate(rgn_tab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,all_rgn))
    !
    do d=1,nmax_dmd
       do i=1,rgnlen
          do j=1,rgnlen
             !
             l=(rgnlen*rgnlen)*(d-1)+rgnlen*(j-1)+i
             !
             do k=ADM_SW,ADM_SE
                select case(k)
                case(ADM_SW)
                   if(j==1) then
                      if(d<=XTMS_K) then
                         i_nb=i
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_NE
                      elseif(d<=XTMS_K*2) then
                         i_nb=i
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_NE
                      elseif(d<=XTMS_K*3) then
                         i_nb=rgnlen
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_SE
                      elseif(d<=XTMS_K*4) then
                         i_nb=rgnlen
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_SE
                      endif
                   else
                      i_nb=i
                      j_nb=j-1
                      d_nb=d
                      edgid_nb=ADM_NE
                   endif
                case(ADM_NW)
                   if(i==1) then
                      if(d<=XTMS_K) then
                         i_nb=rgnlen+1-j
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_NE
                      elseif(d<=XTMS_K*2) then
                         i_nb=rgnlen+1-j
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_NE
                      elseif(d<=XTMS_K*3) then
                         i_nb=rgnlen
                         j_nb=j
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_SE
                      elseif(d<=XTMS_K*4) then
                         i_nb=rgnlen
                         j_nb=j
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_SE
                      endif
                   else
                      i_nb=i-1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_SE
                   endif
                case(ADM_NE)
                   if(j==rgnlen) then
                      if(d<=XTMS_K) then
                         i_nb=1
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_NW
                      elseif(d<=XTMS_K*2) then
                         i_nb=i
                         j_nb=1
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_SW
                      elseif(d<=XTMS_K*3) then
                        i_nb=1
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_NW
                      elseif(d<=XTMS_K*4) then
                         i_nb=i
                         j_nb=1
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_SW
                      endif
                   else
                      i_nb=i
                      j_nb=j+1
                      d_nb=d
                      edgid_nb=ADM_SW
                   endif
                case(ADM_SE)
                   if(i==rgnlen) then
                      if(d<=XTMS_K) then
                         i_nb=1
                         j_nb=j
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_NW
                      elseif(d<=XTMS_K*2) then
                         i_nb=rgnlen+1-j
                         j_nb=1
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_SW
                      elseif(d<=XTMS_K*3) then
                         i_nb=1
                         j_nb=j
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_NW
                      elseif(d<=XTMS_K*4) then
                         i_nb=rgnlen+1-j
                         j_nb=1
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_SW
                      endif
                   else
                      i_nb=i+1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_NW
                   endif
                end select
                !
                l_nb=(rgnlen*rgnlen)*(d_nb-1)+rgnlen*(j_nb-1)+i_nb
                rgn_tab(ADM_RID,k,l)=l_nb
                rgn_tab(ADM_DIR,k,l)=edgid_nb
                !
             enddo
          enddo
       enddo
    enddo
    !
!    nmax_prc=all_rgn
    !
    allocate(mngrgn(nmax_prc))
    allocate(prc_tab(nmax_mng,nmax_prc))
    do m=1,nmax_prc
       if(mod(all_rgn,nmax_prc)/=0) then
          write(*,*) 'Invalid number of process!'
          stop
       else
          mngrgn(m)=all_rgn/nmax_prc
       endif
       prc_tab(1:nmax_mng,m)=-1
       do p=1,mngrgn(m)
          prc_tab(p,m)=(m-1)*(all_rgn/nmax_prc)+p
       enddo
    enddo
    !
    open(fid,file=trim(fname),form='formatted')
    !
    num_of_rgn=all_rgn
    write(fid,nml=rgn_info)
    !
    do l=1,all_rgn
       rgnid=l
       sw=rgn_tab(:,ADM_SW,l)
       nw=rgn_tab(:,ADM_NW,l)
       ne=rgn_tab(:,ADM_NE,l)
       se=rgn_tab(:,ADM_SE,l)
       write(fid,nml=rgn_link_info)
    enddo
    num_of_proc=nmax_prc
    write(fid,nml=proc_info)
    do m=1,nmax_prc
       peid=m
       num_of_mng=mngrgn(m)
       mng_rgnid=prc_tab(:,m)
       write(fid,nml=rgn_mng_info)
    enddo
    !
    close(fid)
    !
  end subroutine generate_mngtab_lcp
  !-----------------------------------------------------------------------------!
  !---------------- for MLCP
  subroutine generate_mngtab_mlcp( rl, nmax_prc, fname )
    use mod_adm, only: &
         nmax_mng => ADM_l_limit
    implicit none

    integer, intent(in) :: rl
    integer, intent(in) :: nmax_prc
    character(len=*), intent(in) :: fname
    !
    integer :: i,j,d
    integer :: i_nb,j_nb,d_nb,edgid_nb
    integer :: l,l_nb
    integer :: k,m,p
    integer :: rgnlen

    integer,save :: nmax_dmd=-1
    !
    integer :: all_rgn
    !
    integer, allocatable :: rgn_tab(:,:,:)
    integer, allocatable :: mngrgn(:)
    integer, allocatable :: prc_tab(:,:)
    !
    integer,parameter :: fid=20
    !
    integer :: num_of_rgn
    namelist / rgn_info / num_of_rgn
    !
    integer :: rgnid
    integer :: &
         sw(ADM_RID:ADM_DIR),&
         nw(ADM_RID:ADM_DIR),&
         ne(ADM_RID:ADM_DIR),&
         se(ADM_RID:ADM_DIR)
    namelist / rgn_link_info / rgnid, sw, nw, ne, se
    !
    integer :: num_of_proc
    namelist /proc_info/ num_of_proc
    !
    integer :: peid
    integer :: num_of_mng
    integer :: mng_rgnid(nmax_mng)
    namelist /rgn_mng_info/ peid, num_of_mng,mng_rgnid

    integer,allocatable :: dmd_data(:,:)

    integer::s
    !---------------------------------------------------------------------------

    nmax_dmd = XTMS_K * (1+XTMS_MLCP_S)
    write(*,*) nmax_dmd, XTMS_K
    allocate(dmd_data(ADM_SW:ADM_SE,nmax_dmd))
    !
    !

    do k=1,XTMS_K
       do s=1,1
          dmd_data(ADM_SW:ADM_SE,(k-1)*(XTMS_MLCP_S+1)+s)=(/ &
               (k-2) *(XTMS_MLCP_S+1)+s+1, &
               (k-2) *(XTMS_MLCP_S+1)+s, &
               (k) *(XTMS_MLCP_S+1)+s, &
               (k-1) *(XTMS_MLCP_S+1)+s+1 &
               /)
       enddo
       do s=2,XTMS_MLCP_S
          dmd_data(ADM_SW:ADM_SE,(k-1)*(XTMS_MLCP_S+1)+s)=(/ &
               (k-2) *(XTMS_MLCP_S+1)+s+1, &
               (k-1) *(XTMS_MLCP_S+1)+s-1, &
               (k) *(XTMS_MLCP_S+1)+s-1, &
               (k-1) *(XTMS_MLCP_S+1)+s+1 &
               /)
       enddo
       do s=XTMS_MLCP_S+1,XTMS_MLCP_S+1
          dmd_data(ADM_SW:ADM_SE,(k-1)*(XTMS_MLCP_S+1)+s)=(/ &
               (k-2) *(XTMS_MLCP_S+1)+s, &
               (k-1) *(XTMS_MLCP_S+1)+s-1, &
               (k) *(XTMS_MLCP_S+1)+s-1, &
               (k) *(XTMS_MLCP_S+1)+s &
               /)
       enddo
    enddo

    do i=ADM_SW,ADM_SE
       do d=1,nmax_dmd
          if (dmd_data(i,d)<1) then
             dmd_data(i,d)=dmd_data(i,d)+nmax_dmd
          elseif(dmd_data(i,d)>nmax_dmd) then
             dmd_data(i,d)=dmd_data(i,d)-nmax_dmd
          endif
       enddo
    enddo

    rgnlen=2**rl
    all_rgn=nmax_dmd*rgnlen*rgnlen
    !
    allocate(rgn_tab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,all_rgn))
    !
    do d=1,nmax_dmd
       do i=1,rgnlen
          do j=1,rgnlen
             !
             l=(rgnlen*rgnlen)*(d-1)+rgnlen*(j-1)+i
             !
             do k=ADM_SW,ADM_SE
                select case(k)
                case(ADM_SW)
                   if(j==1) then
                      if(mod(d,(XTMS_MLCP_S+1))==0) then
                         i_nb=rgnlen
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_SE
                      else
                         i_nb=i
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_NE
                      endif
                   else
                      i_nb=i
                      j_nb=j-1
                      d_nb=d
                      edgid_nb=ADM_NE
                   endif
                case(ADM_NW)
                   if(i==1) then
                      if(mod(d,(XTMS_MLCP_S+1))==1) then
                         i_nb=rgnlen+1-j
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_NE
                      else
                         i_nb=rgnlen
                         j_nb=j
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_SE
                      endif
                   else
                      i_nb=i-1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_SE
                   endif
                case(ADM_NE)
                   if(j==rgnlen) then
                      if(mod(d,(XTMS_MLCP_S+1))==1) then
                         i_nb=1
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_NW
                      else
                         i_nb=i
                         j_nb=1
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_SW
                      endif
                   else
                      i_nb=i
                      j_nb=j+1
                      d_nb=d
                      edgid_nb=ADM_SW
                   endif


                case(ADM_SE)
                   if(i==rgnlen) then
                      if(mod(d,(XTMS_MLCP_S+1))==0) then
                         i_nb=rgnlen+1-j
                         j_nb=1
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_SW
                      else
                         i_nb=1
                         j_nb=j
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_NW
                      endif
                   else
                      i_nb=i+1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_NW
                   endif
                end select
                !
                l_nb=(rgnlen*rgnlen)*(d_nb-1)+rgnlen*(j_nb-1)+i_nb
                rgn_tab(ADM_RID,k,l)=l_nb
                rgn_tab(ADM_DIR,k,l)=edgid_nb
                !
             enddo
          enddo
       enddo
    enddo
    !
!    nmax_prc=all_rgn
    !
    allocate(mngrgn(nmax_prc))
    allocate(prc_tab(nmax_mng,nmax_prc))
    do m=1,nmax_prc
       if(mod(all_rgn,nmax_prc)/=0) then
          write(*,*) 'Invalid number of process!'
          stop
       else
          mngrgn(m)=all_rgn/nmax_prc
       endif
       prc_tab(1:nmax_mng,m)=-1
       do p=1,mngrgn(m)
          prc_tab(p,m)=(m-1)*(all_rgn/nmax_prc)+p
       enddo
    enddo
    !
    open(fid,file=trim(fname),form='formatted')
    !
    num_of_rgn=all_rgn
    write(fid,nml=rgn_info)
    !
    do l=1,all_rgn
       rgnid=l
       sw=rgn_tab(:,ADM_SW,l)
       nw=rgn_tab(:,ADM_NW,l)
       ne=rgn_tab(:,ADM_NE,l)
       se=rgn_tab(:,ADM_SE,l)
       write(fid,nml=rgn_link_info)
    enddo
    num_of_proc=nmax_prc
    write(fid,nml=proc_info)
    do m=1,nmax_prc
       peid=m
       num_of_mng=mngrgn(m)
       mng_rgnid=prc_tab(:,m)
       write(fid,nml=rgn_mng_info)
    enddo
    !
    close(fid)
    !
  end subroutine generate_mngtab_mlcp
  !-----------------------------------------------------------------------------
  !---------------- for MLCP
  subroutine generate_mngtab_mlcp_old( rl, nmax_prc, fname )
    use mod_adm, only: &
         nmax_mng => ADM_l_limit
    implicit none

    integer, intent(in) :: rl
    integer, intent(in) :: nmax_prc
    character(len=*), intent(in) :: fname
    !
    integer :: i,j,d
    integer :: i_nb,j_nb,d_nb,edgid_nb
    integer :: l,l_nb
    integer :: k,m,p
    integer :: rgnlen

    integer,save :: nmax_dmd=-1
    !
    integer :: all_rgn
    !
    integer, allocatable :: rgn_tab(:,:,:)
    integer, allocatable :: mngrgn(:)
    integer, allocatable :: prc_tab(:,:)
    !
    integer,parameter :: fid=20
    !
    integer :: num_of_rgn
    namelist / rgn_info / num_of_rgn
    !
    integer :: rgnid
    integer :: &
         sw(ADM_RID:ADM_DIR),&
         nw(ADM_RID:ADM_DIR),&
         ne(ADM_RID:ADM_DIR),&
         se(ADM_RID:ADM_DIR)
    namelist / rgn_link_info / rgnid, sw, nw, ne, se
    !
    integer :: num_of_proc
    namelist /proc_info/ num_of_proc
    !
    integer :: peid
    integer :: num_of_mng
    integer :: mng_rgnid(nmax_mng)
    namelist /rgn_mng_info/ peid, num_of_mng,mng_rgnid

    integer,allocatable :: dmd_data(:,:)
    !---------------------------------------------------------------------------


    nmax_dmd = XTMS_K * 2
    write(*,*) nmax_dmd, XTMS_K
    allocate(dmd_data(ADM_SW:ADM_SE,nmax_dmd))
    !
    !

    do d=1,XTMS_K
       dmd_data(ADM_SW:ADM_SE, d)=(/ &
            mod(XTMS_K-d+1,XTMS_K)+1+XTMS_K, &
            mod(d-2+XTMS_K,XTMS_K)+1  , &
            mod(d,XTMS_K)+1  , &
            XTMS_K*2+1-d  /)
    enddo
    do d=XTMS_K+1,XTMS_K*2
       dmd_data(ADM_SW:ADM_SE, d)=(/ &
            mod(d,XTMS_K)+XTMS_K+1   , &
            mod(XTMS_K*2-d,XTMS_K)+1  , &
            mod(XTMS_K*2-d+1,XTMS_K)+1  , &
            mod(d-2,XTMS_K)+XTMS_K+1   /)
    enddo

    rgnlen=2**rl
    all_rgn=nmax_dmd*rgnlen*rgnlen
    !
    allocate(rgn_tab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,all_rgn))
    !
    do d=1,nmax_dmd
       do i=1,rgnlen
          do j=1,rgnlen
             !
             l=(rgnlen*rgnlen)*(d-1)+rgnlen*(j-1)+i
             !
             do k=ADM_SW,ADM_SE
                select case(k)
                case(ADM_SW)
                   if(j==1) then
                      if(d<=XTMS_K) then
                         i_nb=i
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_NE
                      else
                         i_nb=rgnlen
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_SE
                      endif
                   else
                      i_nb=i
                      j_nb=j-1
                      d_nb=d
                      edgid_nb=ADM_NE
                   endif
                case(ADM_NW)
                   if(i==1) then
                      if(d<=XTMS_K) then
                         i_nb=rgnlen+1-j
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_NE
                      else
                         i_nb=rgnlen
                         j_nb=j
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_SE
                      endif
                   else
                      i_nb=i-1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_SE
                   endif
                case(ADM_NE)
                   if(j==rgnlen) then
                      if(d<=XTMS_K) then
                         i_nb=1
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_NW
                      else
                         i_nb=i
                         j_nb=1
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_SW
                      endif
                   else
                      i_nb=i
                      j_nb=j+1
                      d_nb=d
                      edgid_nb=ADM_SW
                   endif
                case(ADM_SE)
                   if(i==rgnlen) then
                      if(d<=XTMS_K) then
                         i_nb=1
                         j_nb=j
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_NW
                      else
                         i_nb=rgnlen+1-j
                         j_nb=1
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_SW
                      endif
                   else
                      i_nb=i+1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_NW
                   endif
                end select
                !
                l_nb=(rgnlen*rgnlen)*(d_nb-1)+rgnlen*(j_nb-1)+i_nb
                rgn_tab(ADM_RID,k,l)=l_nb
                rgn_tab(ADM_DIR,k,l)=edgid_nb
                !
             enddo
          enddo
       enddo
    enddo
    !
!    nmax_prc=all_rgn
    !
    allocate(mngrgn(nmax_prc))
    allocate(prc_tab(nmax_mng,nmax_prc))
    do m=1,nmax_prc
       if(mod(all_rgn,nmax_prc)/=0) then
          write(*,*) 'Invalid number of process!'
          stop
       else
          mngrgn(m)=all_rgn/nmax_prc
       endif
       prc_tab(1:nmax_mng,m)=-1
       do p=1,mngrgn(m)
          prc_tab(p,m)=(m-1)*(all_rgn/nmax_prc)+p
       enddo
    enddo
    !
    open(fid,file=trim(fname),form='formatted')
    !
    num_of_rgn=all_rgn
    write(fid,nml=rgn_info)
    !
    do l=1,all_rgn
       rgnid=l
       sw=rgn_tab(:,ADM_SW,l)
       nw=rgn_tab(:,ADM_NW,l)
       ne=rgn_tab(:,ADM_NE,l)
       se=rgn_tab(:,ADM_SE,l)
       write(fid,nml=rgn_link_info)
    enddo
    num_of_proc=nmax_prc
    write(fid,nml=proc_info)
    do m=1,nmax_prc
       peid=m
       num_of_mng=mngrgn(m)
       mng_rgnid=prc_tab(:,m)
       write(fid,nml=rgn_mng_info)
    enddo
    !
    close(fid)
    !
  end subroutine generate_mngtab_mlcp_old

  !-----------------------------------------------------------------------------
  subroutine generate_mngtab_periodic_1dmd( rl, nmax_prc, fname )
    use mod_adm, only: &
         nmax_mng => ADM_l_limit
    implicit none

    integer, intent(in)          :: rl
    integer, intent(in)          :: nmax_prc
    character(len=*), intent(in) :: fname
    !
    integer :: i,j,d
    integer :: i_nb,j_nb,d_nb,edgid_nb
    integer :: l,l_nb
    integer :: k,m,p
    integer :: rgnlen
    integer, parameter :: nmax_dmd=1
    integer :: all_rgn
    !
    integer, allocatable :: rgn_tab(:,:,:)
    integer, allocatable :: mngrgn(:)
    integer, allocatable :: prc_tab(:,:)
    !
    integer,parameter :: fid=20
    !
    integer :: num_of_rgn
    namelist / rgn_info / num_of_rgn
    !
    integer :: rgnid
    integer :: &
         sw(ADM_RID:ADM_DIR),&
         nw(ADM_RID:ADM_DIR),&
         ne(ADM_RID:ADM_DIR),&
         se(ADM_RID:ADM_DIR)
    namelist / rgn_link_info / rgnid, sw, nw, ne, se
    !
    integer :: num_of_proc
    namelist /proc_info/ num_of_proc
    !
    integer :: peid
    integer :: num_of_mng
    integer :: mng_rgnid(nmax_mng)
    namelist /rgn_mng_info/ peid, num_of_mng, mng_rgnid

    integer :: dmd_data(ADM_SW:ADM_SE,nmax_dmd)
    !---------------------------------------------------------------------------

    dmd_data(ADM_SW:ADM_SE, 1)=(/1,1,1,1/)
    !
    rgnlen  = 2**rl
    all_rgn = nmax_dmd*rgnlen**2
    !
    allocate(rgn_tab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,all_rgn))
    !
    do d=1,nmax_dmd
       do i=1,rgnlen
          do j=1,rgnlen
             !
             l=(rgnlen*rgnlen)*(d-1)+rgnlen*(j-1)+i
             !
             do k=ADM_SW,ADM_SE
                select case(k)
                case(ADM_SW)
                   if(j==1) then
                    i_nb=i
                    j_nb=rgnlen
                    d_nb=dmd_data(ADM_SW,d)
                    edgid_nb=ADM_NE
                   else
                      i_nb=i
                      j_nb=j-1
                      d_nb=d
                      edgid_nb=ADM_NE
                   endif
                case(ADM_NW)
                   if(i==1) then
                    i_nb=rgnlen
                    j_nb=j
                    d_nb=dmd_data(ADM_NW,d)
                    edgid_nb=ADM_SE
                   else
                      i_nb=i-1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_SE
                   endif
                case(ADM_NE)
                   if(j==rgnlen) then
                    i_nb=i
                    j_nb=1
                    d_nb=dmd_data(ADM_NE,d)
                    edgid_nb=ADM_SW
                   else
                      i_nb=i
                      j_nb=j+1
                      d_nb=d
                      edgid_nb=ADM_SW
                   endif
                case(ADM_SE)
                   if(i==rgnlen) then
                    i_nb=1
                    j_nb=j
                    d_nb=dmd_data(ADM_SE,d)
                    edgid_nb=ADM_NW
                   else
                      i_nb=i+1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_NW
                   endif
                end select
                !
                l_nb=(rgnlen*rgnlen)*(d_nb-1)+rgnlen*(j_nb-1)+i_nb
                rgn_tab(ADM_RID,k,l)=l_nb
                rgn_tab(ADM_DIR,k,l)=edgid_nb
                !
             enddo
          enddo
       enddo
    enddo
    !
    allocate(mngrgn(nmax_prc))
    allocate(prc_tab(nmax_mng,nmax_prc))
    do m=1,nmax_prc
       if(mod(all_rgn,nmax_prc)/=0) then
          write(*,*) 'Invalid number of process!'
          stop
       else
          mngrgn(m)=all_rgn/nmax_prc
       endif
       prc_tab(1:nmax_mng,m)=-1
       do p=1,mngrgn(m)
          prc_tab(p,m)=(m-1)*(all_rgn/nmax_prc)+p
       enddo
    enddo
    !
    open(fid,file=trim(fname),form='formatted')
    !
    num_of_rgn=all_rgn
    write(fid,nml=rgn_info)
    !
    do l=1,all_rgn
       rgnid = l
       sw    = rgn_tab(:,ADM_SW,l)
       nw    = rgn_tab(:,ADM_NW,l)
       ne    = rgn_tab(:,ADM_NE,l)
       se    = rgn_tab(:,ADM_SE,l)
       write(fid,nml=rgn_link_info)
    enddo
    num_of_proc=nmax_prc
    write(fid,nml=proc_info)
    do m=1,nmax_prc
       peid=m
       num_of_mng=mngrgn(m)
       mng_rgnid=prc_tab(:,m)
       write(fid,nml=rgn_mng_info)
    enddo
    !
    close(fid)
    !
  end subroutine generate_mngtab_periodic_1dmd

  !-----------------------------------------------------------------------------
  subroutine generate_mngtab_1dmd_on_sphere( rl, nmax_prc, fname )
    use mod_adm, only: &
         nmax_mng => ADM_l_limit
    implicit none

    integer, intent(in) :: rl
    integer, intent(in) :: nmax_prc
    character(len=*), intent(in) :: fname
    !
    integer :: i,j,d
    integer :: i_nb,j_nb,d_nb,edgid_nb
    integer :: l,l_nb
    integer :: k,m,p
    integer :: rgnlen
    integer, parameter :: nmax_dmd=10
    !
    integer :: all_rgn
    !
    integer, allocatable :: rgn_tab(:,:,:)
    integer, allocatable :: mngrgn(:)
    integer, allocatable :: prc_tab(:,:)
    !
    integer, parameter :: fid=20
    !
    integer :: num_of_rgn
    namelist / rgn_info / num_of_rgn
    !
    integer :: rgnid
    integer :: &
         sw(ADM_RID:ADM_DIR),&
         nw(ADM_RID:ADM_DIR),&
         ne(ADM_RID:ADM_DIR),&
         se(ADM_RID:ADM_DIR)
    namelist / rgn_link_info / rgnid, sw, nw, ne, se
    !
    integer :: num_of_proc
    namelist /proc_info/ num_of_proc
    !
    integer :: peid
    integer :: num_of_mng
    integer :: mng_rgnid(nmax_mng)
    namelist /rgn_mng_info/ peid, num_of_mng,mng_rgnid

    integer :: dmd_data(ADM_SW:ADM_SE,nmax_dmd)
    !---------------------------------------------------------------------------

    dmd_data(ADM_SW:ADM_SE, 1)=(/ 6, 5, 2,10/)

    rgnlen=2**rl
    all_rgn=nmax_dmd*rgnlen*rgnlen
    !
    allocate(rgn_tab(ADM_RID:ADM_DIR,ADM_SW:ADM_SE,all_rgn))
    !
    do d=1,nmax_dmd
       do i=1,rgnlen
          do j=1,rgnlen
             !
             l=(rgnlen*rgnlen)*(d-1)+rgnlen*(j-1)+i
             !
             do k=ADM_SW,ADM_SE
                select case(k)
                case(ADM_SW)
                   if(j==1) then
                      if(d<=5) then
                         i_nb=i
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_NE
                      else
                         i_nb=rgnlen
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_SW,d)
                         edgid_nb=ADM_SE
                      endif
                   else
                      i_nb=i
                      j_nb=j-1
                      d_nb=d
                      edgid_nb=ADM_NE
                   endif
                case(ADM_NW)
                   if(i==1) then
                      if(d<=5) then
                         i_nb=rgnlen+1-j
                         j_nb=rgnlen
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_NE
                      else
                         i_nb=rgnlen
                         j_nb=j
                         d_nb=dmd_data(ADM_NW,d)
                         edgid_nb=ADM_SE
                      endif
                   else
                      i_nb=i-1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_SE
                   endif
                case(ADM_NE)
                   if(j==rgnlen) then
                      if(d<=5) then
                         i_nb=1
                         j_nb=rgnlen+1-i
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_NW
                      else
                         i_nb=i
                         j_nb=1
                         d_nb=dmd_data(ADM_NE,d)
                         edgid_nb=ADM_SW
                      endif
                   else
                      i_nb=i
                      j_nb=j+1
                      d_nb=d
                      edgid_nb=ADM_SW
                   endif
                case(ADM_SE)
                   if(i==rgnlen) then
                      if(d<=5) then
                         i_nb=1
                         j_nb=j
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_NW
                      else
                         i_nb=rgnlen+1-j
                         j_nb=1
                         d_nb=dmd_data(ADM_SE,d)
                         edgid_nb=ADM_SW
                      endif
                   else
                      i_nb=i+1
                      j_nb=j
                      d_nb=d
                      edgid_nb=ADM_NW
                   endif
                end select
                !
                l_nb=(rgnlen*rgnlen)*(d_nb-1)+rgnlen*(j_nb-1)+i_nb
                rgn_tab(ADM_RID,k,l)=l_nb
                rgn_tab(ADM_DIR,k,l)=edgid_nb
                !
             enddo
          enddo
       enddo
    enddo
    !
!    nmax_prc=all_rgn
    !
    allocate(mngrgn(nmax_prc))
    allocate(prc_tab(nmax_mng,nmax_prc))
    do m=1,nmax_prc
       if(mod(all_rgn,nmax_prc)/=0) then
          write(*,*) 'Invalid number of process!'
          write(*,*) all_rgn, nmax_prc
          stop
       else
          mngrgn(m)=all_rgn/nmax_prc
       endif
       prc_tab(1:nmax_mng,m)=-1
       do p=1,mngrgn(m)
          prc_tab(p,m)=(m-1)*(all_rgn/nmax_prc)+p
       enddo
    enddo
    !
    open(fid,file=trim(fname),form='formatted')

    num_of_rgn=all_rgn

    num_of_rgn = num_of_rgn/10
    do l=1,num_of_rgn
       do k=ADM_SW,ADM_SE
          if (rgn_tab(ADM_RID,k,l) > num_of_rgn) then
             rgn_tab(ADM_RID,k,l) = l
             rgn_tab(ADM_DIR,k,l) = k
          endif
       enddo
    enddo
    write(fid,nml=rgn_info)
    write(6,nml=rgn_info)

    do l=1,num_of_rgn ! M.Hara110604
       rgnid=l
       sw=rgn_tab(:,ADM_SW,l)
       nw=rgn_tab(:,ADM_NW,l)
       ne=rgn_tab(:,ADM_NE,l)
       se=rgn_tab(:,ADM_SE,l)
       write(fid,nml=rgn_link_info)
       write(6,nml=rgn_link_info)
    enddo
    num_of_proc=nmax_prc
    print *, nmax_prc, num_of_proc
    num_of_proc = num_of_proc/10
    write(fid,nml=proc_info)
    write(6,nml=proc_info)
    do m=1,num_of_proc ! M.Hara110604
       peid=m
       num_of_mng=mngrgn(m)
       mng_rgnid=prc_tab(:,m)
       write(fid,nml=rgn_mng_info)
       write(6,nml=rgn_mng_info)
    enddo

    close(fid)

    return
  end subroutine generate_mngtab_1dmd_on_sphere

end program mkmnginfo

