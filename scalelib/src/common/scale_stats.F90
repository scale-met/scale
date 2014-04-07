!-------------------------------------------------------------------------------
!> module Statistics
!!
!! @par Description
!!          global statistics module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-11-21 (H.Yashiro)  [mod] Spin-off from scale_STAT
!!
!<
#include "inc_openmp.h"
module scale_stats
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: STAT_setup
  public :: STAT_total
  public :: STAT_detail

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public, save :: STAT_checktotal = .false. !< calc&report variable totals to logfile?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private, save :: STAT_use_globalcomm = .false. !< calculate total with global communication?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine STAT_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_STATS / &
       STAT_checktotal, &
       STAT_use_globalcomm

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[STAT]/Categ[COMMON]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_STATS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_STATS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_STATS)

    return
  end subroutine STAT_setup

  !-----------------------------------------------------------------------------
  !> Calc volume/area-weighted global sum
  subroutine STAT_total( allstatval, var, varname )
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_datatype
    use scale_grid_real, only: &
       area => REAL_AREA, &
       vol  => REAL_VOL
    implicit none

    real(RP),         intent(out) :: allstatval !< volume/area-weighted total
    real(RP),         intent(in)  :: var(:,:,:) !< 3D value
    character(len=*), intent(in)  :: varname    !< name of item

    real(RP) :: statval
    integer  :: ksize

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    ksize = size(var(:,:,:),1)

    statval = 0.0_RP
    if ( ksize == KA ) then ! 3D
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2) reduction(+:statval)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          statval = statval + var(k,i,j) * vol(k,i,j)
       enddo
       enddo
       enddo
    elseif( ksize == 1 ) then ! 2D
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2) reduction(+:statval)
       do j = JS, JE
       do i = IS, IE
          statval = statval + var(1,i,j) * area(i,j)
       enddo
       enddo
    endif

    if ( .not. ( statval > -1.0_RP .or. statval < 1.0_RP ) ) then ! must be NaN
       write(*,*) 'xxx [STAT_total] NaN is detected for ', trim(varname), ' in rank ', PRC_myrank
       call PRC_MPIstop
    endif

    if ( STAT_use_globalcomm ) then
       call PROF_rapstart('COMM Allreduce MPI')
       ! All reduce
       call MPI_Allreduce( statval,              &
                           allstatval,           &
                           1,                    &
                           COMM_datatype,        &
                           MPI_SUM,              &
                           MPI_COMM_WORLD,       &
                           ierr                  )

       call PROF_rapend  ('COMM Allreduce MPI')

       ! statistics over the all node
       if ( varname /= "" ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,1PE24.17)') &
                     '[', varname, '] SUM(global) =', allstatval
       endif
    else
       allstatval = statval

       ! statistics on each node
       if ( varname /= "" ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,1PE24.17)') &
                     '[', varname, '] SUM(local)  =', statval
       endif
    endif

    return
  end subroutine STAT_total

  !-----------------------------------------------------------------------------
  !> Search global maximum & minimum value
  subroutine STAT_detail(var, varname)
    use scale_process, only: &
       PRC_nmax,   &
       PRC_myrank
    use scale_const, only: &
       CONST_UNDEF8, &
       CONST_UNDEF2
    use scale_comm, only: &
       COMM_datatype
    implicit none

    real(RP),         intent(inout) :: var(:,:,:,:) !< values
    character(len=*), intent(in)    :: varname(:)   !< name of item

    logical :: halomask(KA,IA,JA)

    real(RP), allocatable :: statval   (:,:,:)
    integer,  allocatable :: statidx   (:,:,:,:)
    real(RP), allocatable :: allstatval(:,:)
    integer,  allocatable :: allstatidx(:,:,:)
    integer               :: vsize

    integer :: ierr
    integer :: v, p
    !---------------------------------------------------------------------------

    vsize = size(var(:,:,:,:),4)

    halomask(:,:,:) = .false.
    halomask(KS:KE,IS:IE,JS:JE) = .true.

    allocate( statval(  vsize,2,0:PRC_nmax-1) ); statval(:,:,:)   = CONST_UNDEF8
    allocate( statidx(3,vsize,2,0:PRC_nmax-1) ); statidx(:,:,:,:) = CONST_UNDEF2

    allocate( allstatval(  vsize,2) ); allstatval(:,:)   = CONST_UNDEF8
    allocate( allstatidx(1,vsize,2) ); allstatidx(:,:,:) = CONST_UNDEF2

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Variable Statistics ***'
    do v = 1, vsize
       statval(  v,1,PRC_myrank) = maxval(var(:,:,:,v),mask=halomask)
       statval(  v,2,PRC_myrank) = minval(var(:,:,:,v),mask=halomask)
       statidx(:,v,1,PRC_myrank) = maxloc(var(:,:,:,v),mask=halomask)
       statidx(:,v,2,PRC_myrank) = minloc(var(:,:,:,v),mask=halomask)

! statistics on each node
!       if( IO_L ) write(IO_FID_LOG,*) '*** [', trim(varname(v)), ']'
!       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,3(I5,A))') '*** MAX = ', &
!                                             statval(  v,1,PRC_myrank),'(', &
!                                             statidx(1,v,1,PRC_myrank),',', &
!                                             statidx(2,v,1,PRC_myrank),',', &
!                                             statidx(3,v,1,PRC_myrank),')'
!       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,3(I5,A))') '*** MIN = ', &
!                                             statval(  v,2,PRC_myrank),'(', &
!                                             statidx(1,v,2,PRC_myrank),',', &
!                                             statidx(2,v,2,PRC_myrank),',', &
!                                             statidx(3,v,2,PRC_myrank),')'
    enddo

    ! MPI broadcast
    call PROF_rapstart('COMM Bcast MPI')
    do p = 0, PRC_nmax-1
       call MPI_Bcast( statval(1,1,p),   &
                       vsize*2,          &
                       COMM_datatype,    &
                       p,                &
                       MPI_COMM_WORLD,   &
                       ierr              )
       call MPI_Bcast( statidx(1,1,1,p), &
                       3*vsize*2,        &
                       MPI_INTEGER,      &
                       p,                &
                       MPI_COMM_WORLD,   &
                       ierr              )
    enddo
    call PROF_rapend  ('COMM Bcast MPI')

    do v = 1, vsize
       allstatval(v,1)   = maxval(statval(v,1,:))
       allstatval(v,2)   = minval(statval(v,2,:))
       allstatidx(:,v,1) = maxloc(statval(v,1,:))-1
       allstatidx(:,v,2) = minloc(statval(v,2,:))-1
       if( IO_L ) write(IO_FID_LOG,*) '[', trim(varname(v)), ']'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,4(I5,A))') '  MAX =', &
                                                    allstatval(  v,1), '(', &
                                                    allstatidx(1,v,1), ',', &
                                      statidx(1,v,1,allstatidx(1,v,1)),',', &
                                      statidx(2,v,1,allstatidx(1,v,1)),',', &
                                      statidx(3,v,1,allstatidx(1,v,1)),')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,4(I5,A))') '  MIN =', &
                                                    allstatval(  v,2), '(', &
                                                    allstatidx(1,v,2), ',', &
                                      statidx(1,v,2,allstatidx(1,v,2)),',', &
                                      statidx(2,v,2,allstatidx(1,v,2)),',', &
                                      statidx(3,v,2,allstatidx(1,v,2)),')'
    enddo

    deallocate( statval )
    deallocate( statidx )

    deallocate( allstatval )
    deallocate( allstatidx )

    return
  end subroutine STAT_detail

end module scale_stats
!-------------------------------------------------------------------------------
