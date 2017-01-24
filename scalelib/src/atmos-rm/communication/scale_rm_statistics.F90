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
module scale_rm_statistics
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_process, only: &
     PRC_LOCAL_COMM_WORLD
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

  interface STAT_total
     module procedure STAT_total_2D
     module procedure STAT_total_3D
  end interface STAT_total

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: STATISTICS_checktotal = .false. !< calc&report variable totals to logfile?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: STATISTICS_use_globalcomm = .false. !< calculate total with global communication?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine STAT_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_STATISTICS / &
       STATISTICS_checktotal, &
       STATISTICS_use_globalcomm

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[STATISTICS] / Categ[ATMOS-RM COMM] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_STATISTICS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_STATISTICS. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_STATISTICS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Caluculate statistics?                     : ', STATISTICS_checktotal
    if( IO_L ) write(IO_FID_LOG,*) '*** Allow global communication for statistics? : ', STATISTICS_use_globalcomm
    if ( STATISTICS_use_globalcomm ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** => Global total is calculated using MPI_ALLreduce.'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** => Local total is calculated in each process.'
    endif

    return
  end subroutine STAT_setup

  !-----------------------------------------------------------------------------
  !> Calc volume/area-weighted global sum
  subroutine STAT_total_2D( allstatval, var, varname )
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_datatype
    use scale_grid_real, only: &
       area => REAL_AREA
    implicit none

    real(RP),         intent(out) :: allstatval !< volume/area-weighted total
    real(RP),         intent(in)  :: var(IA,JA) !< 3D value
    character(len=*), intent(in)  :: varname    !< name of item

    character(len=24) :: varname_trim
    real(RP) :: statval

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    varname_trim = trim(varname)

    statval = 0.0_RP
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2) reduction(+:statval)
    do j = JS, JE
    do i = IS, IE
       statval = statval + var(i,j) * area(i,j)
    enddo
    enddo

    if ( .NOT. ( statval > -1.0_RP .OR. statval < 1.0_RP ) ) then ! must be NaN
       write(*,*) 'xxx [STAT_total] NaN is detected for ', varname_trim, ' in rank ', PRC_myrank
       call PRC_MPIstop
    endif

    if ( STATISTICS_use_globalcomm ) then
       call PROF_rapstart('COMM_Allreduce', 2)
       ! All reduce
       call MPI_Allreduce( statval,              &
                           allstatval,           &
                           1,                    &
                           COMM_datatype,        &
                           MPI_SUM,              &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr                  )

       call PROF_rapend  ('COMM_Allreduce', 2)

       ! statistics over the all node
       if ( varname_trim /= "" ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,ES24.17)') &
                     '[', varname_trim, '] SUM(global) = ', allstatval
       endif
    else
       allstatval = statval

       ! statistics on each node
       if ( varname_trim /= "" ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,ES24.17)') &
                     '[', varname_trim, '] SUM(local)  = ', statval
       endif
    endif

    return
  end subroutine STAT_total_2D

  !-----------------------------------------------------------------------------
  !> Calc volume/area-weighted global sum
  subroutine STAT_total_3D( allstatval, var, varname )
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_datatype
    use scale_grid_real, only: &
       vol  => REAL_VOL
    implicit none

    real(RP),         intent(out) :: allstatval    !< volume/area-weighted total
    real(RP),         intent(in)  :: var(KA,IA,JA) !< 3D value
    character(len=*), intent(in)  :: varname       !< name of item

    character(len=24) :: varname_trim
    real(RP) :: statval

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    varname_trim = trim(varname)

    statval = 0.0_RP
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2) reduction(+:statval)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       statval = statval + var(k,i,j) * vol(k,i,j)
    enddo
    enddo
    enddo

    if ( .NOT. ( statval > -1.0_RP .OR. statval < 1.0_RP ) ) then ! must be NaN
       write(*,*) 'xxx [STAT_total] NaN is detected for ', varname_trim, ' in rank ', PRC_myrank
       call PRC_MPIstop
    endif

    if ( STATISTICS_use_globalcomm ) then
       call PROF_rapstart('COMM_Allreduce', 2)
       ! All reduce
       call MPI_Allreduce( statval,              &
                           allstatval,           &
                           1,                    &
                           COMM_datatype,        &
                           MPI_SUM,              &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr                  )

       call PROF_rapend  ('COMM_Allreduce', 2)

       ! statistics over the all node
       if ( varname_trim /= "" ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,ES24.17)') &
                     '[', varname_trim, '] SUM(global) = ', allstatval
       endif
    else
       allstatval = statval

       ! statistics on each node
       if ( varname_trim /= "" ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,ES24.17)') &
                     '[', varname_trim, '] SUM(local)  = ', statval
       endif
    endif

    return
  end subroutine STAT_total_3D

  !-----------------------------------------------------------------------------
  !> Search global maximum & minimum value
  subroutine STAT_detail(var, varname, supress_globalcomm)
    use scale_process, only: &
       PRC_nprocs,   &
       PRC_myrank
    use scale_const, only: &
       CONST_UNDEF8, &
       CONST_UNDEF2
    use scale_comm, only: &
       COMM_datatype
    implicit none

    real(RP),         intent(inout) :: var(:,:,:,:) !< values
    character(len=*), intent(in)    :: varname(:)   !< name of item
    logical,          intent(in), optional :: supress_globalcomm !< supress global comm.?

    logical , allocatable :: halomask  (:,:,:)
    real(RP), allocatable :: statval   (:,:,:)
    integer,  allocatable :: statidx   (:,:,:,:)
    real(RP), allocatable :: allstatval(:,:)
    integer,  allocatable :: allstatidx(:,:,:)
    integer               :: ksize, vsize
    logical               :: do_globalcomm

    integer :: ierr
    integer :: v, p
    !---------------------------------------------------------------------------

    do_globalcomm = STATISTICS_use_globalcomm
    if ( present(supress_globalcomm) ) then
       if ( supress_globalcomm ) then
          do_globalcomm = .false.
       endif
    endif

    ksize = size(var(:,:,:,:),1)
    vsize = size(var(:,:,:,:),4)

    allocate( halomask(ksize,IA,JA) ); halomask(:,:,:) = .false.

    if ( ksize == KA ) then
       halomask(KS:KE,IS:IE,JS:JE) = .true.
    else
       halomask(:,IS:IE,JS:JE) = .true.
    endif

    allocate( statval(  vsize,2,0:PRC_nprocs-1) ); statval(:,:,:)   = CONST_UNDEF8
    allocate( statidx(3,vsize,2,0:PRC_nprocs-1) ); statidx(:,:,:,:) = CONST_UNDEF2

    allocate( allstatval(  vsize,2) ); allstatval(:,:)   = CONST_UNDEF8
    allocate( allstatidx(1,vsize,2) ); allstatidx(:,:,:) = CONST_UNDEF2

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Variable Statistics ***'
    do v = 1, vsize
       statval(  v,1,PRC_myrank) = maxval(var(:,:,:,v),mask=halomask)
       statval(  v,2,PRC_myrank) = minval(var(:,:,:,v),mask=halomask)
       statidx(:,v,1,PRC_myrank) = maxloc(var(:,:,:,v),mask=halomask)
       statidx(:,v,2,PRC_myrank) = minloc(var(:,:,:,v),mask=halomask)
    enddo

    if ( do_globalcomm ) then
       call PROF_rapstart('COMM_Bcast', 2)
       do p = 0, PRC_nprocs-1

          call MPI_Bcast( statval(1,1,p),       &
                          vsize*2,              &
                          COMM_datatype,        &
                          p,                    &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )

          call MPI_Bcast( statidx(1,1,1,p),     &
                          3*vsize*2,            &
                          MPI_INTEGER,          &
                          p,                    &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr                  )

       enddo
       call PROF_rapend  ('COMM_Bcast', 2)

       do v = 1, vsize
          allstatval(v,1)   = maxval(statval(v,1,:))
          allstatval(v,2)   = minval(statval(v,2,:))
          allstatidx(:,v,1) = maxloc(statval(v,1,:))-1
          allstatidx(:,v,2) = minloc(statval(v,2,:))-1
          if( IO_L ) write(IO_FID_LOG,*) '[', trim(varname(v)), ']'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,4(I5,A))') '  MAX =', &
                                                       allstatval(  v,1), '(', &
                                                       allstatidx(1,v,1), ',', &
                                         statidx(1,v,1,allstatidx(1,v,1)),',', &
                                         statidx(2,v,1,allstatidx(1,v,1)),',', &
                                         statidx(3,v,1,allstatidx(1,v,1)),')'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,4(I5,A))') '  MIN =', &
                                                       allstatval(  v,2), '(', &
                                                       allstatidx(1,v,2), ',', &
                                         statidx(1,v,2,allstatidx(1,v,2)),',', &
                                         statidx(2,v,2,allstatidx(1,v,2)),',', &
                                         statidx(3,v,2,allstatidx(1,v,2)),')'
       enddo
    else
       ! statistics on each node
       do v = 1, vsize
          if( IO_L ) write(IO_FID_LOG,*) '*** [', trim(varname(v)), ']'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,3(I5,A))') '*** MAX = ', &
                                                statval(  v,1,PRC_myrank),'(', &
                                                statidx(1,v,1,PRC_myrank),',', &
                                                statidx(2,v,1,PRC_myrank),',', &
                                                statidx(3,v,1,PRC_myrank),')'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,3(I5,A))') '*** MIN = ', &
                                                statval(  v,2,PRC_myrank),'(', &
                                                statidx(1,v,2,PRC_myrank),',', &
                                                statidx(2,v,2,PRC_myrank),',', &
                                                statidx(3,v,2,PRC_myrank),')'
       enddo
    endif

    if( IO_L ) write(IO_FID_LOG,*)

    deallocate( halomask )

    deallocate( statval )
    deallocate( statidx )

    deallocate( allstatval )
    deallocate( allstatidx )

    return
  end subroutine STAT_detail

end module scale_rm_statistics
