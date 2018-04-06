!-------------------------------------------------------------------------------
!> module Statistics
!!
!! @par Description
!!          global statistics module
!!
!! @author Team SCALE
!!
!<
#include "inc_openmp.h"
module scale_statistics
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_prc, only: &
     PRC_LOCAL_COMM_WORLD
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: STATISTICS_setup
  public :: STATISTICS_total
  public :: STATISTICS_detail

  interface STATISTICS_total
     module procedure STATISTICS_total_2D
     module procedure STATISTICS_total_3D
  end interface STATISTICS_total

  interface STATISTICS_detail
     module procedure STATISTICS_detail_2D
     module procedure STATISTICS_detail_3D
  end interface STATISTICS_detail
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
  subroutine STATISTICS_setup
    use scale_prc, only: &
       PRC_abort
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
       call PRC_abort
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
  end subroutine STATISTICS_setup

  !-----------------------------------------------------------------------------
  !> Calc domain sum and area-weighted mean
  subroutine STATISTICS_total_2D( &
       IA, IS, IE, JA, JS, JE, &
       var, varname, &
       area, total,  &
       log_suppress, &
       mean, sum     )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm, only: &
       COMM_datatype
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP),         intent(in) :: var(IA,JA)  !< 3D value
    character(len=*), intent(in) :: varname     !< name of item
    real(RP),         intent(in) :: area(IA,JA) !< area of the grid cell
    real(RP),         intent(in) :: total       !< total area

    logical,  intent(in),  optional :: log_suppress !< suppress log output
    real(RP), intent(out), optional :: mean !< area-weighted mean
    real(RP), intent(out), optional :: sum  !< domain sum

    real(RP) :: statval
    real(RP) :: sendbuf(2), recvbuf(2)
    real(RP) :: sum_, mean_

    logical :: suppress_
    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    statval = 0.0_RP
    if ( var(IS,JS) /= UNDEF ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2) reduction(+:statval)
       do j = JS, JE
       do i = IS, IE
          statval = statval + var(i,j) * area(i,j)
       end do
       end do
    end if

    if ( .NOT. ( statval > -1.0_RP .OR. statval < 1.0_RP ) ) then ! must be NaN
       write(*,*) 'xxx [STATISTICS_total] NaN is detected for ', trim(varname), ' in rank ', PRC_myrank
       call PRC_abort
    endif

    if ( present(log_suppress) ) then
       suppress_ = log_suppress
    else
       suppress_ = .false.
    end if

    if ( STATISTICS_use_globalcomm ) then
       call PROF_rapstart('COMM_Allreduce', 2)
       sendbuf(1) = statval
       sendbuf(2) = total
       ! All reduce
       call MPI_Allreduce( sendbuf(:), recvbuf(:), &
                           2,                      &
                           COMM_datatype,          &
                           MPI_SUM,                &
                           PRC_LOCAL_COMM_WORLD,   &
                           ierr                    )
       call PROF_rapend  ('COMM_Allreduce', 2)

       sum_  = recvbuf(1)
       mean_ = recvbuf(1) / recvbuf(2)
       ! statistics over the all node
       if ( .not. suppress_ ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A24,A,ES24.17)') &
                     '[', trim(varname), '] MEAN(global) = ', mean_
       endif
    else
       sum_ = statval
       mean_ = statval / total

       ! statistics on each node
       if ( .not. suppress_ ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A24,A,ES24.17)') &
                     '[', trim(varname), '] MEAN(local)  = ', mean_
       endif
    endif

    if ( present(mean) ) mean = mean_
    if ( present(sum ) ) sum  = sum_

    return
  end subroutine STATISTICS_total_2D

  !-----------------------------------------------------------------------------
  !> Calc domain sum and volume-weighted mean
  subroutine STATISTICS_total_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       var, varname, &
       vol, total,   &
       log_suppress, &
       mean, sum     )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm, only: &
       COMM_datatype
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP),         intent(in) :: var(KA,IA,JA) !< 3D value
    character(len=*), intent(in) :: varname       !< name of item
    real(RP),         intent(in) :: vol(KA,IA,JA) !< volume of the grid cell
    real(RP),         intent(in) :: total         !< total volume

    logical,  intent(in),  optional :: log_suppress !< suppress log output
    real(RP), intent(out), optional :: mean !< volume/area-weighted total
    real(RP), intent(out), optional :: sum  !< domain sum

    real(RP) :: statval
    real(RP) :: sendbuf(2), recvbuf(2)
    real(RP) :: mean_, sum_

    logical :: suppress_
    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    statval = 0.0_RP
    if ( var(KS,IS,JS) /= UNDEF ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2) reduction(+:statval)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          statval = statval + var(k,i,j) * vol(k,i,j)
       enddo
       enddo
       enddo
    end if

    if ( .NOT. ( statval > -1.0_RP .OR. statval < 1.0_RP ) ) then ! must be NaN
       write(*,*) 'xxx [STATISTICS_total] NaN is detected for ', trim(varname), ' in rank ', PRC_myrank
       call PRC_abort
    endif

    if ( present(log_suppress) ) then
       suppress_ = log_suppress
    else
       suppress_ = .false.
    end if

    if ( STATISTICS_use_globalcomm ) then
       call PROF_rapstart('COMM_Allreduce', 2)
       sendbuf(1) = statval
       sendbuf(2) = total
       ! All reduce
       call MPI_Allreduce( sendbuf(:), recvbuf(:), &
                           2,                      &
                           COMM_datatype,          &
                           MPI_SUM,                &
                           PRC_LOCAL_COMM_WORLD,   &
                           ierr                    )
       call PROF_rapend  ('COMM_Allreduce', 2)

       sum_  = recvbuf(1)
       mean_ = recvbuf(1) / recvbuf(2)
       ! statistics over the all node
       if ( .not. suppress_ ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A24,A,ES24.17)') &
                     '[', trim(varname), '] MEAN(global) = ', mean_
       endif
    else
       sum_  = statval
       mean_ = statval / total

       ! statistics on each node
       if ( .not. suppress_ ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A24,A,ES24.17)') &
                     '[', trim(varname), '] MEAN(local)  = ', mean_
       endif
    endif

    if ( present(mean) ) mean = mean_
    if ( present(sum ) ) sum  = sum_

    return
  end subroutine STATISTICS_total_3D

  !-----------------------------------------------------------------------------
  !> Search global maximum & minimum value
  subroutine STATISTICS_detail_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, VA, &
       varname, var, &
       local         )
    use scale_prc, only: &
       PRC_nprocs,   &
       PRC_myrank
    use scale_const, only: &
       CONST_UNDEF, &
       CONST_UNDEF2
    use scale_comm, only: &
       COMM_datatype
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: VA

    character(len=*), intent(in) :: varname(VA)      !< name of item
    real(RP),         intent(in) :: var(KA,IA,JA,VA) !< values

    logical,          intent(in), optional :: local  !< calc in local node

    real(RP) :: statval_l (  VA,2)
    integer  :: statidx_l (3,VA,2)
    real(RP) :: statval   (  VA,2,0:PRC_nprocs-1)
    integer  :: statidx   (3,VA,2,0:PRC_nprocs-1)
    real(RP) :: allstatval(VA,2)
    integer  :: allstatidx(VA,2)
    logical :: do_globalcomm

    integer :: k, i, j
    integer :: ierr
    integer :: v, p
    !---------------------------------------------------------------------------

    do_globalcomm = STATISTICS_use_globalcomm
    if ( present(local) ) do_globalcomm = ( .not. local )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Variable Statistics ***'
    do v = 1, VA
       statval_l(  v,:) = var(KS,IS,JS,v)
       statidx_l(1,v,:) = KS
       statidx_l(2,v,:) = IS
       statidx_l(3,v,:) = JS
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if ( var(k,i,j,v) > statval_l(v,1) ) then
             statval_l(  v,1) = var(k,i,j,v)
             statidx_l(1,v,1) = k
             statidx_l(2,v,1) = i
             statidx_l(3,v,1) = j
          end if
          if ( var(k,i,j,v) < statval_l(v,2) ) then
             statval_l(  v,2) = var(k,i,j,v)
             statidx_l(1,v,2) = k
             statidx_l(2,v,2) = i
             statidx_l(3,v,2) = j
          end if
       end do
       end do
       end do
    enddo

    if ( do_globalcomm ) then
       call PROF_rapstart('COMM_Bcast', 2)

       call MPI_AllGather( statval_l(:,:),       &
                           VA*2,                 &
                           COMM_datatype,        &
                           statval(:,:,:),       &
                           VA*2,                 &
                           COMM_datatype,        &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr )

       call MPI_AllGather( statidx_l(:,:,:),     &
                           3*VA*2,               &
                           MPI_INTEGER,          &
                           statidx(:,:,:,:),     &
                           3*VA*2,               &
                           MPI_INTEGER,          &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr )

       call PROF_rapend  ('COMM_Bcast', 2)

       do v = 1, VA
          allstatval(v,1) = statval(v,1,0)
          allstatval(v,2) = statval(v,2,0)
          allstatidx(v,:) = 0
          do p = 1, PRC_nprocs-1
             if ( statval(v,1,p) > allstatval(v,1) ) then
                allstatval(v,1) = statval(v,1,p)
                allstatidx(v,1) = p
             end if
             if ( statval(v,2,p) < allstatval(v,2) ) then
                allstatval(v,2) = statval(v,2,p)
                allstatidx(v,2) = p
             end if
          end do
          if( IO_L ) write(IO_FID_LOG,*) '[', trim(varname(v)), ']'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,4(I5,A))') '  MAX =', &
                                                       allstatval(v,1), ' (rank=', &
                                                       allstatidx(v,1), '; ', &
                                         statidx(1,v,1,allstatidx(v,1)),',', &
                                         statidx(2,v,1,allstatidx(v,1)),',', &
                                         statidx(3,v,1,allstatidx(v,1)),')'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,4(I5,A))') '  MIN =', &
                                                       allstatval(v,2), ' (rank=', &
                                                       allstatidx(v,2), '; ', &
                                         statidx(1,v,2,allstatidx(v,2)),',', &
                                         statidx(2,v,2,allstatidx(v,2)),',', &
                                         statidx(3,v,2,allstatidx(v,2)),')'
       enddo
    else
       ! statistics on each node
       do v = 1, VA
          if( IO_L ) write(IO_FID_LOG,*) '*** [', trim(varname(v)), ']'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,3(I5,A))') '*** MAX = ', &
                                                statval_l(  v,1),' (', &
                                                statidx_l(1,v,1),',', &
                                                statidx_l(2,v,1),',', &
                                                statidx_l(3,v,1),')'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,3(I5,A))') '*** MIN = ', &
                                                statval_l(  v,2),' (', &
                                                statidx_l(1,v,2),',', &
                                                statidx_l(2,v,2),',', &
                                                statidx_l(3,v,2),')'
       enddo
    endif

    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine STATISTICS_detail_3D

  subroutine STATISTICS_detail_2D( &
       IA, IS, IE, JA, JS, JE, VA, &
       varname, var, &
       local         )
    use scale_prc, only: &
       PRC_nprocs,   &
       PRC_myrank
    use scale_const, only: &
       CONST_UNDEF, &
       CONST_UNDEF2
    use scale_comm, only: &
       COMM_datatype
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: VA

    character(len=*), intent(in) :: varname(VA)   !< name of item
    real(RP),         intent(in) :: var(IA,JA,VA) !< values

    logical,          intent(in), optional :: local ! calc in local node

    real(RP) :: statval_l (  VA,2)
    integer  :: statidx_l (2,VA,2)
    real(RP) :: statval   (  VA,2,0:PRC_nprocs-1)
    integer  :: statidx   (2,VA,2,0:PRC_nprocs-1)
    real(RP) :: allstatval(VA,2)
    integer  :: allstatidx(VA,2)
    logical :: do_globalcomm

    integer :: i, j
    integer :: ierr
    integer :: v, p
    !---------------------------------------------------------------------------

    do_globalcomm = STATISTICS_use_globalcomm
    if ( present(local) ) do_globalcomm = ( .not. local )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Variable Statistics ***'
    do v = 1, VA
       statval_l(  v,:) = var(IS,JS,v)
       statidx_l(1,v,:) = IS
       statidx_l(2,v,:) = JS
       do j = JS, JE
       do i = IS, IE
          if ( var(i,j,v) > statval_l(v,1) ) then
             statval_l(  v,1) = var(i,j,v)
             statidx_l(1,v,1) = i
             statidx_l(2,v,1) = j
          end if
          if ( var(i,j,v) < statval_l(v,2) ) then
             statval_l(  v,2) = var(i,j,v)
             statidx_l(1,v,2) = i
             statidx_l(2,v,2) = j
          end if
       end do
       end do
    enddo

    if ( do_globalcomm ) then
       call PROF_rapstart('COMM_Bcast', 2)

       call MPI_AllGather( statval_l(:,:),       &
                           VA*2,                 &
                           COMM_datatype,        &
                           statval(:,:,:),       &
                           VA*2,                 &
                           COMM_datatype,        &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr )

       call MPI_AllGather( statidx_l(:,:,:),     &
                           2*VA*2,               &
                           MPI_INTEGER,          &
                           statidx(:,:,:,:),     &
                           2*VA*2,               &
                           MPI_INTEGER,          &
                           PRC_LOCAL_COMM_WORLD, &
                           ierr )

       call PROF_rapend  ('COMM_Bcast', 2)

       do v = 1, VA
          allstatval(v,1) = statval(v,1,0)
          allstatval(v,2) = statval(v,2,0)
          allstatidx(v,:) = 0
          do p = 1, PRC_nprocs-1
             if ( statval(v,1,p) > allstatval(v,1) ) then
                allstatval(v,1) = statval(v,1,p)
                allstatidx(v,1) = p
             end if
             if ( statval(v,2,p) < allstatval(v,2) ) then
                allstatval(v,2) = statval(v,2,p)
                allstatidx(v,2) = p
             end if
          end do
          if( IO_L ) write(IO_FID_LOG,*) '[', trim(varname(v)), ']'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,3(I5,A))') '  MAX =', &
                                                       allstatval(v,1), ' (rank=', &
                                                       allstatidx(v,1), '; ', &
                                         statidx(1,v,1,allstatidx(v,1)),',', &
                                         statidx(2,v,1,allstatidx(v,1)),')'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,3(I5,A))') '  MIN =', &
                                                       allstatval(v,2), ' (rank=', &
                                                       allstatidx(v,2), '; ', &
                                         statidx(1,v,2,allstatidx(v,2)),',', &
                                         statidx(2,v,2,allstatidx(v,2)),')'
       enddo
    else
       ! statistics on each node
       do v = 1, VA
          if( IO_L ) write(IO_FID_LOG,*) '*** [', trim(varname(v)), ']'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,2(I5,A))') '*** MAX = ', &
                                                statval_l(  v,1),' (', &
                                                statidx_l(1,v,1),',', &
                                                statidx_l(2,v,1),')'
          if( IO_L ) write(IO_FID_LOG,'(1x,A,ES17.10,A,2(I5,A))') '*** MIN = ', &
                                                statval_l(  v,2),' (', &
                                                statidx_l(1,v,2),',', &
                                                statidx_l(2,v,2),')'
       enddo
    endif

    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine STATISTICS_detail_2D

end module scale_statistics
