!-------------------------------------------------------------------------------
!> module Statistics
!!
!! @par Description
!!          global statistics module
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_statistics
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi ! TODO: to replace functions in scale_comm module
  use scale_precision
  use scale_io
  use scale_prof
  use scale_const, only: &
     EPS => CONST_EPS
  use scale_prc, only: &
     PRC_LOCAL_COMM_WORLD
#ifdef _OPENACC
  use openacc
#endif
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
  public :: STATISTICS_horizontal_mean
  public :: STATISTICS_horizontal_min
  public :: STATISTICS_horizontal_max

  public :: STATISTICS_summation
  public :: STATISTICS_average
  public :: STATISTICS_variance
  public :: STATISTICS_stddev

  public :: STATISTICS_covariance
  public :: STATISTICS_correlation
  public :: STATISTICS_regression
  public :: STATISTICS_lag_correlation
  public :: STATISTICS_partial_correlation

  public :: STATISTICS_undef_replace
  public :: STATISTICS_undef_embed
  public :: STATISTICS_undef_arraysize

  interface STATISTICS_total
     module procedure STATISTICS_total_2D
     module procedure STATISTICS_total_3D
  end interface STATISTICS_total

  interface STATISTICS_detail
     module procedure STATISTICS_detail_2D
     module procedure STATISTICS_detail_3D
  end interface STATISTICS_detail

  interface STATISTICS_horizontal_mean
     module procedure STATISTICS_horizontal_mean_2D
     module procedure STATISTICS_horizontal_mean_3D
  end interface STATISTICS_horizontal_mean

  interface STATISTICS_horizontal_max
     module procedure STATISTICS_horizontal_max_2D
     module procedure STATISTICS_horizontal_max_3D
  end interface STATISTICS_horizontal_max

  interface STATISTICS_horizontal_min
     module procedure STATISTICS_horizontal_min_2D
     module procedure STATISTICS_horizontal_min_3D
  end interface STATISTICS_horizontal_min

  interface STATISTICS_summation
    module procedure STATISTICS_summation_1D
    module procedure STATISTICS_summation_2D
    module procedure STATISTICS_summation_3D
  end interface STATISTICS_summation

  interface STATISTICS_average
    module procedure STATISTICS_average_1D
    module procedure STATISTICS_average_2D
    module procedure STATISTICS_average_3D
  end interface STATISTICS_average

  interface STATISTICS_variance
    module procedure STATISTICS_variance_1D
    module procedure STATISTICS_variance_2D
    module procedure STATISTICS_variance_3D
  end interface STATISTICS_variance

  interface STATISTICS_stddev
    module procedure STATISTICS_stddev_1D
    module procedure STATISTICS_stddev_2D
    module procedure STATISTICS_stddev_3D
  end interface STATISTICS_stddev

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
    use scale_comm_cartesC, only: &
       COMM_setup
    implicit none

    namelist / PARAM_STATISTICS / &
       STATISTICS_checktotal, &
       STATISTICS_use_globalcomm

    integer :: ierr
    !---------------------------------------------------------------------------

    call COMM_setup

    LOG_NEWLINE
    LOG_INFO("STATISTICS_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_STATISTICS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("STATISTICS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("STATISTICS_setup",*) 'Not appropriate names in namelist PARAM_STATISTICS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_STATISTICS)

    LOG_NEWLINE
    LOG_INFO("STATISTICS_setup",*) 'Caluculate total statistics for monitoring? : ', STATISTICS_checktotal
    if ( STATISTICS_use_globalcomm ) then
       LOG_INFO_CONT(*) '=> The total is calculated for the global domain.'
    else
       LOG_INFO_CONT(*) '=> The total is calculated for the local domain.'
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
       global,       &
       mean, sum     )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       EPS   => CONST_EPS, &
       UNDEF => CONST_UNDEF
    implicit none

    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP),         intent(in) :: var(IA,JA)  !< 3D value
    character(len=*), intent(in) :: varname     !< name of item
    real(RP),         intent(in) :: area(IA,JA) !< area of the grid cell
    real(RP),         intent(in) :: total       !< total area

    logical,  intent(in),  optional :: log_suppress !< suppress log output
    logical,  intent(in),  optional :: global       !< global or local sum
    real(RP), intent(out), optional :: mean !< area-weighted mean
    real(DP), intent(out), optional :: sum  !< domain sum

    real(DP) :: statval
    real(DP) :: sendbuf(2), recvbuf(2)
    real(DP) :: sum_, mean_

    logical :: suppress_, global_
    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    statval = 0.0_DP
    !$acc update host(var(IS,JS)) if(acc_is_present(var))
    if ( var(IS,JS) /= UNDEF ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2) reduction(+:statval)
       !$acc kernels copyin(var, area) if(acc_is_present(var))
       !$acc loop reduction(+:statval)
       do j = JS, JE
       !$acc loop reduction(+:statval)
       do i = IS, IE
          statval = statval + var(i,j) * area(i,j)
       end do
       end do
       !$acc end kernels
    end if

    if ( .NOT. ( statval > -1.0_DP .OR. statval < 1.0_DP ) ) then ! must be NaN
       LOG_ERROR("STATISTICS_total_2D",*) 'NaN is detected for ', trim(varname), ' in rank ', PRC_myrank
       call PRC_abort
    endif

    if ( present(log_suppress) ) then
       suppress_ = log_suppress
    else
       suppress_ = .false.
    end if

    if ( present(global) ) then
       global_ = global
    else
       global_ = STATISTICS_use_globalcomm
    end if

    if ( global_ ) then
       call PROF_rapstart('COMM_Allreduce', 2)
       sendbuf(1) = statval
       sendbuf(2) = total
       ! All reduce
       call MPI_Allreduce( sendbuf(:), recvbuf(:), &
                           2,                      &
                           MPI_DOUBLE_PRECISION,   &
                           MPI_SUM,                &
                           PRC_LOCAL_COMM_WORLD,   &
                           ierr                    )
       call PROF_rapend  ('COMM_Allreduce', 2)

       if ( recvbuf(2) < EPS ) then
          sum_  = UNDEF
          mean_ = UNDEF
       else
          sum_  = recvbuf(1)
          mean_ = recvbuf(1) / recvbuf(2)
       end if
       ! statistics over the all node
       if ( .not. suppress_ ) then ! if varname is empty, suppress output
          LOG_INFO("STATISTICS_total_2D",'(1x,A,A24,A,ES24.17)') &
                     '[', trim(varname), '] MEAN(global) = ', mean_
       endif
    else
       if ( total < EPS ) then
          sum_  = UNDEF
          mean_ = UNDEF
       else
          sum_  = statval
          mean_ = statval / total
       end if

       ! statistics on each node
       if ( .not. suppress_ ) then ! if varname is empty, suppress output
          LOG_INFO("STATISTICS_total_2D",'(1x,A,A24,A,ES24.17)') &
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
       global,       &
       mean, sum     )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       EPS   => CONST_EPS , &
       UNDEF => CONST_UNDEF
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP),         intent(in) :: var(KA,IA,JA) !< 3D value
    character(len=*), intent(in) :: varname       !< name of item
    real(RP),         intent(in) :: vol(KA,IA,JA) !< volume of the grid cell
    real(RP),         intent(in) :: total         !< total volume

    logical,  intent(in),  optional :: log_suppress !< suppress log output
    logical,  intent(in),  optional :: global       !< global or local sum
    real(RP), intent(out), optional :: mean !< volume/area-weighted total
    real(DP), intent(out), optional :: sum  !< domain sum

    real(DP) :: statval
    real(DP) :: sendbuf(2), recvbuf(2)
    real(DP) :: mean_, sum_

    logical :: suppress_, global_

    real(DP) :: work
    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    statval = 0.0_DP
    !$acc update host(var(KS,IS,JS)) if(acc_is_present(var))
    if ( var(KS,IS,JS) /= UNDEF ) then
       !$omp parallel do OMP_SCHEDULE_ reduction(+:statval) &
       !$omp private(work)
       !$acc kernels copyin(var, vol)  if(acc_is_present(var))
       !$acc loop reduction(statval)
       do j = JE, JS, -1
       !$acc loop reduction(statval)
       do i = IE, IS, -1
#ifdef _OPENACC
          !$acc loop reduction(statval)
          do k = KE, KS, -1
             statval = statval + var(k,i,j) * vol(k,i,j)
          enddo
#else
          work = 0.0_RP
          do k = KE, KS, -1
             work = work + var(k,i,j) * vol(k,i,j)
          enddo
          statval = statval + work
#endif
       enddo
       enddo
       !$acc end kernels
    end if

    if ( .NOT. ( statval > -1.0_DP .OR. statval < 1.0_DP ) ) then ! must be NaN
       LOG_ERROR("STATISTICS_total_3D",*) 'NaN is detected for ', trim(varname), ' in rank ', PRC_myrank
       call PRC_abort
    endif

    if ( present(log_suppress) ) then
       suppress_ = log_suppress
    else
       suppress_ = .false.
    end if

    if ( present(global) ) then
       global_ = global
    else
       global_ = STATISTICS_use_globalcomm
    end if

    if ( global_ ) then
       call PROF_rapstart('COMM_Allreduce', 2)
       sendbuf(1) = statval
       sendbuf(2) = total
       ! All reduce
       call MPI_Allreduce( sendbuf(:), recvbuf(:), &
                           2,                      &
                           MPI_DOUBLE_PRECISION,   &
                           MPI_SUM,                &
                           PRC_LOCAL_COMM_WORLD,   &
                           ierr                    )
       call PROF_rapend  ('COMM_Allreduce', 2)

       if ( recvbuf(2) < EPS ) then
          sum_  = UNDEF
          mean_ = UNDEF
       else
          sum_  = recvbuf(1)
          mean_ = recvbuf(1) / recvbuf(2)
       end if
       ! statistics over the all node
       if ( .not. suppress_ ) then ! if varname is empty, suppress output
          LOG_INFO("STATISTICS_total_3D",'(1x,A,A24,A,ES24.17)') &
                     '[', trim(varname), '] MEAN(global) = ', mean_
       endif
    else
       if ( total < EPS ) then
          sum_  = UNDEF
          mean_ = UNDEF
       else
          sum_  = statval
          mean_ = statval / total
       end if

       ! statistics on each node
       if ( .not. suppress_ ) then ! if varname is empty, suppress output
          LOG_INFO("STATISTICS_total_3D",'(1x,A,A24,A,ES24.17)') &
                     '[', trim(varname), '] MEAN(local)  = ', mean_
       endif
    endif

    if ( present(mean) ) mean = mean_
    if ( present(sum ) ) sum  = sum_

    return
  end subroutine STATISTICS_total_3D

  !-----------------------------------------------------------------------------
  !> Calc horizontal mean value
  subroutine STATISTICS_horizontal_mean_2D( &
       IA, IS, IE, JA, JS, JE, &
       var, area, &
       varmean    )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in)  :: var (IA,JA)
    real(RP), intent(in)  :: area(IA,JA)
    real(RP), intent(out) :: varmean

    real(DP) :: statval   (2)
    real(DP) :: allstatval(2)
    real(DP) :: s1, s2

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    s1 = 0.0_DP
    s2 = 0.0_DP
    !$omp parallel do reduction(+:s1,s2)
    !$acc kernels copyin(area, var) if(acc_is_present(var))
    !$acc loop reduction(+:s1,s2)
    do j = JS, JE
    !$acc loop reduction(+:s1,s2)
    do i = IS, IE
       if ( var(i,j) /= UNDEF ) then
          s1 = s1 + area(i,j) * var(i,j)
          s2 = s2 + area(i,j)
       endif
    enddo
    enddo
    !$acc end kernels
    statval(1) = s1
    statval(2) = s2

    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    call MPI_Allreduce( statval(:),           &
                        allstatval(:),        &
                        2,                    &
                        MPI_DOUBLE_PRECISION, &
                        MPI_SUM,              &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    call PROF_rapend  ('COMM_Allreduce', 2)

    if ( allstatval(2) > 0.0_DP ) then
       varmean = allstatval(1) / allstatval(2)
    else
       varmean = UNDEF
    end if

    return
  end subroutine STATISTICS_horizontal_mean_2D

  subroutine STATISTICS_horizontal_mean_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       var, area, &
       varmean    )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: var (KA,IA,JA)
    real(RP), intent(in)  :: area(   IA,JA)
    real(RP), intent(out) :: varmean(KA)

    real(DP) :: statval   (KS:KE,2)
    real(DP) :: allstatval(KS:KE,2)

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$acc data copyin(var, area) copyout(varmean) create(statval, allstatval) if(acc_is_present(var))

    !$acc kernels if(acc_is_present(var))
    statval(:,:) = 0.0_DP
    !$acc end kernels

!    !$omp parallel do reduction(+:statval)
    !$acc kernels if(acc_is_present(var))
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent
    do i = IS, IE
    do k = KS, KE
       if ( var(k,i,j) /= UNDEF ) then
          !$acc atomic update
          statval(k,1) = statval(k,1) + area(i,j) * var(k,i,j)
          !$acc end atomic
          !$acc atomic update
          statval(k,2) = statval(k,2) + area(i,j)
          !$acc end atomic
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    !$acc host_data use_device(statval, allstatval) if(acc_is_present(var))
    call MPI_Allreduce( statval   (:,:),      &
                        allstatval(:,:),      &
                        (KE-KS+1)*2,          &
                        MPI_DOUBLE_PRECISION, &
                        MPI_SUM,              &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    !$acc end host_data
    call PROF_rapend  ('COMM_Allreduce', 2)

    !$acc kernels if(acc_is_present(var))
    do k = KS, KE
       if ( allstatval(k,2) > 0.0_DP ) then
          varmean(k) = allstatval(k,1) / allstatval(k,2)
       else
          varmean(k) = UNDEF
       end if
    enddo
    !$acc end kernels
    !$acc kernels if(acc_is_present(var))
    do k = 1, KS-1
       varmean(k) = UNDEF
    end do
    !$acc end kernels
    !$acc kernels if(acc_is_present(var))
    do k = KE+1, KA
       varmean(k) = UNDEF
    end do
    !$acc end kernels

    !$acc end data

    return
  end subroutine STATISTICS_horizontal_mean_3D

  !-----------------------------------------------------------------------------
  !> Calc horizontal minimum value
  subroutine STATISTICS_horizontal_min_2D( &
       IA, IS, IE, JA, JS, JE, &
       var, &
       varmin )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       HUGE  => CONST_HUGE
    use scale_comm_cartesC, only: &
       COMM_datatype
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: var(IA,JA)
    real(RP), intent(out) :: varmin

    real(RP) :: statval
    real(RP) :: allstatval

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    statval = HUGE
    !$omp parallel do reduction(min:statval)
    !$acc kernels copyin(var) if(acc_is_present(var))
    !$acc loop reduction(min:statval)
    do j = JS, JE
    !$acc loop reduction(min:statval)
    do i = IS, IE
       if ( var(i,j) /= UNDEF .and. var(i,j) < statval ) then
          statval = var(i,j)
       endif
    enddo
    enddo
    !$acc end kernels

    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    call MPI_Allreduce( statval,              &
                        allstatval,           &
                        1,                    &
                        COMM_datatype,        &
                        MPI_MIN,              &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    call PROF_rapend  ('COMM_Allreduce', 2)

    if ( allstatval < HUGE ) then
       varmin = allstatval
    else
       varmin = UNDEF
    end if

    return
  end subroutine STATISTICS_horizontal_min_2D

  subroutine STATISTICS_horizontal_min_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       var, &
       varmin )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       HUGE  => CONST_HUGE
    use scale_comm_cartesC, only: &
       COMM_datatype
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: var(KA,IA,JA)
    real(RP), intent(out) :: varmin(KA)

    real(RP) :: statval   (KA)
    real(RP) :: allstatval(KA)

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$acc data copyin(var) copyout(varmin) create(statval, allstatval) if(acc_is_present(var))

    !$acc kernels if(acc_is_present(var))
    statval(:) = HUGE
    !$acc end kernels

!    !$omp parallel do reduction(min:statval)
    !$acc kernels if(acc_is_present(var))
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent
    do i = IS, IE
    do k = KS, KE
#ifdef _OPENACC
       if ( var(k,i,j) /= UNDEF )then
          !$acc atomic update
          statval(k) = min( var(k,i,j), statval(k) )
          !$acc end atomic
       endif
#else
       if ( var(k,i,j) /= UNDEF .and. var(k,i,j) < statval(k) ) then
          statval(k) = var(k,i,j)
       endif
#endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    !$acc host_data use_device(statval, allstatval) if(acc_is_present(var))
    call MPI_Allreduce( statval   (KS:KE),    &
                        allstatval(KS:KE),    &
                        KE-KS+1,              &
                        COMM_datatype,        &
                        MPI_MIN,              &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    !$acc end host_data
    call PROF_rapend  ('COMM_Allreduce', 2)

    !$acc kernels if(acc_is_present(var))
    do k = KS, KE
       if ( allstatval(k) < HUGE ) then
          varmin(k) = allstatval(k)
       else
          varmin(k) = UNDEF
       end if
    enddo
    !$acc end kernels
    !$acc kernels if(acc_is_present(var))
    do k = 1, KS-1
       varmin(k) = UNDEF
    end do
    !$acc end kernels
    !$acc kernels if(acc_is_present(var))
    do k = KE+1, KA
       varmin(k) = UNDEF
    end do
    !$acc end kernels

    !$acc end data

    return
  end subroutine STATISTICS_horizontal_min_3D

  !-----------------------------------------------------------------------------
  !> Calc horizontal maximum value
  subroutine STATISTICS_horizontal_max_2D( &
       IA, IS, IE, JA, JS, JE, &
       var, &
       varmax )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       HUGE  => CONST_HUGE
    use scale_comm_cartesC, only: &
       COMM_datatype
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: var(IA,JA)
    real(RP), intent(out) :: varmax

    real(RP) :: statval
    real(RP) :: allstatval

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    statval = - HUGE
    !$omp parallel do reduction(max:statval)
    !$acc kernels copyin(var) if(acc_is_present(var))
    !$acc loop reduction(max:statval)
    do j = JS, JE
    !$acc loop reduction(max:statval)
    do i = IS, IE
       if ( var(i,j) /= UNDEF .and. var(i,j) > statval ) then
          statval = var(i,j)
       endif
    enddo
    enddo
    !$acc end kernels

    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    call MPI_Allreduce( statval,              &
                        allstatval,           &
                        1,                    &
                        COMM_datatype,        &
                        MPI_MAX,              &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    call PROF_rapend  ('COMM_Allreduce', 2)

    if ( allstatval > - HUGE ) then
       varmax = allstatval
    else
       varmax = UNDEF
    end if

    return
  end subroutine STATISTICS_horizontal_max_2D

  subroutine STATISTICS_horizontal_max_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       var, &
       varmax )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       HUGE  => CONST_HUGE
    use scale_comm_cartesC, only: &
       COMM_datatype
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: var(KA,IA,JA)
    real(RP), intent(out) :: varmax(KA)

    real(RP) :: statval   (KA)
    real(RP) :: allstatval(KA)

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$acc data copyin(var) copyout(varmax) create(statval, allstatval) if(acc_is_present(var))

    !$acc kernels if(acc_is_present(var))
    statval(:) = - HUGE
    !$acc end kernels

!    !$omp parallel do reduction(max:statval)
    !$acc kernels if(acc_is_present(var))
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent
    do i = IS, IE
    do k = KS, KE
#ifdef _OPENACC
       if ( var(k,i,j) /= UNDEF ) then
          !$acc atomic update
          statval(k) = max( var(k,i,j), statval(k) )
          !$acc end atomic
       endif
#else
       if ( var(k,i,j) /= UNDEF .and. var(k,i,j) > statval(k) ) then
          statval(k) = var(k,i,j)
       endif
#endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    !$acc host_data use_device(statval, allstatval) if(acc_is_present(var))
    call MPI_Allreduce( statval   (KS:KE),    &
                        allstatval(KS:KE),    &
                        KE-KS+1,              &
                        COMM_datatype,        &
                        MPI_MAX,              &
                        PRC_LOCAL_COMM_WORLD, &
                        ierr                  )
    !$acc end host_data
    call PROF_rapend  ('COMM_Allreduce', 2)

    !$acc kernels if(acc_is_present(var))
    do k = KS, KE
       if ( allstatval(k) > - HUGE ) then
          varmax(k) = allstatval(k)
       else
          varmax(k) = UNDEF
       end if
    enddo
    !$acc end kernels
    !$acc kernels if(acc_is_present(var))
    do k = 1, KS-1
       varmax(k) = UNDEF
    end do
    !$acc end kernels
    !$acc kernels if(acc_is_present(var))
    do k = KE+1, KA
       varmax(k) = UNDEF
    end do
    !$acc end kernels

    !$acc end data

    return
  end subroutine STATISTICS_horizontal_max_3D

  !-----------------------------------------------------------------------------
  !> Search global maximum & minimum value
  subroutine STATISTICS_detail_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, VA, &
       varname, var, &
       local         )
    use scale_prc, only: &
       PRC_nprocs
    use scale_comm_cartesC, only: &
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

    !$acc update host(var) if ( acc_is_present(var) )

    LOG_NEWLINE
    LOG_INFO("STATISTICS_detail_3D",*) 'Variable Statistics '
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
          LOG_INFO_CONT(*) '[', trim(varname(v)), ']'
          LOG_INFO_CONT('(1x,A,ES17.10,A,4(I5,A))') '  MAX =', &
                                                       allstatval(v,1), ' (rank=', &
                                                       allstatidx(v,1), '; ', &
                                         statidx(1,v,1,allstatidx(v,1)),',', &
                                         statidx(2,v,1,allstatidx(v,1)),',', &
                                         statidx(3,v,1,allstatidx(v,1)),')'
          LOG_INFO_CONT('(1x,A,ES17.10,A,4(I5,A))') '  MIN =', &
                                                       allstatval(v,2), ' (rank=', &
                                                       allstatidx(v,2), '; ', &
                                         statidx(1,v,2,allstatidx(v,2)),',', &
                                         statidx(2,v,2,allstatidx(v,2)),',', &
                                         statidx(3,v,2,allstatidx(v,2)),')'
       enddo
    else
       ! statistics on each node
       do v = 1, VA
          LOG_INFO_CONT(*) '[', trim(varname(v)), ']'
          LOG_INFO_CONT('(1x,A,ES17.10,A,3(I5,A))') 'MAX = ', &
                                                statval_l(  v,1),' (', &
                                                statidx_l(1,v,1),',', &
                                                statidx_l(2,v,1),',', &
                                                statidx_l(3,v,1),')'
          LOG_INFO_CONT('(1x,A,ES17.10,A,3(I5,A))') 'MIN = ', &
                                                statval_l(  v,2),' (', &
                                                statidx_l(1,v,2),',', &
                                                statidx_l(2,v,2),',', &
                                                statidx_l(3,v,2),')'
       enddo
    endif

    LOG_NEWLINE

    return
  end subroutine STATISTICS_detail_3D

  subroutine STATISTICS_detail_2D( &
       IA, IS, IE, JA, JS, JE, VA, &
       varname, var, &
       local         )
    use scale_prc, only: &
       PRC_nprocs
    use scale_comm_cartesC, only: &
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

    !$acc update host(var) if ( acc_is_present(var) )

    LOG_NEWLINE
    LOG_INFO("STATISTICS_detail_2D",*) 'Variable Statistics '
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
          LOG_INFO_CONT(*) '[', trim(varname(v)), ']'
          LOG_INFO_CONT('(1x,A,ES17.10,A,3(I5,A))') '  MAX =', &
                                                       allstatval(v,1), ' (rank=', &
                                                       allstatidx(v,1), '; ', &
                                         statidx(1,v,1,allstatidx(v,1)),',', &
                                         statidx(2,v,1,allstatidx(v,1)),')'
          LOG_INFO_CONT('(1x,A,ES17.10,A,3(I5,A))') '  MIN =', &
                                                       allstatval(v,2), ' (rank=', &
                                                       allstatidx(v,2), '; ', &
                                         statidx(1,v,2,allstatidx(v,2)),',', &
                                         statidx(2,v,2,allstatidx(v,2)),')'
       enddo
    else
       ! statistics on each node
       do v = 1, VA
          LOG_INFO_CONT(*) '[', trim(varname(v)), ']'
          LOG_INFO_CONT('(1x,A,ES17.10,A,2(I5,A))') 'MAX = ', &
                                                statval_l(  v,1),' (', &
                                                statidx_l(1,v,1),',', &
                                                statidx_l(2,v,1),')'
          LOG_INFO_CONT('(1x,A,ES17.10,A,2(I5,A))') 'MIN = ', &
                                                statval_l(  v,2),' (', &
                                                statidx_l(1,v,2),',', &
                                                statidx_l(2,v,2),')'
       enddo
    endif

    LOG_NEWLINE

    return
  end subroutine STATISTICS_detail_2D

  !-----------------------------------------------------------------------------
  !
  ! Summation (1D):
  !   Compenstaed summation of ARRAY (1D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_summation_1D( &
      IA,    & ! (in)
      ARRAY, & ! (in)
      UNDEF  ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA
    real(RP), intent(in) :: ARRAY(IA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_summation_1D

    ! work
    integer :: i

    real(RP) :: tmp_ARRAY(IA)
    real(RP) :: tmp
    real(RP) :: res
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( any( abs( ARRAY(:) - UNDEF ) > EPS ) ) then
        where( abs( ARRAY(:) - UNDEF ) > EPS )
          tmp_ARRAY(:) = ARRAY(:)
        else where
          tmp_ARRAY(:) = 0.0_RP
        end where
      else
        STATISTICS_summation_1D = UNDEF
        return ! end function
      end if
    else
      tmp_ARRAY(:) = ARRAY(:)
    end if

    STATISTICS_summation_1D = 0.0_RP

    tmp = 0.0_RP
    res = 0.0_RP

    ! Kahan's Compensated Summation (Kahan 1965)
    do i = 1, IA
      tmp = STATISTICS_summation_1D + ( tmp_ARRAY(i) + res )
      res = ( tmp_ARRAY(i) + res ) - ( tmp - STATISTICS_summation_1D )
      STATISTICS_summation_1D = tmp
    end do

    return
  end function STATISTICS_summation_1D

  !-----------------------------------------------------------------------------
  !
  ! Summation (2D):
  !   Compenstaed summation of ARRAY (2D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_summation_2D( &
      IA, JA, & ! (in)
      ARRAY,  & ! (in)
      UNDEF   ) ! (in)
    implicit none

    ! arguments
    integer , intent(in) :: IA
    integer , intent(in) :: JA
    real(RP), intent(in) :: ARRAY(IA,JA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_summation_2D

    ! work
    integer :: i, j

    real(RP) :: tmp_ARRAY(IA,JA)
    real(RP) :: tmp
    real(RP) :: res
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( any( abs( ARRAY(:,:) - UNDEF ) > EPS ) ) then
        where( abs( ARRAY(:,:) - UNDEF ) > EPS )
          tmp_ARRAY(:,:) = ARRAY(:,:)
        else where
          tmp_ARRAY(:,:) = 0.0_RP
        end where
      else
        STATISTICS_summation_2D = UNDEF
        return ! end function
      end if
    else
      tmp_ARRAY(:,:) = ARRAY(:,:)
    end if

    STATISTICS_summation_2D = 0.0_RP

    tmp = 0.0_RP
    res = 0.0_RP

    ! Kahan's Compensated Summation (Kahan 1965)
    do j = 1, JA
    do i = 1, IA
      tmp = STATISTICS_summation_2D + ( tmp_ARRAY(i,j) + res )
      res = ( tmp_ARRAY(i,j) + res ) - ( tmp - STATISTICS_summation_2D )
      STATISTICS_summation_2D = tmp
    end do
    end do

    return
  end function STATISTICS_summation_2D

  !-----------------------------------------------------------------------------
  !
  ! Summation (3D):
  !   Compenstaed summation of ARRAY (3D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_summation_3D( &
      KA, IA, JA, & ! (in)
      ARRAY,      & ! (in)
      UNDEF       ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: KA
    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    real(RP), intent(in) :: ARRAY(KA,IA,JA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_summation_3D

    ! work
    integer :: i, j, k

    real(RP) :: tmp_ARRAY(KA,IA,JA)
    real(RP) :: tmp
    real(RP) :: res
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( any( abs( ARRAY(:,:,:) - UNDEF ) > EPS ) ) then
        where( abs( ARRAY(:,:,:) - UNDEF ) > EPS )
          tmp_ARRAY(:,:,:) = ARRAY(:,:,:)
        else where
          tmp_ARRAY(:,:,:) = 0.0_RP
        end where
      else
        STATISTICS_summation_3D = UNDEF
        return ! end function
      end if
    else
      tmp_ARRAY(:,:,:) = ARRAY(:,:,:)
    end if

    STATISTICS_summation_3D = 0.0_RP

    tmp = 0.0_RP
    res = 0.0_RP

    ! Kahan's Compensated Summation (Kahan 1965)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
      tmp = STATISTICS_summation_3D + ( tmp_ARRAY(i,j,k) + res )
      res = ( tmp_ARRAY(i,j,k) + res ) - ( tmp - STATISTICS_summation_3D )
      STATISTICS_summation_3D = tmp
    end do
    end do
    end do

    return
  end function STATISTICS_summation_3D

  !-----------------------------------------------------------------------------
  !
  ! Average (1D):
  !   average of ARRAY (1D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_average_1D( &
      IA,    & ! (in)
      ARRAY, & ! (in)
      UNDEF  ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA
    real(RP), intent(in) :: ARRAY(IA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_average_1D
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( any( abs( ARRAY(:) - UNDEF ) > EPS ) ) then
        STATISTICS_average_1D = STATISTICS_summation( IA, ARRAY(:), UNDEF ) &
                              / real( count( abs( ARRAY(:) - UNDEF ) > EPS ), kind=RP )
      else
        STATISTICS_average_1D = UNDEF
      end if
    else
      STATISTICS_average_1D = STATISTICS_summation( IA, ARRAY(:) ) &
                            / real( IA, kind=RP )
    end if

    return
  end function STATISTICS_average_1D

  !-----------------------------------------------------------------------------
  !
  ! Average (2D):
  !   average of ARRAY (2D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_average_2D( &
      IA, JA, & ! (in)
      ARRAY,  & ! (in)
      UNDEF   ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    real(RP), intent(in) :: ARRAY(IA,JA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_average_2D
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( any( abs( ARRAY(:,:) - UNDEF ) > EPS ) ) then
        STATISTICS_average_2D = STATISTICS_summation( IA, JA, ARRAY(:,:), UNDEF ) &
                              / real( count( abs( ARRAY(:,:) - UNDEF ) > EPS ), kind=RP )
      else
        STATISTICS_average_2D = UNDEF
      end if
    else
      STATISTICS_average_2D = STATISTICS_summation( IA, JA, ARRAY(:,:) ) &
                            / real( IA*JA, kind=RP )
    end if

    return
  end function STATISTICS_average_2D

  !-----------------------------------------------------------------------------
  !
  ! Average (3D):
  !   average of ARRAY (3D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_average_3D( &
      KA, IA, JA, & ! (in)
      ARRAY,      & ! (in)
      UNDEF       ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: KA
    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    real(RP), intent(in) :: ARRAY(KA,IA,JA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_average_3D
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( any( abs( ARRAY(:,:,:) - UNDEF ) > EPS ) ) then
        STATISTICS_average_3D = STATISTICS_summation( KA, IA, JA, ARRAY(:,:,:), UNDEF ) &
                              / real( count( abs( ARRAY(:,:,:) - UNDEF ) > EPS ), kind=RP )
      else
        STATISTICS_average_3D = UNDEF
      end if
    else
      STATISTICS_average_3D = STATISTICS_summation( KA, IA, JA, ARRAY(:,:,:) ) &
                            / real( KA*IA*JA, kind=RP )
    end if

    return
  end function STATISTICS_average_3D

  !-----------------------------------------------------------------------------
  !
  ! Variance (1D):
  !   variance of ARRAY (1D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_variance_1D( &
      IA,    & ! (in)
      ARRAY, & ! (in)
      UNDEF  ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA
    real(RP), intent(in) :: ARRAY(IA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_variance_1D

    ! work
    real(RP) :: tmp_ARRAY(IA)
    real(RP) :: tmp
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( count( abs( ARRAY(:) - UNDEF ) > EPS ) > 1 ) then
        tmp = STATISTICS_average( IA, ARRAY(:), UNDEF )

        where( abs( ARRAY(:) - UNDEF ) > EPS )
          tmp_ARRAY(:) = ( ARRAY(:) - tmp )**2
        else where
          tmp_ARRAY(:) = UNDEF
        end where

        STATISTICS_variance_1D = STATISTICS_summation( IA, tmp_ARRAY(:), UNDEF ) &
                               / real( count( abs( ARRAY(:) - UNDEF ) > EPS ) - 1, kind=RP )
      else
        STATISTICS_variance_1D = UNDEF
      end if
    else
      STATISTICS_variance_1D = STATISTICS_summation( IA, ( ARRAY(:) - STATISTICS_average( IA, ARRAY(:) ) )**2 ) &
                             / real( IA - 1, kind=RP )
    end if

    return
  end function STATISTICS_variance_1D

  !-----------------------------------------------------------------------------
  !
  ! Variance (2D):
  !   variance of ARRAY (2D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_variance_2D( &
      IA, JA, & ! (in)
      ARRAY,  & ! (in)
      UNDEF   ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    real(RP), intent(in) :: ARRAY(IA,JA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_variance_2D

    ! work
    real(RP) :: tmp_ARRAY(IA,JA)
    real(RP) :: tmp
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( count( abs( ARRAY(:,:) - UNDEF ) > EPS ) > 1 ) then
        tmp = STATISTICS_average( IA, JA, ARRAY(:,:), UNDEF )

        where( abs( ARRAY(:,:) - UNDEF ) > EPS )
          tmp_ARRAY(:,:) = ( ARRAY(:,:) - tmp )**2
        else where
          tmp_ARRAY(:,:) = UNDEF
        end where

        STATISTICS_variance_2D = STATISTICS_summation( IA, JA, tmp_ARRAY(:,:), UNDEF ) &
                               / real( count( abs( ARRAY(:,:) - UNDEF ) > EPS ) - 1, kind=RP )
      else
        STATISTICS_variance_2D = UNDEF
      end if
    else
      STATISTICS_variance_2D = STATISTICS_summation( IA, JA, ( ARRAY(:,:) - STATISTICS_average( IA, JA, ARRAY(:,:) ) )**2 ) &
                             / real( IA*JA - 1, kind=RP )
    end if

    return
  end function STATISTICS_variance_2D

  !-----------------------------------------------------------------------------
  !
  ! Variance (3D):
  !   variance of ARRAY (3D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_variance_3D( &
      KA, IA, JA, & ! (in)
      ARRAY,      & ! (in)
      UNDEF       ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: KA
    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    real(RP), intent(in) :: ARRAY(KA,IA,JA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_variance_3D

    ! work
    real(RP) :: tmp_ARRAY(KA,IA,JA)
    real(RP) :: tmp
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( count( abs( ARRAY(:,:,:) - UNDEF ) > EPS ) > 1 ) then
        tmp = STATISTICS_average( KA, IA, JA, ARRAY(:,:,:), UNDEF )

        where( abs( ARRAY(:,:,:) - UNDEF ) > EPS )
          tmp_ARRAY(:,:,:) = ( ARRAY(:,:,:) - tmp )**2
        else where
          tmp_ARRAY(:,:,:) = UNDEF
        end where

        STATISTICS_variance_3D = STATISTICS_summation(  KA, IA, JA, tmp_ARRAY(:,:,:), UNDEF ) &
                               / real( count( abs( ARRAY(:,:,:) - UNDEF ) > EPS ) - 1, kind=RP )
      else
        STATISTICS_variance_3D = UNDEF
      end if
    else
      STATISTICS_variance_3D = STATISTICS_summation( KA, IA, JA, ( ARRAY(:,:,:) - STATISTICS_average( KA, IA, JA, ARRAY(:,:,:) ) )**2 ) &
                             / real( KA*IA*JA - 1, kind=RP )
    end if

    return
  end function STATISTICS_variance_3D

  !-----------------------------------------------------------------------------
  !
  ! Standard deviation (1D):
  !   standard deviation of ARRAY (1D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_stddev_1D( &
      IA,    & ! (in)
      ARRAY, & ! (in)
      UNDEF  ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA
    real(RP), intent(in) :: ARRAY(IA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_stddev_1D
    !---------------------------------------------------------------------------

    STATISTICS_stddev_1D = sqrt( STATISTICS_variance( IA, ARRAY(:), UNDEF ) )

    return
  end function STATISTICS_stddev_1D

  !-----------------------------------------------------------------------------
  !
  ! Standard deviation (2D):
  !   standard deviation of ARRAY (2D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_stddev_2D( &
      IA, JA, & ! in)
      ARRAY,  & ! (in)
      UNDEF   ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    real(RP), intent(in) :: ARRAY(IA,JA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_stddev_2D
    !---------------------------------------------------------------------------

    STATISTICS_stddev_2D = sqrt( STATISTICS_variance( IA, JA, ARRAY(:,:), UNDEF ) )

    return
  end function STATISTICS_stddev_2D

  !-----------------------------------------------------------------------------
  !
  ! Standard deviation (3D):
  !   standard deviation of ARRAY (3D)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_stddev_3D( &
      KA, IA, JA, & ! (in)
      ARRAY,      & ! (in)
      UNDEF       ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: KA
    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    real(RP), intent(in) :: ARRAY(KA,IA,JA)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_stddev_3D
    !---------------------------------------------------------------------------

    STATISTICS_stddev_3D = sqrt( STATISTICS_variance( KA, IA, JA, ARRAY(:,:,:), UNDEF ) )

    return
  end function STATISTICS_stddev_3D

  !-----------------------------------------------------------------------------
  !
  ! Covariance:
  !   covariance between ARRAY1 and ARRAY2
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_covariance( &
      IA1, IA2, & ! (in)
      ARRAY1,   & ! (in)
      ARRAY2,   & ! (in)
      UNDEF     ) ! (in)
    implicit none

    ! argument
    integer,  intent(in) :: IA1, IA2
    real(RP), intent(in) :: ARRAY1(IA1)
    real(RP), intent(in) :: ARRAY2(IA2)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_covariance

    ! works
    real(RP) :: tmp_ARRAY1(IA1)
    real(RP) :: tmp_ARRAY2(IA2)

    real(RP), allocatable :: fixed_ARRAY1(:)
    real(RP), allocatable :: fixed_ARRAY2(:)

    integer :: n
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      n = STATISTICS_undef_arraysize( IA1, IA2, ARRAY1(:), ARRAY2(:), UNDEF )

      if( n > 1 ) then
        ! embed undefined value each other
        tmp_ARRAY1 = STATISTICS_undef_embed( IA1, IA2, ARRAY1(:), ARRAY2(:), UNDEF )
        tmp_ARRAY2 = STATISTICS_undef_embed( IA2, IA1, ARRAY2(:), ARRAY1(:), UNDEF )

        allocate( fixed_ARRAY1(n) )
        allocate( fixed_ARRAY2(n) )

        ! generate arrays without undefined value
        fixed_ARRAY1(:) = pack( tmp_ARRAY1(:), abs( tmp_ARRAY1(:) - UNDEF ) > EPS )
        fixed_ARRAY2(:) = pack( tmp_ARRAY2(:), abs( tmp_ARRAY2(:) - UNDEF ) > EPS )

        STATISTICS_covariance = STATISTICS_summation( n, ( fixed_ARRAY1(:) - STATISTICS_average( n, ARRAY1(:), UNDEF ) ) &
                                                       * ( fixed_ARRAY2(:) - STATISTICS_average( n, ARRAY2(:), UNDEF ) ) ) &
                              / real( n-1, kind=RP )

        deallocate( fixed_ARRAY1 )
        deallocate( fixed_ARRAY2 )
      else
        STATISTICS_covariance = UNDEF
      end if
    else
      n = STATISTICS_undef_arraysize( IA1, IA2, ARRAY1(:), ARRAY2(:) )
      ! simple covariance
      STATISTICS_covariance = STATISTICS_summation( n, ( ARRAY1(:) - STATISTICS_average( n, ARRAY1(:) ) ) &
                                                     * ( ARRAY2(:) - STATISTICS_average( n, ARRAY2(:) ) ) ) &
                            / real( n-1, kind=RP )
    end if

    return
  end function STATISTICS_covariance

  !-----------------------------------------------------------------------------
  !
  ! Correlation coefficient:
  !   correlation coffiecient between ARRAY1 and ARRAY2
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_correlation( &
      IA1, IA2, & ! (in)
      ARRAY1,   & ! (in)
      ARRAY2,   & ! (in)
      UNDEF     ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA1, IA2
    real(RP), intent(in) :: ARRAY1(IA1)
    real(RP), intent(in) :: ARRAY2(IA2)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_correlation
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( abs( STATISTICS_covariance( IA1, IA2, ARRAY1(:), ARRAY2(:), UNDEF ) - UNDEF ) > EPS .and. &
          abs( STATISTICS_stddev( IA1, ARRAY1(:), UNDEF ) - UNDEF ) > EPS                     .and. &
          abs( STATISTICS_stddev( IA2, ARRAY2(:), UNDEF ) - UNDEF ) > EPS                     .and. &
          abs( STATISTICS_stddev( IA1, ARRAY1(:), UNDEF )         ) > 0.0_RP                  .and. &
          abs( STATISTICS_stddev( IA2, ARRAY2(:), UNDEF )         ) > 0.0_RP                        ) then
        ! correlation coefficient without undefined value
        STATISTICS_correlation = STATISTICS_covariance( IA1, IA2, ARRAY1(:), ARRAY2(:), UNDEF ) &
                               / ( STATISTICS_stddev( IA1, ARRAY1(:), UNDEF ) &
                                 * STATISTICS_stddev( IA2, ARRAY2(:), UNDEF ) )
      else
        STATISTICS_correlation = UNDEF
      end if
    else
      ! simple correlation coefficient
      STATISTICS_correlation = STATISTICS_covariance( IA1, IA2, ARRAY1(:), ARRAY2(:) ) &
                             / ( STATISTICS_stddev( IA1, ARRAY1(:) ) &
                               * STATISTICS_stddev( IA2, ARRAY2(:) ) )
    end if

    return
  end function STATISTICS_correlation

  !-----------------------------------------------------------------------------
  !
  ! Regression coefficient:
  !   regression coffiecient of ARRAY1 to ARRAY2
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_regression( &
      IA1, IA2, & ! (in)
      ARRAY1,   & ! (in)
      ARRAY2,   & ! (in)
      UNDEF     ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA1, IA2
    real(RP), intent(in) :: ARRAY1(IA1)
    real(RP), intent(in) :: ARRAY2(IA2)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real :: STATISTICS_regression
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      if( abs( STATISTICS_correlation( IA1, IA2, ARRAY1(:), ARRAY2(:), UNDEF ) - UNDEF ) > EPS .and. &
          abs( STATISTICS_stddev( IA1, ARRAY1(:), UNDEF ) - UNDEF ) > EPS                      .and. &
          abs( STATISTICS_stddev( IA2, ARRAY2(:), UNDEF ) - UNDEF ) > EPS                      .and. &
          abs( STATISTICS_stddev( IA2, ARRAY2(:), UNDEF )         ) > 0.0_RP                         ) then
        ! regression coefficient without undefined value
        STATISTICS_regression = STATISTICS_correlation( IA1, IA2, ARRAY1(:), ARRAY2(:), UNDEF ) &
                              * STATISTICS_stddev( IA1, ARRAY1(:), UNDEF ) &
                              / STATISTICS_stddev( IA2, ARRAY2(:), UNDEF )
      else
        STATISTICS_regression = UNDEF
      end if
    else
      ! simple regression coefficient
      STATISTICS_regression = STATISTICS_correlation( IA1, IA2, ARRAY1(:), ARRAY2(:) ) &
                            * STATISTICS_stddev( IA1, ARRAY1(:) ) &
                            / STATISTICS_stddev( IA2, ARRAY2(:) )
    end if

    return
  end function STATISTICS_regression

  !-----------------------------------------------------------------------------
  !
  ! Lag correlation coefficient:
  !   correlation between ARRAY1(t) and ARRAY2(t+LAG)
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_lag_correlation( &
      IA1, IA2, & ! (in)
      ARRAY1,   & ! (in)
      ARRAY2,   & ! (in)
      LAG,      & ! (in)
      UNDEF     ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA1, IA2
    real(RP), intent(in) :: ARRAY1(IA1)
    real(RP), intent(in) :: ARRAY2(IA2)

    integer,  intent(in) :: LAG

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_lag_correlation

    ! works
    real(RP) :: tmp(IA2)
    !---------------------------------------------------------------------------

    tmp(:) = eoshift( ARRAY2(:), LAG, UNDEF )

    STATISTICS_lag_correlation = STATISTICS_correlation( IA1, IA2, ARRAY1(:), tmp(:), UNDEF )

    return
  end function STATISTICS_lag_correlation

  !-----------------------------------------------------------------------------
  !
  ! Partial correlation coefficient:
  !   correlation between ARRAY1 and ARRAY2 excepting influence of ARRAY3
  !
  !-----------------------------------------------------------------------------
  function STATISTICS_partial_correlation( &
      IA1, IA2, IA3, & ! (in)
      ARRAY1,        & ! (in)
      ARRAY2,        & ! (in)
      ARRAY3,        & ! (in)
      UNDEF          ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA1, IA2, IA3
    real(RP), intent(in) :: ARRAY1(IA1)
    real(RP), intent(in) :: ARRAY2(IA2)
    real(RP), intent(in) :: ARRAY3(IA3)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_partial_correlation

    ! works
    real(RP) :: corr12
    real(RP) :: corr13
    real(RP) :: corr23
    !---------------------------------------------------------------------------

    if( present( UNDEF ) ) then
      corr12 = STATISTICS_correlation( IA1, IA2, ARRAY1(:), ARRAY2(:), UNDEF )
      corr13 = STATISTICS_correlation( IA1, IA3, ARRAY1(:), ARRAY3(:), UNDEF )
      corr23 = STATISTICS_correlation( IA2, IA3, ARRAY2(:), ARRAY3(:), UNDEF )

      if( corr12 == UNDEF .or. &
          corr13 == UNDEF .or. &
          corr23 == UNDEF      ) then
        ! return undefined value if correlation can't calculate
        STATISTICS_partial_correlation = UNDEF
      else if( corr12 >= 1.0_RP .or. &
               corr13 >= 1.0_RP .or. &
               corr23 >= 1.0_RP      ) then
        ! return undefined value if correlation is 1
        STATISTICS_partial_correlation = UNDEF
      else
        ! partial correlation without undefined value
        STATISTICS_partial_correlation = ( corr12 - ( corr13 * corr23 ) ) &
                                       / sqrt( ( 1.0_RP - corr13**2 ) * ( 1.0_RP - corr23**2 ) )
      endif
    else
      corr12 = STATISTICS_correlation( IA1, IA2, ARRAY1(:), ARRAY2(:) )
      corr13 = STATISTICS_correlation( IA1, IA3, ARRAY1(:), ARRAY3(:) )
      corr23 = STATISTICS_correlation( IA2, IA3, ARRAY2(:), ARRAY3(:) )

      ! simple partial correlation
      STATISTICS_partial_correlation = ( corr12 - ( corr13 * corr23 ) ) &
                                     / sqrt( ( 1.0_RP - corr13**2 ) * ( 1.0_RP - corr23**2 ) )
    end if

    return
  end function STATISTICS_partial_correlation

  !-----------------------------------------------------------------------------
  function STATISTICS_undef_replace( &
      IA,     & ! (in)
      ARRAY,  & ! (in)
      UNDEF1, & ! (in)
      UNDEF2  ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA
    real(RP), intent(in) :: ARRAY(IA)
    real(RP), intent(in) :: UNDEF1
    real(RP), intent(in) :: UNDEF2

    ! function result
    real(RP) :: STATISTICS_undef_replace(IA)

    ! works
    integer :: i
    !---------------------------------------------------------------------------

    do i = 1, IA
      if( abs( ARRAY(i) - UNDEF2 ) > EPS ) then
        STATISTICS_undef_replace(i) = ARRAY(i)
      else
        STATISTICS_undef_replace(i) = UNDEF1
      end if
    end do

    return
  end function STATISTICS_undef_replace

  !-----------------------------------------------------------------------------
  function STATISTICS_undef_embed( &
      IA1, IA2, & ! (in)
      ARRAY1,   & ! (in)
      ARRAY2,   & ! (in)
      UNDEF     ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA1, IA2
    real(RP), intent(in) :: ARRAY1(IA1)
    real(RP), intent(in) :: ARRAY2(IA2)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    real(RP) :: STATISTICS_undef_embed(IA1)

    ! works
    integer :: i
    !---------------------------------------------------------------------------

    ! check the length of input arrays
    if( IA1 /= IA2 ) then
      write(*,*) "Error: different array size !!"
      stop
    end if

    do i = 1, IA1
      if( abs( ARRAY2(i) - UNDEF ) > EPS ) then
        STATISTICS_undef_embed(i) = ARRAY1(i)
      else
        STATISTICS_undef_embed(i) = UNDEF
      end if
    end do
  end function STATISTICS_undef_embed

  !-----------------------------------------------------------------------------
  function STATISTICS_undef_arraysize( &
      IA1, IA2, & ! (in)
      ARRAY1,   & ! (in)
      ARRAY2,   & ! (in)
      UNDEF     ) ! (in)
    implicit none

    ! arguments
    integer,  intent(in) :: IA1, IA2
    real(RP), intent(in) :: ARRAY1(IA1)
    real(RP), intent(in) :: ARRAY2(IA2)

    real(RP), intent(in), optional :: UNDEF

    ! function result
    integer :: STATISTICS_undef_arraysize

    ! works
    integer :: i, cnt
    !---------------------------------------------------------------------------

    ! check the length of input arrays
    if( IA1 /= IA2 ) then
      write(*,*) "Error: different array size !!"
      stop
    end if

    if( present( UNDEF ) ) then
      cnt = 0
      do i = 1, IA1
        if( abs( ARRAY1(i) - UNDEF ) > EPS .and. &
            abs( ARRAY2(i) - UNDEF ) > EPS       ) then
          cnt = cnt + 1
        end if
      end do
      ! array size excepting number of undefined values
      STATISTICS_undef_arraysize = cnt
    else
      ! simple array size
      STATISTICS_undef_arraysize = IA1
    end if

    return
  end function STATISTICS_undef_arraysize

end module scale_statistics
