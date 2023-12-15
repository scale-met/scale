!-------------------------------------------------------------------------------
!> module Parameter Estimation
!!
!! @par Description
!!          Parameter Estimation
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_da_param_estimation
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: DA_param_estimation_setup
  public :: DA_param_estimation_finalize
  public :: DA_param_estimation_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: estimation_getvar
  private :: estimation_putvar
  private :: estimation_log
  private :: phys2func
  private :: func2phys

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: MAXNUM_PARAM_ESTIMATION = 100

  integer :: NUM_PARAM_ESTIMATION
  logical :: DO_PARAM_ESTIMATION

  character(len=H_SHORT) :: ESTIMATION_TARGET(MAXNUM_PARAM_ESTIMATION)

  logical :: PEST_TRANSLATION(MAXNUM_PARAM_ESTIMATION)

  real(RP) :: PEST_ULIMIT(MAXNUM_PARAM_ESTIMATION)
  real(RP) :: PEST_LLIMIT(MAXNUM_PARAM_ESTIMATION)

  real(RP), allocatable :: param_array(:,:)

  integer :: datatype

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine DA_param_estimation_setup
    use mpi
    use scale_prc, only: &
       PRC_abort
    use scale_comm_ensemble, only: &
       ENSEMBLE_nprocs => COMM_ENSEMBLE_nprocs, &
       ENSEMBLE_myrank => COMM_ENSEMBLE_myrank
    implicit none

    namelist / PARAM_DA_PARAM_ESTIMATION / &
       ESTIMATION_TARGET, &
       PEST_TRANSLATION,  &
       PEST_ULIMIT,       &
       PEST_LLIMIT

    integer :: v
    integer :: ierr
    !---------------------------------------------------------------------------

    NUM_PARAM_ESTIMATION = 0
    DO_PARAM_ESTIMATION = .false.

    ESTIMATION_TARGET(:) = ''

    PEST_TRANSLATION(:) = .false.

    PEST_ULIMIT(:) = huge(0.0_RP)
    PEST_LLIMIT(:) = -huge(0.0_RP)

    !--- gather parameters from ensemble communication
    if     ( RP == SP ) then
       datatype = MPI_REAL
    else if( RP == DP ) then
       datatype = MPI_DOUBLE_PRECISION
    else
       LOG_ERROR("DA_param_estimation_setup",*) 'The precision has not been implemented yet:', RP
       call PRC_abort
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_DA_PARAM_ESTIMATION,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("DA_param_estimation_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("DA_param_estimation_setup",*) 'Not appropriate names in namelist PARAM_DA_PARAM_ESTIMATION. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_DA_PARAM_ESTIMATION)

    do v = 1, MAXNUM_PARAM_ESTIMATION
      if( trim(ESTIMATION_TARGET(v)) /= '' ) then
        NUM_PARAM_ESTIMATION = NUM_PARAM_ESTIMATION + 1
      end if
    end do

    if( NUM_PARAM_ESTIMATION > MAXNUM_PARAM_ESTIMATION ) then
       LOG_ERROR("DA_param_estimation_setup",*) 'NUM_PARAM_ESTIMATION has exceeded the limit! ( MAX = ', MAXNUM_PARAM_ESTIMATION, ')'
       call PRC_abort
    end if

    if( NUM_PARAM_ESTIMATION > 0 ) then
      LOG_INFO("DA_param_estimation_setup",*) 'Set up Parameter Estimation module. The number of estimation target:', NUM_PARAM_ESTIMATION
      DO_PARAM_ESTIMATION = .true.

      allocate( param_array( ENSEMBLE_nprocs, NUM_PARAM_ESTIMATION ) )

      call estimation_getvar
      call estimation_log

    end if

    return
  end subroutine DA_param_estimation_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine DA_param_estimation_finalize
    implicit none
    !---------------------------------------------------------------------------

    if( allocated( param_array ) ) deallocate( param_array )

    return
  end subroutine DA_param_estimation_finalize

  !-----------------------------------------------------------------------------
  !> Data Assimilation
  subroutine DA_param_estimation_update
    use scale_comm_ensemble, only: &
      ENSEMBLE_nprocs => COMM_ENSEMBLE_nprocs
    use scale_letkf, only: &
      LETKF_param_estimation_system
    implicit none

    integer :: n, v
    !---------------------------------------------------------------------------

    if( .NOT. DO_PARAM_ESTIMATION ) return

    !--- get parameters for estimation
    call estimation_getvar

    do v = 1, NUM_PARAM_ESTIMATION
    do n = 1, ENSEMBLE_nprocs
       if( PEST_TRANSLATION(v) ) then ! encoding function by Kotsuki (2015)
          param_array(n,v) = phys2func( v, param_array(n,v) )
       end if
    end do
    end do

    !--- execute parameter estimation
    call LETKF_param_estimation_system( NUM_PARAM_ESTIMATION, &
                                        param_array(:,:)      )

    do v = 1, NUM_PARAM_ESTIMATION
    do n = 1, ENSEMBLE_nprocs
       if( PEST_TRANSLATION(v) ) then ! decoding function by Kotsuki (2015)
          param_array(n,v) = func2phys( v, param_array(n,v) )
       end if
    end do
    end do

    !--- overwrite parameters
    call estimation_putvar
    call estimation_log

    return
  end subroutine DA_param_estimation_update

  subroutine estimation_getvar
    use mpi
    use scale_prc, only: &
      PRC_abort
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_comm_ensemble, only: &
      ENSEMBLE_world  => COMM_ENSEMBLE_world,  &
      ENSEMBLE_nprocs => COMM_ENSEMBLE_nprocs, &
      ENSEMBLE_myrank => COMM_ENSEMBLE_myrank
    use scale_atmos_phy_mp_tomita08, only: &
      ATMOS_PHY_MP_tomita08_getvar
    use scale_atmos_phy_cp_kf, only: &
      ATMOS_PHY_CP_kf_getvar
    implicit none

    real(RP), allocatable :: send(:)
    real(RP), allocatable :: recv(:)

    integer :: n, v
    integer :: ierr
    !---------------------------------------------------------------------------

    do v = 1, NUM_PARAM_ESTIMATION
      select case( trim(ESTIMATION_TARGET(v)) )
      case('MP_tomita08_drag_g')
         call ATMOS_PHY_MP_tomita08_getvar( out_drag_g=param_array(1+ENSEMBLE_myrank,v) )
      case('MP_tomita08_re_qc')
         call ATMOS_PHY_MP_tomita08_getvar( out_re_qc=param_array(1+ENSEMBLE_myrank,v) )
      case('MP_tomita08_re_qi')
         call ATMOS_PHY_MP_tomita08_getvar( out_re_qi=param_array(1+ENSEMBLE_myrank,v) )
      case('MP_tomita08_Cr')
         call ATMOS_PHY_MP_tomita08_getvar( out_Cr=param_array(1+ENSEMBLE_myrank,v) )
      case('MP_tomita08_Cs')
         call ATMOS_PHY_MP_tomita08_getvar( out_Cs=param_array(1+ENSEMBLE_myrank,v) )
      case('CP_kf_delcape')
         call ATMOS_PHY_CP_kf_getvar( out_delcape=param_array(1+ENSEMBLE_myrank,v) )
      case('CP_kf_deeplifetime')
         call ATMOS_PHY_CP_kf_getvar( out_deeplifetime=param_array(1+ENSEMBLE_myrank,v) )
      case('CP_kf_shallowlifetime')
         call ATMOS_PHY_CP_kf_getvar( out_shallowlifetime=param_array(1+ENSEMBLE_myrank,v) )
      case default
         LOG_ERROR("DA_param_estimation_update",*) 'Invalid ESTIMATION_TARGET:', trim(ESTIMATION_TARGET(v))
         call PRC_abort
      end select
    end do

    allocate( send(NUM_PARAM_ESTIMATION*ENSEMBLE_nprocs) )
    allocate( recv(NUM_PARAM_ESTIMATION*ENSEMBLE_nprocs) )

    send(:) = UNDEF
    recv(:) = UNDEF

    do v = 1, NUM_PARAM_ESTIMATION
       send( v + ENSEMBLE_myrank*NUM_PARAM_ESTIMATION ) = param_array(1+ENSEMBLE_myrank,v)
    end do

    n = ENSEMBLE_myrank * NUM_PARAM_ESTIMATION + 1

    call MPI_ALLGATHER( send(n),              &
                        NUM_PARAM_ESTIMATION, &
                        datatype,             &
                        recv(1),              &
                        NUM_PARAM_ESTIMATION, &
                        datatype,             &
                        ENSEMBLE_world,       &
                        ierr                  )

    do v = 1, NUM_PARAM_ESTIMATION
    do n = 1, ENSEMBLE_nprocs
       param_array(n,v) = recv( v + (n-1)*NUM_PARAM_ESTIMATION )
    end do
    end do

    deallocate( send )
    deallocate( recv )

    return
  end subroutine estimation_getvar

  subroutine estimation_putvar
    use scale_prc, only: &
      PRC_abort
    use scale_comm_ensemble, only: &
      ENSEMBLE_myrank => COMM_ENSEMBLE_myrank
    use scale_atmos_phy_mp_tomita08, only: &
      ATMOS_PHY_MP_tomita08_putvar
    use scale_atmos_phy_cp_kf, only: &
      ATMOS_PHY_CP_kf_putvar
    implicit none

    integer :: v
    !---------------------------------------------------------------------------

    do v = 1, NUM_PARAM_ESTIMATION
      select case( trim(ESTIMATION_TARGET(v)) )
      case('MP_tomita08_drag_g')
         call ATMOS_PHY_MP_tomita08_putvar( in_drag_g=param_array(1+ENSEMBLE_myrank,v) )
      case('MP_tomita08_re_qc')
         call ATMOS_PHY_MP_tomita08_putvar( in_re_qc=param_array(1+ENSEMBLE_myrank,v) )
      case('MP_tomita08_re_qi')
         call ATMOS_PHY_MP_tomita08_putvar( in_re_qi=param_array(1+ENSEMBLE_myrank,v) )
      case('MP_tomita08_Cr')
         call ATMOS_PHY_MP_tomita08_putvar( in_Cr=param_array(1+ENSEMBLE_myrank,v) )
      case('MP_tomita08_Cs')
         call ATMOS_PHY_MP_tomita08_putvar( in_Cs=param_array(1+ENSEMBLE_myrank,v) )
      case('CP_kf_delcape')
         call ATMOS_PHY_CP_kf_putvar( in_delcape=param_array(1+ENSEMBLE_myrank,v) )
      case('CP_kf_deeplifetime')
         call ATMOS_PHY_CP_kf_putvar( in_deeplifetime=param_array(1+ENSEMBLE_myrank,v) )
      case('CP_kf_shallowlifetime')
         call ATMOS_PHY_CP_kf_putvar( in_shallowlifetime=param_array(1+ENSEMBLE_myrank,v) )
      case default
         LOG_ERROR("DA_param_estimation_update",*) 'Invalid ESTIMATION_TARGET:', trim(ESTIMATION_TARGET(v))
         call PRC_abort
      end select
    end do

    return
  end subroutine estimation_putvar

  subroutine estimation_log
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_comm_ensemble, only: &
      ENSEMBLE_myrank => COMM_ENSEMBLE_myrank
    use scale_time, only: &
      TIME_NOWDATE
    use scale_statistics, only: &
      AVERAGE => STATISTICS_AVERAGE, &
      STDDEV  => STATISTICS_STDDEV
    implicit none

    character(len=20) :: nowdate

    integer :: i, v
    !---------------------------------------------------------------------------

    do v = 1, NUM_PARAM_ESTIMATION
      nowdate = '****/**/**-**:**:** '
      write( nowdate( 1: 4), '(i4)' ) TIME_NOWDATE(1)
      write( nowdate( 6: 7), '(i2)' ) TIME_NOWDATE(2)
      write( nowdate( 9:10), '(i2)' ) TIME_NOWDATE(3)
      write( nowdate(12:13), '(i2)' ) TIME_NOWDATE(4)
      write( nowdate(15:16), '(i2)' ) TIME_NOWDATE(5)
      write( nowdate(18:19), '(i2)' ) TIME_NOWDATE(6)
      do i = 1, 19
        if( nowdate(i:i) == ' ' ) nowdate(i:i) = '0'
      end do

      LOG_INFO("DA_param_estimation_update",*) 'result [time, name, rank-var, mean, std]: ', &
        nowdate, trim(ESTIMATION_TARGET(v)), param_array(1+ENSEMBLE_myrank,v), AVERAGE( param_array(:,v), UNDEF ), STDDEV( param_array(:,v), UNDEF )
    end do

    return
  end subroutine estimation_log

  !---------------------------------------------------------------------
  !  Based on the NICAM-LETKF parameter estimation code
  !  created by by Shunji Kotsuki, RIKEN AICS (Dec 2015)
  !---------------------------------------------------------------------
  function phys2func(pidx,phys) ! phys --> func
    implicit none
  
    real(RP) :: phys2func

    integer,  intent(in) :: pidx
    real(RP), intent(in) :: phys

    real(RP) :: pmax, pmin, pmean
    real(RP) :: tmp_phys
  
    pmax = PEST_ULIMIT(pidx)
    pmin = PEST_LLIMIT(pidx)
  
    pmean = 0.5_RP * (pmax + pmin)
    tmp_phys = (phys - pmean) / (pmax - pmean)
  
    phys2func = 0.5_RP * log( (1.0_RP + tmp_phys) / (1.0_RP - tmp_phys) )
    !phys2func = tanh((phys - pmean) / (pmax - pmean))

    if (phys <= pmin .or. pmax <= phys) then
      LOG_INFO('DA_param_estimation (phys2func)',*) 'Waring: Estimated parameter value is outside of the limit! (pmin/pmax/phys)', pmin, pmax, phys
    endif
  
    return
  end function phys2func

  function func2phys(pidx,func) ! func --> phys
    implicit none
  
    real(RP) :: func2phys

    integer,  intent(in) :: pidx
    real(RP), intent(in) :: func

    real(RP) :: pmax, pmin, pmean
  
    pmax = PEST_ULIMIT(pidx)
    pmin = PEST_LLIMIT(pidx)
  
    pmean = 0.5_RP * (pmax + pmin)
  
    func2phys = tanh(func) * (pmax - pmean) + pmean
    !func2phys = 0.5_RP * log( (1.0_RP + func) / (1.0_RP - func) ) * (pmax - pmean) + pmean
  
    if (func2phys <= pmin .or. pmax <= func2phys) then
      LOG_INFO('DA_param_estimation (func2phys)',*) 'Waring: Estimated parameter value is outside of the limit! (pmin/pmax/phys)', pmin, pmax, func2phys
    endif

    return
  end function func2phys

end module mod_da_param_estimation
