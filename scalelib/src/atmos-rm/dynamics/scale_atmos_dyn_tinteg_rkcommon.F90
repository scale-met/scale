!-------------------------------------------------------------------------------
!> module scale_atmos_dyn_tinteg_rkcommon
!!
!! @par Description
!!      A modoule providing utility for some runge-kutta schemes
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_rkcommon
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_const, only:    &
     UNDEF  => CONST_UNDEF, &
     EPS => CONST_EPS

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public :: RKInfo
    real(RP), allocatable :: rkcoef_a(:,:)
    real(RP), allocatable :: rkcoef_b(:)
    real(RP), allocatable :: work0(:,:,:,:)
    real(RP), allocatable :: work(:,:,:,:,:)
    real(RP), allocatable :: buf(:,:,:,:)
    integer :: var_num
    integer :: register_num
    integer :: stage_num

    integer, allocatable :: comm_ind(:)
    logical :: flux_flag
  end type


  public :: ATMOS_DYN_Tinteg_RKCommon_setup
  public :: ATMOS_DYN_Tinteg_RKCommon_rkwork_alloc
  public :: ATMOS_DYN_Tinteg_RKCommon_rkwork_dealloc
  public :: ATMOS_DYN_Tinteg_RKCommon_comm
  public :: ATMOS_DYN_Tinteg_RKCommon_comm_wait

  public :: ATMOS_DYN_Tinteg_RKCommon_nextstage
  public :: ATMOS_DYN_Tinteg_RKCommon_updateVar
  public :: ATMOS_DYN_Tinteg_RKCommon_updateFlux

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------

  subroutine ATMOS_DYN_Tinteg_RKCommon_setup( this,                  &
      rk_stage_num, rk_register_num, rkcoef_a, rkcoef_b,             &
      varname_list, is_type_flux, alloc_rkwork_flag, comm_id_offset )

    use scale_comm_cartesC, only: COMM_vars8_init      
    implicit none
    
    type(RKInfo), intent(inout) :: this
    integer, intent(in) :: rk_stage_num
    integer, intent(in) :: rk_register_num
    real(RP), intent(in) :: rkcoef_a(rk_stage_num,rk_stage_num)
    real(RP), intent(in) :: rkcoef_b(rk_stage_num)
    character(*), intent(in) :: varname_list(:)
    logical, intent(in), optional :: is_type_flux
    logical, intent(in), optional :: alloc_rkwork_flag
    integer, intent(in), optional :: comm_id_offset

    integer :: reg_id
    integer :: var_id
    character(H_MID) :: rktag
    !------------------------------------------

    this%var_num = size(varname_list)
    this%stage_num = rk_stage_num
    this%register_num = rk_register_num
    if (present(is_type_flux)) then
      this%flux_flag = is_type_flux      
    else
      this%flux_flag = .false.
    end if


    allocate( this%rkcoef_a(this%stage_num,this%stage_num) )
    allocate( this%rkcoef_b(this%stage_num) )
    this%rkcoef_a(:,:) = rkcoef_a(:,:)
    this%rkcoef_b(:) = rkcoef_b(:)

    ! If necessary, allocate working arrays to register the tendencies at some RK stagaes. 

    if ( present(alloc_rkwork_flag) ) then
      if (alloc_rkwork_flag) call ATMOS_DYN_Tinteg_RKCommon_rkwork_alloc( this )
    end if

    ! Allocate an array to store variables at next RK stage

    allocate( this%buf(KA,IA,JA,this%var_num) )
    
    !$omp parallel workshare
    this%buf(:,:,:,:)   = UNDEF
    !$omp end parallel workshare

    ! Prepair some information for communication with MPI

    if ( .not. this%flux_flag ) then
      allocate( this%comm_ind(this%var_num) )
      do var_id = 1, this%var_num
        this%comm_ind(var_id) = var_id
      end do
      if (present(comm_id_offset)) &
        this%comm_ind(:) = comm_id_offset + this%comm_ind(:) 

      do var_id = 1, this%var_num
        write(rktag,'(a,a)') trim(varname_list(var_id)), 'RK'
        call COMM_vars8_init( trim(rktag), this%buf(:,:,:,var_id), this%comm_ind(var_id) )
      end do
    end if

    return
  end subroutine ATMOS_DYN_Tinteg_RKCommon_setup

  subroutine ATMOS_DYN_Tinteg_RKCommon_rkwork_alloc( this )
    implicit none

    type(RKInfo), intent(inout) :: this
    !------------------------------------------

    if (this%register_num > 0 .and. (.not. this%flux_flag) ) then

      allocate( this%work0(KA,IA,JA,this%var_num) )
      allocate( this%work(KA,IA,JA,this%var_num,this%register_num) )

#ifdef DEBUG      
      !$omp parallel workshare
      this%work0 (:,:,:,:)   = UNDEF
      this%work  (:,:,:,:,:) = UNDEF
      !$omp end parallel workshare
#endif
    end if    

    return
  end subroutine ATMOS_DYN_Tinteg_RKCommon_rkwork_alloc

  subroutine ATMOS_DYN_Tinteg_RKCommon_rkwork_dealloc( this )
    implicit none

    type(RKInfo), intent(inout) :: this
    !------------------------------------------

    if ( allocated(this%work) ) then
      deallocate( this%work0 )
      deallocate( this%work )
    end if    

    return
  end subroutine ATMOS_DYN_Tinteg_RKCommon_rkwork_dealloc

  subroutine ATMOS_DYN_Tinteg_RKCommon_comm( this )
    use scale_comm_cartesC, only: COMM_vars8  
    implicit none

    type(RKInfo), intent(inout) :: this

    integer :: var_id
    !------------------------------------------

    do var_id = 1, this%var_num
      call COMM_vars8( this%buf(:,:,:,var_id), this%comm_ind(var_id) )
    end do 

    return
  end subroutine ATMOS_DYN_Tinteg_RKCommon_comm

  subroutine ATMOS_DYN_Tinteg_RKCommon_comm_wait( this )
    use scale_comm_cartesC, only: COMM_wait
    implicit none

    type(RKInfo), intent(inout) :: this

    integer :: var_id
    !------------------------------------------

    do var_id = 1, this%var_num
      call COMM_wait( this%buf(:,:,:,var_id), this%comm_ind(var_id), .false. )
    end do 
    
    return
  end subroutine ATMOS_DYN_Tinteg_RKCommon_comm_wait

  subroutine ATMOS_DYN_TInteg_RKCommon_nextstage( this, nowstage, io, jo, ko, dt )
    implicit none
    type(RKInfo), intent(inout) :: this
    integer, intent(in) :: nowstage
    integer, intent(in) :: io(this%var_num)
    integer, intent(in) :: jo(this%var_num)
    integer, intent(in) :: ko(this%var_num)   
    real(RP), intent(in) :: dt
 
    integer :: i, j, k, iv, rks
    real(RP) :: var0
 
    real(RP) :: a_(this%stage_num)
    !--------------------------------------

    a_(:) = dt * this%RKCoef_a(nowstage+1,:)

    !$omp parallel private(j,i,k,var0,iv,rks)

    do iv=1, this%var_num
      !$omp do collapse(2)
      do j=JS-jo(iv), JE
      do i=IS-io(iv), IE
      do k=KS-ko(iv), KE
          var0 = this%work0(k,i,j,iv)
          this%work(k,i,j,iv,nowstage) =  this%work(k,i,j,iv,nowstage) - var0
          this%buf(k,i,j,iv) = var0 + a_(nowstage) * this%work(k,i,j,iv,nowstage)
      end do
      end do
      end do
      !$omp end do
    end do

    do rks=1, nowstage-1
      if ( abs(this%RKCoef_a(nowstage+1,rks)) < EPS ) cycle
      do iv=1, this%var_num
        !$omp do collapse(2)
        do j=JS-jo(iv), JE
        do i=IS-io(iv), IE
        do k=KS-ko(iv), KE
            this%buf(k,i,j,iv) = this%buf(k,i,j,iv) + a_(rks) * this%work(k,i,j,iv,rks)
        end do
        end do
        end do
        !$omp end do
      end do        
    end do

    !$omp end parallel

    return
  end subroutine ATMOS_DYN_TInteg_RKCommon_nextstage

  subroutine ATMOS_DYN_Tinteg_RKCommon_updateVar( this, io, jo, ko, vs, ve, dt, var ) 
    implicit none
 
    type(RKInfo), intent(inout) :: this
    integer, intent(in) :: io(this%var_num)
    integer, intent(in) :: jo(this%var_num)
    integer, intent(in) :: ko(this%var_num)
    real(RP), intent(in) :: dt
    integer, intent(in) :: vs, ve
    real(RP), intent(inout) :: var(KA,IA,JA,vs:ve)
 
    integer :: i, j, k, iv, rks
    real(RP) :: var0
 
    real(RP) :: b_(this%stage_num)
    !--------------------------------------
 
    b_(:) = dt * this%rkcoef_b(:)
 
    !$omp parallel private(k,i,j,var0)

    rks = this%stage_num
    do iv=vs, ve
      !$omp do private(k,i,j,var0) collapse(2)
      do j=JS-jo(iv), JE
      do i=IS-io(iv), IE
      do k=KS-ko(iv), KE      
          var0 = this%work0(k,i,j,iv)
          var(k,i,j,iv) = var0  + b_(rks) * ( this%work(k,i,j,iv,rks) - var0 )
      end do
      end do
      end do
    end do
    do rks=1, this%stage_num-1
      do iv=vs, ve
        !$omp do private(k,i,j,var0) collapse(2)
        do j=JS-jo(iv), JE
        do i=IS-io(iv), IE
        do k=KS-ko(iv), KE      
          var(k,i,j,iv) = var(k,i,j,iv) + b_(rks) * this%work(k,i,j,iv,rks)
        end do
        end do
        end do
      end do
    end do
    !$omp end parallel
    
    return
  end subroutine ATMOS_DYN_Tinteg_RKCommon_updateVar

  subroutine ATMOS_DYN_Tinteg_RKCommon_updateFlux( this, nowstage, io, jo, ko, va_, flux ) 
    implicit none

    type(RKInfo), intent(inout) :: this
    integer, intent(in) :: nowstage
    integer, intent(in) :: va_
    integer, intent(in) :: io, jo, ko
    real(RP), intent(inout) :: flux(KA,IA,JA,va_)

    integer :: i, j, k, iv
    !--------------------------------------
 
    if ( nowstage == 1) then
      !$omp parallel do private(iv,k,j,i) collapse(3)
      do iv=1, va_
      do j=JS-jo, JE
      do i=IS-io, IE
      do k=KS-ko, KE
          flux(k,i,j,iv) = this%rkcoef_b(nowstage) * this%buf(k,i,j,iv)
      end do
      end do
      end do
      end do  
    else
      !$omp parallel do private(iv,k,j,i) collapse(3)
      do iv=1, va_
      do j=JS-jo, JE
      do i=IS-io, IE
      do k=KS-ko, KE
        flux(k,i,j,iv) = flux(k,i,j,iv) + this%rkcoef_b(nowstage) * this%buf(k,i,j,iv)
      end do
      end do
      end do
      end do
    end if
 
    return
  end subroutine ATMOS_DYN_Tinteg_RKCommon_updateFlux

end module scale_atmos_dyn_tinteg_rkcommon