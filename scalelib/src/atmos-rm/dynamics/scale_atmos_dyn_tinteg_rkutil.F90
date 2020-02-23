!-------------------------------------------------------------------------------
!> module scale_atmos_dyn_tinteg_rkutil
!!
!! @par Description
!!      Utility for some runge-kutta schemes
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_rkutil
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
#if defined DEBUG || defined QUICKDEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public :: RKUtil
    integer, allocatable :: comm_ind(:)
    real(RP), allocatable :: work0(:,:,:,:)
    real(RP), allocatable :: work(:,:,:,:,:)
    real(RP), allocatable :: buf(:,:,:,:)
    integer :: var_num
    integer :: register_num
    integer :: stage_num
  end type

  public :: ATMOS_DYN_Tinteg_RKUtil_setup
  public :: ATMOS_DYN_Tinteg_RKUtil_rkwork_alloc
  public :: ATMOS_DYN_Tinteg_RKUtil_rkwork_dealloc
  public :: ATMOS_DYN_Tinteg_RKUtil_comm
  public :: ATMOS_DYN_Tinteg_RKUtil_comm_wait

  !- 7 stage RK with 6th order --------------

  public :: ATMOS_DYN_Tinteg_RKUtil_7s6o_nextstage
  public :: ATMOS_DYN_Tinteg_RKUtil_7s6o_updateVar
  public :: ATMOS_DYN_Tinteg_RKUtil_7s6o_updateFlux

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !----------------------------------------------------------------------------

  !- coeffecients for 7 stage RK with 6th order by Butcher (1964) --------------

  real(RP), parameter, public :: RKCoef_a_7s6o_Butcher1964(7,7) = reshape( &
  (/ 0.0_RP, 1.0_RP/3.0_RP,        0.0_RP,  1.0_RP/12.0_RP, -1.0_RP/16.0_RP,         0.0_RP,   9.0_RP/44.0_RP,    &
     0.0_RP,        0.0_RP, 2.0_RP/3.0_RP,   1.0_RP/3.0_RP,   9.0_RP/8.0_RP,  9.0_RP/8.0_RP,  -9.0_RP/11.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP, -1.0_RP/12.0_RP, -3.0_RP/16.0_RP,   -3_RP/8.0_RP,  63.0_RP/44.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,  -3.0_RP/8.0_RP, -3.0_RP/4.0_RP,  18.0_RP/11.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,          0.0_RP,  1.0_RP/2.0_RP,           0.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,          0.0_RP,         0.0_RP, -16.0_RP/11.0_RP,    &
     0.0_RP,        0.0_RP,        0.0_RP,          0.0_RP,          0.0_RP,         0.0_RP,          0.0_RP /),  &
  shape(RKCoef_a_7s6o_Butcher1964) )

 real(RP), parameter, public :: RKCoef_b_7s6o_Butcher1964(7) = &
     1.0_RP/120.0_RP * (/ 11.0_RP, 0.0_RP, 81.0_RP, 81.0_RP, -32.0_RP, -32.0_RP, 11.0_RP /)

  
  !- coeffecients for 7 stage RK with 6th order and extended region of stability by Lawson (1967) --------------

  real(RP), parameter, public :: RKCoef_a_7s6o_Lawson1967(7,7) = reshape( &
  (/ 0.0_RP, 3.0_RP/19.0_RP,  9.0_RP/152.0_RP,   94474764.0_RP/318611987.0_RP,       -76607525678.0_RP/925997907411.0_RP,        -113193410749715476.0_RP/1376008387821185625.0_RP,                510341547912673.0_RP/1709758911034368.0_RP,    &
     0.0_RP,         0.0_RP, 27.0_RP/152.0_RP, -310753854.0_RP/318611987.0_RP,             309768324.0_RP/200562683.0_RP,                              68309142.0_RP/42280325.0_RP,                               -3074637.0_RP/21410624.0_RP,    &
     0.0_RP,         0.0_RP,           0.0_RP,  375818328.0_RP/318611987.0_RP,  -57882086555344.0_RP/37088653028409.0_RP, -9901869473098663108168.0_RP/5940196722617929711875.0_RP,          205532548800199165.0_RP/6225256605226855824.0_RP,    &
     0.0_RP,         0.0_RP,           0.0_RP,                         0.0_RP, 643400862141470.0_RP/704684407539771.0_RP,  8947230518934447694268.0_RP/9333225588784524496875.0_RP, 32370527990426718666299.0_RP/90521226376106372167680.0_RP,    &
     0.0_RP,         0.0_RP,           0.0_RP,                         0.0_RP,                                    0.0_RP,          -8377112295767292.0_RP/1089624335851065625.0_RP,          2610287999955961017.0_RP/236243323046620160.0_RP,    &
     0.0_RP,         0.0_RP,           0.0_RP,                         0.0_RP,                                    0.0_RP,                                                   0.0_RP,         -2690946369187951875.0_RP/253991013039290368.0_RP,    &
     0.0_RP,         0.0_RP,           0.0_RP,                         0.0_RP,                                    0.0_RP,                                                   0.0_RP,                                                    0.0_RP /), &
  shape(RKCoef_a_7s6o_Lawson1967) )

  real(RP), parameter, public :: RKCoef_b_7s6o_Lawson1967(7) = &
    (/                     119490041.0_RP/1597112640.0_RP,                                     0.0_RP,  55710603179056.0_RP/168638187800205.0_RP, &
       5739605598843081731.0_RP/28834038834414422400.0_RP, 1477688286853979.0_RP/291957783566400.0_RP, -298030839900625.0_RP/62778200252544.0_RP, &
                                                                                                                      5352656.0_RP/65415735.0_RP /)


  real(RP), parameter, public :: RKCoef_a_7s6o(7,7) =  RKCoef_a_7s6o_Lawson1967                                                                                                          
  real(RP), parameter, public :: RKCoef_b_7s6o(7)   =  RKCoef_b_7s6o_Lawson1967                                                                                                             

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

  subroutine ATMOS_DYN_Tinteg_RKUtil_setup( this, rk_stage_num, rk_register_num, varname_list, &
      comm_id_offset, alloc_rkwork_flag )
    use scale_const, only: &
      UNDEF  => CONST_UNDEF
    use scale_comm_cartesC, only: &
      COMM_vars8_init      
    implicit none
    
    type(RKUtil), intent(inout) :: this
    integer, intent(in) :: rk_stage_num
    integer, intent(in) :: rk_register_num
    character(*), intent(in) :: varname_list(:)
    integer, intent(in) :: comm_id_offset
    logical, intent(in) :: alloc_rkwork_flag

    integer :: reg_id
    integer :: var_id
    character(H_MID) :: rktag
    !------------------------------------------

    this%var_num = size(varname_list)
    this%stage_num = rk_stage_num
    this%register_num = rk_register_num

    if ( alloc_rkwork_flag ) call ATMOS_DYN_Tinteg_RKUtil_rkwork_alloc( this )

    !-
    allocate( this%buf(KA,IA,JA,this%var_num) )
    allocate( this%comm_ind(this%var_num) )
    
    !$omp parallel workshare
    this%buf(:,:,:,:)   = UNDEF
    !$omp end parallel workshare

    do var_id = 1, this%var_num
      write(rktag,'(a,a,i2.2)') trim(varname_list(var_id)), 'RK'
      this%comm_ind(var_id) = comm_id_offset + var_id
      call COMM_vars8_init( trim(rktag), this%buf(:,:,:,var_id), this%comm_ind(var_id) )
    end do

    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_setup

  subroutine ATMOS_DYN_Tinteg_RKUtil_rkwork_alloc( this )
    implicit none

    type(RKUtil), intent(inout) :: this
    !------------------------------------------

    if (this%register_num > 0) then

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
  end subroutine ATMOS_DYN_Tinteg_RKUtil_rkwork_alloc

  subroutine ATMOS_DYN_Tinteg_RKUtil_rkwork_dealloc( this )
    implicit none

    type(RKUtil), intent(inout) :: this
    !------------------------------------------

    if ( allocated(this%work) ) then
      deallocate( this%work0 )
      deallocate( this%work )
    end if    

    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_rkwork_dealloc

  subroutine ATMOS_DYN_Tinteg_RKUtil_comm( this )
    use scale_comm_cartesC, only: COMM_vars8  
    implicit none

    type(RKUtil), intent(inout) :: this

    integer :: var_id
    !------------------------------------------

    do var_id = 1, this%var_num
      call COMM_vars8( this%buf(:,:,:,var_id), this%comm_ind(var_id) )
    end do 

    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_comm

  subroutine ATMOS_DYN_Tinteg_RKUtil_comm_wait( this )
    use scale_comm_cartesC, only: COMM_wait
    implicit none

    type(RKUtil), intent(inout) :: this

    integer :: var_id
    !------------------------------------------

    do var_id = 1, this%var_num
      call COMM_wait( this%buf(:,:,:,var_id), this%comm_ind(var_id), .false. )
    end do 
    
    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_comm_wait

  !*************************************************************

  !* RK7s6o ------------------------------------------------

  subroutine ATMOS_DYN_Tinteg_RKUtil_7s6o_nextstage( rk, nowstage, io, jo, ko, dt )

    implicit none
    type(RKUtil), intent(inout) :: rk
    integer, intent(in) :: nowstage
    integer, intent(in) :: io(rk%var_num)
    integer, intent(in) :: jo(rk%var_num)
    integer, intent(in) :: ko(rk%var_num)   
    real(RP), intent(in) :: dt
 
    integer :: i, j, k, iv
    real(RP) :: var0
 
    real(RP) :: a_(rk%stage_num)
    !--------------------------------------
 
    a_(:) = dt * RKCoef_a_7s6o(nowstage+1,:)
 
    select case(nowstage)
    case(1)
       !$omp parallel do private(iv,k,j,i,var0) collapse(3)
       do iv=1, rk%var_num
       do k=KS-ko(iv), KE
       do j=JS-jo(iv), JE
       do i=IS-io(iv), IE
          var0 = rk%work0(k,i,j,iv)
          rk%work(k,i,j,iv,1) = rk%work(k,i,j,iv,1) - var0

          rk%buf(k,i,j,iv) = var0 + a_(1) * rk%work(k,i,j,iv,1)
       end do
       end do
       end do
       end do 
      case(2)
        !$omp parallel do private(iv,k,j,i,var0) collapse(3)
        do iv=1, rk%var_num
        do k=KS-ko(iv), KE
        do j=JS-jo(iv), JE
        do i=IS-io(iv), IE
           var0 = rk%work0(k,i,j,iv)
           rk%work(k,i,j,iv,2) =  rk%work(k,i,j,iv,2) - var0
     
           rk%buf(k,i,j,iv) = var0 &
              + a_(1) * rk%work(k,i,j,iv,1) + a_(2) * rk%work(k,i,j,iv,2)
        end do
        end do
        end do
        end do       
    case(3)
       !$omp parallel do private(iv,k,j,i,var0) collapse(3)
       do iv=1, rk%var_num
       do k=KS-ko(iv), KE
       do j=JS-jo(iv), JE
       do i=IS-io(iv), IE
          var0 = rk%work0(k,i,j,iv)
          rk%work(k,i,j,iv,3) =  rk%work(k,i,j,iv,3) - var0
    
          rk%buf(k,i,j,iv) = var0 &
             + a_(1) * rk%work(k,i,j,iv,1) + a_(2) * rk%work(k,i,j,iv,2) &
             + a_(3) * rk%work(k,i,j,iv,3) 
       end do
       end do
       end do
       end do
    case(4)
       !$omp parallel do private(iv,k,j,i,var0) collapse(3)
       do iv=1, rk%var_num
       do k=KS-ko(iv), KE
       do j=JS-jo(iv), JE
       do i=IS-io(iv), IE
          var0 = rk%work0(k,i,j,iv)
          rk%work(k,i,j,iv,4) = rk%work(k,i,j,iv,4) - var0
 
          rk%buf(k,i,j,iv) = var0 &
             + a_(1) * rk%work(k,i,j,iv,1) + a_(2) * rk%work(k,i,j,iv,2) &
             + a_(3) * rk%work(k,i,j,iv,3) + a_(4) * rk%work(k,i,j,iv,4) 
       end do
       end do
       end do
       end do
    case(5)
       !$omp parallel do private(iv,k,j,i,var0) collapse(3)
       do iv=1, rk%var_num
       do k=KS-ko(iv), KE
       do j=JS-jo(iv), JE
       do i=IS-io(iv), IE
          var0 = rk%work0(k,i,j,iv)
          rk%work(k,i,j,iv,5) =  rk%work(k,i,j,iv,5) - var0
    
          rk%buf(k,i,j,iv) = var0 &
             + a_(1) * rk%work(k,i,j,iv,1) + a_(2) * rk%work(k,i,j,iv,2) &
             + a_(3) * rk%work(k,i,j,iv,3) + a_(4) * rk%work(k,i,j,iv,4) &
             + a_(5) * rk%work(k,i,j,iv,5) 
       end do
       end do
       end do
       end do
    case(6)
       !$omp parallel do private(iv,k,j,i,var0) collapse(3)
       do iv=1, rk%var_num
       do k=KS-ko(iv), KE
       do j=JS-jo(iv), JE
       do i=IS-io(iv), IE
          var0 = rk%work0(k,i,j,iv)
          rk%work(k,i,j,iv,6) = rk%work(k,i,j,iv,6) - var0
    
          rk%buf(k,i,j,iv) = var0 &
             + a_(1) * rk%work(k,i,j,iv,1) + a_(2) * rk%work(k,i,j,iv,2) &
             + a_(3) * rk%work(k,i,j,iv,3) + a_(4) * rk%work(k,i,j,iv,4) &
             + a_(5) * rk%work(k,i,j,iv,5) + a_(6) * rk%work(k,i,j,iv,6) 
       end do
       end do
       end do
       end do
    end select
 
    return
   end subroutine ATMOS_DYN_Tinteg_RKUtil_7s6o_nextstage
 
  subroutine ATMOS_DYN_Tinteg_RKUtil_7s6o_updateVar( var, rk, io, jo, ko, vs, ve, dt ) 
    implicit none
 
    integer, intent(in) :: vs, ve
    real(RP), intent(inout) :: var(KA,IA,JA,vs:ve)
    type(RKUtil), intent(inout) :: rk
    integer, intent(in) :: io(rk%var_num)
    integer, intent(in) :: jo(rk%var_num)
    integer, intent(in) :: ko(rk%var_num)
    real(RP), intent(in) :: dt
 
    integer :: i, j, k, iv
    real(RP) :: var0
 
    real(RP) :: b_(rk%stage_num)
    real(RP) :: a77
    !--------------------------------------
 
    b_(:) = dt * RKCoef_b_7s6o(:)
 
    !$omp parallel do private(iv,k,j,i,var0) collapse(3)
    do iv=vs, ve
    do k=KS-ko(iv), KE
    do j=JS-jo(iv), JE
    do i=IS-io(iv), IE
      var0 = rk%work0(k,i,j,iv)
      var(k,i,j,iv) = var0  &
               + b_(1) *  rk%work(k,i,j,iv,1) + b_(3) *  rk%work(k,i,j,iv,3)        &
               + b_(4) *  rk%work(k,i,j,iv,4) + b_(5) *  rk%work(k,i,j,iv,5)        &
               + b_(6) *  rk%work(k,i,j,iv,6) + b_(7) * (rk%work(k,i,j,iv,7) - var0 ) 
    end do
    end do
    end do
    end do
 
    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_7s6o_updateVar
 
  subroutine ATMOS_DYN_Tinteg_RKUtil_7s6o_updateFlux( flux, rkwork, io, jo, ko, va_ ) 
    implicit none
 
    integer, intent(in) :: va_
    real(RP), intent(inout) :: flux(KA,IA,JA,va_)
    real(RP), intent(inout) :: rkwork(KA,IA,JA,va_,7)
    integer, intent(in) :: io, jo, ko
 
    real(RP) :: b_(7)
    integer :: i, j, k, iv
    !--------------------------------------
 
    b_(:) = RKCoef_b_7s6o(:)

    !$omp parallel do private(iv,k,j,i) collapse(3)
    do iv=1, va_
    do k=KS-ko, KE
    do j=JS-jo, JE
    do i=IS-io, IE
      flux(k,i,j,iv) = &
         b_(1) * rkwork(k,i,j,iv,1) + b_(3) * rkwork(k,i,j,iv,3)  &
       + b_(4) * rkwork(k,i,j,iv,4) + b_(5) * rkwork(k,i,j,iv,5)  &
       + b_(6) * rkwork(k,i,j,iv,6) + b_(7) * rkwork(k,i,j,iv,7)
    end do
    end do
    end do
    end do
 
    return
  end subroutine ATMOS_DYN_Tinteg_RKUtil_7s6o_updateFlux

end module scale_atmos_dyn_tinteg_rkutil