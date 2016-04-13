!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          calc perturbation
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-09-05 (H.Yashiro)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use scale_grid, only : &
       CZ => GRID_CZ, CX => GRID_CX, FX => GRID_FX, &
       CDX => GRID_CDX, &
       RCDX => GRID_RCDX, RCDY => GRID_RCDY, RCDZ => GRID_RCDZ

  use scale_atmos_numeric_fdm_def, only: &
       & FLXEVALTYPE_CD2, FLXEVALTYPE_UD1,  &
       & FLXEVALTYPE_CD4, FLXEVALTYPE_UD3,  &
       & FLXEVALTYPE_CD6, FLXEVALTYPE_UD5,  &
       & VL_ZXY

  use scale_atmos_numeric_fdm, only: &
      & ATMOS_NUMERIC_FDM_setup,           &
      & ATMOS_NUMERIC_FDM_EvalFlux,        &
      & ATMOS_NUMERIC_FDM_EvolveVar
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private, save :: USER_do = .false. !< do user step?

  integer, parameter :: I_Qadv = 1
  real(RP), allocatable :: QadvIni(:,:,:)
  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    use mod_atmos_dyn_vars, only: &
         & PROG
    use mod_atmos_vars, only: &
         & DENS, QTRC

    implicit none
    !---------------------------------------------------------------------------

    write(*,*) "Set value of PROG(:,:,:,I_Qadv).."
    allocate(QadvIni(KA,IA,JA))
    
    PROG(:,:,:,I_Qadv) = DENS(:,:,:) * QTRC(:,:,:,I_NC)
    QadvIni = QTRC(:,:,:,I_NC)

!    call convCheck()
    
    ! calculate diagnostic value and input to history buffer
    call USER_step

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_grid, only : &
         CZ => GRID_CZ, CX => GRID_CX, &
         CDX => GRID_CDX
    
    
    use scale_time, only: &
         DTSEC => TIME_DTSEC, &
         NOWTSEC => TIME_NOWSEC
    
    use mod_atmos_vars, only: &
       DENS, &
       RHOT, QTRC
    
    use mod_atmos_dyn_vars, only: &
         & PROG
    use scale_history, only: &
       HIST_in
    implicit none

    real(RP) :: PT_diff(KA,IA,JA), ExactSol(KA,IA,JA)
    real(RP), parameter :: &
         WaveNumCOS = 2d0,     &         
         RxCOSBELL = 3.0e3_RP, &
         RxRECT    = 1.5e3_RP, &
         XIni = 10.0e3_RP,     &
         Lx   = 20.0e3_RP,     &
         PI = acos(-1.0_RP)
    
    real(RP) :: r, xshift, x_, intwork1, intwork2, dist2
    real(RP) :: l2error, linf
    integer :: k, i, j

    
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       xshift = 40.0_RP * (NOWTSEC - DTSEC)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE

          x_ = CX(i) - xshift
          if(x_ < 0.0_RP) then
             x_ = Lx - mod(abs(x_),Lx)
!             if(k==KS.and.j==JS) write(*,*) "shift:", CX(i), x_
          else if(x_ > Lx) then
             x_ = mod(x_,Lx)
          end if

          !* Cos 
          ! ExactSol(k,i,j) = cos( WaveNumCOS * 2.0_RP * PI / Lx *  x_ )

          !* CosBell
          !dist2 = (x_ - XIni)**2/RxCOSBELL**2
          !ExactSol(k,i,j) = cos( 0.5_RP * PI * sqrt(min(dist2, 1.0_RP)) )**2

          !* Rect
          dist2 = (x_ - XIni)**2/RxRECT**2
          if(dist2 <= 1d0) then
             ExactSol(k,i,j) = 1d0
          else
             ExactSol(k,i,j) = 0d0
          end if
       enddo
       enddo
       enddo
       !
!       l2error = sum( sqrt(  (PROG(KS,IS:IE,JS,I_Qadv)/DENS(KS,IS:IE,JS) - ExactSol(KS,IS:IE,JS))**2 ) * CDX(IS:IE) )  !/ sum( ExactSol(KS,IS:IE,JS)**2 * CDX(IS:IE) ) )
       l2error = sqrt( sum( (PROG(KS,IS:IE,JS,I_Qadv)/DENS(KS,IS:IE,JS) - ExactSol(KS,IS:IE,JS))**2 * CDX(IS:IE) )  / sum( ExactSol(KS,IS:IE,JS)**2 * CDX(IS:IE) ) )
       linf = maxval( abs( PROG(KS,IS:IE,JS,I_Qadv)/DENS(KS,IS:IE,JS) - ExactSol(KS,IS:IE,JS) ) ) / maxval( abs(ExactSol(KS,IS:IE,JS)) )


!!$       intwork1 = num_int( &
!!$            & (PROG(KS,IS:IE,JS,I_Qadv)/DENS(KS,IS:IE,JS) - ExactSol(KS,IS:IE,JS))**2, &
!!$            & CX(IS), CX(IE) )
!!$       intwork2 = num_int( &
!!$            & ExactSol(KS,IS:IE,JS)**2, &
!!$            & CX(IS), CX(IE) )
!!$       l2error = sqrt( intwork1 / intwork2 )
!!$       
!       l2error = sqrt( sum( (QTRC(KS,IS:IE,JS,I_NC) - ExactSol(KS,IS:IE,JS))**2 * CDX(IS:IE) ) ) / sum( CDX(IS:IE) ) 
!       PT_diff = PROG(:,:,:,I_Qadv)/DENS - ExactSol
!       PT_diff = QTRC(:,:,:,I_NC) - ExactSol
!       if(mod(NOWTSEC,1.0_RP)==0.0_RP)then
       write(*,*) "t=", NOWTSEC, "l2error", l2error
!       end if
       call HIST_in( PROG(:,:,:,1) / DENS, 'Qadv', 'mass concentration of tracer in advection test', '1' )
       PT_diff = l2error       
       call HIST_in( PT_diff, 'l2error', 'l2error', '1' )
       PT_diff = linf
       call HIST_in( PT_diff, 'linf', 'linf', '1' )

       call HIST_in( (PROG(:,:,:,I_Qadv)/DENS - ExactSol)**2, 'NC_rk', 'NC_rk', '1' )

       !       call HIST_in( ExactSol, 'l2error', 'l2error', '1' )
    endif

    return
  end subroutine USER_step

  subroutine ConvCheck()

    use mod_atmos_dyn_vars, only: &
         & PROG
    use mod_atmos_vars, only: &
         & DENS, QTRC
    
    integer :: IIS, IIE, JJS, JJE, i, j
    real(RP), dimension(KA,IA,JA,3) :: mflx_hi
    real(RP), dimension(KA,IA,JA) :: one, zero, VARTMP, exactRHS, RHS
    real(RP) :: Lx, l2error
    integer :: FlxEvalTypeID
    integer, parameter :: WaveNum = 2
    real(RP), parameter :: PI = acos(-1.0_RP)
    
    Lx = FX(IE) - FX(IS-1)
    mflx_hi = 0.0_RP
    
    call ATMOS_NUMERIC_FDM_setup(1, 1)
    
    do j = 1, JA
    do i = 1, IA
      mflx_hi(:,i,j,XDIR) = 1.0_RP
      VARTMP(:,i,j) = cos( WaveNum * 2.0_RP * PI / Lx *  CX(i) )
      exactRHS(:,i,j) = - WaveNum * 2.0_RP * PI / Lx * sin( WaveNum * 2.0_RP * PI / Lx *  CX(i) )
    enddo
    enddo

    write(*,*) "Lx=", Lx
    write(*,*) FX(1:IA)
    PROG(:,:,:,I_Qadv) = DENS(:,:,:) * VARTMP(:,:,:)
    return
    
    call eval_RHS(VARTMP, mflx_hi, exactRHS, FLXEVALTYPE_UD1, "UD1")
    call eval_RHS(VARTMP, mflx_hi, exactRHS, FLXEVALTYPE_CD2, "CD2")
    call eval_RHS(VARTMP, mflx_hi, exactRHS, FLXEVALTYPE_UD3, "UD3")
    call eval_RHS(VARTMP, mflx_hi, exactRHS, FLXEVALTYPE_CD4, "CD4")
    call eval_RHS(VARTMP, mflx_hi, exactRHS, FLXEVALTYPE_UD5, "UD5")
    call eval_RHS(VARTMP, mflx_hi, exactRHS, FLXEVALTYPE_CD6, "CD6")
    stop
    
  end subroutine ConvCheck
 
  subroutine eval_RHS(var, mflx_hi, exactRHS, FlxEvalTypeID, label)

    real(RP), intent(in) :: var(KA,IA,JA)
    real(RP), intent(in) :: mflx_hi(KA,IA,JA,3), exactRHS(KA,IA,JA)
    integer, intent(in) :: FlxEvalTypeID
    character(*), intent(in) :: label

    integer :: IIS, IIE, JJS, JJE, i, j
    real(RP), dimension(KA,IA,JA,3) :: qflx_hi
    real(RP), dimensioN(KA,IA,JA) :: RHS, one, zero
    real(RP) :: l2error, linf
    
    one = 1.0_RP; zero = 0.0_RP
    
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
       call ATMOS_NUMERIC_FDM_EvalFlux( qflx_hi,                                                &  ! (inout)
        & FlxEvalTypeID, VL_ZXY,                                                               &  ! (in)
        & var, one, mflx_hi(:,:,:,XDIR), mflx_hi(:,:,:,YDIR), mflx_hi(:,:,:,ZDIR), .false., &  ! (in)
        & IIS, IIE, JJS, JJE, KS, KE  )                          ! (in)      
    enddo
    enddo

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
      do j = JJS, JJE
      do i = IIS, IIE 
         RHS(:,i,j) = (qflx_hi(:,i,j,XDIR) - qflx_hi(:,i-1,j,XDIR)) * RCDX(i)
      enddo
      enddo
    enddo
    enddo
    
    l2error = sqrt(   sum( (RHS(KS,IS:IE,JS) - exactRHS(KS,IS:IE,JS))**2 * CDX(IS:IE) )     &
         &          / sum( exactRHS(KS,IS:IE,JS)**2 * CDX(IS:IE) )                                 )
    
    linf = maxval( abs( RHS(KS,IS:IE,JS) - exactRHS(KS,IS:IE,JS) ) ) / maxval( abs(exactRHS(KS,IS:IE,JS) ) )

!!$    write(*,*) RHS(KS,IS:IE,JS)
!!$    write(*,*) exactRHS(KS,IS:IE,JS)    
    write(*,*) "FlxEvalTypeID=", trim(label), ", l2error=", l2error, " linf=", linf
    
  end subroutine eval_RHS
  
  function num_int(fx, xa, xb) result(intval)
    real(RP), intent(in) :: fx(IS:IE)
    real(RP), intent(in) :: xa, xb
    real(RP) :: intval
    
    real(RP) :: dh
    integer :: n, i

    n = size(fx)
    dh = (xb - xa)/dble(n)

    intval = fx(IS) + fx(IE)
    do i=1, n-1
       if(mod(i,2) /= 0) then
          intval = intval + 4.0_RP * fx(IS+i)
       else
          intval = intval + 2.0_RP * fx(IS+i)
       end if
    end do

    intval = dh * intval / 3.0_RP
!    intval = intval + 0.5_RP * dh * (fx(IS) + fx(IE))
  end function num_int
  
end module mod_user
