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
  use scale_gridtrans  
  use scale_index
  use scale_tracer
  use scale_const, only: &
       PI => CONST_PI
  use scale_grid, only: &
       CZ   => GRID_CZ, &
       CX   => GRID_CX, &
       FX   => GRID_FX, &
       CDX  => GRID_CDX, &
       RCDX => GRID_RCDX
  use mpi
  use scale_comm, only: &
       COMM_datatype, &
       COMM_world
  use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_flux_setup, &
       ATMOS_DYN_FVM_fluxX_XYZ
  use scale_time, only: &
       DTSEC => TIME_DTSEC, &
       NOWTSEC => TIME_NOWSEC
  use mod_atmos_vars, only: &
       DENS, &
       RHOT, &
       QTRC
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
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

  real(RP), parameter :: WaveNumCOS = 2.0_RP
  real(RP), parameter :: RxCOSBELL  = 3.0E3_RP
  real(RP), parameter :: RxRECT     = 1.5E3_RP
  real(RP), parameter :: XIni       = 10.0E3_RP

  character(len=H_SHORT), private, save :: InitShape 
  real(RP),               private, save :: Lx
  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only : &
       FXG   => GRID_FXG, &
    implicit none

    namelist / PARAM_USER / &
       USER_do,             &
       InitShape
    

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

    Lx = FXG(IEG) - FXG(ISG-1)

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

    implicit none

    integer :: k, i, j
    
    !---------------------------------------------------------------------------

    if ( NOWTSEC == 0.0_RP ) then
       select case(InitShape)
       case('COS')
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,I_NC) = cos( WaveNumCOS * 2.0_RP * PI / Lx *  CX(i) )
          enddo
          enddo
          enddo
       end select
    end if
    
    ! calculate diagnostic value and input to history buffer    
    call USER_step

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    use scale_history, only: &
       HIST_in

    implicit none

    real(RP) :: PT_diff(KA,IA,JA), ExactSol(KA,IA,JA)
    
    real(RP) :: r, xshift, x_, dist2
    real(RP) :: l2_error
    real(RP) :: linf_error
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

          select case(InitShape)
          case("COS")
             ExactSol(k,i,j) = cos( WaveNumCOS * 2.0_RP * PI / Lx *  x_ )
          case("BUBBLE") ! COSBELL
             dist2 = (x_ - XIni)**2/RxCOSBELL**2
             ExactSol(k,i,j) = cos( 0.5_RP * PI * sqrt(min(dist2, 1.0_RP)) )**2
          case("RECT")
             dist2 = (x_ - XIni)**2/RxRECT**2
             if(dist2 <= 1d0) then
                ExactSol(k,i,j) = 1d0
             else
                ExactSol(k,i,j) = 0d0
             end if
          end select
          
       enddo
       enddo
       enddo

       l2_error = calc_l2error( QTRC(:,:,:,I_NC), ExactSol )
       linf_error = calc_linferror(QTRC(:,:,:,I_NC), ExactSol)
       if ( mod(NOWTSEC, 10.0_RP) == 0 ) then
          write(*,*) "t=", NOWTSEC, "l2=", l2_error, "linf=", linf_error
       end if
       call HIST_in( l2_error, 'l2error', 'l2error', '1' )
       call HIST_in( linf_error, 'linferror', 'linferror', '1' )
       call HIST_in( (QTRC(:,:,:,I_NC) - ExactSol)**2, 'NC_diff', 'NC_diff', '1' )
    endif

    return
  end subroutine USER_step

  subroutine ConvCheck()
    use mod_atmos_dyn_vars, only: &
         PROG
    use mod_atmos_vars, only: &
         DENS, &
         QTRC
    integer, parameter :: WaveNum = 2
    real(RP) :: l2error

    integer :: IIS, IIE, JJS, JJE, i, j
    real(RP) :: mflx_hi(KA,IA,JA,3)
    real(RP) :: VARTMP(KA,IA,JA)
    real(RP) :: exactRHS(KA,IA,JA)
    real(RP) :: RHS(KA,IA,JA)

    
    mflx_hi = 0.0_RP
    
!    call ATMOS_NUMERIC_FDM_setup(1, 1)
    
    do j = 1, JA
    do i = 1, IA
      mflx_hi(:,i,j,XDIR) = 1.0_RP
      VARTMP(:,i,j) = cos( WaveNum * 2.0_RP * PI / Lx *  CX(i) )
      exactRHS(:,i,j) = - WaveNum * 2.0_RP * PI / Lx * sin( WaveNum * 2.0_RP * PI / Lx *  CX(i) )
    enddo
    enddo

    write(*,*) "Lx=", Lx
!!$    write(*,*) FX(1:IA)
!    return

    !*********************************************
    
!!$    call eval_RHS(VARTMP, mflx_hi, exactRHS, "UD1")
!!$    call eval_RHS(VARTMP, mflx_hi, exactRHS, "CD2")
!!$    call eval_RHS(VARTMP, mflx_hi, exactRHS, "UD3")
!!$    call eval_RHS(VARTMP, mflx_hi, exactRHS, "CD4")
!!$    call eval_RHS(VARTMP, mflx_hi, exactRHS, "UD5")
!!$    call eval_RHS(VARTMP, mflx_hi, exactRHS, "CD6")
!!$
!!$    stop

    !*********************************************
    
  end subroutine ConvCheck

  
!!$  subroutine eval_RHS(var, mflx_hi, exactRHS, lblFluxScheme)
!!$
!!$    use scale_grid, only : &
!!$         CDZ => GRID_CDZ
!!$    
!!$    real(RP), intent(in) :: var(KA,IA,JA)
!!$    real(RP), intent(in) :: mflx_hi(KA,IA,JA,3), exactRHS(KA,IA,JA)
!!$    character(*), intent(in) :: lblFluxScheme
!!$
!!$    integer :: IIS, IIE, JJS, JJE, i, j
!!$    real(RP) :: qflx_hi(KA,IA,JA,3)
!!$    real(RP) :: GSQRT(KA,IA,JA,7)
!!$    real(RP) :: num_diff(KA,IA,JA,5,3)
!!$    real(RP) :: RHS(KA,IA,JA)
!!$    real(RP) :: l2_error
!!$    real(RP) :: linf_error
!!$
!!$
!!$    call ATMOS_DYN_FVM_flux_setup(lblFluxScheme)
!!$    
!!$    num_diff = 0.0_RP
!!$    GSQRT = 1.0_RP
!!$    
!!$    do JJS = JS, JE, JBLOCK
!!$    JJE = JJS+JBLOCK-1
!!$    do IIS = IS, IE, IBLOCK
!!$    IIE = IIS+IBLOCK-1
!!$
!!$       call ATMOS_DYN_FVM_fluxX_XYZ( qflx_hi(:,:,:,XDIR), & ! (out)
!!$            mflx_hi(:,:,:,XDIR), var, GSQRT(:,:,:,I_UYZ), & ! (in)
!!$            num_diff(:,:,:,I_RHOT,XDIR), & ! (in)
!!$            CDZ, & ! (in)
!!$            IIS, IIE, JJS, JJE ) ! (in)
!!$       
!!$    enddo
!!$    enddo
!!$
!!$    do JJS = JS, JE, JBLOCK
!!$    JJE = JJS+JBLOCK-1
!!$    do IIS = IS, IE, IBLOCK
!!$    IIE = IIS+IBLOCK-1
!!$      do j = JJS, JJE
!!$      do i = IIS, IIE 
!!$         RHS(:,i,j) = (qflx_hi(:,i,j,XDIR) - qflx_hi(:,i-1,j,XDIR)) * RCDX(i)
!!$      enddo
!!$      enddo
!!$    enddo
!!$    enddo
!!$
!!$    l2_error = calc_l2error(RHS, exactRHS)
!!$    linf_error = calc_linferror(RHS, exactRHS)
!!$    write(*,*) "FluxScheme=", trim(lblFluxScheme), ", l2error=", l2_error, " linf=", linf_error
!!$    
!!$  end subroutine eval_RHS

  !----------------------------------------------------------
  
  function calc_l2error(NumSol, ExactSol) result(l2error)
    
    real(RP), intent(in) :: NumSol(KA,IA,JA)
    real(RP), intent(in) :: ExactSol(KA,IA,JA)
    real(RP) :: l2error

    real(RP) :: l2_error_lctmp1, l2_error_lctmp2
    real(RP) :: l2_error_gltmp1, l2_error_gltmp2
    integer :: ierr
    
    l2_error_lctmp1 = sum( (NumSol(KS,IS:IE,JS) - ExactSol(KS,IS:IE,JS))**2  * CDX(IS:IE) )
    l2_error_lctmp2 = sum( ExactSol(KS,IS:IE,JS)**2 * CDX(IS:IE) )

    call MPI_Allreduce( l2_error_lctmp1, l2_error_gltmp1, 1, &
         COMM_datatype, MPI_SUM, COMM_world, ierr )

    call MPI_Allreduce( l2_error_lctmp2, l2_error_gltmp2, 1, &
         COMM_datatype, MPI_SUM, COMM_world, ierr )

    l2error = sqrt(l2_error_gltmp1 / l2_error_gltmp2)
    
  end function calc_l2error

  function calc_linferror(NumSol, ExactSol) result(linf_error)
    
    real(RP), intent(in) :: NumSol(KA,IA,JA)
    real(RP), intent(in) :: ExactSol(KA,IA,JA)
    real(RP) :: linf_error

    real(RP) :: linf_error_lctmp1, linf_error_lctmp2
    real(RP) :: linf_error_gltmp1, linf_error_gltmp2
    integer :: ierr

    linf_error_lctmp1 = maxval( abs(NumSol(KS,IS:IE,JS) - ExactSol(KS,IS:IE,JS)) )
    linf_error_lctmp2 = maxval( abs(ExactSol(KS,IS:IE,JS)) )

    call MPI_Allreduce( linf_error_lctmp1, linf_error_gltmp1, 1, &
         COMM_datatype, MPI_MAX, COMM_world, ierr )

    call MPI_Allreduce( linf_error_lctmp2, linf_error_gltmp2, 1, &
         COMM_datatype, MPI_MAX, COMM_world, ierr )

    linf_error = linf_error_gltmp1 / linf_error_gltmp2
    
  end function calc_linferror
  
  
!!$  function num_int(fx, xa, xb) result(intval)
!!$    real(RP), intent(in) :: fx(IS:IE)
!!$    real(RP), intent(in) :: xa, xb
!!$    real(RP) :: intval
!!$    
!!$    real(RP) :: dh
!!$    integer :: n, i
!!$
!!$    n = size(fx)
!!$    dh = (xb - xa)/dble(n)
!!$
!!$    intval = fx(IS) + fx(IE)
!!$    do i=1, n-1
!!$       if(mod(i,2) /= 0) then
!!$          intval = intval + 4.0_RP * fx(IS+i)
!!$       else
!!$          intval = intval + 2.0_RP * fx(IS+i)
!!$       end if
!!$    end do
!!$
!!$    intval = dh * intval / 3.0_RP
!!$!    intval = intval + 0.5_RP * dh * (fx(IS) + fx(IE))
!!$  end function num_int
  
end module mod_user
