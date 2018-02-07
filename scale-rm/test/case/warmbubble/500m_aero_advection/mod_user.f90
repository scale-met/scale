!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-20 (S.Nishizawa)   [new] split from dynamical core
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
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
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
  real(RP) :: user_h2so4dt = 0.0_RP
  real(RP) :: user_ocgasdt = 0.0_RP
  real(RP) :: user_emitdt = 0.0_RP
  real(RP) :: user_emitpoint(3)
  real(RP) :: user_m0_sulf = 0.0_RP
  real(RP) :: user_dg_sulf = 80.e-9_RP
  real(RP) :: user_sg_sulf = 1.6_RP
  integer  :: user_n_kap = 1
  integer  :: emit_indx(3)
  data user_emitpoint / 0.0_RP, 0.0_RP, 0.0_RP /
  data emit_indx / 0, 0, 0 /
  real(RP),allocatable :: user_aerosol_procs(:,:,:,:) !(n_atr,n_siz_max,n_kap_max,n_ctg)
  real(RP), parameter :: d_min_def = 1.e-9_RP ! default lower bound of 1st size bin
  real(RP), parameter :: d_max_def = 1.e-5_RP ! upper bound of last size bin
  integer, parameter  :: n_kap_def = 1
  integer :: n_siz_max, n_kap_max
  real(RP),allocatable :: d_lw(:,:), d_up(:,:)  !diameter [m]
  real(RP),allocatable :: d_min(:), d_max(:)    !lower and upper bound of 1st size bin (n_ctg)
  integer, allocatable :: n_kap(:)              !number of kappa bins (n_ctg)
  integer, parameter ::  ia_m0  = 1             !1. number conc        [#/m3]
  integer, parameter ::  ia_m2  = 2             !2. 2nd mom conc       [m2/m3]
  integer, parameter ::  ia_m3  = 3             !3. 3rd mom conc       [m3/m3]
  integer, parameter ::  ia_ms  = 4             !4. mass conc          [ug/m3]
  integer, parameter ::  ia_kp  = 5
  real(RP), parameter  :: rhod_ae    = 1.83_RP              ! particle density [g/cm3] sulfate assumed
  real(RP), parameter  :: conv_vl_ms = rhod_ae/1.e-12_RP    ! M3(volume)[m3/m3] to mass[m3/m3]
  !
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
       PRC_MPIstop, &
       PRC_myrank
    use scale_atmos_grid_cartesC, only: &
       CZ => GRID_CARTESC_CZ, &
       CY => GRID_CARTESC_CY, &
       CX => GRID_CARTESC_CX

    implicit none

    namelist / PARAM_USER /    &
               user_h2so4dt,   &
               user_ocgasdt,   &
               user_emitdt,    &
               user_emitpoint, &
               user_m0_sulf,   &
               user_dg_sulf,   &
               user_sg_sulf

    integer :: k, i, j
    integer :: ic, is0
    integer :: ierr
    real(RP) :: dlogd
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
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    !--- Determine the emission point
    emit_indx(:) = 0

    do i = IS, IE
       if( CX(i) <= user_emitpoint(1) .and. CX(i+1) > user_emitpoint(1) ) then
         emit_indx(1) = i
       endif
    enddo

    do j = JS, JE
       if( CY(j) <= user_emitpoint(2) .and. CY(j+1) > user_emitpoint(2) ) then
         emit_indx(2) = j
       endif
    enddo

    do k = KS, KE
       if( CZ(k) <= user_emitpoint(3) .and. CZ(k+1) > user_emitpoint(3) ) then
         emit_indx(3) = k
       endif
    enddo

    if( emit_indx(1) /= 0 .and. emit_indx(2) /= 0 &
                          .and. emit_indx(3) /= 0 ) then
      if( IO_L ) write(IO_FID_LOG,*) '*** Emission point is .'
      if( IO_L ) write(IO_FID_LOG,*) '(', CX(emit_indx(1)), &
                                          CY(emit_indx(2)), &
                                          CZ(emit_indx(3)), ')'
      if( IO_L ) write(IO_FID_LOG,*) '*** Emssion Rank is .', PRC_myrank
    endif

    !--- set up aerosol for emission
    n_siz_max = 0
    n_kap_max = 0
    do ic = 1, AE_CTG
      n_siz_max = max(n_siz_max, NSIZ(ic))
      n_kap_max = max(n_kap_max, NKAP(ic))
    enddo

    allocate( user_aerosol_procs (N_ATR,n_siz_max,n_kap_max,AE_CTG)  )
    allocate( d_lw(n_siz_max,AE_CTG) )
    allocate( d_up(n_siz_max,AE_CTG) )
    allocate( d_min(AE_CTG) )
    allocate( d_max(AE_CTG) )
    allocate( n_kap(AE_CTG) )

    user_aerosol_procs(:,:,:,:) = 0._RP
    d_lw(:,:) = 0._RP
    d_up(:,:) = 0._RP

    d_min(:) = d_min_def
    d_max(:) = d_max_def
    n_kap(:) = n_kap_def
    do ic = 1, AE_CTG
      dlogd = (log(d_max(ic)) - log(d_min(ic)))/float(NSIZ(ic))
      do is0 = 1, NSIZ(ic)  !size bin
        d_lw(is0,ic) = exp(log(d_min(ic))+dlogd* float(is0-1)      )
        d_up(is0,ic) = exp(log(d_min(ic))+dlogd* float(is0)        )
      enddo !is (1:n_siz(ic))
    enddo !ic (1:AE_CTG)

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
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_step
    use scale_const, only: &
       CONST_PI
    use scale_time, only: &
       TIME_NOWSEC
    use scale_atmos_phy_ae, only: &
       QS_AE, &
       QE_AE
    use mod_atmos_phy_ae_vars, only: &
       AE_EMIT => ATMOS_PHY_AE_EMIT
    use mod_atmos_vars, only: &
       DENS
    implicit none

    real(RP) :: m0t, dgt, sgt, m2t, m3t, mst
    real(RP) :: pi6

    integer  :: ic, ik, is0, ia0
    integer  :: kpnt, ipnt, jpnt
    integer  :: k, iq
    !---------------------------------------------------------------------------

    pi6 = CONST_PI / 6._RP

    kpnt = emit_indx(3)
    ipnt = emit_indx(1)
    jpnt = emit_indx(2)

    !--- Add tendency of box model

    ! update conc_h2so4
    if ( i /= 0 .AND. j /= 0 .AND. k /= 0 .AND.  TIME_NOWSEC <= user_emitdt ) then

       !--- Emission from Aerosol
       m0t = user_m0_sulf                               ! total M0 [#/m3]
       dgt = user_dg_sulf                               ! [m]
       sgt = user_sg_sulf                               ! [-]
       m2t = m0t * dgt**2 * exp( 2.0_RP * log(sgt)**2 ) ! total M2 [m2/m3]
       m3t = m0t * dgt**3 * exp( 4.5_RP * log(sgt)**3 ) ! total M3 [m3/m3]
       mst = m3t * pi6 * conv_vl_ms                     ! total Ms [ug/m3]

       do ic  = 1, AE_CTG       ! category
          do ik  = 1, n_kap(ic) ! kappa bin
          do is0 = 1, NSIZ(ic)  ! size bin

             if (    dgt >= d_lw(is0,ic) .AND. dgt < d_up(is0,ic) ) then

                user_aerosol_procs(ia_m0,is0     ,ik,ic) = user_aerosol_procs(ia_m0,is0     ,ik,ic) + m0t          ! [#/m3]
                user_aerosol_procs(ia_m2,is0     ,ik,ic) = user_aerosol_procs(ia_m2,is0     ,ik,ic) + m2t          ! [m2/m3]
                user_aerosol_procs(ia_m3,is0     ,ik,ic) = user_aerosol_procs(ia_m3,is0     ,ik,ic) + m3t          ! [m3/m3]
                user_aerosol_procs(ia_ms,is0     ,ik,ic) = user_aerosol_procs(ia_ms,is0     ,ik,ic) + mst*1.E-9_RP ! [kg/m3]

             elseif( dgt < d_lw(1,ic) ) then

                user_aerosol_procs(ia_m0,1       ,ik,ic) = user_aerosol_procs(ia_m0,1       ,ik,ic) + m0t          ! [#/m3]
                user_aerosol_procs(ia_m2,1       ,ik,ic) = user_aerosol_procs(ia_m2,1       ,ik,ic) + m2t          ! [m2/m3]
                user_aerosol_procs(ia_m3,1       ,ik,ic) = user_aerosol_procs(ia_m3,1       ,ik,ic) + m3t          ! [m3/m3]
                user_aerosol_procs(ia_ms,1       ,ik,ic) = user_aerosol_procs(ia_ms,1       ,ik,ic) + mst*1.E-9_RP ! [kg/m3]

             elseif( dgt >= d_up(NSIZ(ic),ic) ) then

                user_aerosol_procs(ia_m0,NSIZ(ic),ik,ic) = user_aerosol_procs(ia_m0,NSIZ(ic),ik,ic) + m0t          ! [#/m3]
                user_aerosol_procs(ia_m2,NSIZ(ic),ik,ic) = user_aerosol_procs(ia_m2,NSIZ(ic),ik,ic) + m2t          ! [m2/m3]
                user_aerosol_procs(ia_m3,NSIZ(ic),ik,ic) = user_aerosol_procs(ia_m3,NSIZ(ic),ik,ic) + m3t          ! [m3/m3]
                user_aerosol_procs(ia_ms,NSIZ(ic),ik,ic) = user_aerosol_procs(ia_ms,NSIZ(ic),ik,ic) + mst*1.E-9_RP ! [kg/m3]

             endif

          enddo
          enddo
       enddo

       iq = QS_AE
       do ic  = 1, AE_CTG      ! category
          do ik  = 1, NKAP(ic) ! kappa bin
          do is0 = 1, NSIZ(ic) ! size bin
          do ia0 = 1, N_ATR    ! attributes

             do k = 1, kpnt
                AE_EMIT(k,ipnt,jpnt,iq) = user_aerosol_procs(ia0,is0,ik,ic) / DENS(k,ipnt,jpnt) !#,m2,m3,kg/m3 -> #,m2,m3,kg/kg
             enddo

             iq = iq + 1
          enddo
          enddo
          enddo
       enddo

       !--- Emission from Gas
       do k = 1, kpnt
          AE_EMIT(k,ipnt,jpnt,QE_AE-GAS_CTG+IG_H2SO4) = user_h2so4dt
          AE_EMIT(k,ipnt,jpnt,QE_AE-GAS_CTG+IG_CGAS ) = user_ocgasdt
       enddo

    else

       AE_EMIT(:,:,:,:) = 0.0_RP

    endif

    return
  end subroutine USER_step

end module mod_user
