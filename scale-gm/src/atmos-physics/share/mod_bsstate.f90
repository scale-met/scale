!-------------------------------------------------------------------------------
!> Module basic state
!!
!! @par Description
!!          This module is for the set of basic state for non-hydrostatic model
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_bsstate
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: bsstate_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: rho_bs   (:,:,:)
  real(RP), public, allocatable :: rho_bs_pl(:,:,:)
  real(RP), public, allocatable :: pre_bs   (:,:,:)
  real(RP), public, allocatable :: pre_bs_pl(:,:,:)
  real(RP), public, allocatable :: tem_bs   (:,:,:)
  real(RP), public, allocatable :: tem_bs_pl(:,:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: bsstate_input_ref
  private :: bsstate_output_ref
  private :: bsstate_generate
  private :: set_basicstate

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: ref_type  = 'NOBASE' !--- Basic state type
                                               ! 'NOBASE' : no basic state
                                               ! 'INPUT'  : input
                                               ! 'TEM'    : temperature is given.
                                               ! 'TH'     : potential temperature is given.

  character(len=H_LONG),  private :: ref_fname = 'ref.dat'

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine bsstate_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       PRE00 => CONST_PRE00
    use mod_adm, only: &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_gall_pl, &
       ADM_gall,    &
       ADM_kmax,    &
       ADM_kmin
    implicit none

    namelist / BSSTATEPARAM / &
       ref_type,  &
       ref_fname

    real(RP) :: pre_ref(ADM_kall) ! reference pressure
    real(RP) :: tem_ref(ADM_kall) ! reference temperature
    real(RP) :: qv_ref (ADM_kall) ! water vapor
    real(RP) :: rho_ref(ADM_kall) ! density
    real(RP) :: th_ref (ADM_kall) ! potentical temperature (dry)

    integer  :: k
    integer  :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[basic state]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=BSSTATEPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** BSSTATEPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist BSSTATEPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=BSSTATEPARAM)

    !--- allocation of reference variables
    allocate( rho_bs   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( rho_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    rho_bs   (:,:,:) = 0.0_RP
    rho_bs_pl(:,:,:) = 0.0_RP

    allocate( pre_bs   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( pre_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    pre_bs   (:,:,:) = 0.0_RP
    pre_bs_pl(:,:,:) = 0.0_RP

    allocate( tem_bs   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( tem_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    tem_bs   (:,:,:) = 0.0_RP
    tem_bs_pl(:,:,:) = 0.0_RP

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Basic state information ***'
    if( IO_L ) write(IO_FID_LOG,*) '--- Basic state type : ', trim(ref_type)

    if    ( ref_type == 'INPUT' ) then

       call bsstate_input_ref ( ref_fname,  & ! [IN]
                                pre_ref(:), & ! [OUT]
                                tem_ref(:), & ! [OUT]
                                qv_ref (:)  ) ! [OUT]

    elseif( ref_type == 'INIT' ) then

       call bsstate_generate  ( pre_ref(:), & ! [OUT]
                                tem_ref(:), & ! [OUT]
                                qv_ref (:)  ) ! [OUT]

       call bsstate_output_ref( ref_fname,  & ! [IN]
                                pre_ref(:), & ! [IN]
                                tem_ref(:), & ! [IN]
                                qv_ref (:)  ) ! [IN]

    endif

    if ( ref_type == 'INPUT' .OR. ref_type == 'INIT' ) then
       ! set 3-D basic state
       call set_basicstate( pre_ref(:), & ! [IN]
                            tem_ref(:), & ! [IN]
                            qv_ref (:)  ) ! [IN]

       if( IO_L ) write(IO_FID_LOG,*) '-------------------------------------------------------'
       if( IO_L ) write(IO_FID_LOG,*) 'Level   Density  Pressure     Temp. Pot. Tem.        qv'

       do k = ADM_kall, 1, -1
          th_ref (k) = tem_ref(k) * ( PRE00 / pre_ref(k) )**(Rdry/CPdry)
          rho_ref(k) = pre_ref(k) / tem_ref(k) / ( ( 1.0_RP - qv_ref(k) ) * Rdry &
                                              + (          qv_ref(k) ) * Rvap )

          if ( IO_L ) then
             if( k == ADM_kmax ) write(IO_FID_LOG,*) '-------------------------------------------------------'
             write(IO_FID_LOG,'(I4,F12.4,3F10.2,F10.7)') k,rho_ref(k),pre_ref(k),tem_ref(k),th_ref(k),qv_ref(k)
             if( k == ADM_kmin ) write(IO_FID_LOG,*) '-------------------------------------------------------'
          endif
       enddo
    endif

    return
  end subroutine bsstate_setup

  !-----------------------------------------------------------------------------
  subroutine bsstate_input_ref( &
       fname,   &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall
    implicit none

    character(len=*), intent(in)  :: fname
    real(RP),         intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP),         intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP),         intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(DP) :: pre_ref_DP(ADM_kall)
    real(DP) :: tem_ref_DP(ADM_kall)
    real(DP) :: qv_ref_DP (ADM_kall)

    integer  :: fid
    !---------------------------------------------------------------------------

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'old'          )

       read(fid) pre_ref_DP(:)
       read(fid) tem_ref_DP(:)
       read(fid) qv_ref_DP (:)

    close(fid)

    pre_ref(:) = real(pre_ref_DP(:),kind=RP)
    tem_ref(:) = real(tem_ref_DP(:),kind=RP)
    qv_ref (:) = real(qv_ref_DP (:),kind=RP)

    return
  end subroutine bsstate_input_ref

  !-----------------------------------------------------------------------------
  subroutine bsstate_output_ref( &
       fname,   &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_kall
    implicit none

    character(len=*), intent(in)  :: fname
    real(RP),         intent(in)  :: pre_ref(ADM_kall) ! reference pressure
    real(RP),         intent(in)  :: tem_ref(ADM_kall) ! reference temperature
    real(RP),         intent(in)  :: qv_ref (ADM_kall) ! water vapor

    real(DP) :: pre_ref_DP(ADM_kall)
    real(DP) :: tem_ref_DP(ADM_kall)
    real(DP) :: qv_ref_DP (ADM_kall)

    integer :: fid
    !---------------------------------------------------------------------------

    pre_ref_DP(:) = real(pre_ref(:),kind=DP)
    tem_ref_DP(:) = real(tem_ref(:),kind=DP)
    qv_ref_DP (:) = real(qv_ref (:),kind=DP)

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'new'          )

       write(fid) pre_ref_DP(:)
       write(fid) tem_ref_DP(:)
       write(fid) qv_ref_DP (:)

    close(fid)

    return
  end subroutine bsstate_output_ref

  !-----------------------------------------------------------------------------
  subroutine bsstate_generate( &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_gm_statistics, only: &
       GTL_global_mean_eachlayer
    use mod_runconf, only: &
       TRC_vmax, &
       I_QV
    use mod_prgvar, only: &
       prgvar_get_withdiag
    implicit none

    real(RP), intent(out) :: pre_ref(ADM_kall) ! reference pressure
    real(RP), intent(out) :: tem_ref(ADM_kall) ! reference temperature
    real(RP), intent(out) :: qv_ref (ADM_kall) ! water vapor

    real(RP) :: rhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhoge    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogq    (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    real(RP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: pre      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: pre_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: tem      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: tem_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vx       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vy       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vz       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: w        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q        (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    !---------------------------------------------------------------------------

    call prgvar_get_withdiag( rhog,   rhog_pl,   & ! [OUT]
                              rhogvx, rhogvx_pl, & ! [OUT]
                              rhogvy, rhogvy_pl, & ! [OUT]
                              rhogvz, rhogvz_pl, & ! [OUT]
                              rhogw,  rhogw_pl,  & ! [OUT]
                              rhoge,  rhoge_pl,  & ! [OUT]
                              rhogq,  rhogq_pl,  & ! [OUT]
                              rho,    rho_pl,    & ! [OUT]
                              pre,    pre_pl,    & ! [OUT]
                              tem,    tem_pl,    & ! [OUT]
                              vx,     vx_pl,     & ! [OUT]
                              vy,     vy_pl,     & ! [OUT]
                              vz,     vz_pl,     & ! [OUT]
                              w,      w_pl,      & ! [OUT]
                              q,      q_pl       ) ! [OUT]

    call GTL_global_mean_eachlayer( pre(:,:,:),      pre_pl(:,:,:)     , pre_ref(:) )
    call GTL_global_mean_eachlayer( tem(:,:,:),      tem_pl(:,:,:)     , tem_ref(:) )
    call GTL_global_mean_eachlayer( q  (:,:,:,I_QV), q_pl  (:,:,:,I_QV), qv_ref (:) )

    return
  end subroutine bsstate_generate

  !-----------------------------------------------------------------------------
  !> generation of basic state from reference state
  subroutine set_basicstate( &
       pre_ref, &
       tem_ref, &
       qv_ref   )
    use scale_const, only: &
       Rdry => CONST_Rdry, &
       Rvap => CONST_Rvap
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_vmtr, only : &
       VMTR_getIJ_PHI
    use mod_vintrpl, only: &
       VINTRPL_Z2Xi
    use mod_bndcnd, only: &
       BNDCND_thermo
    implicit none

    real(RP), intent(in) :: pre_ref(ADM_kall) ! reference pressure
    real(RP), intent(in) :: tem_ref(ADM_kall) ! reference temperature
    real(RP), intent(in) :: qv_ref (ADM_kall) ! water vapor

    real(RP) :: qv_bs   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qv_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: VMTR_PHI   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_PHI_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: k, l
    !---------------------------------------------------------------------------

    if ( ref_type == 'NOBASE' ) return

    call VMTR_getIJ_PHI( VMTR_PHI, VMTR_PHI_pl )

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       pre_bs(:,k,l) = pre_ref(k)
       tem_bs(:,k,l) = tem_ref(k)
       qv_bs (:,k,l) = qv_ref (k)
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          pre_bs_pl(:,k,l) = pre_ref(k)
          tem_bs_pl(:,k,l) = tem_ref(k)
          qv_bs_pl (:,k,l) = qv_ref (k)
       enddo
       enddo
    endif

    !-- from z-level to zstar-level
    call VINTRPL_Z2Xi( pre_bs(:,:,:), pre_bs_pl(:,:,:) )
    call VINTRPL_Z2Xi( tem_bs(:,:,:), tem_bs_pl(:,:,:) )
    call VINTRPL_Z2Xi( qv_bs (:,:,:), qv_bs_pl (:,:,:) )

    rho_bs(:,:,l) = pre_bs(:,:,l) / tem_bs(:,:,l) / ( ( 1.0_RP-qv_bs(:,:,l) ) * Rdry &
                                                    + (        qv_bs(:,:,l) ) * Rvap )

    if ( ADM_have_pl ) then
       rho_bs_pl(:,:,l) = pre_bs_pl(:,:,l) / tem_bs_pl(:,:,l) / ( ( 1.0_RP-qv_bs_pl(:,:,l) ) * Rdry &
                                                                + (        qv_bs_pl(:,:,l) ) * Rvap )
    endif

    !--- set boundary conditions of basic state
    do l = 1, ADM_lall
       call BNDCND_thermo( ADM_gall,        & ! [IN]
                           rho_bs  (:,:,l), & ! [INOUT]
                           pre_bs  (:,:,l), & ! [INOUT]
                           tem_bs  (:,:,l), & ! [INOUT]
                           VMTR_PHI(:,:,l)  ) ! [IN]
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          call BNDCND_thermo( ADM_gall_pl,        & ! [IN]
                              rho_bs_pl  (:,:,l), & ! [INOUT]
                              pre_bs_pl  (:,:,l), & ! [INOUT]
                              tem_bs_pl  (:,:,l), & ! [INOUT]
                              VMTR_PHI_pl(:,:,l)  ) ! [IN]
       enddo
    endif

    return
  end subroutine set_basicstate

end module mod_bsstate
