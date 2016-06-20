!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Hydrostatic barance
!!
!! @par Description
!!          make hydrostatic profile in the model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_hydrostatic
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_const, only: &
     GRAV    => CONST_GRAV,    &
     Rdry    => CONST_Rdry,    &
     Rvap    => CONST_Rvap,    &
     CVdry   => CONST_CVdry,   &
     CVvap   => CONST_CVvap,   &
     CL      => CONST_CL,      &
     LAPS    => CONST_LAPS,    &
     LAPSdry => CONST_LAPSdry, &
     P00     => CONST_PRE00,   &
     THERMODYN_TYPE => CONST_THERMODYN_TYPE
  use scale_grid, only: &
     GRID_CZ, &
     GRID_FDZ
  use scale_grid_real, only: &
     REAL_CZ, &
     REAL_FZ
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_HYDROSTATIC_setup
  public :: ATMOS_HYDROSTATIC_buildrho
  public :: ATMOS_HYDROSTATIC_buildrho_real
  public :: ATMOS_HYDROSTATIC_buildrho_atmos
  public :: ATMOS_HYDROSTATIC_buildrho_bytemp
  public :: ATMOS_HYDROSTATIC_buildrho_bytemp_atmos

  public :: ATMOS_HYDROSTATIC_buildrho_atmos_0D
  public :: ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D
  public :: ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D

  public :: ATMOS_HYDROSTATIC_barometric_law_mslp
  public :: ATMOS_HYDROSTATIC_barometric_law_pres

  interface ATMOS_HYDROSTATIC_buildrho
     module procedure ATMOS_HYDROSTATIC_buildrho_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_3D
  end interface ATMOS_HYDROSTATIC_buildrho

  interface ATMOS_HYDROSTATIC_buildrho_real
!     module procedure ATMOS_HYDROSTATIC_buildrho_1D ! buildrho_real_1D is completely equal to buildrho_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_real_3D
  end interface ATMOS_HYDROSTATIC_buildrho_real

  interface ATMOS_HYDROSTATIC_buildrho_atmos
     module procedure ATMOS_HYDROSTATIC_buildrho_atmos_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_atmos_3D
  end interface ATMOS_HYDROSTATIC_buildrho_atmos

  interface ATMOS_HYDROSTATIC_buildrho_bytemp
     module procedure ATMOS_HYDROSTATIC_buildrho_bytemp_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_bytemp_3D
  end interface ATMOS_HYDROSTATIC_buildrho_bytemp

  interface ATMOS_HYDROSTATIC_buildrho_bytemp_atmos
     module procedure ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D
  end interface ATMOS_HYDROSTATIC_buildrho_bytemp_atmos

  interface ATMOS_HYDROSTATIC_barometric_law_mslp
     module procedure ATMOS_HYDROSTATIC_barometric_law_mslp_0D
     module procedure ATMOS_HYDROSTATIC_barometric_law_mslp_2D
  end interface ATMOS_HYDROSTATIC_barometric_law_mslp

  interface ATMOS_HYDROSTATIC_barometric_law_pres
     module procedure ATMOS_HYDROSTATIC_barometric_law_pres_0D
     module procedure ATMOS_HYDROSTATIC_barometric_law_pres_2D
  end interface ATMOS_HYDROSTATIC_barometric_law_pres

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_HYDROSTATIC_buildrho_atmos_2D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter :: itelim = 100 !< itelation number limit
  real(RP), private              :: criteria                            !< convergence judgement criteria
  logical,  private              :: HYDROSTATIC_uselapserate  = .false. !< use lapse rate?
  integer,  private              :: HYDROSTATIC_buildrho_real_kref = 1


  real(RP), private :: CV_qv
  real(RP), private :: CV_qc

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_HYDROSTATIC_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_EPS
    implicit none

    NAMELIST / PARAM_ATMOS_HYDROSTATIC / &
       HYDROSTATIC_uselapserate, &
       HYDROSTATIC_buildrho_real_kref

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[HYDROSTATIC] / Categ[ATMOS SHARE] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_HYDROSTATIC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_HYDROSTATIC. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_HYDROSTATIC)

    criteria = sqrt( CONST_EPS )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** use lapse rate for estimation of surface temperature? : ', HYDROSTATIC_uselapserate
    if( IO_L ) write(IO_FID_LOG,*) '*** buildrho conversion criteria : ', criteria

    if ( THERMODYN_TYPE == 'EXACT' ) then
       CV_qv = CVvap
       CV_qc = CL
    elseif( THERMODYN_TYPE == 'SIMPLE' ) then
       CV_qv = CVdry
       CV_qc = CVdry
    endif

    return
  end subroutine ATMOS_HYDROSTATIC_setup

  !-----------------------------------------------------------------------------
  !> Build up density from surface (1D)
  subroutine ATMOS_HYDROSTATIC_buildrho_1D( &
       dens,     &
       temp,     &
       pres,     &
       pott,     &
       qv,       &
       qc,       &
       temp_sfc, &
       pres_sfc, &
       pott_sfc, &
       qv_sfc,   &
       qc_sfc    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens(KA) !< density               [kg/m3]
    real(RP), intent(out) :: temp(KA) !< temperature           [K]
    real(RP), intent(out) :: pres(KA) !< pressure              [Pa]
    real(RP), intent(in)  :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)  :: qc  (KA) !< liquid water          [kg/kg]

    real(RP), intent(out) :: temp_sfc !< surface temperature           [K]
    real(RP), intent(in)  :: pres_sfc !< surface pressure              [Pa]
    real(RP), intent(in)  :: pott_sfc !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc   !< surface water vapor           [kg/kg]
    real(RP), intent(in)  :: qc_sfc   !< surface liquid water          [kg/kg]

    real(RP) :: dens_sfc

    real(RP) :: Rtot_sfc
    real(RP) :: CVtot_sfc
    real(RP) :: CPovCV_sfc
    real(RP) :: Rtot
    real(RP) :: CVtot
    real(RP) :: CPovCV

    real(RP) :: CVovCP_sfc, CPovR, CVovCP, RovCV
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere

    Rtot_sfc   = Rdry  * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + Rvap  * qv_sfc
    CVtot_sfc  = CVdry * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + CV_qv * qv_sfc                       &
               + CV_qc * qc_sfc
    CPovCV_sfc = ( CVtot_sfc + Rtot_sfc ) / CVtot_sfc

    Rtot   = Rdry  * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + Rvap  * qv(KS)
    CVtot  = CVdry * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + CV_qv * qv(KS)                       &
           + CV_qc * qc(KS)
    CPovCV = ( CVtot + Rtot ) / CVtot

    ! density at surface
    CVovCP_sfc = 1.0_RP / CPovCV_sfc
    dens_sfc   = P00 / Rtot_sfc / pott_sfc * ( pres_sfc/P00 )**CVovCP_sfc
    temp_sfc   = pres_sfc / ( dens_sfc * Rtot_sfc )

    ! make density at lowermost cell center
    if ( HYDROSTATIC_uselapserate ) then

       CPovR  = ( CVtot + Rtot ) / Rtot
       CVovCP = 1.0_RP / CPovCV

       temp(KS) = pott_sfc - LAPSdry * GRID_CZ(KS) ! use dry lapse rate
       pres(KS) = P00 * ( temp(KS)/pott(KS) )**CPovR
       dens(KS) = P00 / Rtot / pott(KS) * ( pres(KS)/P00 )**CVovCP

    else ! use itelation

       RovCV = Rtot / CVtot

       dens_s   = 0.0_RP
       dens(KS) = dens_sfc ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(KS)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(KS)

          dhyd = + ( P00 * ( dens_sfc * Rtot_sfc * pott_sfc / P00 )**CPovCV_sfc &
                   - P00 * ( dens_s   * Rtot     * pott(KS) / P00 )**CPovCV     ) / GRID_CZ(KS) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens_sfc + dens_s )                                     ! rho*g

          dgrd = - P00 * ( Rtot * pott(KS) / P00 )**CPovCV / GRID_CZ(KS) &
                 * CPovCV * dens_s**RovCV                           &
                 - 0.5_RP * GRAV

          dens(KS) = dens_s - dhyd/dgrd

          if( dens(KS)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [buildrho 1D sfc] iteration not converged!', &
                     dens(KS),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif

    endif

    !--- from lowermost atmosphere to top of atmosphere
    call ATMOS_HYDROSTATIC_buildrho_atmos_1D( dens(:), & ! [INOUT]
                                              temp(:), & ! [OUT]
                                              pres(:), & ! [OUT]
                                              pott(:), & ! [IN]
                                              qv  (:), & ! [IN]
                                              qc  (:)  ) ! [IN]

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_1D

  !-----------------------------------------------------------------------------
  !> Build up density from surface (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_3D( &
       dens,     &
       temp,     &
       pres,     &
       pott,     &
       qv,       &
       qc,       &
       temp_sfc, &
       pres_sfc, &
       pott_sfc, &
       qv_sfc,   &
       qc_sfc    )
    use scale_process, only: &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_horizontal_mean
    implicit none

    real(RP), intent(out) :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out) :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(out) :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)  :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)  :: qc  (KA,IA,JA) !< liquid water          [kg/kg]

    real(RP), intent(out) :: temp_sfc(1,IA,JA) !< surface temperature           [K]
    real(RP), intent(in)  :: pres_sfc(1,IA,JA) !< surface pressure              [Pa]
    real(RP), intent(in)  :: pott_sfc(1,IA,JA) !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc  (1,IA,JA) !< surface water vapor           [kg/kg]
    real(RP), intent(in)  :: qc_sfc  (1,IA,JA) !< surface liquid water          [kg/kg]

    real(RP) :: dz(KA,IA,JA)

    real(RP) :: dens_sfc  (1,IA,JA)
    real(RP) :: pott_toa  (IA,JA) !< surface potential temperature [K]
    real(RP) :: qv_toa    (IA,JA) !< surface water vapor           [kg/kg]
    real(RP) :: qc_toa    (IA,JA) !< surface liquid water          [kg/kg]
    real(RP) :: dens_1D   (KA)

    real(RP) :: Rtot_sfc  (IA,JA)
    real(RP) :: CVtot_sfc (IA,JA)
    real(RP) :: CPovCV_sfc(IA,JA)
    real(RP) :: Rtot      (IA,JA)
    real(RP) :: CVtot     (IA,JA)
    real(RP) :: CPovCV    (IA,JA)

    real(RP) :: CVovCP_sfc, CPovR, CVovCP

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere

    do j = JSB, JEB
    do i = ISB, IEB
       dz(KS,i,j) = REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j) ! distance from surface to cell center
       do k = KS+1, KE
          dz(k,i,j) = REAL_CZ(k,i,j) - REAL_CZ(k-1,i,j) ! distance from cell center to cell center
       enddo
       dz(KE+1,i,j) = REAL_FZ(KE,i,j) - REAL_CZ(KE,i,j) ! distance from cell center to TOA
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       pott_toa(i,j) = pott(KE,i,j)
       qv_toa  (i,j) = qv  (KE,i,j)
       qc_toa  (i,j) = qc  (KE,i,j)
    enddo
    enddo

    ! density at surface
    do j = JSB, JEB
    do i = ISB, IEB
       Rtot_sfc  (i,j) = Rdry  * ( 1.0_RP - qv_sfc(1,i,j) - qc_sfc(1,i,j) ) &
                       + Rvap  * qv_sfc(1,i,j)
       CVtot_sfc (i,j) = CVdry * ( 1.0_RP - qv_sfc(1,i,j) - qc_sfc(1,i,j) ) &
                       + CV_qv * qv_sfc(1,i,j)                              &
                       + CV_qc * qc_sfc(1,i,j)
       CPovCV_sfc(i,j) = ( CVtot_sfc(i,j) + Rtot_sfc(i,j) ) / CVtot_sfc(i,j)
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       Rtot  (i,j) = Rdry  * ( 1.0_RP - qv(KS,i,j) - qc(KS,i,j) ) &
                   + Rvap  * qv(KS,i,j)
       CVtot (i,j) = CVdry * ( 1.0_RP - qv(KS,i,j) - qc(KS,i,j) ) &
                   + CV_qv * qv(KS,i,j)                           &
                   + CV_qc * qc(KS,i,j)
       CPovCV(i,j) = ( CVtot(i,j) + Rtot(i,j) ) / CVtot(i,j)
    enddo
    enddo

    ! density at surface
    do j = JSB, JEB
    do i = ISB, IEB
       CVovCP_sfc      = 1.0_RP / CPovCV_sfc(i,j)
       dens_sfc(1,i,j) = P00 / Rtot_sfc(i,j) / pott_sfc(1,i,j) * ( pres_sfc(1,i,j)/P00 )**CVovCP_sfc
       temp_sfc(1,i,j) = pres_sfc(1,i,j) / ( dens_sfc(1,i,j) * Rtot_sfc(i,j) )
    enddo
    enddo

    ! make density at lowermost cell center
    if ( HYDROSTATIC_uselapserate ) then

       do j = JSB, JEB
       do i = ISB, IEB
          CPovR  = ( CVtot(i,j) + Rtot(i,j) ) / Rtot(i,j)
          CVovCP = 1.0_RP / CPovCV(i,j)

          temp(KS,i,j) = pott_sfc(1,i,j) - LAPSdry * ( REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j) ) ! use dry lapse rate
          pres(KS,i,j) = P00 * ( temp(KS,i,j)/pott(KS,i,j) )**CPovR
          dens(KS,i,j) = P00 / Rtot(i,j) / pott(KS,i,j) * ( pres(KS,i,j)/P00 )**CVovCP
       enddo
       enddo

    else ! use itelation

       call ATMOS_HYDROSTATIC_buildrho_atmos_2D( dens    (KS,:,:), & ! [OUT]
                                                 temp    (KS,:,:), & ! [OUT]->not used
                                                 pres    (KS,:,:), & ! [OUT]->not used
                                                 pott    (KS,:,:), & ! [IN]
                                                 qv      (KS,:,:), & ! [IN]
                                                 qc      (KS,:,:), & ! [IN]
                                                 dens_sfc(1,:,:),  & ! [IN]
                                                 pott_sfc(1,:,:),  & ! [IN]
                                                 qv_sfc  (1,:,:),  & ! [IN]
                                                 qc_sfc  (1,:,:),  & ! [IN]
                                                 dz      (KS,:,:), & ! [IN]
                                                 KS                ) ! [IN]

    endif

    !--- from lowermost atmosphere to top of atmosphere
    call ATMOS_HYDROSTATIC_buildrho_atmos_3D( dens(:,:,:), & ! [INOUT]
                                              temp(:,:,:), & ! [OUT]
                                              pres(:,:,:), & ! [OUT]
                                              pott(:,:,:), & ! [IN]
                                              qv  (:,:,:), & ! [IN]
                                              qc  (:,:,:), & ! [IN]
                                              dz  (:,:,:)  ) ! [IN]

    call ATMOS_HYDROSTATIC_buildrho_atmos_2D( dens    (KE+1,:,:), & ! [OUT]
                                              temp    (KE+1,:,:), & ! [OUT]->not used
                                              pres    (KE+1,:,:), & ! [OUT]->not used
                                              pott_toa(:,:),      & ! [IN]
                                              qv_toa  (:,:),      & ! [IN]
                                              qc_toa  (:,:),      & ! [IN]
                                              dens    (KE  ,:,:), & ! [IN]
                                              pott    (KE  ,:,:), & ! [IN]
                                              qv      (KE  ,:,:), & ! [IN]
                                              qc      (KE  ,:,:), & ! [IN]
                                              dz      (KE+1,:,:), & ! [IN]
                                              KE+1                ) ! [IN]

    ! density at TOA
    dens(   1:KS-1,:,:) = 0.D0 ! fill dummy
    dens(KE+2:KA  ,:,:) = 0.D0 ! fill dummy

    call COMM_horizontal_mean( dens_1D(:), dens(:,:,:) )
    do j = JSB, JEB
    do i = ISB, IEB
       dens(:,i,j) = dens_1D(:)
    enddo
    enddo

    call ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D( dens    (KE  ,:,:), & ! [OUT]
                                                  temp    (KE  ,:,:), & ! [OUT]->not used
                                                  pres    (KE  ,:,:), & ! [OUT]->not used
                                                  pott    (KE  ,:,:), & ! [IN]
                                                  qv      (KE  ,:,:), & ! [IN]
                                                  qc      (KE  ,:,:), & ! [IN]
                                                  dens    (KE+1,:,:), & ! [IN]
                                                  pott_toa(:,:),      & ! [IN]
                                                  qv_toa  (:,:),      & ! [IN]
                                                  qc_toa  (:,:),      & ! [IN]
                                                  dz      (KE+1,:,:), & ! [IN]
                                                  KE+1                ) ! [IN]

    !--- from top of atmosphere to lowermost atmosphere
    call ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D( dens(:,:,:), & ! [INOUT]
                                                  temp(:,:,:), & ! [OUT]
                                                  pres(:,:,:), & ! [OUT]
                                                  pott(:,:,:), & ! [IN]
                                                  qv  (:,:,:), & ! [IN]
                                                  qc  (:,:,:), & ! [IN]
                                                  dz  (:,:,:)  ) ! [IN]

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_3D

  !-----------------------------------------------------------------------------
  !> Build up density from surface (3D), not to reverse from TOA
  subroutine ATMOS_HYDROSTATIC_buildrho_real_3D( &
       dens, &
       temp, &
       pres, &
       pott, &
       qv,   &
       qc    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out)   :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out)   :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(inout) :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)    :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]

    real(RP) :: dz(KA,IA,JA)

    real(RP) :: pott_toa(IA,JA) !< surface potential temperature [K]
    real(RP) :: qv_toa  (IA,JA) !< surface water vapor           [kg/kg]
    real(RP) :: qc_toa  (IA,JA) !< surface liquid water          [kg/kg]

    real(RP) :: Rtot
    real(RP) :: CVtot
    real(RP) :: CVovCP

    integer  :: kref
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere

    do j = JSB, JEB
    do i = ISB, IEB
       dz(KS,i,j) = REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j) ! distance from surface to cell center
       do k = KS+1, KE
          dz(k,i,j) = REAL_CZ(k,i,j) - REAL_CZ(k-1,i,j) ! distance from cell center to cell center
       enddo
       dz(KE+1,i,j) = REAL_FZ(KE,i,j) - REAL_CZ(KE,i,j) ! distance from cell center to TOA
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       pott_toa(i,j) = pott(KE,i,j)
       qv_toa  (i,j) = qv  (KE,i,j)
       qc_toa  (i,j) = qc  (KE,i,j)
    enddo
    enddo

    kref = HYDROSTATIC_buildrho_real_kref + KS - 1

    ! calc density at reference level
    do j = JSB, JEB
    do i = ISB, IEB
       Rtot = Rdry  * ( 1.0_RP - qv(kref,i,j) - qc(kref,i,j) ) &
            + Rvap  * qv(kref,i,j)
       CVtot = CVdry * ( 1.0_RP - qv(kref,i,j) - qc(kref,i,j) ) &
             + CV_qv * qv(kref,i,j)                           &
             + CV_qc * qc(kref,i,j)
       CVovCP = CVtot / ( CVtot + Rtot )
       dens(kref,i,j) = P00 / ( Rtot * pott(kref,i,j) ) * ( pres(kref,i,j)/P00 )**CVovCP
    enddo
    enddo

    !--- from lowermost atmosphere to top of atmosphere
    call ATMOS_HYDROSTATIC_buildrho_atmos_3D( dens(:,:,:), & ! [INOUT]
                                              temp(:,:,:), & ! [OUT]
                                              pres(:,:,:), & ! [OUT]
                                              pott(:,:,:), & ! [IN]
                                              qv  (:,:,:), & ! [IN]
                                              qc  (:,:,:), & ! [IN]
                                              dz  (:,:,:), & ! [IN]
                                              kref         ) ! [IN]
    call ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D( dens(:,:,:), & ! [INOUT]
                                                  temp(:,:,:), & ! [OUT]
                                                  pres(:,:,:), & ! [OUT]
                                                  pott(:,:,:), & ! [IN]
                                                  qv  (:,:,:), & ! [IN]
                                                  qc  (:,:,:), & ! [IN]
                                                  dz  (:,:,:), & ! [IN]
                                                  kref         ) ! [IN]

    call ATMOS_HYDROSTATIC_buildrho_atmos_2D( dens    (KE+1,:,:), & ! [OUT]
                                              temp    (KE+1,:,:), & ! [OUT]->not used
                                              pres    (KE+1,:,:), & ! [OUT]->not used
                                              pott_toa(:,:),      & ! [IN]
                                              qv_toa  (:,:),      & ! [IN]
                                              qc_toa  (:,:),      & ! [IN]
                                              dens    (KE  ,:,:), & ! [IN]
                                              pott    (KE  ,:,:), & ! [IN]
                                              qv      (KE  ,:,:), & ! [IN]
                                              qc      (KE  ,:,:), & ! [IN]
                                              dz      (KE+1,:,:), & ! [IN]
                                              KE+1                ) ! [IN]

    ! density at TOA
    dens(   1:KS-1,:,:) = 0.D0 ! fill dummy
    dens(KE+2:KA  ,:,:) = 0.D0 ! fill dummy

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_real_3D

  !-----------------------------------------------------------------------------
  !> Build up density (0D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_0D( &
       dens_L2, &
       temp_L2, &
       pres_L2, &
       pott_L2, &
       qv_L2,   &
       qc_L2,   &
       dens_L1, &
       pott_L1, &
       qv_L1,   &
       qc_L1,   &
       dz,      &
       k        )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens_L2 !< density               at layer 2 [kg/m3]
    real(RP), intent(out) :: temp_L2 !< temperature           at layer 2 [K]
    real(RP), intent(out) :: pres_L2 !< pressure              at layer 2 [Pa]
    real(RP), intent(in)  :: pott_L2 !< potential temperature at layer 2 [K]
    real(RP), intent(in)  :: qv_L2   !< water vapor           at layer 2 [kg/kg]
    real(RP), intent(in)  :: qc_L2   !< liquid water          at layer 2 [kg/kg]
    real(RP), intent(in)  :: dens_L1 !< density               at layer 1 [Pa]
    real(RP), intent(in)  :: pott_L1 !< potential temperature at layer 1 [K]
    real(RP), intent(in)  :: qv_L1   !< water vapor           at layer 1 [kg/kg]
    real(RP), intent(in)  :: qc_L1   !< liquid water          at layer 1 [kg/kg]
    real(RP), intent(in)  :: dz      !< distance from layer 1 to layer 2 [m]
    integer,  intent(in)  :: k       !< for monitor

    real(RP) :: Rtot_L1  , Rtot_L2
    real(RP) :: CVtot_L1 , CVtot_L2
    real(RP) :: CPovCV_L1, CPovCV_L2

    real(RP) :: RovCV
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged
    !---------------------------------------------------------------------------

    Rtot_L1   = Rdry  * ( 1.0_RP - qv_L1 - qc_L1 ) &
              + Rvap  * qv_L1
    CVtot_L1  = CVdry * ( 1.0_RP - qv_L1 - qc_L1 ) &
              + CV_qv * qv_L1                      &
              + CV_qc * qc_L1
    CPovCV_L1 = ( CVtot_L1 + Rtot_L1 ) / CVtot_L1

    Rtot_L2   = Rdry  * ( 1.0_RP - qv_L2 - qc_L2 ) &
              + Rvap  * qv_L2
    CVtot_L2  = CVdry * ( 1.0_RP - qv_L2 - qc_L2 ) &
              + CV_qv * qv_L2                      &
              + CV_qc * qc_L2
    CPovCV_L2 = ( CVtot_L2 + Rtot_L2 ) / CVtot_L2

    RovCV = Rtot_L2 / CVtot_L2

    dens_s  = 0.0_RP
    dens_L2 = dens_L1 ! first guess

    converged = .false.
    do ite = 1, itelim
       if ( abs(dens_L2-dens_s) <= criteria ) then
          converged = .true.
          exit
       endif

       dens_s = dens_L2

       dhyd = + ( P00 * ( dens_L1 * Rtot_L1 * pott_L1 / P00 )**CPovCV_L1 &
                - P00 * ( dens_s  * Rtot_L2 * pott_L2 / P00 )**CPovCV_L2 ) / dz & ! dp/dz
              - GRAV * 0.5_RP * ( dens_L1 + dens_s )                              ! rho*g

       dgrd = - P00 * ( Rtot_L2 * pott_L2 / P00 )**CPovCV_L2 / dz &
              * CPovCV_L2 * dens_s**RovCV                         &
              - 0.5_RP * GRAV

       dens_L2 = dens_s - dhyd/dgrd

       if( dens_L2*0.0_RP /= 0.0_RP) exit
    enddo

    if ( .NOT. converged ) then
       write(*,*) 'xxx [buildrho 0D atmos] iteration not converged!', &
                  k,dens_L2,ite,dens_s,dhyd,dgrd
       call PRC_MPIstop
    endif

    pres_L2 = P00 * ( dens_L2 * Rtot_L2 * pott_L2 / P00 )**CPovCV_L2
    temp_L2 = pres_L2 / ( dens_L2 * Rtot_L2 )

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_0D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (1D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_1D( &
       dens, &
       temp, &
       pres, &
       pott, &
       qv,   &
       qc    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: dens(KA) !< density               [kg/m3]
    real(RP), intent(out)   :: temp(KA) !< temperature           [K]
    real(RP), intent(out)   :: pres(KA) !< pressure              [Pa]
    real(RP), intent(in)    :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA) !< liquid water          [kg/kg]

    real(RP) :: Rtot  (KA)
    real(RP) :: CVtot (KA)
    real(RP) :: CPovCV(KA)

    real(RP) :: RovCV
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       Rtot  (k) = Rdry  * ( 1.0_RP - qv(k) - qc(k) ) &
                 + Rvap  * qv(k)
       CVtot (k) = CVdry * ( 1.0_RP - qv(k) - qc(k) ) &
                 + CV_qv * qv(k)                      &
                 + CV_qc * qc(k)
       CPovCV(k) = ( CVtot(k) + Rtot(k) ) / CVtot(k)
    enddo

    do k = KS+1, KE
       RovCV = Rtot(k) / CVtot(k)

       dens_s  = 0.0_RP
       dens(k) = dens(k-1) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(k)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(k)

          dhyd = + ( P00 * ( dens(k-1) * Rtot(k-1) * pott(k-1) / P00 )**CPovCV(k-1) &
                   - P00 * ( dens_s    * Rtot(k  ) * pott(k  ) / P00 )**CPovCV(k  ) ) / GRID_FDZ(k-1) & ! dpdz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s )                                          ! rho*g

          dgrd = - P00 * ( Rtot(k) * pott(k) / P00 )**CPovCV(k) / GRID_FDZ(k-1) &
                 * CPovCV(k) * dens_s**RovCV                               &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

          if( dens(k)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [buildrho 1D atmos] iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo

    do k = KS, KE
       pres(k) = P00 * ( dens(k) * Rtot(k) * pott(k) / P00 )**CPovCV(k)
       temp(k) = pres(k) / ( dens(k) * Rtot(k) )
    enddo

    dens(   1:KS-1) = dens(KS)
    dens(KE+1:KA  ) = dens(KE)
    pres(   1:KS-1) = pres(KS)
    pres(KE+1:KA  ) = pres(KE)
    temp(   1:KS-1) = temp(KS)
    temp(KE+1:KA  ) = temp(KE)

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_1D

  !-----------------------------------------------------------------------------
  !> Build up density (2D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_2D( &
       dens_L2, &
       temp_L2, &
       pres_L2, &
       pott_L2, &
       qv_L2,   &
       qc_L2,   &
       dens_L1, &
       pott_L1, &
       qv_L1,   &
       qc_L1,   &
       dz,      &
       k        )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens_L2(IA,JA) !< density               at layer 2 [kg/m3]
    real(RP), intent(out) :: temp_L2(IA,JA) !< temperature           at layer 2 [K]
    real(RP), intent(out) :: pres_L2(IA,JA) !< pressure              at layer 2 [Pa]
    real(RP), intent(in)  :: pott_L2(IA,JA) !< potential temperature at layer 2 [K]
    real(RP), intent(in)  :: qv_L2  (IA,JA) !< water vapor           at layer 2 [kg/kg]
    real(RP), intent(in)  :: qc_L2  (IA,JA) !< liquid water          at layer 2 [kg/kg]
    real(RP), intent(in)  :: dens_L1(IA,JA) !< density               at layer 1 [Pa]
    real(RP), intent(in)  :: pott_L1(IA,JA) !< potential temperature at layer 1 [K]
    real(RP), intent(in)  :: qv_L1  (IA,JA) !< water vapor           at layer 1 [kg/kg]
    real(RP), intent(in)  :: qc_L1  (IA,JA) !< liquid water          at layer 1 [kg/kg]
    real(RP), intent(in)  :: dz     (IA,JA) !< distance from layer 1 to layer 2 [m]
    integer,  intent(in)  :: k              !< for monitor

    real(RP) :: Rtot_L1  (IA,JA), Rtot_L2  (IA,JA)
    real(RP) :: CVtot_L1 (IA,JA), CVtot_L2 (IA,JA)
    real(RP) :: CPovCV_L1(IA,JA), CPovCV_L2(IA,JA)

    real(RP) :: RovCV
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       Rtot_L1  (i,j) = Rdry  * ( 1.0_RP - qv_L1(i,j) - qc_L1(i,j) ) &
                      + Rvap  * qv_L1(i,j)
       CVtot_L1 (i,j) = CVdry * ( 1.0_RP - qv_L1(i,j) - qc_L1(i,j) ) &
                      + CV_qv * qv_L1(i,j)                           &
                      + CV_qc * qc_L1(i,j)
       CPovCV_L1(i,j) = ( CVtot_L1(i,j) + Rtot_L1(i,j) ) / CVtot_L1(i,j)

       Rtot_L2  (i,j) = Rdry  * ( 1.0_RP - qv_L2(i,j) - qc_L2(i,j) ) &
                      + Rvap  * qv_L2(i,j)
       CVtot_L2 (i,j) = CVdry * ( 1.0_RP - qv_L2(i,j) - qc_L2(i,j) ) &
                      + CV_qv * qv_L2(i,j)                           &
                      + CV_qc * qc_L2(i,j)
       CPovCV_L2(i,j) = ( CVtot_L2(i,j) + Rtot_L2(i,j) ) / CVtot_L2(i,j)
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       RovCV = Rtot_L2(i,j) / CVtot_L2(i,j)

       dens_s       = 0.0_RP
       dens_L2(i,j) = dens_L1(i,j) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens_L2(i,j)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens_L2(i,j)

          dhyd = + ( P00 * ( dens_L1(i,j) * Rtot_L1(i,j) * pott_L1(i,j) / P00 )**CPovCV_L1(i,j) &
                   - P00 * ( dens_s       * Rtot_L2(i,j) * pott_L2(i,j) / P00 )**CPovCV_L2(i,j) ) / dz(i,j) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens_L1(i,j) + dens_s )                                                  ! rho*g

          dgrd = - P00 * ( Rtot_L2(i,j) * pott_L2(i,j) / P00 )**CPovCV_L2(i,j) / dz(i,j) &
                 * CPovCV_L2(i,j) * dens_s**RovCV                                        &
                 - 0.5_RP * GRAV

          dens_L2(i,j) = dens_s - dhyd/dgrd

          if( dens_L2(i,j)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [buildrho 2D atmos] iteration not converged!', &
                     i,j,k,dens_L2(i,j),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       pres_L2(i,j) = P00 * ( dens_L2(i,j) * Rtot_L2(i,j) * pott_L2(i,j) / P00 )**CPovCV_L2(i,j)
       temp_L2(i,j) = pres_L2(i,j) / ( dens_L2(i,j) * Rtot_L2(i,j) )
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_2D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_3D( &
       dens, &
       temp, &
       pres, &
       pott, &
       qv,   &
       qc,   &
       dz,   &
       kref_in )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out)   :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(out)   :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)    :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA,IA,JA) !< distance between the layer (center) [m]
    integer, intent(in), optional :: kref_in

    real(RP) :: Rtot  (KA,IA,JA)
    real(RP) :: CVtot (KA,IA,JA)
    real(RP) :: CPovCV(KA,IA,JA)

    real(RP) :: RovCV
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: kref
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if ( present(kref_in) ) then
       kref = kref_in
    else
       kref = KS
    end if

    do j = JSB, JEB
    do i = ISB, IEB
    do k = kref, KE
       Rtot  (k,i,j) = Rdry  * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                     + Rvap  * qv(k,i,j)
       CVtot (k,i,j) = CVdry * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                     + CV_qv * qv(k,i,j)                          &
                     + CV_qc * qc(k,i,j)
       CPovCV(k,i,j) = ( CVtot(k,i,j) + Rtot(k,i,j) ) / CVtot(k,i,j)
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = kref+1, KE
       RovCV = Rtot(k,i,j) / CVtot(k,i,j)

       dens_s      = 0.0_RP
       dens(k,i,j) = dens(k-1,i,j) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(k,i,j)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(k,i,j)

          dhyd = + ( P00 * ( dens(k-1,i,j) * Rtot(k-1,i,j) * pott(k-1,i,j) / P00 )**CPovCV(k-1,i,j) &
                   - P00 * ( dens_s        * Rtot(k  ,i,j) * pott(k  ,i,j) / P00 )**CPovCV(k  ,i,j) ) / dz(k,i,j) & ! dpdz
                 - GRAV * 0.5_RP * ( dens(k-1,i,j) + dens_s )                                                       ! rho*g

          dgrd = - P00 * ( Rtot(k,i,j) * pott(k,i,j) / P00 )**CPovCV(k,i,j) / dz(k,i,j) &
                 * CPovCV(k,i,j) * dens_s**RovCV                                        &
                 - 0.5_RP * GRAV

          dens(k,i,j) = dens_s - dhyd/dgrd

          if( dens(k,i,j)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [buildrho 3D atmos] iteration not converged!', &
                     k,i,j,dens(k,i,j),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = kref, KE
       pres(k,i,j) = P00 * ( dens(k,i,j) * Rtot(k,i,j) * pott(k,i,j) / P00 )**CPovCV(k,i,j)
       temp(k,i,j) = pres(k,i,j) / ( dens(k,i,j) * Rtot(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_3D

  !-----------------------------------------------------------------------------
  !> Build up density (2D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D( &
       dens_L1, &
       temp_L1, &
       pres_L1, &
       pott_L1, &
       qv_L1,   &
       qc_L1,   &
       dens_L2, &
       pott_L2, &
       qv_L2,   &
       qc_L2,   &
       dz,      &
       k        )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens_L1(IA,JA) !< density               at layer 1 [kg/m3]
    real(RP), intent(out) :: temp_L1(IA,JA) !< temperature           at layer 1 [K]
    real(RP), intent(out) :: pres_L1(IA,JA) !< pressure              at layer 1 [Pa]
    real(RP), intent(in)  :: pott_L1(IA,JA) !< potential temperature at layer 1 [K]
    real(RP), intent(in)  :: qv_L1  (IA,JA) !< water vapor           at layer 1 [kg/kg]
    real(RP), intent(in)  :: qc_L1  (IA,JA) !< liquid water          at layer 1 [kg/kg]
    real(RP), intent(in)  :: dens_L2(IA,JA) !< density               at layer 2 [Pa]
    real(RP), intent(in)  :: pott_L2(IA,JA) !< potential temperature at layer 2 [K]
    real(RP), intent(in)  :: qv_L2  (IA,JA) !< water vapor           at layer 2 [kg/kg]
    real(RP), intent(in)  :: qc_L2  (IA,JA) !< liquid water          at layer 2 [kg/kg]
    real(RP), intent(in)  :: dz     (IA,JA) !< distance from layer 1 to layer 2 [m]
    integer,  intent(in)  :: k              !< for monitor

    real(RP) :: Rtot_L1  (IA,JA), Rtot_L2  (IA,JA)
    real(RP) :: CVtot_L1 (IA,JA), CVtot_L2 (IA,JA)
    real(RP) :: CPovCV_L1(IA,JA), CPovCV_L2(IA,JA)

    real(RP) :: RovCV
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       Rtot_L1  (i,j) = Rdry  * ( 1.0_RP - qv_L1(i,j) - qc_L1(i,j) ) &
                      + Rvap  * qv_L1(i,j)
       CVtot_L1 (i,j) = CVdry * ( 1.0_RP - qv_L1(i,j) - qc_L1(i,j) ) &
                      + CV_qv * qv_L1(i,j)                           &
                      + CV_qc * qc_L1(i,j)
       CPovCV_L1(i,j) = ( CVtot_L1(i,j) + Rtot_L1(i,j) ) / CVtot_L1(i,j)

       Rtot_L2  (i,j) = Rdry  * ( 1.0_RP - qv_L2(i,j) - qc_L2(i,j) ) &
                      + Rvap  * qv_L2(i,j)
       CVtot_L2 (i,j) = CVdry * ( 1.0_RP - qv_L2(i,j) - qc_L2(i,j) ) &
                      + CV_qv * qv_L2(i,j)                           &
                      + CV_qc * qc_L2(i,j)
       CPovCV_L2(i,j) = ( CVtot_L2(i,j) + Rtot_L2(i,j) ) / CVtot_L2(i,j)
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       RovCV = Rtot_L1(i,j) / CVtot_L1(i,j)

       dens_s       = 0.0_RP
       dens_L1(i,j) = dens_L2(i,j) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens_L1(i,j)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens_L1(i,j)

          dhyd = + ( P00 * ( dens_s       * Rtot_L1(i,j) * pott_L1(i,j) / P00 )**CPovCV_L1(i,j) &
                   - P00 * ( dens_L2(i,j) * Rtot_L2(i,j) * pott_L2(i,j) / P00 )**CPovCV_L2(i,j) ) / dz(i,j) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens_s + dens_L2(i,j) )                                                  ! rho*g

          dgrd = + P00 * ( Rtot_L1(i,j) * pott_L1(i,j) / P00 )**CPovCV_L1(i,j) / dz(i,j) &
                 * CPovCV_L1(i,j) * dens_s**RovCV                                        &
                 - 0.5_RP * GRAV

          dens_L1(i,j) = dens_s - dhyd/dgrd

          if( dens_L1(i,j)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [buildrho 2D rev atmos] iteration not converged!', &
                     i,j,k,dens_L1(i,j),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       pres_L1(i,j) = P00 * ( dens_L1(i,j) * Rtot_L1(i,j) * pott_L1(i,j) / P00 )**CPovCV_L1(i,j)
       temp_L1(i,j) = pres_L1(i,j) / ( dens_L1(i,j) * Rtot_L1(i,j) )
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_2D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D( &
       dens, &
       temp, &
       pres, &
       pott, &
       qv,   &
       qc,   &
       dz,   &
       kref_in )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out)   :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(out)   :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)    :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]
    real(RP), intent(in)    :: dz  (KA,IA,JA) !< distance between the layer (center) [m]
    integer, intent(in), optional :: kref_in

    real(RP) :: Rtot  (KA,IA,JA)
    real(RP) :: CVtot (KA,IA,JA)
    real(RP) :: CPovCV(KA,IA,JA)

    real(RP) :: RovCV
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: kref
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if ( present(kref_in) ) then
       kref = kref_in
    else
       kref = KE
    end if

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, kref
       Rtot  (k,i,j) = Rdry  * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                     + Rvap  * qv(k,i,j)
       CVtot (k,i,j) = CVdry * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                     + CV_qv * qv(k,i,j)                          &
                     + CV_qc * qc(k,i,j)
       CPovCV(k,i,j) = ( CVtot(k,i,j) + Rtot(k,i,j) ) / CVtot(k,i,j)
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = kref-1, KS, -1
       RovCV = Rtot(k,i,j) / CVtot(k,i,j)

       dens_s      = 0.0_RP
       dens(k,i,j) = dens(k+1,i,j) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(k,i,j)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(k,i,j)

          dhyd = + ( P00 * ( dens_s        * Rtot(k  ,i,j) * pott(k  ,i,j) / P00 )**CPovCV(k  ,i,j) &
                   - P00 * ( dens(k+1,i,j) * Rtot(k+1,i,j) * pott(k+1,i,j) / P00 )**CPovCV(k+1,i,j) ) / dz(k+1,i,j) & ! dpdz
                 - GRAV * 0.5_RP * ( dens_s + dens(k+1,i,j) )                                                         ! rho*g

          dgrd = + P00 * ( Rtot(k,i,j) * pott(k,i,j) / P00 )**CPovCV(k,i,j) / dz(k+1,i,j) &
                 * CPovCV(k,i,j) * dens_s**RovCV                                          &
                 - 0.5_RP * GRAV

          dens(k,i,j) = dens_s - dhyd/dgrd

          if( dens(k,i,j)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [buildrho 3D rev atmos] iteration not converged!', &
                     k,i,j,dens(k,i,j),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, kref
       pres(k,i,j) = P00 * ( dens(k,i,j) * Rtot(k,i,j) * pott(k,i,j) / P00 )**CPovCV(k,i,j)
       temp(k,i,j) = pres(k,i,j) / ( dens(k,i,j) * Rtot(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_rev_3D

  !-----------------------------------------------------------------------------
  !> Build up density from surface (1D)
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_1D( &
       dens,     &
       pott,     &
       pres,     &
       temp,     &
       qv,       &
       qc,       &
       pott_sfc, &
       pres_sfc, &
       temp_sfc, &
       qv_sfc,   &
       qc_sfc    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens(KA) !< density               [kg/m3]
    real(RP), intent(out) :: pott(KA) !< potential temperature [K]
    real(RP), intent(out) :: pres(KA) !< pressure              [Pa]
    real(RP), intent(in)  :: temp(KA) !< temperature           [K]
    real(RP), intent(in)  :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)  :: qc  (KA) !< liquid water          [kg/kg]

    real(RP), intent(out) :: pott_sfc !< surface potential temperature [K]
    real(RP), intent(in)  :: pres_sfc !< surface pressure              [Pa]
    real(RP), intent(in)  :: temp_sfc !< surface temperature           [K]
    real(RP), intent(in)  :: qv_sfc   !< surface water vapor           [kg/kg]
    real(RP), intent(in)  :: qc_sfc   !< surface liquid water          [kg/kg]

    real(RP) :: dens_sfc

    real(RP) :: Rtot_sfc
    real(RP) :: CVtot_sfc
    real(RP) :: Rtot
    real(RP) :: CVtot

    real(RP) :: RovCP_sfc
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere

    Rtot_sfc   = Rdry  * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + Rvap  * qv_sfc
    CVtot_sfc  = CVdry * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + CV_qv * qv_sfc                       &
               + CV_qc * qc_sfc

    Rtot   = Rdry  * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + Rvap  * qv(KS)
    CVtot  = CVdry * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + CV_qv * qv(KS)                       &
           + CV_qc * qc(KS)

    ! density at surface
    RovCP_sfc = Rtot_sfc / ( CVtot_sfc + Rtot_sfc )
    dens_sfc  = pres_sfc / ( Rtot_sfc * temp_sfc )
    pott_sfc  = temp_sfc * ( P00/pres_sfc )**RovCP_sfc

    ! make density at lowermost cell center
    dens_s   = 0.0_RP
    dens(KS) = dens_sfc ! first guess

    converged = .false.
    do ite = 1, itelim
       if ( abs(dens(KS)-dens_s) <= criteria ) then
          converged = .true.
          exit
       endif

       dens_s = dens(KS)

       dhyd = + ( dens_sfc * Rtot_sfc * temp_sfc &
                - dens_s   * Rtot     * temp(KS) ) / GRID_CZ(KS) & ! dp/dz
              - GRAV * 0.5_RP * ( dens_sfc + dens_s )         ! rho*g

       dgrd = - Rtot * temp(KS) / GRID_CZ(KS) &
              - 0.5_RP * GRAV

       dens(KS) = dens_s - dhyd/dgrd

       if( dens(KS)*0.0_RP /= 0.0_RP) exit
    enddo

    if ( .NOT. converged ) then
       write(*,*) 'xxx [buildrho bytemp 1D sfc] iteration not converged!', &
                  dens(KS),ite,dens_s,dhyd,dgrd
       call PRC_MPIstop
    endif

    !--- from lowermost atmosphere to top of atmosphere
    call ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D( dens(:), & ! [INOUT]
                                                     pott(:), & ! [OUT]
                                                     pres(:), & ! [OUT]
                                                     temp(:), & ! [IN]
                                                     qv  (:), & ! [IN]
                                                     qc  (:)  ) ! [IN]

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_1D

  !-----------------------------------------------------------------------------
  !> Build up density from surface (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_3D( &
       dens,     &
       pott,     &
       pres,     &
       temp,     &
       qv,       &
       qc,       &
       pott_sfc, &
       pres_sfc, &
       temp_sfc, &
       qv_sfc,   &
       qc_sfc    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out) :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(out) :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(in)  :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)  :: qc  (KA,IA,JA) !< liquid water          [kg/kg]

    real(RP), intent(out) :: pott_sfc(1,IA,JA) !< surface potential temperature [K]
    real(RP), intent(in)  :: pres_sfc(1,IA,JA) !< surface pressure              [Pa]
    real(RP), intent(in)  :: temp_sfc(1,IA,JA) !< surface temperature           [K]
    real(RP), intent(in)  :: qv_sfc  (1,IA,JA) !< surface water vapor           [kg/kg]
    real(RP), intent(in)  :: qc_sfc  (1,IA,JA) !< surface liquid water          [kg/kg]

    real(RP) :: dens_sfc  (1,IA,JA)

    real(RP) :: Rtot_sfc  (IA,JA)
    real(RP) :: CVtot_sfc (IA,JA)
    real(RP) :: Rtot      (IA,JA)
    real(RP) :: CVtot     (IA,JA)

    real(RP) :: RovCP_sfc
    real(RP) :: DZ
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: i, j
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere

    do j = JSB, JEB
    do i = ISB, IEB
       Rtot_sfc (i,j) = Rdry  * ( 1.0_RP - qv_sfc(1,i,j) - qc_sfc(1,i,j) ) &
                      + Rvap  * qv_sfc(1,i,j)
       CVtot_sfc(i,j) = CVdry * ( 1.0_RP - qv_sfc(1,i,j) - qc_sfc(1,i,j) ) &
                      + CV_qv * qv_sfc(1,i,j)                              &
                      + CV_qc * qc_sfc(1,i,j)
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       Rtot (i,j) = Rdry  * ( 1.0_RP - qv(KS,i,j) - qc(KS,i,j) ) &
                  + Rvap  * qv(KS,i,j)
       CVtot(i,j) = CVdry * ( 1.0_RP - qv(KS,i,j) - qc(KS,i,j) ) &
                  + CV_qv * qv(KS,i,j)                           &
                  + CV_qc * qc(KS,i,j)
    enddo
    enddo

    ! density at surface
    do j = JSB, JEB
    do i = ISB, IEB
       RovCP_sfc       = Rtot_sfc(i,j) / ( CVtot_sfc(i,j) + Rtot_sfc(i,j) )
       dens_sfc(1,i,j) = pres_sfc(1,i,j) / ( Rtot_sfc(i,j) * temp_sfc(1,i,j) )
       pott_sfc(1,i,j) = temp_sfc(1,i,j) / ( P00/pres_sfc(1,i,j) )**RovCP_sfc
    enddo
    enddo

    ! make density at lowermost cell center
    do j = JSB, JEB
    do i = ISB, IEB
       DZ = REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j)

       dens_s       = 0.0_RP
       dens(KS,i,j) = dens_sfc(1,i,j) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(KS,i,j)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(KS,i,j)

          dhyd = + ( dens_sfc(1,i,j) * Rtot_sfc(i,j) * temp_sfc(1,i,j) &
                   - dens_s          * Rtot    (i,j) * temp   (KS,i,j) ) / DZ & ! dp/dz
                 - GRAV * 0.5_RP * ( dens_sfc(1,i,j) + dens_s )                 ! rho*g

          dgrd = - Rtot(i,j) * temp(KS,i,j) / DZ &
                 - 0.5_RP * GRAV

          dens(KS,i,j) = dens_s - dhyd/dgrd

          if( dens(KS,i,j)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [buildrho bytemp 3D sfc] iteration not converged!', &
                     i,j,dens(KS,i,j),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo
    enddo

    !--- from lowermost atmosphere to top of atmosphere
    call ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D( dens(:,:,:), & ! [INOUT]
                                                     pott(:,:,:), & ! [OUT]
                                                     pres(:,:,:), & ! [OUT]
                                                     temp(:,:,:), & ! [IN]
                                                     qv  (:,:,:), & ! [IN]
                                                     qc  (:,:,:)  ) ! [IN]

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_3D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (1D)
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D( &
       dens, &
       pott, &
       pres, &
       temp, &
       qv,   &
       qc    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: dens(KA) !< density               [kg/m3]
    real(RP), intent(out)   :: pott(KA) !< potential temperature [K]
    real(RP), intent(out)   :: pres(KA) !< pressure              [Pa]
    real(RP), intent(in)    :: temp(KA) !< temperature           [K]
    real(RP), intent(in)    :: qv  (KA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA) !< liquid water          [kg/kg]

    real(RP) :: Rtot  (KA)
    real(RP) :: CVtot (KA)

    real(RP) :: RovCP
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       Rtot (k) = Rdry  * ( 1.0_RP - qv(k) - qc(k) ) &
                + Rvap  * qv(k)
       CVtot(k) = CVdry * ( 1.0_RP - qv(k) - qc(k) ) &
                + CV_qv * qv(k)                      &
                + CV_qc * qc(k)
    enddo

    do k = KS+1, KE

       dens_s  = 0.0_RP
       dens(k) = dens(k-1) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(k)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(k)

          dhyd = + ( dens(k-1) * Rtot(k-1) * temp(k-1)  &
                   - dens_s    * Rtot(k  ) * temp(k  ) ) / GRID_FDZ(k-1) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s )             ! rho*g

          dgrd = - Rtot(k) * temp(k) / GRID_FDZ(k-1) &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

          if( dens(k)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [buildrho bytemp 1D atmos] iteration not converged!', &
                     k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo

    do k = KS, KE
       RovCP   = Rtot(k) / ( CVtot(k) + Rtot(k) )
       pres(k) = dens(k) * Rtot(k) * temp(k)
       pott(k) = temp(k) * ( P00 / pres(k) )**RovCP
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_1D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D( &
       dens, &
       pott, &
       pres, &
       temp, &
       qv,   &
       qc    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out)   :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(out)   :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)    :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]

    real(RP) :: Rtot  (KA,IA,JA)
    real(RP) :: CVtot (KA,IA,JA)

    real(RP) :: RovCP
    real(RP) :: DZ
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       Rtot (k,i,j) = Rdry  * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                    + Rvap  * qv(k,i,j)
       CVtot(k,i,j) = CVdry * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                    + CV_qv * qv(k,i,j)                          &
                    + CV_qc * qc(k,i,j)
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS+1, KE
       DZ = REAL_CZ(k,i,j) - REAL_CZ(k-1,i,j)

       dens_s      = 0.0_RP
       dens(k,i,j) = dens(k-1,i,j) ! first guess

       converged = .false.
       do ite = 1, itelim
          if ( abs(dens(k,i,j)-dens_s) <= criteria ) then
             converged = .true.
             exit
          endif

          dens_s = dens(k,i,j)

          dhyd = + ( dens(k-1,i,j) * Rtot(k-1,i,j) * temp(k-1,i,j) &
                   - dens_s        * Rtot(k  ,i,j) * temp(k  ,i,j) ) / DZ & ! dpdz
                 - GRAV * 0.5_RP * ( dens(k-1,i,j) + dens_s )               ! rho*g

          dgrd = - Rtot(k,i,j) * temp(k,i,j) / DZ &
                 - 0.5_RP * GRAV

          dens(k,i,j) = dens_s - dhyd/dgrd

          if( dens(k,i,j)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [buildrho bytemp 3D atmos] iteration not converged!', &
                     k,i,j,dens(k,i,j),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       RovCP   = Rtot(k,i,j) / ( CVtot(k,i,j) + Rtot(k,i,j) )
       pres(k,i,j) = dens(k,i,j) * Rtot(k,i,j) * temp(k,i,j)
       pott(k,i,j) = temp(k,i,j) * ( P00 / pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D

  !-----------------------------------------------------------------------------
  !> Calculate mean sea-level pressure from barometric law (0D)
  subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_0D( &
       mslp, &
       pres, &
       temp, &
       dz    )
    implicit none

    real(RP), intent(out) :: mslp !< mean sea-level pressure [Pa]
    real(RP), intent(in)  :: pres !< surface pressure        [Pa]
    real(RP), intent(in)  :: temp !< surface air temperature [K]
    real(RP), intent(in)  :: dz   !< surface height from MSL [m]

    ! work
    real(RP) :: TM
    !---------------------------------------------------------------------------

    TM = temp + LAPS * dz * 0.5_RP ! column-mean air temperature

    ! barometric law
    mslp = pres * exp( GRAV * dz / ( Rdry * TM ) )

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_0D

  !-----------------------------------------------------------------------------
  !> Calculate mean sea-level pressure from barometric law (2D)
  subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_2D( &
       mslp, &
       pres, &
       temp, &
       dz    )
    implicit none

    real(RP), intent(out) :: mslp(IA,JA) !< mean sea-level pressure [Pa]
    real(RP), intent(in)  :: pres(IA,JA) !< surface pressure        [Pa]
    real(RP), intent(in)  :: temp(IA,JA) !< surface air temperature [K]
    real(RP), intent(in)  :: dz  (IA,JA) !< surface height from MSL [m]

    ! work
    real(RP) :: TM

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       TM = temp(i,j) + LAPS * dz(i,j) * 0.5_RP ! column-mean air temperature

       ! barometric law
       mslp(i,j) = pres(i,j) * exp( GRAV * dz(i,j) / ( Rdry * TM ) )
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_mslp_2D

  !-----------------------------------------------------------------------------
  !> Calculate surface pressure from barometric law (0D)
  subroutine ATMOS_HYDROSTATIC_barometric_law_pres_0D( &
       pres, &
       mslp, &
       temp, &
       dz    )
    implicit none

    real(RP), intent(out) :: pres !< surface pressure        [Pa]
    real(RP), intent(in)  :: mslp !< mean sea-level pressure [Pa]
    real(RP), intent(in)  :: temp !< surface air temperature [K]
    real(RP), intent(in)  :: dz   !< surface height from MSL [m]

    ! work
    real(RP) :: TM
    !---------------------------------------------------------------------------

    TM = temp + LAPS * dz * 0.5_RP ! column-mean air temperature

    ! barometric law
    pres = mslp / exp( GRAV * dz / ( Rdry * TM ) )

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_pres_0D

  !-----------------------------------------------------------------------------
  !> Calculate surface pressure from barometric law (2D)
  subroutine ATMOS_HYDROSTATIC_barometric_law_pres_2D( &
       pres, &
       mslp, &
       temp, &
       dz    )
    implicit none

    real(RP), intent(out) :: pres(IA,JA) !< surface pressure        [Pa]
    real(RP), intent(in)  :: mslp(IA,JA) !< mean sea-level pressure [Pa]
    real(RP), intent(in)  :: temp(IA,JA) !< surface air temperature [K]
    real(RP), intent(in)  :: dz  (IA,JA) !< surface height from MSL [m]

    ! work
    real(RP) :: TM

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       TM = temp(i,j) + LAPS * dz(i,j) * 0.5_RP ! column-mean air temperature

       ! barometric law
       pres(i,j) = mslp(i,j) / exp( GRAV * dz(i,j) / ( Rdry * TM ) )
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_barometric_law_pres_2D

end module scale_atmos_hydrostatic
!-------------------------------------------------------------------------------
