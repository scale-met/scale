!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_mkinit
  public :: USER_calc_tendency
  public :: USER_update

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
  logical, private :: USER_do = .false. !< do user step?
  real(RP), private :: position_phi_x  = 1.E+3_RP   ![m]
  real(RP), private :: position_phi_y  = 1.E+3_RP   ![m]
  real(RP), private :: position_phi_z  = 1.E+3_RP   ![m]
  real(RP), private :: radius_phi_x = 5.E+1_RP   ![m]
  real(RP), private :: radius_phi_y = 5.E+1_RP   ![m]
  real(RP), private :: radius_phi_z = 5.E+1_RP   ![m]
  real(RP), private :: qvalue = 0.0_RP   ![m]
  !--- Indeces for determining species of cloud particle
  integer, parameter, private :: LTS = 1
  !--- For history output
  integer, private, parameter :: w_nmax = 6
  real(RP), private, allocatable :: w3d(:,:,:,:)
  integer,  private              :: HIST_id(w_nmax)
  character(len=H_SHORT), private :: w_name(w_nmax)
  character(len=H_MID),   private :: w_longname(w_nmax)
  character(len=H_SHORT), private :: w_unit(w_nmax)
  data w_name / 'CRGD_TOT_ana', &
                'Ex_ana', &
                'Ey_ana', &
                'Ez_ana', &
                'Eabs_ana', &
                'Epot_ana' /
  data w_longname / &
                'Charge density', &
                'X component of Electrical Field', &
                'Y component of Electrical Field', &
                'Z component of Electrical Field', &
                'Absolute value of Electrical Field', &
                'Electric Potential' /
  data w_unit / 'nC/m3', &
                'kV/m', &
                'kV/m', &
                'kV/m', &
                'kV/m', &
                'V' /
  integer, parameter :: I_QCRG_C = 1
  integer, parameter :: I_Ex     = 2
  integer, parameter :: I_Ey     = 3
  integer, parameter :: I_Ez     = 4
  integer, parameter :: I_Eabs   = 5
  integer, parameter :: I_Epot   = 6
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config before setup of tracers
  subroutine USER_tracer_setup
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    ! if you want to add tracers, call the TRACER_regist subroutine.
    ! e.g.,
!    integer, parameter     :: NQ = 1
!    integer                :: QS
!    character(len=H_SHORT) :: NAME(NQ)
!    character(len=H_MID)   :: DESC(NQ)
!    character(len=H_SHORT) :: UNIT(NQ)
!
!    data NAME (/ 'name' /)
!    data DESC (/ 'tracer name' /)
!    data UNIT (/ 'kg/kg' /)
    !---------------------------------------------------------------------------

!    call TRACER_regist( QS,   & ! [OUT]
!                        NQ,   & ! [IN]
!                        NAME, & ! [IN]
!                        DESC, & ! [IN]
!                        UNIT  ) ! [IN]

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PI    => CONST_PI
    use scale_atmos_grid_cartesC, only: &
       CDX   => ATMOS_GRID_CARTESC_CDX, &
       CDY   => ATMOS_GRID_CARTESC_CDY, &
       CX    => ATMOS_GRID_CARTESC_CX, &
       CY    => ATMOS_GRID_CARTESC_CY, &
       CZ    => ATMOS_GRID_CARTESC_CZ
    use scale_atmos_phy_lt_sato2019, only: &
       ATMOS_PHY_LT_sato2019_setup
    use scale_file_history, only: &
       FILE_HISTORY_reg
    use mod_atmos_admin, only: &
       MP_TYPE => ATMOS_PHY_MP_TYPE
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
       position_phi_x, &
       position_phi_y, &
       position_phi_z, &
       radius_phi_x, &
       radius_phi_y, &
       radius_phi_z, &
       qvalue

    integer :: ierr, ip
    integer :: dum = 2
    !---------------------------------------------------------------------------

    if( MP_TYPE /= 'TOMITA08' ) then
       LOG_ERROR("USER_setup",*) 'MP_TYPE should be TOMITA08 for this test. Check!'
       call PRC_abort
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    if(  radius_phi_x /= radius_phi_y .or. &
         radius_phi_y /= radius_phi_z .or. &
         radius_phi_x /= radius_phi_z      ) then
       LOG_ERROR("USER_setup",*) 'xxx radius_phi_x, _y, and _z should be same, stop!'
       call PRC_abort
    endif

    !-- for history output
    allocate( w3d(KA,IA,JA,w_nmax) )
    w3d(:,:,:,:) = 0.0_RP

    do ip = 1, w_nmax
       call FILE_HISTORY_reg( w_name(ip), w_longname(ip), w_unit(ip), & ! [IN]
                              HIST_id(ip)                             ) ! [OUT]
    end do

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    implicit none
    !---------------------------------------------------------------------------
    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calculation tendency
  subroutine USER_calc_tendency
    implicit none

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_update
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       EPSair => CONST_EPSair, &
       EPSvac => CONST_EPSvac
    use scale_atmos_grid_cartesC, only: &
       FDX    => ATMOS_GRID_CARTESC_FDX, &
       FDY    => ATMOS_GRID_CARTESC_FDY, &
       FDZ    => ATMOS_GRID_CARTESC_FDZ, &
       CX    => ATMOS_GRID_CARTESC_CX, &
       CY    => ATMOS_GRID_CARTESC_CY, &
       CZ    => ATMOS_GRID_CARTESC_CZ
    use scale_atmos_phy_lt_sato2019, only: &
       ATMOS_PHY_LT_sato2019_tendency
    use scale_file_history, only: &
       FILE_HISTORY_in, &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    use mod_atmos_admin, only: &
       MP_TYPE => ATMOS_PHY_MP_TYPE
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, QS_MP, QE_MP
    use mod_atmos_phy_lt_vars, only: &
       QA_LT, QS_LT, QE_LT
    use mod_atmos_vars, only: &
       DENS, RHOT
    implicit none

    integer  :: k, i, j, m, n, ip
    real(RP) :: distance(3), point
    real(RP) :: dist(2)
    real(RP) :: Efield(KA,IA,JA,4), E_pot(KA,IA,JA), QCRG_out(KA,IA,JA)
    real(RP) :: E_old(KA,IA,JA)
    logical  :: HIST_sw(w_nmax), hist_flag
    real(RP), allocatable :: QCRG_local(:,:,:,:)
    real(RP), allocatable :: QTRC_local(:,:,:,:)
    real(RP), allocatable :: dummy_mp(:,:,:,:)
    real(RP), allocatable :: dummy_lt(:,:,:,:)
    real(RP), allocatable :: dummy_lt2(:,:,:,:)
    real(RP), allocatable :: dummy_sarea(:,:,:,:)
    real(RP), allocatable :: dummy_splt(:,:,:,:)
    !---------------------------------------------------------------------------

    if ( USER_do ) then

       LOG_NEWLINE
       LOG_INFO("USER_update",*) 'USER update'

       allocate( QTRC_local(KA,IA,JA,QS_MP:QE_MP) )
       allocate( QCRG_local(KA,IA,JA,QS_LT:QE_LT) )
       QTRC_local(:,:,:,:) = 0.0_RP
       QCRG_local(:,:,:,:) = 0.0_RP
       allocate( dummy_mp(KA,IA,JA,QS_MP:QE_MP) )
       allocate( dummy_lt(KA,IA,JA,QS_LT:QE_LT) )
       allocate( dummy_lt2(KA,IA,JA,QS_LT:QE_LT) )
       allocate( dummy_splt(KA,IA,JA,3) )
       allocate( dummy_sarea(KA,IA,JA,QA_LT) )
       dummy_mp(:,:,:,:) = 0.0_RP
       dummy_lt(:,:,:,:) = 0.0_RP
       dummy_splt(:,:,:,:) = 0.0_RP
       dummy_sarea(:,:,:,:) = 1.0_RP

       hist_flag = .false.
       do ip = 1, w_nmax
          call FILE_HISTORY_query( HIST_id(ip), HIST_sw(ip) )
          hist_flag = hist_flag .or. HIST_sw(ip)
       end do

       !--- Electric field from analytical solution
       distance(1) = sqrt( (CZ(KA)-position_phi_z)**2 )
       E_pot(:,:,:) = 0.0_RP
       QCRG_out(:,:,:) = 0.0_RP
       do k = 1, KA-1
       do i = 1, IA-1
       do j = 1, JA-1
            distance(3) = sqrt( ( CZ(k)-position_phi_z )**2 &
                              + ( CX(i)-position_phi_x )**2 &
                              + ( CY(j)-position_phi_y )**2 )

            if( distance(3) < radius_phi_z ) then
             QCRG_out(k,i,j) = qvalue*1.E+9_RP*DENS(k,i,j)
             point = 1.0_RP
            else
             QCRG_OUT(k,i,j) = 0.0_RP  !--- outside
             point = 0.0_RP
            endif

            if( point == 0.0_RP ) then
             E_pot(k,i,j) = qvalue*radius_phi_z**3 / (3.0_RP*EPSvac*EPSair*distance(3))
            elseif( point == 1.0_RP ) then
             E_pot(k,i,j) = qvalue/(6.0_RP*EPSvac*EPSair)  &
                          * ( 3.0_RP*radius_phi_z**2-distance(3)**2 )
            endif
       enddo
       enddo
       enddo

       !---- Calculate Electrical Field
       Efield(:,:,:,:) = 0.d0
       do i = IS, IE
       do j = JS, JE
       do k = KS, KE
          Efield(k,i,j,1) = - ( E_pot(k  ,i+1,j  ) - E_pot(k  ,i-1,j  ) ) / FDX(i) * 0.5_RP
          Efield(k,i,j,2) = - ( E_pot(k  ,i  ,j+1) - E_pot(k  ,i  ,j-1) ) / FDY(j) * 0.5_RP
          Efield(k,i,j,3) = - ( E_pot(k+1,i  ,j  ) - E_pot(k-1,i  ,j  ) ) / FDZ(k) * 0.5_RP
          Efield(k,i,j,4) = sqrt( Efield(k,i,j,1)*Efield(k,i,j,1) &
                                + Efield(k,i,j,2)*Efield(k,i,j,2) &
                                + Efield(k,i,j,3)*Efield(k,i,j,3) )

       enddo
       enddo
       enddo

       if ( hist_flag ) then
          if ( HIST_sw(I_QCRG_C  ) ) then
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                w3d(k,i,j,I_QCRG_C) = QCRG_out(k,i,j)
             enddo
             enddo
             enddo
          endif
          if ( HIST_sw(I_Ex  ) ) then
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                w3d(k,i,j,I_Ex  ) = Efield(k,i,j,1  )*1.0E-3_RP ![kV/m]
             enddo
             enddo
             enddo
          endif
          if ( HIST_sw(I_Ey  ) ) then
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                w3d(k,i,j,I_Ey  ) = Efield(k,i,j,2  )*1.0E-3_RP ![kV/m]
             enddo
             enddo
             enddo
          endif
          if ( HIST_sw(I_Ez  ) ) then
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                w3d(k,i,j,I_Ez  ) = Efield(k,i,j,3 )*1.0E-3_RP ![kV/m]
             enddo
             enddo
             enddo
          endif
          if ( HIST_sw(I_Eabs) ) then
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                w3d(k,i,j,I_Eabs) = Efield(k,i,j,4)*1.0E-3_RP ![kV/m]
             enddo
             enddo
             enddo
          endif
          if ( HIST_sw(I_Epot) ) then
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                w3d(k,i,j,I_Epot) = E_pot(k,i,j)
             enddo
             enddo
             enddo
          endif
       end if

       do ip = 1, w_nmax
          if ( HIST_sw(ip) ) call FILE_HISTORY_put( HIST_id(ip), w3d(:,:,:,ip) )
       enddo


       !--- Electric field from Bi-CGSTAB
       QCRG_local(:,:,:,:) = 0.0_RP
       do j = 1, JA-1
       do i = 1, IA-1
       do k = 1, KA-1
         distance(3) = sqrt( ( CZ(k)-position_phi_z )**2 &
                           + ( CX(i)-position_phi_x )**2 &
                           + ( CY(j)-position_phi_y )**2 )

        if( distance(3) < radius_phi_z ) then
          QCRG_local(k,i,j,QS_LT) = qvalue * 1.0E+15_RP * DENS(k,i,j)  ![fC/m3]
        else
          QCRG_local(k,i,j,QS_LT) = 0.0_RP
        endif

       enddo
       enddo
       enddo

       !--- Electric field from QTRC
       Efield(:,:,:,:) = 0.d0
       E_old(:,:,:) = 0.d0
       call ATMOS_PHY_LT_sato2019_tendency( &
            KA, KS, KE, IA, IS, IE, JA, JS, JE, KIJMAX, IMAX, JMAX, MP_TYPE,           & ! [IN]
            QA_MP, QA_LT, DENS(:,:,:),                                                 & ! [IN]
            RHOT(:,:,:), QTRC_local(:,:,:,QS_MP:QE_MP),                                & ! [IN]
            QCRG_local(:,:,:,QS_LT:QE_LT), 1.0_DP,                                     & ! [IN]
            1.0_DP, dummy_splt(:,:,:,:), dummy_sarea(:,:,:,:),                         & ! [IN]
            dummy_mp(:,:,:,QS_MP:QE_MP), dummy_lt(:,:,:,QS_LT:QE_LT),                  & ! [IN]
            E_old(:,:,:), dummy_lt2(:,:,:,QS_LT:QE_LT)                                 ) ! [INOUT]

    endif

    return
  end subroutine USER_update

end module mod_user
