!-------------------------------------------------------------------------------
!> module hydrostatic
!!
!! @par Description
!!          make hydrostatic profile in the model
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-02-20 (H.Yashiro)  [new] Extract from the tool of Y.Miyamoto
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_hydrostatic
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_FID_CONF, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_hydro_buildrho
  public :: ATMOS_hydro_buildrho_fromKS
  public :: ATMOS_hydro_buildrho_1d
  public :: ATMOS_hydro_buildrho_temp
  public :: ATMOS_hydro_buildrho_temp_1d

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: buildrho
  private :: buildrho_temp

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical :: HYDROSTATIC_userapserate = .true.

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Buildup density from surface
  !-----------------------------------------------------------------------------
  subroutine ATMOS_hydro_buildrho( &
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
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens(KA,IA,JA) !< density [kg/m3]
    real(RP), intent(out) :: temp(KA,IA,JA) !< temperature [K]
    real(RP), intent(out) :: pres(KA,IA,JA) !< pressure [Pa]
    real(RP), intent(in)  :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA,IA,JA) !< water vapor [kg/kg]
    real(RP), intent(in)  :: qc  (KA,IA,JA) !< water vapor [kg/kg]

    real(RP), intent(out) :: temp_sfc(1,IA,JA) !< surface temperature [K]
    real(RP), intent(in)  :: pres_sfc(1,IA,JA) !< surface pressure [Pa]
    real(RP), intent(in)  :: pott_sfc(1,IA,JA) !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc  (1,IA,JA) !< surface water vapor [kg/kg]
    real(RP), intent(in)  :: qc_sfc  (1,IA,JA) !< surface water vapor [kg/kg]

    NAMELIST / PARAM_ATMOS_HYDROSTATIC / &
       HYDROSTATIC_userapserate

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[HYDROSTATIC]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_HYDROSTATIC,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used!'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_HYDROSTATIC. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_HYDROSTATIC)

    do j = JS, JE
    do i = IS, IE
       call buildrho( dens(:,i,j),             & ! [OUT]
                      temp(:,i,j),             & ! [OUT]
                      pres(:,i,j),             & ! [OUT]
                      pott(:,i,j),             & ! [IN]
                      qv  (:,i,j),             & ! [IN]
                      qc  (:,i,j),             & ! [IN]
                      temp_sfc(1,i,j),         & ! [OUT]
                      pres_sfc(1,i,j),         & ! [IN]
                      pott_sfc(1,i,j),         & ! [IN]
                      qv_sfc  (1,i,j),         & ! [IN]
                      qc_sfc  (1,i,j),         & ! [IN]
                      HYDROSTATIC_userapserate ) ! [IN]
    enddo
    enddo

    ! fill IHALO & JHALO
    call COMM_vars8( dens(:,:,:), 1 )
    call COMM_wait ( dens(:,:,:), 1 )

    return
  end subroutine ATMOS_hydro_buildrho

  !-----------------------------------------------------------------------------
  subroutine ATMOS_hydro_buildrho_1d( &
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
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens(KA) !< density [kg/m3]
    real(RP), intent(out) :: temp(KA) !< temperature [K]
    real(RP), intent(out) :: pres(KA) !< pressure [Pa]
    real(RP), intent(in)  :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA) !< water vapor [kg/kg]
    real(RP), intent(in)  :: qc  (KA) !< water vapor [kg/kg]

    real(RP), intent(out) :: temp_sfc !< surface temperature [K]
    real(RP), intent(in)  :: pres_sfc !< surface pressure [Pa]
    real(RP), intent(in)  :: pott_sfc !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc   !< surface water vapor [kg/kg]
    real(RP), intent(in)  :: qc_sfc   !< surface water vapor [kg/kg]

    NAMELIST / PARAM_ATMOS_HYDROSTATIC / &
         HYDROSTATIC_userapserate

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[HYDROSTATIC_1d]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_HYDROSTATIC,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used!'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_HYDROSTATIC. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_HYDROSTATIC)

    call buildrho( dens(:),                 & ! [OUT]
                   temp(:),                 & ! [OUT]
                   pres(:),                 & ! [OUT]
                   pott(:),                 & ! [IN]
                   qv  (:),                 & ! [IN]
                   qc  (:),                 & ! [IN]
                   temp_sfc,                & ! [OUT]
                   pres_sfc,                & ! [IN]
                   pott_sfc,                & ! [IN]
                   qv_sfc,                  & ! [IN]
                   qc_sfc,                  & ! [IN]
                   HYDROSTATIC_userapserate ) ! [IN]

    return
  end subroutine ATMOS_hydro_buildrho_1d

  !-----------------------------------------------------------------------------
  subroutine ATMOS_hydro_buildrho_fromKS( &
       dens,     &
       temp,     &
       pres,     &
       pott,     &
       qv,       &
       qc        )
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: dens(KA,IA,JA) !< density [kg/m3]
    real(RP), intent(out)   :: temp(KA,IA,JA) !< temperature [K]
    real(RP), intent(out)   :: pres(KA,IA,JA) !< pressure [Pa]
    real(RP), intent(in)    :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< water vapor [kg/kg]

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[HYDROSTATIC]/Categ[ATMOS]'

    do j = JS, JE
    do i = IS, IE
       call buildrho_fromKS( dens(:,i,j), & ! [INOUT]
                             temp(:,i,j), & ! [OUT]
                             pres(:,i,j), & ! [OUT]
                             pott(:,i,j), & ! [IN]
                             qv  (:,i,j), & ! [IN]
                             qc  (:,i,j)  ) ! [IN]
    enddo
    enddo

    ! fill IHALO & JHALO
    call COMM_vars8( dens(:,:,:), 1 )
    call COMM_wait ( dens(:,:,:), 1 )

    return
  end subroutine ATMOS_hydro_buildrho_fromKS

  !-----------------------------------------------------------------------------
  !> Buildup density from surface with temperature
  !-----------------------------------------------------------------------------
  subroutine ATMOS_hydro_buildrho_temp( &
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
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out) :: dens(KA,IA,JA) !< density [kg/m3]
    real(RP), intent(out) :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(out) :: pres(KA,IA,JA) !< pressure [Pa]
    real(RP), intent(in)  :: temp(KA,IA,JA) !< temperature [K]
    real(RP), intent(in)  :: qv  (KA,IA,JA) !< water vapor [kg/kg]
    real(RP), intent(in)  :: qc  (KA,IA,JA) !< water vapor [kg/kg]

    real(RP), intent(out) :: pott_sfc(1,IA,JA) !< surface potential temperature [K]
    real(RP), intent(in)  :: pres_sfc(1,IA,JA) !< surface pressure [Pa]
    real(RP), intent(in)  :: temp_sfc(1,IA,JA) !< surface temperature [K]
    real(RP), intent(in)  :: qv_sfc  (1,IA,JA) !< surface water vapor [kg/kg]
    real(RP), intent(in)  :: qc_sfc  (1,IA,JA) !< surface water vapor [kg/kg]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       call buildrho_temp( dens(:,i,j),     & ! [OUT]
                           pott(:,i,j),     & ! [OUT]
                           pres(:,i,j),     & ! [OUT]
                           temp(:,i,j),     & ! [IN]
                           qv  (:,i,j),     & ! [IN]
                           qc  (:,i,j),     & ! [IN]
                           pott_sfc(1,i,j), & ! [OUT]
                           pres_sfc(1,i,j), & ! [IN]
                           temp_sfc(1,i,j), & ! [IN]
                           qv_sfc  (1,i,j), & ! [IN]
                           qc_sfc  (1,i,j)  ) ! [IN]
    enddo
    enddo

    ! fill IHALO & JHALO
    call COMM_vars8( dens(:,:,:), 1 )
    call COMM_wait ( dens(:,:,:), 1 )

    return
  end subroutine ATMOS_hydro_buildrho_temp

  !-----------------------------------------------------------------------------
  subroutine ATMOS_hydro_buildrho_temp_1d( &
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
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: dens(KA) !< density [kg/m3]
    real(RP), intent(out) :: pott(KA) !< potential temperature [K]
    real(RP), intent(out) :: pres(KA) !< pressure [Pa]
    real(RP), intent(in)  :: temp(KA) !< temperature [K]
    real(RP), intent(in)  :: qv  (KA) !< water vapor [kg/kg]
    real(RP), intent(in)  :: qc  (KA) !< water vapor [kg/kg]

    real(RP), intent(out) :: pott_sfc !< surface potential temperature [K]
    real(RP), intent(in)  :: pres_sfc !< surface pressure [Pa]
    real(RP), intent(in)  :: temp_sfc !< surface temperature [K]
    real(RP), intent(in)  :: qv_sfc   !< surface water vapor [kg/kg]
    real(RP), intent(in)  :: qc_sfc   !< surface water vapor [kg/kg]
    !---------------------------------------------------------------------------

    call buildrho_temp( dens(:),  & ! [OUT]
                        pott(:),  & ! [OUT]
                        pres(:),  & ! [OUT]
                        temp(:),  & ! [IN]
                        qv  (:),  & ! [IN]
                        qc  (:),  & ! [IN]
                        pott_sfc, & ! [OUT]
                        pres_sfc, & ! [IN]
                        temp_sfc, & ! [IN]
                        qv_sfc  , & ! [IN]
                        qc_sfc    ) ! [IN]

    return
  end subroutine ATMOS_hydro_buildrho_temp_1d

  !-----------------------------------------------------------------------------
  !> Buildup density from surface
  !-----------------------------------------------------------------------------
  subroutine buildrho( &
       dens,        &
       temp,        &
       pres,        &
       pott,        &
       qv,          &
       qc,          &
       temp_sfc,    &
       pres_sfc,    &
       pott_sfc,    &
       qc_sfc,      &
       qv_sfc,      &
       userapserate )
    use mod_const, only: &
       GRAV    => CONST_GRAV,    &
       EPS     => CONST_EPS,     &
       Rdry    => CONST_Rdry,    &
       Rvap    => CONST_Rvap,    &
       CVdry   => CONST_CVdry,   &
       CVvap   => CONST_CVvap,   &
       CL      => CONST_CL,      &
       LASPdry => CONST_LASPdry, &
       P00     => CONST_PRE00
    use mod_grid, only: &
       CZ  => GRID_CZ
    implicit none

    real(RP), intent(out) :: dens(KA) !< density [kg/m3]
    real(RP), intent(out) :: temp(KA) !< temperature [K]
    real(RP), intent(out) :: pres(KA) !< pressure [Pa]
    real(RP), intent(in)  :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)  :: qv  (KA) !< water vapor [kg/kg]
    real(RP), intent(in)  :: qc  (KA) !< water vapor [kg/kg]

    real(RP), intent(out) :: temp_sfc !< surface temperature [K]
    real(RP), intent(in)  :: pres_sfc !< surface pressure [Pa]
    real(RP), intent(in)  :: pott_sfc !< surface potential temperature [K]
    real(RP), intent(in)  :: qv_sfc   !< surface water vapor [kg/kg]
    real(RP), intent(in)  :: qc_sfc   !< surface water vapor [kg/kg]
    logical,  intent(in)  :: userapserate

    real(RP) :: dens_sfc
    real(RP) :: Rtot_sfc
    real(RP) :: CVtot_sfc
    real(RP) :: CPovCV_sfc

    real(RP) :: Rtot
    real(RP) :: CVtot
    real(RP) :: CPovCV
    real(RP) :: CPovR
    real(RP) :: CVovCP
    real(RP) :: RovCV

    real(RP) :: dens_s, dhyd, dgrd
    real(RP) :: criteria

    integer, parameter :: itelim = 100

    integer :: ite
    !---------------------------------------------------------------------------

    criteria = EPS * 5

    ! make density at surface
    Rtot_sfc   = Rdry  * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + Rvap  * qv_sfc
    CVtot_sfc  = CVdry * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + CVvap * qv_sfc                       &
               + CL    * qc_sfc

    CPovCV_sfc = ( CVtot_sfc + Rtot_sfc ) / CVtot_sfc
    CVovCP = CVtot_sfc / ( CVtot_sfc + Rtot_sfc )

    dens_sfc = P00 / Rtot_sfc / pott_sfc * ( pres_sfc/P00 )**CVovCP
    temp_sfc = pres_sfc / ( dens_sfc * Rtot_sfc )

    Rtot  = Rdry  * ( 1.0_RP - qv(KS) - qc(KS) ) &
          + Rvap  * qv(KS)
    CVtot = CVdry * ( 1.0_RP - qv(KS) - qc(KS) ) &
          + CVvap * qv(KS)                      &
          + CL    * qc(KS)
    CPovCV = ( CVtot + Rtot ) / CVtot

    ! make density at lowermost cell center
    if ( userapserate ) then
       CPovR  = ( CVtot + Rtot ) / Rtot
       CVovCP = CVtot / ( CVtot + Rtot )

       temp(KS) = pott_sfc - LASPdry * CZ(KS) ! use dry lapse rate
       pres(KS) = P00 * ( temp(KS)/pott(KS) )**CPovR
       dens(KS) = P00 / Rtot / pott(KS) * ( pres(KS)/P00 )**CVovCP

    else ! use itelation

       RovCV = Rtot / CVtot

       dens_s  = 0.0_RP
       dens(KS) = dens_sfc ! first guess

       do ite = 1, itelim
          if( abs(dens(KS)-dens_s) <= criteria ) exit

          dens_s = dens(KS)

          dhyd = + ( P00 * ( dens_sfc * Rtot_sfc * pott_sfc / P00 )**CPovCV_sfc &
                   - P00 * ( dens_s   * Rtot     * pott(KS) / P00 )**CPovCV  ) / CZ(KS) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens_sfc + dens_s )                                    ! rho*g

          dgrd = - P00 * ( Rtot * pott(KS) / P00 )**CPovCV / CZ(KS) &
                 * CPovCV * dens_s**RovCV                           &
                 - 0.5_RP * GRAV

          dens(KS) = dens_s - dhyd/dgrd

       enddo

       if ( ite > itelim ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', KS, ite, dens(KS), dens_s, dhyd, dgrd
       endif

    endif

    call buildrho_fromKS( dens, temp, pres, &
                          pott, qv, qc      )

    return
  end subroutine buildrho

  !-----------------------------------------------------------------------------
  !> Buildup density from surface
  !-----------------------------------------------------------------------------
  subroutine buildrho_temp( &
       dens,     &
       pott,     &
       pres,     &
       temp,     &
       qv,       &
       qc,       &
       pott_sfc, &
       pres_sfc, &
       temp_sfc, &
       qc_sfc,   &
       qv_sfc    )
    use mod_const, only: &
       GRAV    => CONST_GRAV,    &
       EPS     => CONST_EPS,     &
       Rdry    => CONST_Rdry,    &
       Rvap    => CONST_Rvap,    &
       CVdry   => CONST_CVdry,   &
       CVvap   => CONST_CVvap,   &
       CL      => CONST_CL,      &
       P00     => CONST_PRE00
    use mod_grid, only: &
       CZ  => GRID_CZ, &
       FDZ => GRID_FDZ
    implicit none

    real(RP), intent(out) :: dens(KA) !< density [kg/m3]
    real(RP), intent(out) :: pott(KA) !< potential temperature [K]
    real(RP), intent(out) :: pres(KA) !< pressure [Pa]
    real(RP), intent(in)  :: temp(KA) !< temperature [K]
    real(RP), intent(in)  :: qv  (KA) !< water vapor [kg/kg]
    real(RP), intent(in)  :: qc  (KA) !< water vapor [kg/kg]

    real(RP), intent(out) :: pott_sfc !< surface potential temperature [K]
    real(RP), intent(in)  :: pres_sfc !< surface pressure [Pa]
    real(RP), intent(in)  :: temp_sfc !< surface temperature [K]
    real(RP), intent(in)  :: qv_sfc   !< surface water vapor [kg/kg]
    real(RP), intent(in)  :: qc_sfc   !< surface water vapor [kg/kg]

    real(RP) :: dens_sfc
    real(RP) :: Rtot_sfc
    real(RP) :: CVtot_sfc
    real(RP) :: RovCP_sfc

    real(RP) :: Rtot(KA)
    real(RP) :: CVtot(KA)
    real(RP) :: RovCP(KA)

    real(RP) :: dens_s, dhyd, dgrd
    real(RP) :: criteria

    integer, parameter :: itelim = 100

    integer :: k, ite
    !---------------------------------------------------------------------------

    criteria = EPS * 5

    ! make density at surface
    Rtot_sfc   = Rdry  * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + Rvap  * qv_sfc
    CVtot_sfc  = CVdry * ( 1.0_RP - qv_sfc - qc_sfc ) &
               + CVvap * qv_sfc                       &
               + CL    * qc_sfc

    RovCP_sfc = Rtot_sfc / ( CVtot_sfc + Rtot_sfc )

    dens_sfc = pres_sfc / ( Rtot_sfc * temp_sfc )
    pott_sfc = temp_sfc * ( P00/pres_sfc )**RovCP_sfc

    do k = KS, KE
       Rtot(k)   = Rdry  * ( 1.0_RP - qv(k) - qc(k) ) &
                 + Rvap  * qv(k)
       CVtot(k)  = CVdry * ( 1.0_RP - qv(k) - qc(k) ) &
                 + CVvap * qv(k)                      &
                 + CL    * qc(k)

       RovCP(k) = Rtot(k) / ( CVtot(k) + Rtot(k) )
    enddo

    ! make density at lowermost cell center
    k = KS

    dens_s  = 0.0_RP
    dens(k) = dens_sfc ! first guess

    do ite = 1, itelim
       if( abs(dens(k)-dens_s) <= criteria ) exit

       dens_s = dens(k)

       dhyd = + ( dens_sfc * Rtot_sfc * temp_sfc &
                - dens_s   * Rtot(k)  * temp(k)  ) / CZ(k) & ! dp/dz
              - GRAV * 0.5_RP * ( dens_sfc + dens_s )        ! rho*g

       dgrd = - Rtot(k) * temp(k) / CZ(k) &
              - 0.5_RP * GRAV

       dens(k) = dens_s - dhyd/dgrd

    enddo

    if ( ite > itelim ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k), dens_s, dhyd, dgrd
    endif

    ! make density
    do k = KS+1, KE

       dens_s      = 0.0_RP
       dens(k) = dens(k-1)

       do ite = 1, itelim
          if( abs(dens(k)-dens_s) <= criteria ) exit

          dens_s = dens(k)

          dhyd = + ( dens(k-1) * Rtot(k-1) * temp(k-1)  &
                   - dens_s    * Rtot(k  ) * temp(k  ) ) / FDZ(k-1) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s )             ! rho*g

          dgrd = - Rtot(k) * temp(k) / FDZ(k-1) &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

       enddo

       if ( ite > itelim ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k), dens_s, dhyd, dgrd, Rtot, FDZ(k-1)
       endif
    enddo

    do k = KS, KE
       pres(k) = dens(k) * Rtot(k) * temp(k)
       pott(k) = temp(k) * ( P00 / pres(k) )**RovCP(k)
    enddo

    ! fill KHALO
    dens(   1:KS-1) = dens(KS)
    dens(KE+1:KA  ) = dens(KE)

    return
  end subroutine buildrho_temp

  subroutine buildrho_fromKS( &
       dens,     &
       temp,     &
       pres,     &
       pott,     &
       qv,       &
       qc        )
    use mod_const, only: &
       GRAV    => CONST_GRAV,    &
       EPS     => CONST_EPS,     &
       Rdry    => CONST_Rdry,    &
       Rvap    => CONST_Rvap,    &
       CVdry   => CONST_CVdry,   &
       CVvap   => CONST_CVvap,   &
       CL      => CONST_CL,      &
       P00     => CONST_PRE00
    use mod_grid, only: &
       FDZ => GRID_FDZ
    implicit none

    real(RP), intent(inout) :: dens(KA) !< density [kg/m3]
    real(RP), intent(out)   :: temp(KA) !< temperature [K]
    real(RP), intent(out)   :: pres(KA) !< pressure [Pa]
    real(RP), intent(in)    :: pott(KA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA) !< water vapor [kg/kg]
    real(RP), intent(in)    :: qc  (KA) !< water vapor [kg/kg]

    real(RP) :: Rtot(KA)
    real(RP) :: CVtot(KA)
    real(RP) :: CPovCV(KA)
    real(RP) :: RovCV

    real(RP) :: dens_s, dhyd, dgrd
    real(RP) :: criteria

    integer, parameter :: itelim = 100

    integer :: k, ite
    !---------------------------------------------------------------------------

    criteria = EPS * 5

    do k = KS, KE
       Rtot(k)  = Rdry  * ( 1.0_RP - qv(k) - qc(k) ) &
                + Rvap  * qv(k)
       CVtot(k) = CVdry * ( 1.0_RP - qv(k) - qc(k) ) &
                + CVvap * qv(k)                      &
                + CL    * qc(k)
       CPovCV(k) = ( CVtot(k) + Rtot(k) ) / CVtot(k)
    enddo

    ! make density
    do k = KS+1, KE

       RovCV  = Rtot(k) / CVtot(k)

       dens_s  = 0.0_RP
       dens(k) = dens(k-1)

       do ite = 1, itelim
          if( abs(dens(k)-dens_s) <= criteria ) exit

          dens_s = dens(k)

          dhyd = + ( P00 * ( dens(k-1) * Rtot(k-1) * pott(k-1) / P00 )**CPovCV(k-1) &
                   - P00 * ( dens_s    * Rtot(k  ) * pott(k  ) / P00 )**CPovCV(k  ) ) / FDZ(k-1) & ! dpdz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s )                                         ! rho*g

          dgrd = - P00 * ( Rtot(k) * pott(k) / P00 )**CPovCV(k) / FDZ(k-1) &
                 * CPovCV(k) * dens_s**RovCV                               &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

       enddo

       if ( ite > itelim ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k), dens_s, dhyd, dgrd, Rtot, FDZ(k-1)
       endif
    enddo

    do k = KS, KE
       pres(k) = P00 * ( dens(k) * Rtot(k) * pott(k) / P00 )**CPovCV(k)
       temp(k) = pres(k) / ( dens(k) * Rtot(k) )
    enddo

    ! fill KHALO
    dens(   1:KS-1) = dens(KS)
    dens(KE+1:KA  ) = dens(KE)

    return
  end subroutine buildrho_fromKS

end module mod_atmos_hydrostatic
