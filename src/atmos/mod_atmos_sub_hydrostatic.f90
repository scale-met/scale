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
  public :: ATMOS_hydro_buildrho_1d

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
      qv_sfc,&
      qc_sfc    )
    use mod_const, only : &
       GRAV    => CONST_GRAV,    &
       EPS     => CONST_EPS,     &
       Rdry    => CONST_Rdry,    &
       Rvap    => CONST_Rvap,    &
       CPovR   => CONST_CPovR,   &
       RovCV   => CONST_RovCV,   &
       CVovCP  => CONST_CVovCP,  &
       CPovCV  => CONST_CPovCV,  &
       LASPdry => CONST_LASPdry, &
       P00     => CONST_PRE00
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_grid, only : &
       CZ  => GRID_CZ, &
       FDZ => GRID_FDZ
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

    integer :: i, j
    integer :: ierr

    NAMELIST / PARAM_ATMOS_HYDROSTATIC / &
       HYDROSTATIC_userapserate

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

       call buildrho(dens(:,i,j), temp(:,i,j), pres(:,i,j),         & ! (out)
                     pott(:,i,j), qv(:,i,j), qc(:,i,j),             & ! (in)
                     temp_sfc(1,i,j),                               & ! (out)
                     pres_sfc(1,i,j), pott_sfc(1,i,j),              & ! (in)
                     qv_sfc(1,i,j), qc_sfc(1,i,j),                  & ! (in)
                     HYDROSTATIC_userapserate                       ) ! (in)

    enddo
    enddo

    ! fill IHALO & JHALO
    call COMM_vars8( dens(:,:,:), 1 )
    call COMM_wait ( dens(:,:,:), 1 )

    return
  end subroutine ATMOS_hydro_buildrho

  !-----------------------------------------------------------------------------
  !> Buildup density from surface
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
      qv_sfc,&
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

    integer :: ierr

    NAMELIST / PARAM_ATMOS_HYDROSTATIC / &
         HYDROSTATIC_userapserate

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

    call buildrho(dens, temp, pres,          & ! (out)
         pott, qv, qc,                       & ! (in)
         temp_sfc,                           & ! (out)
         pres_sfc, pott_sfc, qv_sfc, qc_sfc, & ! (in)
         HYDROSTATIC_userapserate            ) ! (in)

    return
  end subroutine ATMOS_hydro_buildrho_1d

  !-----------------------------------------------------------------------------
  !> Buildup density from surface
  !-----------------------------------------------------------------------------
  subroutine buildrho( &
      dens,     &
      temp,     &
      pres,     &
      pott,     &
      qv,       &
      qc,       &
      temp_sfc, &
      pres_sfc, &
      pott_sfc, &
      qc_sfc,   &
      qv_sfc,   &
      userapserate )
    use mod_const, only : &
       GRAV    => CONST_GRAV,    &
       EPS     => CONST_EPS,     &
       Rdry    => CONST_Rdry,    &
       Rvap    => CONST_Rvap,    &
       CVdry   => CONST_CVdry,   &
       CVvap   => CONST_CVvap,   &
       CL      => CONST_CL,      &
       LASPdry => CONST_LASPdry, &
       P00     => CONST_PRE00
    use mod_grid, only : &
       CZ  => GRID_CZ, &
       FDZ => GRID_FDZ
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


    real(RP) :: Qdry_sfc
    real(RP) :: Rmoist_sfc
    real(RP) :: CVtot_sfc
    real(RP) :: dens_sfc

    real(RP) :: Qdry(KA)
    real(RP) :: Rmoist(KA)
    real(RP) :: CVtot(KA)
    real(RP) :: CPtot(KA)

    real(RP) :: CPovCV(KA)

    real(RP) :: dens_s, dhyd, dgrd

    real(RP) :: criteria
    integer, parameter :: itelim = 100

    integer :: k, ite

    !---------------------------------------------------------------------------

    criteria = EPS * 5

    ! make density at surface
    Qdry_sfc = 1.0_RP - qv_sfc - qc_sfc
    Rmoist_sfc = Rdry * Qdry_sfc + Rvap * qv_sfc
    CVtot_sfc = CVdry * Qdry_sfc + CVvap * qv_sfc + CL * qc_sfc

    dens_sfc = P00 / Rmoist_sfc / pott_sfc * ( pres_sfc/P00 )**(CVtot_sfc/(CVtot_sfc+Rmoist_sfc))
    temp_sfc = pres_sfc / ( dens_sfc * Rmoist_sfc )

    do k = KS, KE
       Qdry(k) = 1.0_RP - qv(k) - qc(k)
       Rmoist(k) = Rdry * Qdry(k) + Rvap * qv(k)
       CVtot(k) = CVdry * Qdry(k) + CVvap * qv(k) + CL * qc(k)
       CPtot(k) = CVtot(k) + Rmoist(k)
       CPovCV(k) = CPtot(k) / CVtot(k)
    enddo

    ! make density at lowermost cell center
    k = KS

    if ( userapserate ) then

       temp(k) = pott_sfc - LASPdry * CZ(k) ! use dry lapse rate
       pres(k) = P00 * ( temp(k)/pott(k) )**(CPtot(k)/Rmoist(k))
       dens(k) = P00 / Rmoist(k) / pott(k) * ( pres(k)/P00 )**(CVtot(k)/CPtot(k))

    else ! use itelation

       dens_s      = 0.0_RP
       dens(k) = dens_sfc ! first guess

       do ite = 1, itelim
          if( abs(dens(k)-dens_s) <= criteria ) exit

          dens_s = dens(k)

          dhyd = + ( P00 * ( dens_sfc * Rmoist_sfc * pott_sfc / P00 )**((CVtot_sfc+Rmoist_sfc)/CVtot_sfc) &
                   - P00 * ( dens_s   * Rmoist(k)  * pott(k)  / P00 )**CPovCV(k)) / CZ(k) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens_sfc + dens_s )                                                 ! rho*g

          dgrd = - P00 * ( Rmoist(k) * pott(k) / P00 )**CPovCV(k) / CZ(k) &
                 * CPovCV(k) * dens_s**(Rmoist(k)/CVtot(k))               &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

       enddo

       if ( ite > itelim ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k), dens_s, dhyd, dgrd
       endif

    endif

    ! make density
    do k = KS+1, KE

       dens_s      = 0.0_RP
       dens(k) = dens(k-1)

       do ite = 1, itelim
          if( abs(dens(k)-dens_s) <= criteria ) exit

          dens_s = dens(k)

          dhyd = + ( P00 * ( dens(k-1) * Rmoist(k-1) * pott(k-1) / P00 )**CPovCV(k-1) &
                   - P00 * ( dens_s    * Rmoist(k  ) * pott(k  ) / P00 )**CPovCV(k) ) / FDZ(k-1) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s )                                                ! rho*g

          dgrd = - P00 * ( Rmoist(k) * pott(k) / P00 )**CPovCV(k) / FDZ(k-1) &
                 * CPovCV(k) * dens_s**(Rmoist(k)/CVtot(k))                  &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

       enddo

       if ( ite > itelim ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx iteration not converged!', k, ite, dens(k), dens_s, dhyd, dgrd, Rmoist, FDZ(k-1)
       endif
    enddo

    do k = KS, KE
       pres(k) = P00 * ( dens(k) * Rmoist(k) * pott(k) / P00 )**CPovCV(k)
       temp(k) = pres(k) / ( dens(k) * Rmoist(k) )
    enddo

    ! fill KHALO
    dens(   1:KS-1) = dens(KS)
    dens(KE+1:KA  ) = dens(KE)

    return
  end subroutine buildrho

end module mod_atmos_hydrostatic
