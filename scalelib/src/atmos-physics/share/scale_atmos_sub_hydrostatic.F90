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
     LASPdry => CONST_LASPdry, &
     P00     => CONST_PRE00
  use scale_grid, only: &
     CZ  => GRID_CZ, &
     FDZ => GRID_FDZ
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
  public :: ATMOS_HYDROSTATIC_buildrho_atmos
  public :: ATMOS_HYDROSTATIC_buildrho_bytemp
  public :: ATMOS_HYDROSTATIC_buildrho_bytemp_atmos

  interface ATMOS_HYDROSTATIC_buildrho
     module procedure ATMOS_HYDROSTATIC_buildrho_1D
     module procedure ATMOS_HYDROSTATIC_buildrho_3D
  end interface ATMOS_HYDROSTATIC_buildrho

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
  integer,  private, parameter :: itelim = 100 !< itelation number limit
  real(RP), private,      save :: criteria     !< convergence judgement criteria

  logical,  private,      save :: HYDROSTATIC_uselapserate = .false. !< use lapse rate?

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
       HYDROSTATIC_uselapserate

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[HYDROSTATIC]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_HYDROSTATIC,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_HYDROSTATIC. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_HYDROSTATIC)

    criteria = CONST_EPS * 5

    if( IO_L ) write(IO_FID_LOG,*) '*** buildrho conversion criteria:', criteria

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
               + CVvap * qv_sfc                       &
               + CL    * qc_sfc
    CPovCV_sfc = ( CVtot_sfc + Rtot_sfc ) / CVtot_sfc

    Rtot   = Rdry  * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + Rvap  * qv(KS)
    CVtot  = CVdry * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + CVvap * qv(KS)                       &
           + CL    * qc(KS)
    CPovCV = ( CVtot + Rtot ) / CVtot

    ! density at surface
    CVovCP_sfc = 1.0_RP / CPovCV_sfc
    dens_sfc   = P00 / Rtot_sfc / pott_sfc * ( pres_sfc/P00 )**CVovCP_sfc
    temp_sfc   = pres_sfc / ( dens_sfc * Rtot_sfc )

    ! make density at lowermost cell center
    if ( HYDROSTATIC_uselapserate ) then

       CPovR  = ( CVtot + Rtot ) / Rtot
       CVovCP = 1.0_RP / CPovCV

       temp(KS) = pott_sfc - LASPdry * CZ(KS) ! use dry lapse rate
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
                   - P00 * ( dens_s   * Rtot     * pott(KS) / P00 )**CPovCV     ) / CZ(KS) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens_sfc + dens_s )                                     ! rho*g

          dgrd = - P00 * ( Rtot * pott(KS) / P00 )**CPovCV / CZ(KS) &
                 * CPovCV * dens_s**RovCV                           &
                 - 0.5_RP * GRAV

          dens(KS) = dens_s - dhyd/dgrd

          if( dens(KS)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [buildrho 1D sfc] iteration not converged!', &
                                         dens(KS),ite,dens_s,dhyd,dgrd
          if( IO_L ) write(*         ,*) 'xxx [buildrho 1D sfc] iteration not converged!', &
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

    real(RP) :: dens_sfc  (1,IA,JA)

    real(RP) :: Rtot_sfc  (IA,JA)
    real(RP) :: CVtot_sfc (IA,JA)
    real(RP) :: CPovCV_sfc(IA,JA)
    real(RP) :: Rtot      (IA,JA)
    real(RP) :: CVtot     (IA,JA)
    real(RP) :: CPovCV    (IA,JA)

    real(RP) :: CVovCP_sfc, CPovR, CVovCP, RovCV
    real(RP) :: DZ
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: i, j
    !---------------------------------------------------------------------------

    !--- from surface to lowermost atmosphere

    do j = JS, JE
    do i = IS, IE
       Rtot_sfc  (i,j) = Rdry  * ( 1.0_RP - qv_sfc(1,i,j) - qc_sfc(1,i,j) ) &
                       + Rvap  * qv_sfc(1,i,j)
       CVtot_sfc (i,j) = CVdry * ( 1.0_RP - qv_sfc(1,i,j) - qc_sfc(1,i,j) ) &
                       + CVvap * qv_sfc(1,i,j)                              &
                       + CL    * qc_sfc(1,i,j)
       CPovCV_sfc(i,j) = ( CVtot_sfc(i,j) + Rtot_sfc(i,j) ) / CVtot_sfc(i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Rtot  (i,j) = Rdry  * ( 1.0_RP - qv(KS,i,j) - qc(KS,i,j) ) &
                   + Rvap  * qv(KS,i,j)
       CVtot (i,j) = CVdry * ( 1.0_RP - qv(KS,i,j) - qc(KS,i,j) ) &
                   + CVvap * qv(KS,i,j)                           &
                   + CL    * qc(KS,i,j)
       CPovCV(i,j) = ( CVtot(i,j) + Rtot(i,j) ) / CVtot(i,j)
    enddo
    enddo

    ! density at surface
    do j = JS, JE
    do i = IS, IE
       CVovCP_sfc      = 1.0_RP / CPovCV_sfc(i,j)
       dens_sfc(1,i,j) = P00 / Rtot_sfc(i,j) / pott_sfc(1,i,j) * ( pres_sfc(1,i,j)/P00 )**CVovCP_sfc
       temp_sfc(1,i,j) = pres_sfc(1,i,j) / ( dens_sfc(1,i,j) * Rtot_sfc(i,j) )
    enddo
    enddo

    ! make density at lowermost cell center
    if ( HYDROSTATIC_uselapserate ) then

       do j = JS, JE
       do i = IS, IE
          CPovR  = ( CVtot(i,j) + Rtot(i,j) ) / Rtot(i,j)
          CVovCP = 1.0_RP / CPovCV(i,j)

          temp(KS,i,j) = pott_sfc(1,i,j) - LASPdry * ( REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j) ) ! use dry lapse rate
          pres(KS,i,j) = P00 * ( temp(KS,i,j)/pott(KS,i,j) )**CPovR
          dens(KS,i,j) = P00 / Rtot(i,j) / pott(KS,i,j) * ( pres(KS,i,j)/P00 )**CVovCP
       enddo
       enddo

    else ! use itelation

       do j = JS, JE
       do i = IS, IE
          RovCV = Rtot(i,j) / CVtot(i,j)

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

             dhyd = + ( P00 * ( dens_sfc(1,i,j) * Rtot_sfc(i,j) * pott_sfc(1,i,j) / P00 )**CPovCV_sfc(i,j) &
                      - P00 * ( dens_s          * Rtot    (i,j) * pott   (KS,i,j) / P00 )**CPovCV    (i,j) ) / DZ & ! dp/dz
                    - GRAV * 0.5_RP * ( dens_sfc(1,i,j) + dens_s )                                                  ! rho*g

             dgrd = - P00 * ( Rtot(i,j) * pott(KS,i,j) / P00 )**CPovCV(i,j) / DZ &
                    * CPovCV(i,j) * dens_s**RovCV                                &
                    - 0.5_RP * GRAV

             dens(KS,i,j) = dens_s - dhyd/dgrd

             if( dens(KS,i,j)*0.0_RP /= 0.0_RP) exit
          enddo

          if ( .NOT. converged ) then
             if( IO_L ) write(IO_FID_LOG,*) 'xxx [buildrho 3D sfc] iteration not converged!', &
                                            i,j,dens(KS,i,j),ite,dens_s,dhyd,dgrd
             if( IO_L ) write(*         ,*) 'xxx [buildrho 3D sfc] iteration not converged!', &
                                            i,j,dens(KS,i,j),ite,dens_s,dhyd,dgrd
             call PRC_MPIstop
          endif
       enddo
       enddo

    endif

    !--- from lowermost atmosphere to top of atmosphere
    call ATMOS_HYDROSTATIC_buildrho_atmos_3D( dens(:,:,:), & ! [INOUT]
                                              temp(:,:,:), & ! [OUT]
                                              pres(:,:,:), & ! [OUT]
                                              pott(:,:,:), & ! [IN]
                                              qv  (:,:,:), & ! [IN]
                                              qc  (:,:,:)  ) ! [IN]

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_3D

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
                 + CVvap * qv(k)                      &
                 + CL    * qc(k)
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
                   - P00 * ( dens_s    * Rtot(k  ) * pott(k  ) / P00 )**CPovCV(k  ) ) / FDZ(k-1) & ! dpdz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s )                                          ! rho*g

          dgrd = - P00 * ( Rtot(k) * pott(k) / P00 )**CPovCV(k) / FDZ(k-1) &
                 * CPovCV(k) * dens_s**RovCV                               &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

          if( dens(k)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [buildrho 1D atmos] iteration not converged!', &
                                         k,dens(k),ite,dens_s,dhyd,dgrd
          if( IO_L ) write(*         ,*) 'xxx [buildrho 1D atmos] iteration not converged!', &
                                         k,dens(k),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo

    do k = KS, KE
       pres(k) = P00 * ( dens(k) * Rtot(k) * pott(k) / P00 )**CPovCV(k)
       temp(k) = pres(k) / ( dens(k) * Rtot(k) )
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_1D

  !-----------------------------------------------------------------------------
  !> Build up density from lowermost atmosphere (3D)
  subroutine ATMOS_HYDROSTATIC_buildrho_atmos_3D( &
       dens, &
       temp, &
       pres, &
       pott, &
       qv,   &
       qc    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: dens(KA,IA,JA) !< density               [kg/m3]
    real(RP), intent(out)   :: temp(KA,IA,JA) !< temperature           [K]
    real(RP), intent(out)   :: pres(KA,IA,JA) !< pressure              [Pa]
    real(RP), intent(in)    :: pott(KA,IA,JA) !< potential temperature [K]
    real(RP), intent(in)    :: qv  (KA,IA,JA) !< water vapor           [kg/kg]
    real(RP), intent(in)    :: qc  (KA,IA,JA) !< liquid water          [kg/kg]

    real(RP) :: Rtot  (KA,IA,JA)
    real(RP) :: CVtot (KA,IA,JA)
    real(RP) :: CPovCV(KA,IA,JA)

    real(RP) :: RovCV
    real(RP) :: DZ
    real(RP) :: dens_s, dhyd, dgrd
    integer  :: ite
    logical  :: converged

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Rtot  (k,i,j) = Rdry  * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                     + Rvap  * qv(k,i,j)
       CVtot (k,i,j) = CVdry * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                     + CVvap * qv(k,i,j)                          &
                     + CL    * qc(k,i,j)
       CPovCV(k,i,j) = ( CVtot(k,i,j) + Rtot(k,i,j) ) / CVtot(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS+1, KE
       RovCV = Rtot(k,i,j) / CVtot(k,i,j)

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

          dhyd = + ( P00 * ( dens(k-1,i,j) * Rtot(k-1,i,j) * pott(k-1,i,j) / P00 )**CPovCV(k-1,i,j) &
                   - P00 * ( dens_s        * Rtot(k  ,i,j) * pott(k  ,i,j) / P00 )**CPovCV(k  ,i,j) ) / DZ & ! dpdz
                 - GRAV * 0.5_RP * ( dens(k-1,i,j) + dens_s )                                                ! rho*g

          dgrd = - P00 * ( Rtot(k,i,j) * pott(k,i,j) / P00 )**CPovCV(k,i,j) / DZ &
                 * CPovCV(k,i,j) * dens_s**RovCV                                 &
                 - 0.5_RP * GRAV

          dens(k,i,j) = dens_s - dhyd/dgrd

          if( dens(k,i,j)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [buildrho 3D atmos] iteration not converged!', &
                                         k,i,j,dens(k,i,j),ite,dens_s,dhyd,dgrd
          if( IO_L ) write(*         ,*) 'xxx [buildrho 3D atmos] iteration not converged!', &
                                         k,i,j,dens(k,i,j),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pres(k,i,j) = P00 * ( dens(k,i,j) * Rtot(k,i,j) * pott(k,i,j) / P00 )**CPovCV(k,i,j)
       temp(k,i,j) = pres(k,i,j) / ( dens(k,i,j) * Rtot(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_atmos_3D

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
               + CVvap * qv_sfc                       &
               + CL    * qc_sfc

    Rtot   = Rdry  * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + Rvap  * qv(KS)
    CVtot  = CVdry * ( 1.0_RP - qv(KS) - qc(KS) ) &
           + CVvap * qv(KS)                       &
           + CL    * qc(KS)

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
                - dens_s   * Rtot     * temp(KS) ) / CZ(KS) & ! dp/dz
              - GRAV * 0.5_RP * ( dens_sfc + dens_s )         ! rho*g

       dgrd = - Rtot * temp(KS) / CZ(KS) &
              - 0.5_RP * GRAV

       dens(KS) = dens_s - dhyd/dgrd

       if( dens(KS)*0.0_RP /= 0.0_RP) exit
    enddo

    if ( .NOT. converged ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [buildrho bytemp 1D sfc] iteration not converged!', &
                                      dens(KS),ite,dens_s,dhyd,dgrd
       if( IO_L ) write(*         ,*) 'xxx [buildrho bytemp 1D sfc] iteration not converged!', &
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

    do j = JS, JE
    do i = IS, IE
       Rtot_sfc (i,j) = Rdry  * ( 1.0_RP - qv_sfc(1,i,j) - qc_sfc(1,i,j) ) &
                      + Rvap  * qv_sfc(1,i,j)
       CVtot_sfc(i,j) = CVdry * ( 1.0_RP - qv_sfc(1,i,j) - qc_sfc(1,i,j) ) &
                      + CVvap * qv_sfc(1,i,j)                              &
                      + CL    * qc_sfc(1,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Rtot (i,j) = Rdry  * ( 1.0_RP - qv(KS,i,j) - qc(KS,i,j) ) &
                  + Rvap  * qv(KS,i,j)
       CVtot(i,j) = CVdry * ( 1.0_RP - qv(KS,i,j) - qc(KS,i,j) ) &
                  + CVvap * qv(KS,i,j)                           &
                  + CL    * qc(KS,i,j)
    enddo
    enddo

    ! density at surface
    do j = JS, JE
    do i = IS, IE
       RovCP_sfc       = Rtot_sfc(i,j) / ( CVtot_sfc(i,j) + Rtot_sfc(i,j) )
       dens_sfc(1,i,j) = pres_sfc(1,i,j) / ( Rtot_sfc(i,j) * temp_sfc(1,i,j) )
       pott_sfc(1,i,j) = temp_sfc(1,i,j) / ( P00/pres_sfc(1,i,j) )**RovCP_sfc
    enddo
    enddo

    ! make density at lowermost cell center
    do j = JS, JE
    do i = IS, IE
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
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [buildrho bytemp 3D sfc] iteration not converged!', &
                                         i,j,dens(KS,i,j),ite,dens_s,dhyd,dgrd
          if( IO_L ) write(*         ,*) 'xxx [buildrho bytemp 3D sfc] iteration not converged!', &
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
                + CVvap * qv(k)                      &
                + CL    * qc(k)
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
                   - dens_s    * Rtot(k  ) * temp(k  ) ) / FDZ(k-1) & ! dp/dz
                 - GRAV * 0.5_RP * ( dens(k-1) + dens_s )             ! rho*g

          dgrd = - Rtot(k) * temp(k) / FDZ(k-1) &
                 - 0.5_RP * GRAV

          dens(k) = dens_s - dhyd/dgrd

          if( dens(k)*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [buildrho bytemp 1D atmos] iteration not converged!', &
                                         k,dens(k),ite,dens_s,dhyd,dgrd
          if( IO_L ) write(*         ,*) 'xxx [buildrho bytemp 1D atmos] iteration not converged!', &
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

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Rtot (k,i,j) = Rdry  * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                    + Rvap  * qv(k,i,j)
       CVtot(k,i,j) = CVdry * ( 1.0_RP - qv(k,i,j) - qc(k,i,j) ) &
                    + CVvap * qv(k,i,j)                          &
                    + CL    * qc(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
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
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [buildrho bytemp 3D atmos] iteration not converged!', &
                                         k,i,j,dens(k,i,j),ite,dens_s,dhyd,dgrd
          if( IO_L ) write(*         ,*) 'xxx [buildrho bytemp 3D atmos] iteration not converged!', &
                                         k,i,j,dens(k,i,j),ite,dens_s,dhyd,dgrd
          call PRC_MPIstop
       endif
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RovCP   = Rtot(k,i,j) / ( CVtot(k,i,j) + Rtot(k,i,j) )
       pres(k,i,j) = dens(k,i,j) * Rtot(k,i,j) * temp(k,i,j)
       pott(k,i,j) = temp(k,i,j) * ( P00 / pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_HYDROSTATIC_buildrho_bytemp_atmos_3D

end module scale_atmos_hydrostatic
!-------------------------------------------------------------------------------
