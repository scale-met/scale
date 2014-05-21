!-------------------------------------------------------------------------------
!> module URBAN / Physics Urban Canopy Model (UCM)
!!
!! @par Description
!!          Urban physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_urban_phy_ucm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_PHY_driver_setup
  public :: URBAN_PHY_driver_first
  public :: URBAN_PHY_driver_final

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_PHY_driver_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_urban_vars, only: &
       URBAN_TYPE
    implicit none

    logical :: dummy

    NAMELIST / PARAM_URBAN_UCM / &
       dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[UCM]/Categ[URBAN]'

    if ( URBAN_TYPE /= 'UCM' ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx URBAN_TYPE is not UCM. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_UCM,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN_UCM. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_URBAN_UCM)

    return
  end subroutine URBAN_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for urban submodel
  subroutine URBAN_PHY_driver_first
    use scale_time, only: &
       dt => TIME_DTSEC_URBAN
    use scale_grid_real, only: &
       REAL_lon, &
       REAL_lat
    use mod_urban_vars, only: &
       TR_URB,  &
       TG_URB,  &
       TB_URB,  &
       TC_URB,  &
       QC_URB,  &
       UC_URB,  &
       TRL_URB, &
       TGL_URB, &
       TBL_URB
    use mod_cpl_vars, only: &
       CPL_getUrb
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Urban step: UCM'

    call CPL_getUrb( DZ  (:,:), & ! (out)
                         DENS(:,:), & ! (out)
                         MOMX(:,:), & ! (out)
                         MOMY(:,:), & ! (out)
                         MOMZ(:,:), & ! (out)
                         TEMP(:,:), & ! (out)
                         QV  (:,:), & ! (out)
                         SWD (:,:), & ! (out)
                         LWD (:,:), & ! (out)
                         PREC(:,:)  ) ! (out)

    do j = JS-1, JE+1
    do i = IS-1, IE+1

      TA = TEMP(i,j)                  ! temp at 1st atmospheric level         [K]
      QA = QV(i,j)  /(1-QV(i,j))      ! mixing ratio at 1st atmospheric level [kg/kg]
                                      ! QV specific humidity                  [kg/kg]
      UA = sqrt(          &
            ( MOMZ(i,j)               )**2 &
          + ( MOMX(i-1,j) + MOMX(i,j) )**2 &
          + ( MOMY(i,j-1) + MOMY(i,j) )**2 &
          ) / DENS(i,j) * 0.5_RP
                                  ! wind speed at 1st atmospheric level   [m/s]
      U1 = 0.5_RP * ( MOMX(i-1,j) + MOMX(i,j) ) / DENS(i,j)
                                  ! u at 1st atmospheric level            [m/s]
      V1 = 0.5_RP * ( MOMY(i,j-1) + MOMY(i,j) ) / DENS(i,j)
                                  ! v at 1st atmospheric level            [m/s]
      SSG  = SWD(i,j)             ! downward total short wave radiation   [W/m/m]
      LLG  = LWD(i,j)             ! downward long wave radiation          [W/m/m]
      RAIN = PREC(i,j)            ! precipitation                         [mm/h]
      RHOO = DENS(i,j)            ! air density                           [kg/m^3]
      ZA   = DZ(i,j)              ! first atmospheric level               [m]
      XLON = REAL_lon(i,j)        ! longitude                             [deg]
      XLAT = REAL_lat(i,j)        ! latitude                              [deg]

      TR   = TR_URB(i,j)
      TB   = TB_URB(i,j)
      TG   = TG_URB(i,j)
      TC   = TC_URB(i,j)
      QC   = QC_URB(i,j)
      UC   = UC_URB(i,j)

      do k = UKS, UKE
        TRL(k) = TRL_URB(k,i,j)
        TBL(k) = TBL_URB(k,i,j)
        TGL(k) = TGL_URB(k,i,j)
      enddo

      call urban(LSOLAR,                             & ! (in)
                 TA, QA, UA, U1, V1, ZA,             & ! (in)
                 SSG, LLG, RAIN, RHOO, XLON, XLAT,   & ! (in)
                 TR, TB, TG, TC, QC, UC,             & ! (inout)
                 TRL, TBL, TGL,                      & ! (inout)
                 TS, SH, LH, SW, LW, G               ) ! (out)

      do k = UKS, UKE
        TRL_URB(k,i,j) = TRL(k)
        TBL_URB(k,i,j) = TBL(k)
        TGL_URB(k,i,j) = TGL(k)
      enddo

      TR_URB(i,j) = TR
      TB_URB(i,j) = TB
      TG_URB(i,j) = TG
      TC_URB(i,j) = TC
      QC_URB(i,j) = QC
      UC_URB(i,j) = UC
      TS_URB(i,j) = TS

      SHFLX_URB(i,j)  = SH  ! sensible heat flux               [W/m/m]
      LHFLX_URB(i,j)  = LH  ! latent heat flux                 [W/m/m]
      GHFLX_URB(i,j)  = G   ! heat flux into the ground        [W/m/m]
      SWUFLX_URB(i,j) = SW  ! upward short wave radiation flux [W/m/m]
      LWUFLX_URB(i,j) = LW  ! upward long wave radiation flux  [W/m/m]

    end do
    end do

    return
  end subroutine URBAN_PHY_driver_first

  subroutine URBAN_PHY_driver_final
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine URBAN_PHY_driver_final

end module mod_urban_phy_ucm
