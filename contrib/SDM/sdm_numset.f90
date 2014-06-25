!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Set the number of super-droplets
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-06-14 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_numset

  implicit none
  private
  public :: sdm_numset 

contains
  subroutine sdm_numset(         &
      sdm_extbuf,                &
      sdm_aslset, sdm_sdnmlvol,  &
      sdm_inisdnc,               &
      sdm_aslfmsdnc,             &
      sdm_zlower, sdm_zupper,    &
      minzph, sdininum_s2c,      &
      sdfmnum_s2c, sdnum_s2c,    &
      ni_s2c, nj_s2c, nk_s2c,zph_crs )

    use scale_precision
    use scale_topography, only: &
         TOPO_Zsfc
    use scale_grid_real, only: &
!         REAL_CZ
         REAL_FZ
    use scale_grid_index, only: &
         IE,IS,KE,KS,JE,JS,IA,KA,JA ! S:start, E: end of active grids. A: num of grids including HALO.
    use m_sdm_common, only: &
         xmax_sdm, &    ! Width of the domain. xmax_sdm = FX(IE)-FX(IS-1)
         ymax_sdm       ! Depth of the domain. ymax_sdm = FY(JE)-FY(JS-1)
    
    ! These argument variables are all defined in m_sdm_common. Perhaps we should make it arguments only or common module only?
    integer, intent(in) :: sdm_extbuf    ! Rate of increase for extra buffer of super droplets
    integer, intent(in) :: sdm_aslset    ! Option for aerosol species
    real(RP),intent(in) :: sdm_sdnmlvol  ! Normal volume for number concentration
    real(RP),intent(in) :: sdm_inisdnc   ! Number of super droplets at initial per sdm_sdnmlvol
    real(RP),intent(in) :: sdm_aslfmsdnc ! Number of super droplets at aerosol formation per sdm_sdnmlvol
    real(RP),intent(in) :: sdm_zlower    ! Lowest limitaion of initial super droplets distribution
    real(RP),intent(in) :: sdm_zupper    ! Highest limitaion of initial super droplets distribution
    real(RP),intent(out) :: minzph       ! Lowest surface height of the domain
    real(RP),intent(out) :: sdininum_s2c ! Initial number of SDs used 
    integer,intent(out) :: sdfmnum_s2c   ! Number of SDs added when the aerosol formation func is called
    integer, intent(out) :: sdnum_s2c    ! Total number of SDs allocated
    integer, intent(out) :: ni_s2c, nj_s2c, nk_s2c ! Num of grids used by SDs
    real(RP),intent(out) :: zph_crs(KA,IA,JA) ! This is Real_CZ. Waste of memory. Will be removed in the near future.
    
    ! Work variables
    real(RP) :: bufsiz
    real(RP) :: dtmp          ! Temporary
    integer :: i, j, k        ! index

    ! Num of grids used by SDs
    ni_s2c = IE-IS+1
    nj_s2c = JE-JS+1
    nk_s2c = KE-KS+1

    minzph = TOPO_Zsfc(IS,JS) ! Is TOPO_Zsfc defined on the center? face?
    do j = JS, JE             ! Assuming it's on the center
       do i = IS, IE
          minzph = min( minzph, TOPO_Zsfc(i,j) )
       enddo
    enddo

    ! Waste of memory. zph_crs will be removed in the near future. At least we can use POINTER.
    ! Note also that CReSS uses (i,j,k) but SCALE uses (k,i,j). Be careful..
    do k = 1, KA
       do j = 1, JA
          do i = 1, IA
!             zph_crs(k,i,j) = REAL_CZ(k,i,j)
             zph_crs(k,i,j) = REAL_FZ(k,i,j)
          enddo
       enddo
    enddo
    
    !### setting number of super-droplets at initial ###!
    sdininum_s2c = real(xmax_sdm,kind=RP) * real(ymax_sdm,kind=RP)         &
         * real(sdm_zupper-(sdm_zlower+minzph),kind=RP)       &
         * real(sdm_inisdnc,kind=RP)/real(sdm_sdnmlvol,kind=RP)
    
    !### setting number of super-droplets at aerosol formation ###!
    if( abs(sdm_aslset)>10 ) then
       dtmp = real(xmax_sdm,kind=RP) * real(ymax_sdm,kind=RP)              &
            * real(sdm_zupper-(sdm_zlower+minzph),kind=RP)            &
            * real(sdm_aslfmsdnc,kind=RP)/real(sdm_sdnmlvol,kind=RP)
       sdfmnum_s2c = nint(dtmp)
    else
       sdfmnum_s2c = 1
    end if

    !### setting buffer size of super-droplets ###!
    bufsiz = nint(sdininum_s2c)
    bufsiz = nint( real(bufsiz) * (real(sdm_extbuf)*1.E-2_RP) )
    sdnum_s2c = nint(sdininum_s2c) + bufsiz
    
    return
  end subroutine sdm_numset
end module m_sdm_numset
