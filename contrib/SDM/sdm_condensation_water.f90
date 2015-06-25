!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Condensation of water to the particles
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
!! @li      2014-07-11 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-18 (Y.Sato)  [rev] Change Teten's formula to SATURATION library of SCALE Library 
!! @li      2014-07-24 (Y.Sato)  [mod] Modify a bug relating to the revision on 2014/07/18
!! @li      2014-12-12 (Y.Sato)  [mod] modify epsva, LatHet, and DNS_RL as those used in SCALE Library
!! @li      2014-12-19 (Y.Sato)  [mod] modify the location for defining LatHet
!! @li      2014-12-24 (Y.Sato)  [mod] Modify the Latent Heat for SCALE library 
!! @li      2015-06-25 (S.Shima) [add] fapp_start/stop added for performance monitoring
!! @li      2015-06-26 (S.Shima) [add] OCL added for Auto parallelization on K/FX10
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_condensation_water
  use scale_precision

  implicit none
  private
  public :: sdm_condevp,sdm_condevp_updatefluid

contains
  subroutine sdm_condevp(sdm_aslset,                     &
       sdm_aslmw,sdm_aslion,sdm_dtevl, &
       pres_scale,t_scale,qv_scale, &
       sd_num,sd_numasl,sd_x,sd_y,     &
       sd_r,sd_asl,sd_ri,sd_rj,sd_rk              )
    
    use scale_const, only: &
         cp   => CONST_CPdry, &
         p0   => CONST_PRE00, &
         t0   => CONST_TEM00, &
         es0  => CONST_PSAT0, &
         epsva=> CONST_EPSvap, &
         rd   => CONST_Rdry
    use scale_grid_index, only: &
         IA,JA,KA
    use scale_atmos_saturation, only: &
        ATMOS_SATURATION_pres2qsat_liq, &
        ATMOS_SATURATION_psat_liq
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         mass_amsul, ion_amsul, mass_nacl, ion_nacl, &
         VALID2INVALID, &
         RLRv_D, LatGas, L_RL_K, CurveF, ASL_FF
    ! Input variables
    integer,  intent(in) :: sdm_aslset
    real(RP), intent(in) :: sdm_aslmw(20)
    real(RP), intent(in) :: sdm_aslion(20)
    real(RP), intent(in) :: sdm_dtevl  ! tims step of {condensation/evaporation} process
    real(RP), intent(in) :: pres_scale(KA,IA,JA) ! Pressure
    real(RP), intent(in) :: t_scale(KA,IA,JA) ! Temperature
    real(RP), intent(in) :: qv_scale(KA,IA,JA)   ! Water vapor mixing ratio
    integer,  intent(in) :: sd_num      ! number of super-droplets
    integer,  intent(in) :: sd_numasl   ! number of kind of chemical material contained as water-soluble aerosol in super droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)   ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)   ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)   ! index[k/real] of super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_r(1:sd_num) ! equivalent radius of super-droplets

    ! Internal shared variables
    real(RP):: sd_aslmw(1:22) ! Molecular mass of chemical material contained as water-soluble aerosol in super droplets (default+20)
    real(RP):: sd_aslion(1:22) ! Degree of ion dissociation of chemical material contained as water-soluble aerosol in super droplets (default+20)

    ! Work variables
    real(RP) :: dmask(1:22)  ! mask for vactorization
    real(RP) :: p_sd      ! pressure of the grid contained the SD
    real(RP) :: t_sd      ! temperature of the grid contained the SD
    real(RP) :: qv_sd     ! water-vapor of the grid contained the SD
    real(RP) :: qvs_sd    ! Saturation mixing ratio of the grid contained the SD.
    real(RP) :: es_sd     ! Saturation vapor pressure of the grid contained the SD.
    real(RP) :: ss_sd     ! Degree of super-saturation of the grid contained the SD.
    real(RP) :: Fac_dd    ! Fd in growth EQ.(R.R.Rogers)
    real(RP) :: Fac_kk    ! Fk in growth EQ.(R.R.Rogers)
    real(RP) :: Rc        ! Rc
    real(RP) :: ivt_sd    ! 1.d0 / t_sd
    real(RP) :: dtivfdk   ! dt / ( Fk + Td )
    real(RP) :: crd       ! sd_r before interation
    real(RP) :: rdi       ! sd_r after interation
    real(RP) :: a         ! temporary
    real(RP) :: a3        ! temporary
    real(RP) :: b         ! temporary
    real(RP) :: eq_a      ! temporary for growth EQ.
    real(RP) :: eq_b      ! temporary for growth EQ.
    real(RP) :: eq_c      ! temporary for growth EQ.
    real(RP) :: new_rd    ! temporary for interation
    real(RP) :: rd2       ! temporary for interation
    real(RP) :: rd5       ! temporary for interation
    real(RP) :: rd7       ! temporary for interation
    real(RP) :: dterm1    ! temporary for interation
    real(RP) :: dterm2    ! temporary for interation
    real(RP) :: dterm3    ! temporary for interation
    real(RP) :: dtmp      ! temporary for interation
    real(RP) :: rddvcp    ! rd / cp
    integer :: idx_nasl(1:22)  ! index for vactorization

    integer :: i, j, k, n, s, t, it               ! index
    ! Parameters
    integer, parameter :: itr_max = 25   ! iteration number
!!$    real(RP), parameter :: epsva = 0.622_RP   ! Molecular weight ratio of vapor/air
    !---------------------------------------------------------------------
      
    ! Section specification for fapp profiler
    call fapp_start("sdm_condevp",1,1)

    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

    rddvcp = real(rd,kind=RP)/real(cp,kind=RP)

    !### aerosol type ###!
    
    if( abs(mod(sdm_aslset,10))==1 ) then

       !### numasl=1 @ init+rest : (NH4)2SO4 ###!

       sd_aslmw(1)  = mass_amsul
       sd_aslion(1) = ion_amsul

    else if( abs(mod(sdm_aslset,10))==2 ) then

       if( abs(sdm_aslset)==2 ) then
          
          !### numasl=1 @ init : NaCl ###!

          sd_aslmw(1)  = mass_nacl
          sd_aslion(1) = ion_nacl

       else if( abs(sdm_aslset)==12 ) then

          !### numasl=2 @ init : NaCl, rest : (NH4)2SO4 ###!

          sd_aslmw(1) = mass_amsul
          sd_aslmw(2) = mass_nacl
          sd_aslion(1) = ion_amsul
          sd_aslion(2) = ion_nacl

       end if

    else if( abs(mod(sdm_aslset,10))==3 ) then

       !### numasl>=2 @ init+rest : (NH4)2SO4, NaCl, ... ###!

       sd_aslmw(1) = mass_amsul
       sd_aslmw(2) = mass_nacl

       sd_aslion(1) = ion_amsul
       sd_aslion(2) = ion_nacl

       !! Must be a Bug. This cannot be simply commented out
       do n=1,20
       !            call getrname( id_sdm_aslmw  + (n-1), sd_aslmw(n+2)  )
       !            call getrname( id_sdm_aslion + (n-1), sd_aslion(n+2) )
       end do

    end if

    do n=1,22

       if( n<=sd_numasl ) then
          idx_nasl(n) = n
          dmask(n) = 1.0_RP
       else
          idx_nasl(n) = sd_numasl
          dmask(n) = 0.0_RP
       end if

    end do


    ! Condensation of water
!OCL INDEPENDENT
    do n=1,sd_num

       !### Skip invalid super-droplets ###!

       if( sd_rk(n)<VALID2INVALID ) cycle

       !### Get the location and variables of Super-Droplets ###!
       
       i = floor(sd_ri(n))+1
       j = floor(sd_rj(n))+1
       k = floor(sd_rk(n))+1

       p_sd  = pres_scale(k,i,j)
       t_sd  = t_scale(k,i,j)
       qv_sd = qv_scale(k,i,j)

       !### Calculate degree of super-saturation ###!

       ivt_sd = 1.0_RP / t_sd

       a = 1.0_RP / ( t_sd - 35.860_RP )
       b = a * ( t_sd - real(t0,kind=RP) )

!!       es_sd  = real(es0,kind=RP) * exp(17.269_RP*b)   !! Teten-eq.
!!       qvs_sd = epsva * es_sd / ( p_sd - es_sd )

       !! When ice phase is implemented ss_ice should be also implement
!!$       call ATMOS_SATURATION_pres2qsat_ice( qvsi_sd,t_sd,p_sd )
!!$       call ATMOS_SATURATION_psat_ice( esi_sd,t_sd )
       call ATMOS_SATURATION_pres2qsat_liq( qvs_sd,t_sd,p_sd )
       call ATMOS_SATURATION_psat_liq( es_sd,t_sd )
       ss_sd  = qv_sd / qvs_sd          !! degree of super-saturation

       !### Set parameters for growth EQ. of the radius ###!
       
       !  LatGas = LatHet / GasV_C
       !  L_RL_K = LatHet * DNS_RL / Heat_C
       !  RLRv_D = DNS_RL * GasV_C / Diff_C
       !  ASL_RR = ASL_FF * ION_asl / WGTasl
       
       Fac_dd  = RLRv_D * t_sd / es_sd
       Fac_kk  = ( LatGas*ivt_sd - 1.0_RP ) * L_RL_K * ivt_sd
       dtivFdk = real( sdm_dtevl, kind=RP ) / ( Fac_dd + Fac_kk )

       eq_a  = CurveF * ivt_sd

       eq_b = 0.0_RP

!OCL UNROLL('full'),NOSWP  
       do t=1,22

          s = idx_nasl(t)

          dtmp = sd_asl(n,s) * (real(sd_aslion(s),kind=RP)            &
               / real(sd_aslmw(s),kind=RP))
          eq_b = eq_b + dmask(t) * dtmp

       end do

       eq_b = eq_b * ASL_FF

       Rc    = sqrt(eq_b/eq_a)
       eq_c  = ss_sd - 1.0_RP

       a  = eq_a / eq_c
       b  = eq_b / eq_c
       a3 = a * a * a

       !### advance the particle radius ###!
       
       crd    = sd_r(n)
       new_rd = sd_r(n)

       if( ( ss_sd > 1.0_RP ) .and. ( a3 < b*(27.0_RP/4.0_RP) ) ) then
          new_rd = 1.0E-3_RP
       end if

       if( new_rd<Rc ) then
          new_rd = Rc
       end if

       !== iteration ( newton-raphson method ) ==!

       rdi = new_rd

!OCL FRECIPRO,UNROLL
       do it=1,itr_max

          ! Considering Kohler-curve (R.R.Rogers)
          
          rd2 = rdi * rdi
          rd5 = rd2 * rd2 * rdi
          rd7 = rd5 * rd2

          dterm1 = dtivFdk * ( rd5*eq_c + (eq_b-eq_a*rd2)*rd2 )
          dterm2 = rd5 * crd * crd
          dterm3 = rd5 - dtivFdk*( -3.0_RP*eq_b + eq_a*rd2 )

          dtmp = rd2 - ( rd7 - 2.0_RP*dterm1 - dterm2 ) / dterm3

          if( dtmp<=0.0_RP ) dtmp = rdi * rdi * 1.E-4_RP

          rdi = sqrt(dtmp)

       end do

       sd_r(n) = rdi   !! particle radius at future

    end do

    ! Section specification for fapp profiler
    call fapp_stop("sdm_condevp",1,1)

    return
  end subroutine sdm_condevp
  !-----------------------------------------------------------------------------
  subroutine sdm_condevp_updatefluid(RHOT,QTRC,DENS,rhowp_sdm,rhowc_sdm)
    use scale_tracer, only: &
         I_QV,QAD=>QA
    use scale_atmos_thermodyn, only: &
        ATMOS_THERMODYN_templhv
    use scale_grid_index, only: &
         IA,JA,KA,IS,IE,JS,JE,KS,KE
    use m_sdm_fluidconv, only: &
         sdm_rhot_qtrc2cpexnr, &
         sdm_rhot_qtrc2p_t
    use m_sdm_common, only: &
         LatHet
    ! Input variables
    real(RP), intent(in) :: rhowp_sdm(KA,IA,JA) ! density of liquid water at preveous step
    real(RP), intent(in) :: rhowc_sdm(KA,IA,JA) ! density of liquid water at current step
    ! Output variables
    real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]
    real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    ! Work variables
    real(RP) :: rhov(KA,IA,JA),cpexnr(KA,IA,JA)
    real(RP) :: delta_rhow, dens_old
    integer  :: i, j, k ! index
    real(RP) :: lhv(KA,IA,JA), pre(KA,IA,JA), temp(KA,IA,JA)
    !-------------------------------------------------------------------7--
     
    ! exchange vapor
    do k = 1, KA
    do i = 1, IA
    do j = 1, JA
       rhov(k,i,j)=DENS(k,i,j)*QTRC(k,i,j,I_QV)
    end do
    end do
    end do

    do k=KS,KE
    do j=JS,JE
    do i=IS,IE
       delta_rhow = rhowc_sdm(k,i,j)-rhowp_sdm(k,i,j)
       dens_old = DENS(k,i,j)
       
       DENS(k,i,j) = DENS(k,i,j) - delta_rhow
       rhov(k,i,j) = rhov(k,i,j) - delta_rhow
       QTRC(k,i,j,I_QV) = rhov(k,i,j) / DENS(k,i,j)
       RHOT(k,i,j) = RHOT(k,i,j) * (DENS(k,i,j)/dens_old)  
    end do
    end do
    end do

    call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pre,temp)
    call ATMOS_THERMODYN_templhv(lhv,temp)

    ! exchange heat
    call sdm_rhot_qtrc2cpexnr(RHOT,QTRC,DENS,cpexnr)

    do k=KS,KE
    do j=JS,JE
    do i=IS,IE
       delta_rhow = rhowc_sdm(k,i,j)-rhowp_sdm(k,i,j)
       
!!$       RHOT(k,i,j) = RHOT(k,i,j) + LatHet*delta_rhow / cpexnr(k,i,j)  
       RHOT(k,i,j) = RHOT(k,i,j) + lhv(k,i,j)*delta_rhow / cpexnr(k,i,j)  
    end do
    end do
    end do

    return
  end subroutine sdm_condevp_updatefluid
end module m_sdm_condensation_water
