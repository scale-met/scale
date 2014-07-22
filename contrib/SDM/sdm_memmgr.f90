!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Memory management of the SDM variables
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
!! @li      2014-06-23 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-14 (S.Shima) [rev] Removed unused variables
!! @li      2014-07-18 (Y.Sato)  [add] add QTRC_sdm
!! @li      2014-07-22 (Y.Sato)  [mod] modify the bug at the allocation of QTRC_sdm
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_memmgr

  implicit none
  private
  public :: sdm_allocinit

contains
  !-----------------------------------------------------------------------------
  subroutine sdm_allocinit
    use scale_gridtrans, only: &
       I_XYZ, I_XYW,    &
       GTRANS_GSQRT, &
       GTRANS_J13G,  &
       GTRANS_J23G,  &
       GTRANS_J33G
    use scale_tracer_sdm, only: &
         QA
    use scale_grid_index, only: &
         IE,IS,KE,KS,JE,JS,IA,KA,JA ! S:start, E: end of active grids. A: num of grids including HALO.
    use m_sdm_common
    integer :: i, j, k, n, s
    integer :: bndsdmdim, bufsiz ! For send/receive buffer. bndsdmdim is the number of variables for each sd. bufsiz is the num of sds that the buffer can store.
                                 ! There exists a bug on the multiplicity transfer, and should be fixed in the near future. 
    !------------------------------------------------------
    ! Get numbers of kind of chemical material contained as
    ! water-soluble aerosol in super droplets.
    if( abs(mod(sdm_aslset,10))==1 ) then
       !### init+rest : (NH4)2SO4 ###!
       sdnumasl_s2c = 1
    else if( abs(mod(sdm_aslset,10))==2 ) then
       if( abs(sdm_aslset)==2 ) then
          !### init : NaCl ###!
          sdnumasl_s2c = 1
       else if( abs(sdm_aslset)==12 ) then
          !### init : NaCl, rest : (NH4)2SO4 ###!
          sdnumasl_s2c = 2
       end if
    else if( abs(mod(sdm_aslset,10))==3 ) then
       !### init+rest : (NH4)2SO4, NaCl, ... ###!
       sdnumasl_s2c = 2   !! default
       do n=1,20
          if( sdm_aslmw(n)>0.0_RP ) then
             sdnumasl_s2c = n + 2
          end if
       end do
    end if
    bndsdmdim = 8 + sdnumasl_s2c   !! n,x,y,rk,u,v,w,r + asl(1:sdnumasl)
    
    !--- Allocate arrays for SDM
    allocate(sdn_s2c(1:sdnum_s2c))
    allocate(sdri_s2c(1:sdnum_s2c))
    allocate(sdrj_s2c(1:sdnum_s2c))
    allocate(sdrk_s2c(1:sdnum_s2c))
    allocate(sdx_s2c(1:sdnum_s2c))
    allocate(sdy_s2c(1:sdnum_s2c))
    allocate(sdz_s2c(1:sdnum_s2c))
    allocate(sdr_s2c(1:sdnum_s2c))
    allocate(sdu_s2c(1:sdnum_s2c))
    allocate(sdv_s2c(1:sdnum_s2c))
    allocate(sdvz_s2c(1:sdnum_s2c))
    allocate(sdasl_s2c(1:sdnum_s2c,1:sdnumasl_s2c))
    allocate(sdn_fm(1:sdfmnum_s2c))
    allocate(sdri_fm(1:sdfmnum_s2c))
    allocate(sdrj_fm(1:sdfmnum_s2c))
    allocate(sdrk_fm(1:sdfmnum_s2c))
    allocate(sdx_fm(1:sdfmnum_s2c))
    allocate(sdy_fm(1:sdfmnum_s2c))
    allocate(sdz_fm(1:sdfmnum_s2c))
    allocate(sdr_fm(1:sdfmnum_s2c))
    allocate(sdvz_fm(1:sdfmnum_s2c))
    allocate(sdasl_fm(1:sdnum_s2c,1:sdnumasl_s2c))

    allocate(sdrkl_s2c(IA,JA))
    allocate(sdrku_s2c(IA,JA))

    ! These variables are not needed for SCALE. They are for leap-frog scheme of CReSS.
    allocate(sdn_tmp(1:1))
    allocate(sdrk_tmp(1:1))
    allocate(sdx_tmp(1:1))
    allocate(sdy_tmp(1:1))
    allocate(sdz_tmp(1:1))
    allocate(sdr_tmp(1:1))
    allocate(sdu_tmp(1:1))
    allocate(sdv_tmp(1:1))
    allocate(sdvz_tmp(1:1))
    allocate(sdasl_tmp(1:1,1:1))

    allocate(rand_s2c(1:sdnum_s2c))
    allocate(sortid_s2c(1:sdnum_s2c))
    allocate(sortkey_s2c(1:sdnum_s2c))
    allocate(sortfreq_s2c(1:ni_s2c*nj_s2c*nk_s2c+1))
    allocate(sorttag_s2c(1:ni_s2c*nj_s2c*nk_s2c+2))
    ! These CReSS related variables should not be used for SCALE. Remove these in the near future.
!!$    allocate(rhod_crs(KA,IA,JA))
!!$    allocate(rhoc_sdm(KA,IA,JA))
!!$    allocate(rhor_sdm(KA,IA,JA))
!!$    allocate(rhoa_sdm(KA,IA,JA))

    bufsiz = nint( sdininum_s2c )
    bufsiz = nint( real(bufsiz) * ( real(sdm_extbuf)*1.E-2_RP) )
    allocate(rbuf(1:bufsiz,bndsdmdim,1:2))
    allocate(sbuf(1:bufsiz,bndsdmdim,1:2))
    allocate(sdm_itmp1(1:ni_s2c*nj_s2c*nk_s2c+2))
    allocate(sdm_itmp2(1:ni_s2c*nj_s2c*nk_s2c+2))
    allocate(sdm_itmp3(1:ni_s2c*nj_s2c*nk_s2c+2))

    ! OpenMP not supported for SCLAE-LES-SDM. Rmove nomp in the near future.
    allocate(sd_itmp1(1:sdnum_s2c,1:nomp))
    allocate(sd_itmp2(1:sdnum_s2c,1:nomp))
    allocate(sd_itmp3(1:sdnum_s2c,1:nomp))

    allocate(sd_dtmp1(1:sdnum_s2c))

    allocate( prr_crs(IA,JA,1:2) )
!    allocate( j31(KA,IA,JA) )
!    allocate( j32(KA,IA,JA) )
!    allocate( jcb(KA,IA,JA) )
!    allocate( jcb8w(KA,IA,JA) )
!    allocate( mf(IA,JA) )
    allocate( dxiv_sdm(IA) )
    allocate( dyiv_sdm(JA) )
!    allocate( dziv_sdm(KA) )
    allocate( dx_sdm(IA) )
    allocate( dy_sdm(JA) )
!    allocate( dz_sdm(KA) )

    allocate(QTRC_sdm(KA,IA,JA,QA))

    ! Initialize allocated array
    QTRC_sdm(:,:,:,:) = 0.0_RP
!!$    rhod_crs(:,:,:) = 0.0_RP
!!$    rhoc_sdm(:,:,:) = 0.0_RP
!!$    rhor_sdm(:,:,:) = 0.0_RP
!!$    rhoa_sdm(:,:,:) = 0.0_RP
    do n = 1, sdnum_s2c
       sortid_s2c(n) = 0
       sortkey_s2c(n) = 0
       rand_s2c(n) = 0.0_RP
       sdn_s2c(n) = int(0,kind=DP)
       sdx_s2c(n) = 0.0_RP
       sdy_s2c(n) = 0.0_RP
       sdz_s2c(n) = 0.0_RP
       sdr_s2c(n) = 0.0_RP
       sdu_s2c(n) = 0.0_RP
       sdv_s2c(n) = 0.0_RP
       sdvz_s2c(n) = 0.0_RP

       sd_itmp1(n,1:nomp) = 0
       sd_itmp2(n,1:nomp) = 0
       sd_itmp3(n,1:nomp) = 0

       sd_dtmp1(n) = 0.0_RP
    enddo

    do s = 1, sdnumasl_s2c
       do n = 1, sdnum_s2c
          sdasl_s2c(n,s) = 0.0_RP
       enddo
    enddo

    if( abs(sdm_aslset) > 10 ) then
       do n = 1, sdfmnum_s2c
         sdn_fm(n) = int(0,kind=DP)
         sdx_fm(n) = 0.0_RP
         sdy_fm(n) = 0.0_RP
         sdz_fm(n) = 0.0_RP
         sdr_fm(n) = 0.0_RP
         sdvz_fm(n) = 0.0_RP
       enddo
       do s = 1, sdnumasl_s2c
            do n = 1, sdnum_s2c
              sdasl_fm(n,s) = 0.0_RP
            enddo
       enddo
    endif

    do n = 1, ni_s2c*nj_s2c*nk_s2c+1
       sortfreq_s2c(n) = 0
    enddo

    do n = 1, ni_s2c*nj_s2c*nk_s2c+2
       sorttag_s2c(n) = 0
       sdm_itmp1(n) = 0
       sdm_itmp2(n) = 0
       sdm_itmp3(n) = 0
    enddo

    do k = 1,2
    do j = 1, bndsdmdim
    do i = 1, bufsiz
      rbuf(i,j,k) = 0.0_RP
      sbuf(i,j,k) = 0.0_RP
    enddo
    enddo
    enddo

!!$    ! Check the evaluation of Jacobian later.
!!$    do k = 1, KA
!!$    do i = 1, IA
!!$    do j = 1, JA
!!$      j31(k,i,j) = GTRANS_J13G(k,i,j,I_XYZ) / GTRANS_GSQRT(k,i,j,I_XYZ)
!!$      j32(k,i,j) = GTRANS_J23G(k,i,j,I_XYZ) / GTRANS_GSQRT(k,i,j,I_XYZ)
!!$!      jcb(k,i,j) = GTRANS_GSQRT(k,i,j,I_XYZ)
!!$!      jcb8w(k,i,j) = GTRANS_GSQRT(k,i,j,I_XYW)
!!$    enddo
!!$    enddo
!!$    enddo

    return

   end subroutine sdm_allocinit
end module m_sdm_memmgr
