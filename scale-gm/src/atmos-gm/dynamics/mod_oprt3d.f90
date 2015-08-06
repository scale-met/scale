!-------------------------------------------------------------------------------
!> Module 3D Operator
!!
!! @par Description
!!          This module contains the subroutines for differential operators using vertical metrics.
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_oprt3d
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_adm, only: &
     ADM_LOG_FID,    &
     TI  => ADM_TI,  &
     TJ  => ADM_TJ,  &
     AI  => ADM_AI,  &
     AIJ => ADM_AIJ, &
     AJ  => ADM_AJ,  &
     K0  => ADM_KNONE
  use mod_gmtr, only: &
     P_RAREA => GMTR_P_RAREA, &
     T_RAREA => GMTR_T_RAREA, &
     W1      => GMTR_T_W1,    &
     W2      => GMTR_T_W2,    &
     W3      => GMTR_T_W3,    &
     HNX     => GMTR_A_HNX,   &
     HNY     => GMTR_A_HNY,   &
     HNZ     => GMTR_A_HNZ,   &
     HTX     => GMTR_A_HTX,   &
     HTY     => GMTR_A_HTY,   &
     HTZ     => GMTR_A_HTZ,   &
     TNX     => GMTR_A_TNX,   &
     TNY     => GMTR_A_TNY,   &
     TNZ     => GMTR_A_TNZ,   &
     TN2X    => GMTR_A_TN2X,  &
     TN2Y    => GMTR_A_TN2Y,  &
     TN2Z    => GMTR_A_TN2Z
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OPRT3D_divdamp

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OPRT3D_divdamp( &
       ddivdx, ddivdx_pl, &
       ddivdy, ddivdy_pl, &
       ddivdz, ddivdz_pl, &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl   )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl,  &
       ADM_kmin,     &
       ADM_kmax
    use mod_grd, only: &
       GRD_rdgz
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    use mod_oprt, only: &
       OPRT_nstart, &
       OPRT_nend,   &
       cinterp_TN,  &
       cinterp_HN,  &
       cinterp_TRA, &
       cinterp_PRA
    use mod_vmtr, only: &
       VMTR_RGAM,        &
       VMTR_RGAM_pl,     &
       VMTR_RGAMH,       &
       VMTR_RGAMH_pl,    &
       VMTR_RGSQRTH,     &
       VMTR_RGSQRTH_pl,  &
       VMTR_C2WfactGz,   &
       VMTR_C2WfactGz_pl
    implicit none

    real(RP), intent(out) :: ddivdx   (ADM_gall   ,ADM_kall,ADM_lall   ) ! tendency
    real(RP), intent(out) :: ddivdx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vx { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vy { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vz { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w  { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: sclt         (ADM_gall   ,TI:TJ) ! scalar on the hexagon vertex
    real(RP) :: sclt_pl      (ADM_gall_pl)
    real(RP) :: sclt_rhogw
    real(RP) :: sclt_rhogw_pl

    real(RP) :: rhogvx_vm   (ADM_gall   )          ! rho*vx / vertical metrics
    real(RP) :: rhogvx_vm_pl(ADM_gall_pl)
    real(RP) :: rhogvy_vm   (ADM_gall   )          ! rho*vy / vertical metrics
    real(RP) :: rhogvy_vm_pl(ADM_gall_pl)
    real(RP) :: rhogvz_vm   (ADM_gall   )          ! rho*vz / vertical metrics
    real(RP) :: rhogvz_vm_pl(ADM_gall_pl)
    real(RP) :: rhogw_vm    (ADM_gall,   ADM_kall) ! rho*w  / vertical metrics
    real(RP) :: rhogw_vm_pl (ADM_gall_pl,ADM_kall)

    integer :: nstart, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: gall, gall_1d, gmin, kmin, kmax
    integer :: g, k, l, v, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT3D_divdamp')

    gall    = ADM_gall
    gall_1d = ADM_gall_1d
    gmin    = ADM_gmin
    kmin    = ADM_kmin
    kmax    = ADM_kmax

    nstart = suf(ADM_gmin-1,ADM_gmin-1)
    nend   = suf(ADM_gmax,  ADM_gmax  )

    do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          rhogw_vm(g,k) = ( VMTR_C2WfactGz(g,k,1,l) * rhogvx(g,k  ,l) &
                          + VMTR_C2WfactGz(g,k,2,l) * rhogvx(g,k-1,l) &
                          + VMTR_C2WfactGz(g,k,3,l) * rhogvy(g,k  ,l) &
                          + VMTR_C2WfactGz(g,k,4,l) * rhogvy(g,k-1,l) &
                          + VMTR_C2WfactGz(g,k,5,l) * rhogvz(g,k  ,l) &
                          + VMTR_C2WfactGz(g,k,6,l) * rhogvz(g,k-1,l) &
                          ) * VMTR_RGAMH(g,k,l)                       & ! horizontal contribution
                        + rhogw(g,k,l) * VMTR_RGSQRTH(g,k,l)            ! vertical   contribution
       enddo
       enddo
       do g = 1, ADM_gall
          rhogw_vm(g,ADM_kmin  ) = 0.0_RP
          rhogw_vm(g,ADM_kmax+1) = 0.0_RP
       enddo

       !$omp parallel default(none),private(g,n,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1,sclt_rhogw), &
       !$omp shared(OPRT_nstart,OPRT_nend,gall,gall_1d,gmin,kmin,kmax,nstart,nend,l,ADM_have_sgp, &
       !$omp GRD_rdgz,VMTR_RGAM,rhogvx,rhogvy,rhogvz,ddivdx,ddivdy,ddivdz,sclt, &
       !$omp rhogvx_vm,rhogvy_vm,rhogvz_vm,rhogw_vm,cinterp_TN,cinterp_HN,cinterp_PRA,cinterp_TRA)
       do k = kmin, kmax

          !$omp do
          do g = 1, gall
             rhogvx_vm(g) = rhogvx(g,k,l) * VMTR_RGAM(g,k,l)
             rhogvy_vm(g) = rhogvy(g,k,l) * VMTR_RGAM(g,k,l)
             rhogvz_vm(g) = rhogvz(g,k,l) * VMTR_RGAM(g,k,l)
          enddo
          !$omp end do

          !$omp do
          do n = nstart, nend
             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + gall_1d

             sclt_rhogw = ( ( rhogw_vm(ij,k+1) + rhogw_vm(ip1j,k+1) + rhogw_vm(ip1jp1,k+1) ) &
                          - ( rhogw_vm(ij,k  ) + rhogw_vm(ip1j,k  ) + rhogw_vm(ip1jp1,k  ) ) &
                          ) / 3.0_RP * GRD_rdgz(k)

             sclt(n,TI) = ( - (rhogvx_vm(ij    )+rhogvx_vm(ip1j  )) * cinterp_TN(ij  ,l,AI ,1) &
                            - (rhogvx_vm(ip1j  )+rhogvx_vm(ip1jp1)) * cinterp_TN(ip1j,l,AJ ,1) &
                            + (rhogvx_vm(ip1jp1)+rhogvx_vm(ij    )) * cinterp_TN(ij  ,l,AIJ,1) &
                            - (rhogvy_vm(ij    )+rhogvy_vm(ip1j  )) * cinterp_TN(ij  ,l,AI ,2) &
                            - (rhogvy_vm(ip1j  )+rhogvy_vm(ip1jp1)) * cinterp_TN(ip1j,l,AJ ,2) &
                            + (rhogvy_vm(ip1jp1)+rhogvy_vm(ij    )) * cinterp_TN(ij  ,l,AIJ,2) &
                            - (rhogvz_vm(ij    )+rhogvz_vm(ip1j  )) * cinterp_TN(ij  ,l,AI ,3) &
                            - (rhogvz_vm(ip1j  )+rhogvz_vm(ip1jp1)) * cinterp_TN(ip1j,l,AJ ,3) &
                            + (rhogvz_vm(ip1jp1)+rhogvz_vm(ij    )) * cinterp_TN(ij  ,l,AIJ,3) &
                          ) * 0.5_RP * cinterp_TRA(ij,l,TI) &
                        + sclt_rhogw
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart, nend
             ij     = n
             ijp1   = n     + gall_1d
             ip1jp1 = n + 1 + gall_1d

             sclt_rhogw = ( ( rhogw_vm(ij,k+1) + rhogw_vm(ijp1,k+1) + rhogw_vm(ip1jp1,k+1) ) &
                          - ( rhogw_vm(ij,k  ) + rhogw_vm(ijp1,k  ) + rhogw_vm(ip1jp1,k  ) ) &
                          ) / 3.0_RP * GRD_rdgz(k)

             sclt(n,TJ) = ( - (rhogvx_vm(ij    )+rhogvx_vm(ip1jp1)) * cinterp_TN(ij  ,l,AIJ,1) &
                            + (rhogvx_vm(ip1jp1)+rhogvx_vm(ijp1  )) * cinterp_TN(ijp1,l,AI ,1) &
                            + (rhogvx_vm(ijp1  )+rhogvx_vm(ij    )) * cinterp_TN(ij  ,l,AJ ,1) &
                            - (rhogvy_vm(ij    )+rhogvy_vm(ip1jp1)) * cinterp_TN(ij  ,l,AIJ,2) &
                            + (rhogvy_vm(ip1jp1)+rhogvy_vm(ijp1  )) * cinterp_TN(ijp1,l,AI ,2) &
                            + (rhogvy_vm(ijp1  )+rhogvy_vm(ij    )) * cinterp_TN(ij  ,l,AJ ,2) &
                            - (rhogvz_vm(ij    )+rhogvz_vm(ip1jp1)) * cinterp_TN(ij  ,l,AIJ,3) &
                            + (rhogvz_vm(ip1jp1)+rhogvz_vm(ijp1  )) * cinterp_TN(ijp1,l,AI ,3) &
                            + (rhogvz_vm(ijp1  )+rhogvz_vm(ij    )) * cinterp_TN(ij  ,l,AJ ,3) &
                          ) * 0.5_RP * cinterp_TRA(ij,l,TJ) &
                        + sclt_rhogw
          enddo
          !$omp end do

          !$omp do
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             im1j   = n - 1
             ijm1   = n     - gall_1d
             im1jm1 = n - 1 - gall_1d

             ddivdx(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(ij,    l,AI ,1) &
                               + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(ij,    l,AIJ,1) &
                               + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(ij,    l,AJ ,1) &
                               - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(im1j,  l,AI ,1) &
                               - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(im1jm1,l,AIJ,1) &
                               - ( sclt(ijm1  ,TJ) + sclt(im1jm1,TI) ) * cinterp_HN(ijm1,  l,AJ ,1) &
                             ) * 0.5_RP * cinterp_PRA(ij,l)

             ddivdy(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(ij,    l,AI ,2) &
                               + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(ij,    l,AIJ,2) &
                               + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(ij,    l,AJ ,2) &
                               - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(im1j,  l,AI ,2) &
                               - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(im1jm1,l,AIJ,2) &
                               - ( sclt(ijm1  ,TJ) + sclt(im1jm1,TI) ) * cinterp_HN(ijm1,  l,AJ ,2) &
                             ) * 0.5_RP * cinterp_PRA(ij,l)

             ddivdz(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(ij,    l,AI ,3) &
                               + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(ij,    l,AIJ,3) &
                               + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(ij,    l,AJ ,3) &
                               - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(im1j,  l,AI ,3) &
                               - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(im1jm1,l,AIJ,3) &
                               - ( sclt(ijm1  ,TJ) + sclt(im1jm1,TI) ) * cinterp_HN(ijm1,  l,AJ ,3) &
                             ) * 0.5_RP * cinterp_PRA(ij,l)
          enddo

          if ( ADM_have_sgp(l) ) then
             !$omp master
             n = suf(gmin,gmin)

             ij     = n
             im1j   = n - 1
             ijm1   = n     - gall_1d
             im1jm1 = n - 1 - gall_1d

             sclt(im1jm1,TI) = sclt(ijm1,TJ) ! copy

             ddivdx(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(ij,    l,AI ,1) &
                               + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(ij,    l,AIJ,1) &
                               + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(ij,    l,AJ ,1) &
                               - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(im1j,  l,AI ,1) &
                               - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(im1jm1,l,AIJ,1) &
                             ) * 0.5_RP * cinterp_PRA(ij,l)

             ddivdy(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(ij,    l,AI ,2) &
                               + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(ij,    l,AIJ,2) &
                               + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(ij,    l,AJ ,2) &
                               - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(im1j,  l,AI ,2) &
                               - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(im1jm1,l,AIJ,2) &
                             ) * 0.5_RP * cinterp_PRA(ij,l)

             ddivdz(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(ij,    l,AI ,3) &
                               + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(ij,    l,AIJ,3) &
                               + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(ij,    l,AJ ,3) &
                               - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(im1j,  l,AI ,3) &
                               - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(im1jm1,l,AIJ,3) &
                             ) * 0.5_RP * cinterp_PRA(ij,l)
             !$omp end master
             !$omp barrier
          endif

       enddo
       !$omp end parallel

       do g = 1, ADM_gall
          ddivdx(g,ADM_kmin-1,l) = 0.0_RP
          ddivdx(g,ADM_kmax+1,l) = 0.0_RP
          ddivdy(g,ADM_kmin-1,l) = 0.0_RP
          ddivdy(g,ADM_kmax+1,l) = 0.0_RP
          ddivdz(g,ADM_kmin-1,l) = 0.0_RP
          ddivdz(g,ADM_kmax+1,l) = 0.0_RP
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             rhogw_vm_pl(g,k) = ( VMTR_C2WfactGz_pl(g,k,1,l) * rhogvx_pl(g,k  ,l) &
                                + VMTR_C2WfactGz_pl(g,k,2,l) * rhogvx_pl(g,k-1,l) &
                                + VMTR_C2WfactGz_pl(g,k,3,l) * rhogvy_pl(g,k  ,l) &
                                + VMTR_C2WfactGz_pl(g,k,4,l) * rhogvy_pl(g,k-1,l) &
                                + VMTR_C2WfactGz_pl(g,k,5,l) * rhogvz_pl(g,k  ,l) &
                                + VMTR_C2WfactGz_pl(g,k,6,l) * rhogvz_pl(g,k-1,l) &
                                ) * VMTR_RGAMH_pl(g,k,l)                          & ! horizontal contribution
                              + rhogw_pl(g,k,l) * VMTR_RGSQRTH_pl(g,k,l)            ! vertical   contribution
          enddo
          enddo
          do g = 1, ADM_gall_pl
             rhogw_vm_pl(g,ADM_kmin  ) = 0.0_RP
             rhogw_vm_pl(g,ADM_kmax+1) = 0.0_RP
          enddo

          n = ADM_GSLF_PL

          do k = ADM_kmin, ADM_kmax
             do v = 1, ADM_gall_pl
                rhogvx_vm_pl(v) = rhogvx_pl(v,k,l) * VMTR_RGAM_pl(v,k,l)
                rhogvy_vm_pl(v) = rhogvy_pl(v,k,l) * VMTR_RGAM_pl(v,k,l)
                rhogvz_vm_pl(v) = rhogvz_pl(v,k,l) * VMTR_RGAM_pl(v,k,l)
             enddo

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                sclt_rhogw_pl = ( ( rhogw_vm_pl(n,k+1) + rhogw_vm_pl(ij,k+1) + rhogw_vm_pl(ijp1,k+1) ) &
                                - ( rhogw_vm_pl(n,k  ) + rhogw_vm_pl(ij,k  ) + rhogw_vm_pl(ijp1,k  ) ) &
                                ) / 3.0_RP * GRD_rdgz(k)

                sclt_pl(v) = ( + ( rhogvx_vm_pl(n   ) + rhogvx_vm_pl(ij  ) ) * GMTR_A_var_pl(ij,  k0,l,TNX ) &
                               + ( rhogvy_vm_pl(n   ) + rhogvy_vm_pl(ij  ) ) * GMTR_A_var_pl(ij,  k0,l,TNY ) &
                               + ( rhogvz_vm_pl(n   ) + rhogvz_vm_pl(ij  ) ) * GMTR_A_var_pl(ij,  k0,l,TNZ ) &
                               + ( rhogvx_vm_pl(ij  ) + rhogvx_vm_pl(ijp1) ) * GMTR_A_var_pl(ij,  k0,l,TN2X) &
                               + ( rhogvy_vm_pl(ij  ) + rhogvy_vm_pl(ijp1) ) * GMTR_A_var_pl(ij,  k0,l,TN2Y) &
                               + ( rhogvz_vm_pl(ij  ) + rhogvz_vm_pl(ijp1) ) * GMTR_A_var_pl(ij,  k0,l,TN2Z) &
                               - ( rhogvx_vm_pl(ijp1) + rhogvx_vm_pl(n   ) ) * GMTR_A_var_pl(ijp1,k0,l,TNX ) &
                               - ( rhogvy_vm_pl(ijp1) + rhogvy_vm_pl(n   ) ) * GMTR_A_var_pl(ijp1,k0,l,TNY ) &
                               - ( rhogvz_vm_pl(ijp1) + rhogvz_vm_pl(n   ) ) * GMTR_A_var_pl(ijp1,k0,l,TNZ ) &
                             ) * 0.5_RP * GMTR_T_var_pl(ij,k0,l,T_RAREA) &
                           + sclt_rhogw_pl
             enddo

             ddivdx_pl(n,k,l) = 0.0_RP
             ddivdy_pl(n,k,l) = 0.0_RP
             ddivdz_pl(n,k,l) = 0.0_RP

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijm1 = v - 1
                if( ijm1 < ADM_gmin_pl ) ijm1 = ADM_gmax_pl ! cyclic condition

                ddivdx_pl(n,k,l) = ddivdx_pl(n,k,l) + ( sclt_pl(ijm1) + sclt_pl(ij) ) * GMTR_A_var_pl(ij,k0,l,HNX)
                ddivdy_pl(n,k,l) = ddivdy_pl(n,k,l) + ( sclt_pl(ijm1) + sclt_pl(ij) ) * GMTR_A_var_pl(ij,k0,l,HNY)
                ddivdz_pl(n,k,l) = ddivdz_pl(n,k,l) + ( sclt_pl(ijm1) + sclt_pl(ij) ) * GMTR_A_var_pl(ij,k0,l,HNZ)
             enddo

             ddivdx_pl(n,k,l) = ddivdx_pl(n,k,l) * 0.5_RP * GMTR_P_var_pl(n,k0,l,P_RAREA)
             ddivdy_pl(n,k,l) = ddivdy_pl(n,k,l) * 0.5_RP * GMTR_P_var_pl(n,k0,l,P_RAREA)
             ddivdz_pl(n,k,l) = ddivdz_pl(n,k,l) * 0.5_RP * GMTR_P_var_pl(n,k0,l,P_RAREA)
          enddo

       enddo
    else
       ddivdx_pl(:,:,:) = 0.0_RP
       ddivdy_pl(:,:,:) = 0.0_RP
       ddivdz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT3D_divdamp')

    return
  end subroutine OPRT3D_divdamp

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_oprt3d
