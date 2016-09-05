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
       ddivdx,    ddivdx_pl,    &
       ddivdy,    ddivdy_pl,    &
       ddivdz,    ddivdz_pl,    &
       rhogvx,    rhogvx_pl,    &
       rhogvy,    rhogvy_pl,    &
       rhogvz,    rhogvz_pl,    &
       rhogw,     rhogw_pl,     &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl, &
       RGSQRTH,   RGSQRTH_pl,   &
       RGAM,      RGAM_pl,      &
       RGAMH,     RGAMH_pl,     &
       C2WfactGz, C2WfactGz_pl  )
    use mod_adm, only: &
       TI  => ADM_TI,  &
       TJ  => ADM_TJ,  &
       ADM_nxyz,       &
       ADM_have_pl,    &
       ADM_have_sgp,   &
       ADM_vlink,      &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_kall,       &
       ADM_jall,       &
       ADM_iall,       &
       ADM_gall_pl,    &
       ADM_jmin,       &
       ADM_jmax,       &
       ADM_imin,       &
       ADM_imax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl,    &
       ADM_kmin,       &
       ADM_kmax
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR, &
       GRD_rdgz
    implicit none

    real(RP), intent(out) :: ddivdx      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! tendency
    real(RP), intent(out) :: ddivdx_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdy      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdy_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdz      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdz_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! rho*vx { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvx_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! rho*vy { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvy_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! rho*vz { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvz_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw       (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! rho*w  { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogw_pl    (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
!    real(RP), intent(in)  :: coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   )
!    real(RP), intent(in)  :: coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl)
!    real(RP), intent(in)  :: coef_diff   (ADM_nxyz,ADM_gall,1:6        ,ADM_lall   )
!    real(RP), intent(in)  :: coef_diff_pl(ADM_nxyz,         1:ADM_vlink,ADM_lall_pl)
    real(RP), intent(in)  :: coef_intp   (ADM_iall,ADM_jall,1:3        ,ADM_nxyz,TI:TJ,ADM_lall   )
    real(RP), intent(in)  :: coef_intp_pl(ADM_gall_pl      ,1:3        ,ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(in)  :: coef_diff   (ADM_iall,ADM_jall,1:6        ,ADM_nxyz,      ADM_lall   )
    real(RP), intent(in)  :: coef_diff_pl(                  1:ADM_vlink,ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(in)  :: RGSQRTH     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: RGSQRTH_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: RGAM        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: RGAM_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: RGAMH       (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: RGAMH_pl    (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: C2WfactGz   (ADM_iall,ADM_jall,ADM_kall,6,ADM_lall   )
    real(RP), intent(in)  :: C2WfactGz_pl(ADM_gall_pl      ,ADM_kall,6,ADM_lall_pl)

    real(RP) :: sclt   (ADM_iall,ADM_jall,TI:TJ) ! scalar on the hexagon vertex
    real(RP) :: sclt_pl(ADM_gall_pl      )
    real(RP) :: sclt_rhogw
    real(RP) :: sclt_rhogw_pl

    real(RP) :: rhogvx_vm   (ADM_iall,ADM_jall)          ! rho*vx / vertical metrics
    real(RP) :: rhogvx_vm_pl(ADM_gall_pl      )
    real(RP) :: rhogvy_vm   (ADM_iall,ADM_jall)          ! rho*vy / vertical metrics
    real(RP) :: rhogvy_vm_pl(ADM_gall_pl      )
    real(RP) :: rhogvz_vm   (ADM_iall,ADM_jall)          ! rho*vz / vertical metrics
    real(RP) :: rhogvz_vm_pl(ADM_gall_pl      )
    real(RP) :: rhogw_vm    (ADM_iall,ADM_jall,ADM_kall) ! rho*w  / vertical metrics
    real(RP) :: rhogw_vm_pl (ADM_gall_pl      ,ADM_kall)

    integer  :: imin, imax, jmin, jmax, kall, kmin, kmax, lall
    integer  :: ij, ijp1, ijm1

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT3D_divdamp',2)

    imin = ADM_imin
    imax = ADM_imax
    jmin = ADM_jmin
    jmax = ADM_jmax
    kall = ADM_kall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    !$omp parallel default(none),private(i,j,k,l,sclt_rhogw), &
    !$omp shared(imin,imax,jmin,jmax,kall,kmin,kmax,lall,ADM_have_sgp,GRD_rdgz, &
    !$omp ddivdx,ddivdy,ddivdz,rhogvx,rhogvy,rhogvz,rhogw,sclt,coef_intp,coef_diff, &
    !$omp rhogvx_vm,rhogvy_vm,rhogvz_vm,rhogw_vm,C2WfactGz,RGAMH,RGSQRTH,RGAM)
    do l = 1, lall

       !$omp do schedule(static)
       do k = kmin+1, kmax
       do j = jmin-1, jmax+1
       do i = imin-1, imax+1
          rhogw_vm(i,j,k) = ( C2WfactGz(i,j,k,1,l) * rhogvx(i,j,k  ,l) &
                            + C2WfactGz(i,j,k,2,l) * rhogvx(i,j,k-1,l) &
                            + C2WfactGz(i,j,k,3,l) * rhogvy(i,j,k  ,l) &
                            + C2WfactGz(i,j,k,4,l) * rhogvy(i,j,k-1,l) &
                            + C2WfactGz(i,j,k,5,l) * rhogvz(i,j,k  ,l) &
                            + C2WfactGz(i,j,k,6,l) * rhogvz(i,j,k-1,l) &
                            ) * RGAMH(i,j,k,l)                         & ! horizontal contribution
                          + rhogw(i,j,k,l) * RGSQRTH(i,j,k,l)            ! vertical   contribution
       enddo
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin-1, jmax+1
       do i = imin-1, imax+1
          rhogw_vm(i,j,kmin  ) = 0.0_RP
          rhogw_vm(i,j,kmax+1) = 0.0_RP
       enddo
       enddo
       !$omp end do

       do k = kmin, kmax

          !$omp do schedule(static)
          do j = jmin-1, jmax+1
          do i = imin-1, imax+1
             rhogvx_vm(i,j) = rhogvx(i,j,k,l) * RGAM(i,j,k,l)
             rhogvy_vm(i,j) = rhogvy(i,j,k,l) * RGAM(i,j,k,l)
             rhogvz_vm(i,j) = rhogvz(i,j,k,l) * RGAM(i,j,k,l)
          enddo
          enddo
          !$omp end do

          !$omp do schedule(static)
          do j = jmin-1, jmax
          do i = imin-1, imax
             sclt_rhogw = ( ( rhogw_vm(i,j,k+1) + rhogw_vm(i+1,j,k+1) + rhogw_vm(i+1,j+1,k+1) ) &
                          - ( rhogw_vm(i,j,k  ) + rhogw_vm(i+1,j,k  ) + rhogw_vm(i+1,j+1,k  ) ) &
                          ) / 3.0_RP * GRD_rdgz(k)

             sclt(i,j,TI) = coef_intp(i,j,1,XDIR,TI,l) * rhogvx_vm(i  ,j  ) &
                          + coef_intp(i,j,2,XDIR,TI,l) * rhogvx_vm(i+1,j  ) &
                          + coef_intp(i,j,3,XDIR,TI,l) * rhogvx_vm(i+1,j+1) &
                          + coef_intp(i,j,1,YDIR,TI,l) * rhogvy_vm(i  ,j  ) &
                          + coef_intp(i,j,2,YDIR,TI,l) * rhogvy_vm(i+1,j  ) &
                          + coef_intp(i,j,3,YDIR,TI,l) * rhogvy_vm(i+1,j+1) &
                          + coef_intp(i,j,1,ZDIR,TI,l) * rhogvz_vm(i  ,j  ) &
                          + coef_intp(i,j,2,ZDIR,TI,l) * rhogvz_vm(i+1,j  ) &
                          + coef_intp(i,j,3,ZDIR,TI,l) * rhogvz_vm(i+1,j+1) &
                          + sclt_rhogw
          enddo
          enddo
          !$omp end do nowait

          !$omp do schedule(static)
          do j = jmin-1, jmax
          do i = imin-1, imax
             sclt_rhogw = ( ( rhogw_vm(i,j,k+1) + rhogw_vm(i+1,j+1,k+1) + rhogw_vm(i,j+1,k+1) ) &
                          - ( rhogw_vm(i,j,k  ) + rhogw_vm(i+1,j+1,k  ) + rhogw_vm(i,j+1,k  ) ) &
                          ) / 3.0_RP * GRD_rdgz(k)

             sclt(i,j,TJ) = coef_intp(i,j,1,XDIR,TJ,l) * rhogvx_vm(i  ,j  ) &
                          + coef_intp(i,j,2,XDIR,TJ,l) * rhogvx_vm(i+1,j+1) &
                          + coef_intp(i,j,3,XDIR,TJ,l) * rhogvx_vm(i  ,j+1) &
                          + coef_intp(i,j,1,YDIR,TJ,l) * rhogvy_vm(i  ,j  ) &
                          + coef_intp(i,j,2,YDIR,TJ,l) * rhogvy_vm(i+1,j+1) &
                          + coef_intp(i,j,3,YDIR,TJ,l) * rhogvy_vm(i  ,j+1) &
                          + coef_intp(i,j,1,ZDIR,TJ,l) * rhogvz_vm(i  ,j  ) &
                          + coef_intp(i,j,2,ZDIR,TJ,l) * rhogvz_vm(i+1,j+1) &
                          + coef_intp(i,j,3,ZDIR,TJ,l) * rhogvz_vm(i  ,j+1) &
                          + sclt_rhogw
          enddo
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then ! pentagon
             sclt(imin-1,jmin-1,TI) = sclt(imin,jmin-1,TJ)
          endif

          !$omp do schedule(static)
          do j = jmin, jmax
          do i = imin, imax
             ddivdx(i,j,k,l) = coef_diff(i,j,1,XDIR,l) * ( sclt(i  ,j  ,TI) + sclt(i  ,j  ,TJ) ) &
                             + coef_diff(i,j,2,XDIR,l) * ( sclt(i  ,j  ,TJ) + sclt(i-1,j  ,TI) ) &
                             + coef_diff(i,j,3,XDIR,l) * ( sclt(i-1,j  ,TI) + sclt(i-1,j-1,TJ) ) &
                             + coef_diff(i,j,4,XDIR,l) * ( sclt(i-1,j-1,TJ) + sclt(i-1,j-1,TI) ) &
                             + coef_diff(i,j,5,XDIR,l) * ( sclt(i-1,j-1,TI) + sclt(i  ,j-1,TJ) ) &
                             + coef_diff(i,j,6,XDIR,l) * ( sclt(i  ,j-1,TJ) + sclt(i  ,j  ,TI) )
          enddo
          enddo
          !$omp end do nowait

          !$omp do schedule(static)
          do j = jmin, jmax
          do i = imin, imax
             ddivdy(i,j,k,l) = coef_diff(i,j,1,YDIR,l) * ( sclt(i  ,j  ,TI) + sclt(i  ,j  ,TJ) ) &
                             + coef_diff(i,j,2,YDIR,l) * ( sclt(i  ,j  ,TJ) + sclt(i-1,j  ,TI) ) &
                             + coef_diff(i,j,3,YDIR,l) * ( sclt(i-1,j  ,TI) + sclt(i-1,j-1,TJ) ) &
                             + coef_diff(i,j,4,YDIR,l) * ( sclt(i-1,j-1,TJ) + sclt(i-1,j-1,TI) ) &
                             + coef_diff(i,j,5,YDIR,l) * ( sclt(i-1,j-1,TI) + sclt(i  ,j-1,TJ) ) &
                             + coef_diff(i,j,6,YDIR,l) * ( sclt(i  ,j-1,TJ) + sclt(i  ,j  ,TI) )
          enddo
          enddo
          !$omp end do nowait

          !$omp do schedule(static)
          do j = jmin, jmax
          do i = imin, imax
             ddivdz(i,j,k,l) = coef_diff(i,j,1,ZDIR,l) * ( sclt(i  ,j  ,TI) + sclt(i  ,j  ,TJ) ) &
                             + coef_diff(i,j,2,ZDIR,l) * ( sclt(i  ,j  ,TJ) + sclt(i-1,j  ,TI) ) &
                             + coef_diff(i,j,3,ZDIR,l) * ( sclt(i-1,j  ,TI) + sclt(i-1,j-1,TJ) ) &
                             + coef_diff(i,j,4,ZDIR,l) * ( sclt(i-1,j-1,TJ) + sclt(i-1,j-1,TI) ) &
                             + coef_diff(i,j,5,ZDIR,l) * ( sclt(i-1,j-1,TI) + sclt(i  ,j-1,TJ) ) &
                             + coef_diff(i,j,6,ZDIR,l) * ( sclt(i  ,j-1,TJ) + sclt(i  ,j  ,TI) )
          enddo
          enddo
          !$omp end do nowait

          !$omp workshare
          ddivdx(:,jmin-1,k,l) = 0.0_RP
          ddivdy(:,jmin-1,k,l) = 0.0_RP
          ddivdz(:,jmin-1,k,l) = 0.0_RP
          ddivdx(:,jmax+1,k,l) = 0.0_RP
          ddivdy(:,jmax+1,k,l) = 0.0_RP
          ddivdz(:,jmax+1,k,l) = 0.0_RP
          ddivdx(imin-1,:,k,l) = 0.0_RP
          ddivdy(imin-1,:,k,l) = 0.0_RP
          ddivdz(imin-1,:,k,l) = 0.0_RP
          ddivdx(imax+1,:,k,l) = 0.0_RP
          ddivdy(imax+1,:,k,l) = 0.0_RP
          ddivdz(imax+1,:,k,l) = 0.0_RP
          !$omp end workshare
       enddo

       !$omp workshare
       ddivdx(:,:,kmin-1,l) = 0.0_RP
       ddivdx(:,:,kmax+1,l) = 0.0_RP
       ddivdy(:,:,kmin-1,l) = 0.0_RP
       ddivdy(:,:,kmax+1,l) = 0.0_RP
       ddivdz(:,:,kmin-1,l) = 0.0_RP
       ddivdz(:,:,kmax+1,l) = 0.0_RP
       !$omp end workshare
    enddo ! l loop
    !$omp end parallel

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do ij = 1, ADM_gall_pl
             rhogw_vm_pl(ij,k) = ( C2WfactGz_pl(ij,k,1,l) * rhogvx_pl(ij,k  ,l) &
                                 + C2WfactGz_pl(ij,k,2,l) * rhogvx_pl(ij,k-1,l) &
                                 + C2WfactGz_pl(ij,k,3,l) * rhogvy_pl(ij,k  ,l) &
                                 + C2WfactGz_pl(ij,k,4,l) * rhogvy_pl(ij,k-1,l) &
                                 + C2WfactGz_pl(ij,k,5,l) * rhogvz_pl(ij,k  ,l) &
                                 + C2WfactGz_pl(ij,k,6,l) * rhogvz_pl(ij,k-1,l) &
                                 ) * RGAMH_pl(ij,k,l)                           & ! horizontal contribution
                               + rhogw_pl(ij,k,l) * RGSQRTH_pl(ij,k,l)            ! vertical   contribution
          enddo
          enddo
          do ij = 1, ADM_gall_pl
             rhogw_vm_pl(ij,ADM_kmin  ) = 0.0_RP
             rhogw_vm_pl(ij,ADM_kmax+1) = 0.0_RP
          enddo

          n = ADM_gslf_pl

          do k = ADM_kmin, ADM_kmax
             do v = 1, ADM_gall_pl
                rhogvx_vm_pl(v) = rhogvx_pl(v,k,l) * RGAM_pl(v,k,l)
                rhogvy_vm_pl(v) = rhogvy_pl(v,k,l) * RGAM_pl(v,k,l)
                rhogvz_vm_pl(v) = rhogvz_pl(v,k,l) * RGAM_pl(v,k,l)
             enddo

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

                sclt_rhogw_pl = ( ( rhogw_vm_pl(n,k+1) + rhogw_vm_pl(ij,k+1) + rhogw_vm_pl(ijp1,k+1) ) &
                                - ( rhogw_vm_pl(n,k  ) + rhogw_vm_pl(ij,k  ) + rhogw_vm_pl(ijp1,k  ) ) &
                                ) / 3.0_RP * GRD_rdgz(k)

                sclt_pl(ij) = coef_intp_pl(v,1,XDIR,l) * rhogvx_vm_pl(n   ) &
                            + coef_intp_pl(v,2,XDIR,l) * rhogvx_vm_pl(ij  ) &
                            + coef_intp_pl(v,3,XDIR,l) * rhogvx_vm_pl(ijp1) &
                            + coef_intp_pl(v,1,YDIR,l) * rhogvy_vm_pl(n   ) &
                            + coef_intp_pl(v,2,YDIR,l) * rhogvy_vm_pl(ij  ) &
                            + coef_intp_pl(v,3,YDIR,l) * rhogvy_vm_pl(ijp1) &
                            + coef_intp_pl(v,1,ZDIR,l) * rhogvz_vm_pl(n   ) &
                            + coef_intp_pl(v,2,ZDIR,l) * rhogvz_vm_pl(ij  ) &
                            + coef_intp_pl(v,3,ZDIR,l) * rhogvz_vm_pl(ijp1) &
                            + sclt_rhogw_pl
             enddo

             ddivdx_pl(:,k,l) = 0.0_RP
             ddivdy_pl(:,k,l) = 0.0_RP
             ddivdz_pl(:,k,l) = 0.0_RP

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijm1 = v - 1
                if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl ! cyclic condition

                ddivdx_pl(n,k,l) = ddivdx_pl(n,k,l) + coef_diff_pl(v-1,XDIR,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
                ddivdy_pl(n,k,l) = ddivdy_pl(n,k,l) + coef_diff_pl(v-1,YDIR,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
                ddivdz_pl(n,k,l) = ddivdz_pl(n,k,l) + coef_diff_pl(v-1,ZDIR,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
             enddo
          enddo

          do ij = 1, ADM_gall_pl
             ddivdx_pl(ij,ADM_kmin-1,l) = 0.0_RP
             ddivdx_pl(ij,ADM_kmax+1,l) = 0.0_RP
             ddivdy_pl(ij,ADM_kmin-1,l) = 0.0_RP
             ddivdy_pl(ij,ADM_kmax+1,l) = 0.0_RP
             ddivdz_pl(ij,ADM_kmin-1,l) = 0.0_RP
             ddivdz_pl(ij,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo
    else
       ddivdx_pl(:,:,:) = 0.0_RP
       ddivdy_pl(:,:,:) = 0.0_RP
       ddivdz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT3D_divdamp',2)

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
