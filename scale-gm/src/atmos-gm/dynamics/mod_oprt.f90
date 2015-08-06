!-------------------------------------------------------------------------------
!> Module operator
!!
!! @par Description
!!          This module contains the subroutines for differential operators
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_oprt
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_adm, only: &
     ADM_LOG_FID,      &
     TI  => ADM_TI,    &
     TJ  => ADM_TJ,    &
     AI  => ADM_AI,    &
     AIJ => ADM_AIJ,   &
     AJ  => ADM_AJ,    &
     K0  => ADM_KNONE, &
     ADM_nxyz,         &
     ADM_lall,         &
     ADM_lall_pl,      &
     ADM_gall,         &
     ADM_gall_pl,      &
     ADM_kall
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
  public :: OPRT_setup
  public :: OPRT_divergence
  public :: OPRT_gradient
  public :: OPRT_laplacian
  public :: OPRT_diffusion
  public :: OPRT_horizontalize_vec
  public :: OPRT_vorticity
  public :: OPRT_divdamp

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: OPRT_nstart
  integer, public :: OPRT_nend

#ifdef _FIXEDINDEX_
  real(RP), public              :: cdiv       (ADM_gall   ,ADM_lall   ,0:6,            ADM_nxyz)
  real(RP), public              :: cdiv_pl    (ADM_gall_pl,ADM_lall_pl,0:ADM_gall_pl-1,ADM_nxyz)
  real(RP), public              :: cgrad      (ADM_gall   ,ADM_lall   ,0:6            ,ADM_nxyz)
  real(RP), public              :: cgrad_pl   (ADM_gall_pl,ADM_lall_pl,0:ADM_gall_pl-1,ADM_nxyz)
  real(RP), public              :: clap       (ADM_gall   ,ADM_lall   ,0:6                     )
  real(RP), public              :: clap_pl    (ADM_gall_pl,ADM_lall_pl,0:ADM_gall_pl-1         )

  real(RP), public              :: cinterp_TN (ADM_gall,ADM_lall,AI:AJ,ADM_nxyz)
  real(RP), public              :: cinterp_HN (ADM_gall,ADM_lall,AI:AJ,ADM_nxyz)
  real(RP), public              :: cinterp_TRA(ADM_gall,ADM_lall,TI:TJ         )
  real(RP), public              :: cinterp_PRA(ADM_gall,ADM_lall               )
#else
  real(RP), public, allocatable :: cdiv       (:,:,:,:) ! coefficient for divergence operator
  real(RP), public, allocatable :: cdiv_pl    (:,:,:,:)
  real(RP), public, allocatable :: cgrad      (:,:,:,:) ! coefficient for gradient operator
  real(RP), public, allocatable :: cgrad_pl   (:,:,:,:)
  real(RP), public, allocatable :: clap       (:,:,:)   ! coefficient for laplacian operator
  real(RP), public, allocatable :: clap_pl    (:,:,:)

  real(RP), public, allocatable :: cinterp_TN (:,:,:,:) ! coefficient for diffusion operator
  real(RP), public, allocatable :: cinterp_HN (:,:,:,:)
  real(RP), public, allocatable :: cinterp_TRA(:,:,:)
  real(RP), public, allocatable :: cinterp_PRA(:,:)
#endif

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
  subroutine OPRT_setup
    use mod_adm, only: &
       ADM_have_pl,    &
       ADM_have_sgp,   &
       ADM_gall_1d,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl
    use mod_gmtr, only: &
       GMTR_P_var,    &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var,    &
       GMTR_A_var_pl
    implicit none

    integer :: n0,n1,n2,n3,n4

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, l, m, md, v
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[oprt]/Category[common share]'

    ! dummy call
    call PROF_rapstart('OPRT_divergence')
    call PROF_rapend  ('OPRT_divergence')
    call PROF_rapstart('OPRT_gradient')
    call PROF_rapend  ('OPRT_gradient')
    call PROF_rapstart('OPRT_laplacian')
    call PROF_rapend  ('OPRT_laplacian')
    call PROF_rapstart('OPRT_diffusion')
    call PROF_rapend  ('OPRT_diffusion')
    call PROF_rapstart('OPRT_horizontalize_vec')
    call PROF_rapend  ('OPRT_horizontalize_vec')
    call PROF_rapstart('OPRT_divdamp')
    call PROF_rapend  ('OPRT_divdamp')
    call PROF_rapstart('OPRT3D_divdamp')
    call PROF_rapend  ('OPRT3D_divdamp')

    OPRT_nstart = suf(ADM_gmin,ADM_gmin)
    OPRT_nend   = suf(ADM_gmax,ADM_gmax)

    !---< setup coefficient of divergence operator >
    write(ADM_LOG_FID,*) '*** setup coefficient of divergence operator'

#ifndef _FIXEDINDEX_
    allocate( cdiv       (ADM_gall   ,ADM_lall   ,0:6,            ADM_nxyz) )
    allocate( cdiv_pl    (ADM_gall_pl,ADM_lall_pl,0:ADM_gall_pl-1,ADM_nxyz) )
    allocate( cgrad      (ADM_gall   ,ADM_lall   ,0:6            ,ADM_nxyz) )
    allocate( cgrad_pl   (ADM_gall_pl,ADM_lall_pl,0:ADM_gall_pl-1,ADM_nxyz) )
    allocate( clap       (ADM_gall   ,ADM_lall   ,0:6                     ) )
    allocate( clap_pl    (ADM_gall_pl,ADM_lall_pl,0:ADM_gall_pl-1         ) )

    allocate( cinterp_TN (ADM_gall,ADM_lall,AI:AJ,ADM_nxyz) )
    allocate( cinterp_HN (ADM_gall,ADM_lall,AI:AJ,ADM_nxyz) )
    allocate( cinterp_TRA(ADM_gall,ADM_lall,TI:TJ         ) )
    allocate( cinterp_PRA(ADM_gall,ADM_lall               ) )
#endif

    do l = 1, ADM_lall

       do m = 1, 3
          md = m + HNX - 1
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             ! ij
             cdiv(n,l,0,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                             ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1j
             cdiv(n,l,1,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1jp1
             cdiv(n,l,2,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijp1
             cdiv(n,l,3,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5_RP*GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1j
             cdiv(n,l,4,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1jm1
             cdiv(n,l,5,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                             ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijm1
             cdiv(n,l,6,m) = ( - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                             ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
          enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          do m = 1, 3
             md = m + HNX - 1

             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             ! ij
             cdiv(n,l,0,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5_RP * GMTR_P_var(n,k0,l,P_RAREA)
             ! ip1j
             cdiv(n,l,1,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                               + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                             ) * 0.5_RP * GMTR_P_var(n,k0,l,P_RAREA)
             ! ip1jp1
             cdiv(n,l,2,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                             ) * 0.5_RP * GMTR_P_var(n,k0,l,P_RAREA)
             ! ijp1
             cdiv(n,l,3,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                               + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                             ) * 0.5_RP * GMTR_P_var(n,k0,l,P_RAREA)
             ! im1j
             cdiv(n,l,4,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                               - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5_RP * GMTR_P_var(n,k0,l,P_RAREA)
             ! im1jm1
             cdiv(n,l,5,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                               - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5_RP * GMTR_P_var(n,k0,l,P_RAREA)
             ! ijm1
             cdiv(n,l,6,m) = ( + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                               - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                             ) * 0.5_RP * GMTR_P_var(n,k0,l,P_RAREA)
          enddo
       endif

    enddo ! loop l

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
          do m = 1, 3
             md = m + HNX - 1

             cdiv_pl(n,l,0,m) = 0.0_RP
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                cdiv_pl(n,l,0,m) = cdiv_pl(n,l,0,m) + ( GMTR_T_var_pl(ij,k0,l,W1) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                                      + GMTR_T_var_pl(ij,k0,l,W1) * GMTR_A_var_pl(ijp1,k0,l,md) )
             enddo
             cdiv_pl(n,l,0,m) = cdiv_pl(n,l,0,m) * 0.5_RP * GMTR_P_var_pl(n,k0,l,P_RAREA)

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                cdiv_pl(n,l,v-1,m) = ( + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ijm1,k0,l,md) &
                                       + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                       + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                       + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ijp1,k0,l,md) &
                                     ) * 0.5_RP * GMTR_P_var_pl(n,k0,l,P_RAREA)
             enddo
          enddo ! loop m
       enddo ! loop l
    endif



    !---< setup coefficient of gradient operator >

    write(ADM_LOG_FID,*) '*** setup coefficient of gradient operator'

    do l = 1, ADM_lall

       do m = 1, 3
          md = m + HNX - 1
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             ! ij
             cgrad(n,l,0,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W3) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                - 2.0_RP * GMTR_A_var(ij    ,k0,l,AI ,md)                          &
                                - 2.0_RP * GMTR_A_var(ij    ,k0,l,AIJ,md)                          &
                                - 2.0_RP * GMTR_A_var(ij    ,k0,l,AJ ,md)                          &
                                + 2.0_RP * GMTR_A_var(im1j  ,k0,l,AI ,md)                          &
                                + 2.0_RP * GMTR_A_var(im1jm1,k0,l,AIJ,md)                          &
                                + 2.0_RP * GMTR_A_var(ijm1  ,k0,l,AJ ,md)                          &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1j
             cgrad(n,l,1,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1jp1
             cgrad(n,l,2,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijp1
             cgrad(n,l,3,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1j
             cgrad(n,l,4,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1jm1
             cgrad(n,l,5,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijm1
             cgrad(n,l,6,m) = ( - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TI,W2) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ijm1  ,k0,l,AJ ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
          enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          do m = 1, 3
             md = m + HNX - 1

             ij     = n
             ip1j   = n + 1
             ijp1   = n     + ADM_gall_1d
             ip1jp1 = n + 1 + ADM_gall_1d
             im1j   = n - 1
             ijm1   = n     - ADM_gall_1d
             im1jm1 = n - 1 - ADM_gall_1d

             ! ij
             cgrad(n,l,0,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                - 2.0_RP * GMTR_A_var(ij    ,k0,l,AI ,md)                          &
                                - 2.0_RP * GMTR_A_var(ij    ,k0,l,AIJ,md)                          &
                                - 2.0_RP * GMTR_A_var(ij    ,k0,l,AJ ,md)                          &
                                + 2.0_RP * GMTR_A_var(im1j  ,k0,l,AI ,md)                          &
                                + 2.0_RP * GMTR_A_var(im1jm1,k0,l,AIJ,md)                          &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1j
             cgrad(n,l,1,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ip1jp1
             cgrad(n,l,2,m) = ( + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                                + GMTR_T_var(ij    ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W2) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijp1
             cgrad(n,l,3,m) = ( + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AIJ,md) &
                                + GMTR_T_var(ij    ,k0,l,TJ,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                + GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1j
             cgrad(n,l,4,m) = ( + GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(ij    ,k0,l,AJ ,md) &
                                - GMTR_T_var(im1j  ,k0,l,TI,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W3) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! im1jm1
             cgrad(n,l,5,m) = ( - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1j  ,k0,l,AI ,md) &
                                - GMTR_T_var(im1jm1,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
             ! ijm1
             cgrad(n,l,6,m) = ( - GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(im1jm1,k0,l,AIJ,md) &
                                + GMTR_T_var(ijm1  ,k0,l,TJ,W1) * GMTR_A_var(ij    ,k0,l,AI ,md) &
                              ) * 0.5_RP * GMTR_P_var(ij,k0,l,P_RAREA)
          enddo
       endif
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
          do m = 1, 3
             md = m + HNX - 1

             cgrad_pl(n,l,0,m) = 0.0_RP
             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

                cgrad_pl(n,l,0,m) = cgrad_pl(n,l,0,m) &
                                  + 2.0_RP * ( GMTR_T_var_pl(ij,k0,l,W1) - 1.0_RP ) * GMTR_A_var_pl(ijp1,k0,l,md)
             enddo
             cgrad_pl(n,l,0,m) = cgrad_pl(n,l,0,m) * 0.5_RP * GMTR_P_var_pl(n,k0,l,P_RAREA)

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                cgrad_pl(n,l,v-1,m) = ( + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ijm1,k0,l,md) &
                                        + GMTR_T_var_pl(ijm1,k0,l,W3) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                        + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ij  ,k0,l,md) &
                                        + GMTR_T_var_pl(ij  ,k0,l,W2) * GMTR_A_var_pl(ijp1,k0,l,md) &
                                      ) * 0.5_RP * GMTR_P_var_pl(n,k0,l,P_RAREA)
             enddo
          enddo ! loop m
       enddo ! loop l
    endif

    ! ---- setup coefficient of laplacian operator

    write(ADM_LOG_FID,*) '*** setup coefficient of laplacian operator'

    do l = 1, ADM_lall

       do n = OPRT_nstart, OPRT_nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d
          im1jm1 = n - 1 - ADM_gall_1d

          ! ij
          clap(ij,l,0) = ( - 1.0_RP * GMTR_A_var(ij  ,k0,l,AI ,TNX) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNX) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 2.0_RP * GMTR_A_var(ip1j,k0,l,AJ ,TNX) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ijm1,k0,l,AJ ,TNX) * GMTR_T_var(ijm1,k0,l,TJ,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AI ,TNX) * GMTR_T_var(ijm1,k0,l,TJ,T_RAREA) &
                           + 2.0_RP * GMTR_A_var(ijm1,k0,l,AIJ,TNX) * GMTR_T_var(ijm1,k0,l,TJ,T_RAREA) &
                         ) * GMTR_A_var(ij,k0,l,AI ,HNX)

          clap(ij,l,0) = clap(ij,l,0) &
                       + ( - 1.0_RP * GMTR_A_var(ij  ,k0,l,AI ,TNY) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNY) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 2.0_RP * GMTR_A_var(ip1j,k0,l,AJ ,TNY) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ijm1,k0,l,AJ ,TNY) * GMTR_T_var(ijm1,k0,l,TJ,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AI ,TNY) * GMTR_T_var(ijm1,k0,l,TJ,T_RAREA) &
                           + 2.0_RP * GMTR_A_var(ijm1,k0,l,AIJ,TNY) * GMTR_T_var(ijm1,k0,l,TJ,T_RAREA) &
                         ) * GMTR_A_var(ij,k0,l,AI ,HNY)

          clap(ij,l,0) = clap(ij,l,0) &
                       + ( - 1.0_RP * GMTR_A_var(ij  ,k0,l,AI ,TNZ) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNZ) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 2.0_RP * GMTR_A_var(ip1j,k0,l,AJ ,TNZ) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ijm1,k0,l,AJ ,TNZ) * GMTR_T_var(ijm1,k0,l,TJ,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AI ,TNZ) * GMTR_T_var(ijm1,k0,l,TJ,T_RAREA) &
                           + 2.0_RP * GMTR_A_var(ijm1,k0,l,AIJ,TNZ) * GMTR_T_var(ijm1,k0,l,TJ,T_RAREA) &
                         ) * GMTR_A_var(ij,k0,l,AI ,HNZ)

          clap(ij,l,0) = clap(ij,l,0) &
                       + ( - 1.0_RP * GMTR_A_var(ij  ,k0,l,AI ,TNX) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNX) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 2.0_RP * GMTR_A_var(ip1j,k0,l,AJ ,TNX) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           - 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNX) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AJ ,TNX) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           - 2.0_RP * GMTR_A_var(ijp1,k0,l,AI ,TNX) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                         ) * GMTR_A_var(ij,k0,l,AIJ,HNX)

          clap(ij,l,0) = clap(ij,l,0) &
                       + ( - 1.0_RP * GMTR_A_var(ij  ,k0,l,AI ,TNY) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNY) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 2.0_RP * GMTR_A_var(ip1j,k0,l,AJ ,TNY) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           - 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNY) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AJ ,TNY) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           - 2.0_RP * GMTR_A_var(ijp1,k0,l,AI ,TNY) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                         ) * GMTR_A_var(ij,k0,l,AIJ,HNY)

          clap(ij,l,0) = clap(ij,l,0) &
                       + ( - 1.0_RP * GMTR_A_var(ij  ,k0,l,AI ,TNZ) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNZ) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           + 2.0_RP * GMTR_A_var(ip1j,k0,l,AJ ,TNZ) * GMTR_T_var(ij  ,k0,l,TI,T_RAREA) &
                           - 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNZ) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           - 2.0_RP * GMTR_A_var(ijp1,k0,l,AI ,TNZ) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AJ ,TNZ) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                         ) * GMTR_A_var(ij,k0,l,AIJ,HNZ)

          clap(ij,l,0) = clap(ij,l,0) &
                       + ( - 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNX) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AJ ,TNX) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           - 2.0_RP * GMTR_A_var(ijp1,k0,l,AI ,TNX) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           - 1.0_RP * GMTR_A_var(ij  ,k0,l,AJ ,TNX) * GMTR_T_var(im1j,k0,l,TI,T_RAREA) &
                           - 1.0_RP * GMTR_A_var(im1j,k0,l,AI ,TNX) * GMTR_T_var(im1j,k0,l,TI,T_RAREA) &
                           - 2.0_RP * GMTR_A_var(im1j,k0,l,AIJ,TNX) * GMTR_T_var(im1j,k0,l,TI,T_RAREA) &
                         ) * GMTR_A_var(ij,k0,l,AJ ,HNX)

          clap(ij,l,0) = clap(ij,l,0) &
                       + ( - 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNY) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AJ ,TNY) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           - 2.0_RP * GMTR_A_var(ijp1,k0,l,AI ,TNY) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           - 1.0_RP * GMTR_A_var(ij  ,k0,l,AJ ,TNY) * GMTR_T_var(im1j,k0,l,TI,T_RAREA) &
                           - 1.0_RP * GMTR_A_var(im1j,k0,l,AI ,TNY) * GMTR_T_var(im1j,k0,l,TI,T_RAREA) &
                           - 2.0_RP * GMTR_A_var(im1j,k0,l,AIJ,TNY) * GMTR_T_var(im1j,k0,l,TI,T_RAREA) &
                         ) * GMTR_A_var(ij,k0,l,AJ ,HNY)

          clap(ij,l,0) = clap(ij,l,0) &
                       + ( - 1.0_RP * GMTR_A_var(ij  ,k0,l,AIJ,TNZ) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           + 1.0_RP * GMTR_A_var(ij  ,k0,l,AJ ,TNZ) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           - 2.0_RP * GMTR_A_var(ijp1,k0,l,AI ,TNZ) * GMTR_T_var(ij  ,k0,l,TJ,T_RAREA) &
                           - 1.0_RP * GMTR_A_var(ij  ,k0,l,AJ ,TNZ) * GMTR_T_var(im1j,k0,l,TI,T_RAREA) &
                           - 1.0_RP * GMTR_A_var(im1j,k0,l,AI ,TNZ) * GMTR_T_var(im1j,k0,l,TI,T_RAREA) &
                           - 2.0_RP * GMTR_A_var(im1j,k0,l,AIJ,TNZ) * GMTR_T_var(im1j,k0,l,TI,T_RAREA) &
                         ) * GMTR_A_var(ij,k0,l,AJ ,HNZ)

          clap(ij,l,0) = clap(ij,l,0) * GMTR_P_var(ij,k0,l,P_RAREA) / 12.0_RP

          clap(ij,l,0) = clap(ij,l,0) + ( &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -2*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! ip1j
          clap(ij,l,1) = ( &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! ip1jp1
          clap(ij,l,2) = ( &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! ijp1
          clap(ij,l,3) = ( &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! im1j
          clap(ij,l,4) = ( &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! im1jm1
          clap(ij,l,5) = ( &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      -2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! ijm1
          clap(ij,l,6) = ( &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +2*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +2*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +2*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      +2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TI,T_RAREA)*GMTR_A_var(ijm1,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          n = suf(ADM_gmin,ADM_gmin)

          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d
          im1jm1 = n - 1 - ADM_gall_1d

          ! 0: ij
          clap(ij,l,0) = ( &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP  ! Y.Niwa add 06/08/22

          clap(ij,l,0) = clap(ij,l,0) + ( &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +2*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! ip1j
          clap(ij,l,1) = ( &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -2*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(ij    ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(ijm1  ,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNX)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNY)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(ijm1  ,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1  ,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! ip1jp1
          clap(ij,l,2) = ( &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          -1*GMTR_A_var(ip1j,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! ijp1
          clap(ij,l,3) = ( &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AIJ,HNZ) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNX)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNY)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(ijp1,k0,l,AI ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AIJ,TNZ)*GMTR_T_var(ij  ,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          -1*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! im1j
          clap(ij,l,4) = ( &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNX)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNX) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          +1*GMTR_A_var(im1j,k0,l,AIJ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNY)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNY) &
          -1*GMTR_A_var(im1j,k0,l,AI ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
          +2*GMTR_A_var(ij  ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j,k0,l,TI,T_RAREA)*GMTR_A_var(ij,k0,l,AJ,HNZ) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1j  ,k0,l,AIJ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNX)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNY)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -2*GMTR_A_var(ij    ,k0,l,AJ ,TNZ)*GMTR_T_var(im1j  ,k0,l,TI,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -2*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! im1jm1
          clap(ij,l,5) = ( &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNX) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNY) &
       +2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1j,k0,l,AI,HNZ) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(im1jm1,k0,l,AJ ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(im1jm1,k0,l,AIJ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNX)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNY)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(im1j  ,k0,l,AI ,TNZ)*GMTR_T_var(im1jm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

          ! ijm1
          clap(ij,l,6) = ( &
      -1*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      -1*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      -1*GMTR_A_var(ijm1,k0,l,AJ ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      +2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNX) &
      +2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNY) &
      +2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(im1jm1,k0,l,AIJ,HNZ) &
      -1*GMTR_A_var(ijm1,k0,l,AIJ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
      +1*GMTR_A_var(ijm1,k0,l,AJ ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
      -2*GMTR_A_var(ij  ,k0,l,AI ,TNX)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNX) &
      -1*GMTR_A_var(ijm1,k0,l,AIJ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
      +1*GMTR_A_var(ijm1,k0,l,AJ ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
      -1*GMTR_A_var(ijm1,k0,l,AIJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
      +1*GMTR_A_var(ijm1,k0,l,AJ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
      -2*GMTR_A_var(ij  ,k0,l,AI ,TNZ)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNZ) &
      -2*GMTR_A_var(ij  ,k0,l,AI ,TNY)*GMTR_T_var(ijm1,k0,l,TJ,T_RAREA)*GMTR_A_var(ij,k0,l,AI,HNY) &
         )*GMTR_P_var(ij,k0,l,P_RAREA)/12.0_RP

      endif
    enddo

    if ( ADM_have_pl ) then

      n =ADM_gslf_pl
      n0=ADM_gmin_pl
      n1=ADM_gmin_pl+1
      n2=ADM_gmin_pl+2
      n3=ADM_gmin_pl+3
      n4=ADM_gmin_pl+4

      do l = 1,ADM_lall_pl
          clap_pl(n,l,0)=( &
                       +1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0_RP   ! Y.Niwa add 060822

          clap_pl(n,l,0)= clap_pl(n,l,0) + ( &                        ! Y.Niwa add 060822
                       +1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0_RP
          !
          ! n0
          clap_pl(n,l,1)=( &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0_RP
          !
          ! n1
          clap_pl(n,l,2)=( &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2X)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Y)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n0,k0,l,TN2Z)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n0,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0_RP
          !
          ! n2
          clap_pl(n,l,3)=( &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n1,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2X)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Y)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n1,k0,l,TN2Z)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n1,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n1,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0_RP
          !
          ! n3
          clap_pl(n,l,4)=( &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n2,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2X)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Y)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n2,k0,l,TN2Z)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n2,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n2,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0_RP
          !
          ! n4
          clap_pl(n,l,5)=( &
                       +1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n0,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n3,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2X)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2X)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNX)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNX) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Y)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Y)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNY)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNY) &
                       +1*GMTR_A_var_pl(n3,k0,l,TN2Z)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -2*GMTR_A_var_pl(n3,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       -1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n3,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +1*GMTR_A_var_pl(n4,k0,l,TN2Z)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                       +2*GMTR_A_var_pl(n0,k0,l,TNZ)*GMTR_T_var_pl(n4,k0,l,T_RAREA)*GMTR_A_var_pl(n4,k0,l,HNZ) &
                      )*GMTR_P_var_pl(n,k0,l,P_RAREA)/12.0_RP
      enddo
    endif




    ! ---- setup coefficient of diffusion operator

    write(ADM_LOG_FID,*) '*** setup coefficient of diffusion operator'

    do m = AI, AJ
    do l = 1, ADM_lall
    do n = 1, ADM_gall
       cinterp_TN(n,l,m,1) = GMTR_A_var(n,k0,l,m,TNX)
       cinterp_TN(n,l,m,2) = GMTR_A_var(n,k0,l,m,TNY)
       cinterp_TN(n,l,m,3) = GMTR_A_var(n,k0,l,m,TNZ)

       cinterp_HN(n,l,m,1) = GMTR_A_var(n,k0,l,m,HNX)
       cinterp_HN(n,l,m,2) = GMTR_A_var(n,k0,l,m,HNY)
       cinterp_HN(n,l,m,3) = GMTR_A_var(n,k0,l,m,HNZ)
    enddo
    enddo
    enddo

    do l = 1, ADM_lall
    do n = 1, ADM_gall
       cinterp_TRA(n,l,TI) = GMTR_T_var(n,k0,l,TI,T_RAREA)
       cinterp_TRA(n,l,TJ) = GMTR_T_var(n,k0,l,TJ,T_RAREA)

       cinterp_PRA(n,l)    = GMTR_P_var(n,k0,l,P_RAREA)
    enddo
    enddo

    return
  end subroutine OPRT_setup

  !-----------------------------------------------------------------------------
  !> horizontal divergence operator
  subroutine OPRT_divergence( &
       scl, scl_pl, &
       vx,  vx_pl,  &
       vy,  vy_pl,  &
       vz,  vz_pl   )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gslf_pl, &
       ADM_gmax_pl
    implicit none

    real(RP), intent(out) :: scl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: sclx(ADM_gall)
    real(RP) :: scly(ADM_gall)
    real(RP) :: sclz(ADM_gall)

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: gall, gall_1d, kall
    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_divergence')

    gall    = ADM_gall
    gall_1d = ADM_gall_1d
    kall    = ADM_kall

    do l = 1, ADM_lall
       !$omp parallel default(none),private(n,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
       !$omp shared(OPRT_nstart,OPRT_nend,gall,gall_1d,kall,l,scl,sclx,scly,sclz,cdiv,vx,vy,vz)
       do k = 1, kall

          !$omp do
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + gall_1d
             ip1jp1 = n + 1 + gall_1d
             im1j   = n - 1
             ijm1   = n     - gall_1d
             im1jm1 = n - 1 - gall_1d

             sclx(n) = cdiv(n,l,0,1) * vx(ij    ,k,l) &
                     + cdiv(n,l,1,1) * vx(ip1j  ,k,l) &
                     + cdiv(n,l,2,1) * vx(ip1jp1,k,l) &
                     + cdiv(n,l,3,1) * vx(ijp1  ,k,l) &
                     + cdiv(n,l,4,1) * vx(im1j  ,k,l) &
                     + cdiv(n,l,5,1) * vx(im1jm1,k,l) &
                     + cdiv(n,l,6,1) * vx(ijm1  ,k,l)
          enddo
          !$omp end do nowait

          !$omp do
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + gall_1d
             ip1jp1 = n + 1 + gall_1d
             im1j   = n - 1
             ijm1   = n     - gall_1d
             im1jm1 = n - 1 - gall_1d

             scly(n) = cdiv(n,l,0,2) * vy(ij    ,k,l) &
                     + cdiv(n,l,1,2) * vy(ip1j  ,k,l) &
                     + cdiv(n,l,2,2) * vy(ip1jp1,k,l) &
                     + cdiv(n,l,3,2) * vy(ijp1  ,k,l) &
                     + cdiv(n,l,4,2) * vy(im1j  ,k,l) &
                     + cdiv(n,l,5,2) * vy(im1jm1,k,l) &
                     + cdiv(n,l,6,2) * vy(ijm1  ,k,l)
          enddo
          !$omp end do nowait

          !$omp do
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + gall_1d
             ip1jp1 = n + 1 + gall_1d
             im1j   = n - 1
             ijm1   = n     - gall_1d
             im1jm1 = n - 1 - gall_1d

             sclz(n) = cdiv(n,l,0,3) * vz(ij    ,k,l) &
                     + cdiv(n,l,1,3) * vz(ip1j  ,k,l) &
                     + cdiv(n,l,2,3) * vz(ip1jp1,k,l) &
                     + cdiv(n,l,3,3) * vz(ijp1  ,k,l) &
                     + cdiv(n,l,4,3) * vz(im1j  ,k,l) &
                     + cdiv(n,l,5,3) * vz(im1jm1,k,l) &
                     + cdiv(n,l,6,3) * vz(ijm1  ,k,l)
          enddo
          !$omp end do nowait

          !$omp do
          do n = 1, OPRT_nstart-1
             scl(n,k,l) = 0.0_RP
          enddo
          !$omp end do nowait

          !$omp do
          do n = OPRT_nend+1, gall
             scl(n,k,l) = 0.0_RP
          enddo
          !$omp end do

          !$omp do
          do n = OPRT_nstart, OPRT_nend
             scl(n,k,l) = sclx(n) + scly(n) + sclz(n)
          enddo
          !$omp end do

       enddo
       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          scl_pl(:,k,l) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + ( cdiv_pl(n,l,v-1,1) * vx_pl(v,k,l) &
                                             + cdiv_pl(n,l,v-1,2) * vy_pl(v,k,l) &
                                             + cdiv_pl(n,l,v-1,3) * vz_pl(v,k,l) )
          enddo
       enddo
       enddo
    else
       scl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_divergence')

    return
  end subroutine OPRT_divergence

  !-----------------------------------------------------------------------------
  !> horizontal gradient operator
  subroutine OPRT_gradient( &
       scl,  scl_pl, &
       grad, grad_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gslf_pl, &
       ADM_gmax_pl, &
       ADM_nxyz
    implicit none

    real(RP), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: grad   (ADM_gall   ,ADM_kall,ADM_lall   ,ADM_nxyz)
    real(RP), intent(out) :: grad_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,ADM_nxyz)

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v, d
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_gradient')

!OCL SERIAL
    do d = 1, ADM_nxyz
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
       do n = OPRT_nstart, OPRT_nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d
          im1jm1 = n - 1 - ADM_gall_1d

          grad(n,k,l,d) = cgrad(n,l,0,d) * scl(ij    ,k,l) &
                        + cgrad(n,l,1,d) * scl(ip1j  ,k,l) &
                        + cgrad(n,l,2,d) * scl(ip1jp1,k,l) &
                        + cgrad(n,l,3,d) * scl(ijp1  ,k,l) &
                        + cgrad(n,l,4,d) * scl(im1j  ,k,l) &
                        + cgrad(n,l,5,d) * scl(im1jm1,k,l) &
                        + cgrad(n,l,6,d) * scl(ijm1  ,k,l)
       enddo
       grad(          1:OPRT_nstart-1,k,l,d) = 0.0_RP
       grad(OPRT_nend+1:ADM_gall     ,k,l,d) = 0.0_RP
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do d = 1, ADM_nxyz
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          grad_pl(:,k,l,d) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             grad_pl(n,k,l,d) = grad_pl(n,k,l,d) + cgrad_pl(n,l,v-1,d) * scl_pl(v,k,l)
          enddo
       enddo
       enddo
       enddo
    else
       grad_pl(:,:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_gradient')

    return
  end subroutine OPRT_gradient

  !-----------------------------------------------------------------------------
  !> horizontal laplacian operator
  subroutine OPRT_laplacian( &
       dscl, dscl_pl, &
       scl,  scl_pl   )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gslf_pl, &
       ADM_gmax_pl
    implicit none

    real(RP), intent(out) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_laplacian')

!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
       do n = OPRT_nstart, OPRT_nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d
          im1jm1 = n - 1 - ADM_gall_1d

          dscl(n,k,l) = clap(n,l,0) * scl(ij    ,k,l) &
                      + clap(n,l,1) * scl(ip1j  ,k,l) &
                      + clap(n,l,2) * scl(ip1jp1,k,l) &
                      + clap(n,l,3) * scl(ijp1  ,k,l) &
                      + clap(n,l,4) * scl(im1j  ,k,l) &
                      + clap(n,l,5) * scl(im1jm1,k,l) &
                      + clap(n,l,6) * scl(ijm1  ,k,l)
       enddo
       dscl(          1:OPRT_nstart-1,k,l) = 0.0_RP
       dscl(OPRT_nend+1:ADM_gall     ,k,l) = 0.0_RP
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          dscl_pl(:,k,l) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             dscl_pl(n,k,l) = dscl_pl(n,k,l) + clap_pl(n,l,v-1) * scl_pl(v,k,l)
          enddo
       enddo
       enddo
    else
       dscl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_laplacian')

    return
  end subroutine OPRT_laplacian

  !-----------------------------------------------------------------------------
  !> horizontal diffusion operator
  subroutine OPRT_diffusion( &
       dscl, dscl_pl, &
       scl,  scl_pl,  &
       kh,   kh_pl    )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    implicit none

    real(RP), intent(out) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP)  :: vxt    (ADM_gall   ,TI:TJ)
    real(RP)  :: vyt    (ADM_gall   ,TI:TJ)
    real(RP)  :: vzt    (ADM_gall   ,TI:TJ)
    real(RP)  :: flux   (ADM_gall   ,AI:AJ)
    real(RP)  :: vxt_pl (ADM_gall_pl)
    real(RP)  :: vyt_pl (ADM_gall_pl)
    real(RP)  :: vzt_pl (ADM_gall_pl)
    real(RP)  :: flux_pl(ADM_gall_pl)

    real(RP) :: u1, u2, u3, smean

    integer :: nstart1, nstart2, nstart3, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: gall, gall_1d, gmin, kall
    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_diffusion')

    gall    = ADM_gall
    gall_1d = ADM_gall_1d
    gmin    = ADM_gmin
    kall    = ADM_kall

    nstart1 = suf(ADM_gmin-1,ADM_gmin-1)
    nstart2 = suf(ADM_gmin-1,ADM_gmin  )
    nstart3 = suf(ADM_gmin  ,ADM_gmin-1)
    nend    = suf(ADM_gmax  ,ADM_gmax  )

    do l = 1, ADM_lall
       !$omp parallel default(none),private(n,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1,smean,u1,u2,u3), &
       !$omp shared(OPRT_nstart,OPRT_nend,gall,gall_1d,gmin,kall,nstart1,nstart2,nstart3,nend,l,ADM_have_sgp, &
       !$omp dscl,scl,kh,vxt,vyt,vzt,flux,cinterp_TN,cinterp_HN,cinterp_PRA,cinterp_TRA)
       do k = 1, kall

          !$omp do
          do n = nstart1, nend
             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + gall_1d

             smean = ( scl(ij,k,l) + scl(ip1j,k,l) + scl(ip1jp1,k,l) ) / 3.0_RP

             u1 = 0.5_RP * (scl(ij    ,k,l)+scl(ip1j  ,k,l)) - smean
             u2 = 0.5_RP * (scl(ip1j  ,k,l)+scl(ip1jp1,k,l)) - smean
             u3 = 0.5_RP * (scl(ip1jp1,k,l)+scl(ij    ,k,l)) - smean

             vxt(n,TI) = ( - u1 * cinterp_TN(ij  ,l,AI ,1) &
                           - u2 * cinterp_TN(ip1j,l,AJ ,1) &
                           + u3 * cinterp_TN(ij  ,l,AIJ,1) ) * cinterp_TRA(ij,l,TI)
             vyt(n,TI) = ( - u1 * cinterp_TN(ij  ,l,AI ,2) &
                           - u2 * cinterp_TN(ip1j,l,AJ ,2) &
                           + u3 * cinterp_TN(ij  ,l,AIJ,2) ) * cinterp_TRA(ij,l,TI)
             vzt(n,TI) = ( - u1 * cinterp_TN(ij  ,l,AI ,3) &
                           - u2 * cinterp_TN(ip1j,l,AJ ,3) &
                           + u3 * cinterp_TN(ij  ,l,AIJ,3) ) * cinterp_TRA(ij,l,TI)
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             ij     = n
             ijp1   = n     + gall_1d
             ip1jp1 = n + 1 + gall_1d

             smean = ( scl(ij,k,l) + scl(ip1jp1,k,l) + scl(ijp1,k,l) ) / 3.0_RP

             u1 = 0.5_RP * (scl(ij    ,k,l)+scl(ip1jp1,k,l)) - smean
             u2 = 0.5_RP * (scl(ip1jp1,k,l)+scl(ijp1  ,k,l)) - smean
             u3 = 0.5_RP * (scl(ijp1  ,k,l)+scl(ij    ,k,l)) - smean

             vxt(n,TJ) = ( - u1 * cinterp_TN(ij  ,l,AIJ,1) &
                           + u2 * cinterp_TN(ijp1,l,AI ,1) &
                           + u3 * cinterp_TN(ij  ,l,AJ ,1) ) * cinterp_TRA(ij,l,TJ)
             vyt(n,TJ) = ( - u1 * cinterp_TN(ij  ,l,AIJ,2) &
                           + u2 * cinterp_TN(ijp1,l,AI ,2) &
                           + u3 * cinterp_TN(ij  ,l,AJ ,2) ) * cinterp_TRA(ij,l,TJ)
             vzt(n,TJ) = ( - u1 * cinterp_TN(ij  ,l,AIJ,3) &
                           + u2 * cinterp_TN(ijp1,l,AI ,3) &
                           + u3 * cinterp_TN(ij  ,l,AJ ,3) ) * cinterp_TRA(ij,l,TJ)
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then ! pentagon
             !$omp master
             vxt(suf(gmin-1,gmin-1),TI) = vxt(suf(gmin,gmin-1),TJ)
             vyt(suf(gmin-1,gmin-1),TI) = vyt(suf(gmin,gmin-1),TJ)
             vzt(suf(gmin-1,gmin-1),TI) = vzt(suf(gmin,gmin-1),TJ)
             !$omp end master
             !$omp barrier
          endif

          !$omp do
          do n = nstart2, nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + gall_1d
             ip1jp1 = n + 1 + gall_1d
             im1j   = n - 1
             ijm1   = n     - gall_1d

             flux(n,AI ) = 0.25_RP * ( (vxt(ijm1,TJ)+vxt(ij,TI)) * cinterp_HN(ij,l,AI ,1) &
                                    + (vyt(ijm1,TJ)+vyt(ij,TI)) * cinterp_HN(ij,l,AI ,2) &
                                    + (vzt(ijm1,TJ)+vzt(ij,TI)) * cinterp_HN(ij,l,AI ,3) &
                                    ) * (kh(ij,k,l)+kh(ip1j,k,l))
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             ij     = n
             ip1jp1 = n + 1 + gall_1d

             flux(n,AIJ) = 0.25_RP * ( (vxt(ij ,TI)+vxt(ij ,TJ)) * cinterp_HN(ij,l,AIJ,1) &
                                    + (vyt(ij ,TI)+vyt(ij ,TJ)) * cinterp_HN(ij,l,AIJ,2) &
                                    + (vzt(ij ,TI)+vzt(ij ,TJ)) * cinterp_HN(ij,l,AIJ,3) &
                                    ) * (kh(ij,k,l)+kh(ip1jp1,k,l))
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart3, nend
             ij     = n
             ijp1   = n     + gall_1d
             im1j   = n - 1

             flux(n,AJ ) = 0.25_RP * ( (vxt(ij,TJ)+vxt(im1j,TI)) * cinterp_HN(ij,l,AJ ,1) &
                                    + (vyt(ij,TJ)+vyt(im1j,TI)) * cinterp_HN(ij,l,AJ ,2) &
                                    + (vzt(ij,TJ)+vzt(im1j,TI)) * cinterp_HN(ij,l,AJ ,3) &
                                    ) * (kh(ij,k,l)+kh(ijp1,k,l))
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then ! pentagon
             !$omp master
             flux(suf(gmin,gmin-1),AJ) = 0.0_RP
             !$omp end master
             !$omp barrier
          endif

          !$omp do
          do n = OPRT_nstart, OPRT_nend
             ij     = n
             im1j   = n - 1
             im1jm1 = n - 1 - gall_1d
             ijm1   = n     - gall_1d

             dscl(n,k,l) = ( flux(ij,AI ) - flux(im1j  ,AI ) &
                           + flux(ij,AIJ) - flux(im1jm1,AIJ) &
                           + flux(ij,AJ ) - flux(ijm1  ,AJ ) ) * cinterp_PRA(ij,l)
          enddo
          !$omp end do nowait

          !$omp do
          do n = 1, OPRT_nstart-1
             dscl(n,k,l) = 0.0_RP
          enddo
          !$omp end do nowait

          !$omp do
          do n = OPRT_nend+1, gall
             dscl(n,k,l) = 0.0_RP
          enddo
          !$omp end do

       enddo
       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

             smean = ( scl_pl(n,k,l) + scl_pl(ij,k,l) + scl_pl(ijp1,k,l) ) / 3.0_RP

             u1 = 0.5_RP * (scl_pl(n   ,k,l)+scl_pl(ij  ,k,l)) - smean
             u2 = 0.5_RP * (scl_pl(ij  ,k,l)+scl_pl(ijp1,k,l)) - smean
             u3 = 0.5_RP * (scl_pl(ijp1,k,l)+scl_pl(n   ,k,l)) - smean

             vxt_pl(v) = ( + u1 * GMTR_A_var_pl(ij  ,k0,l,TNX ) &
                           + u2 * GMTR_A_var_pl(ij  ,k0,l,TN2X) &
                           - u3 * GMTR_A_var_pl(ijp1,k0,l,TNX ) ) * GMTR_T_var_pl(ij,k0,l,T_RAREA)
             vyt_pl(v) = ( + u1 * GMTR_A_var_pl(ij  ,k0,l,TNY ) &
                           + u2 * GMTR_A_var_pl(ij  ,k0,l,TN2Y) &
                           - u3 * GMTR_A_var_pl(ijp1,k0,l,TNY ) ) * GMTR_T_var_pl(ij,k0,l,T_RAREA)
             vzt_pl(v) = ( + u1 * GMTR_A_var_pl(ij  ,k0,l,TNZ ) &
                           + u2 * GMTR_A_var_pl(ij  ,k0,l,TN2Z) &
                           - u3 * GMTR_A_var_pl(ijp1,k0,l,TNZ ) ) * GMTR_T_var_pl(ij,k0,l,T_RAREA)
          enddo

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             flux_pl(v) = 0.25_RP * ( (vxt_pl(ij)+vxt_pl(ijm1)) * GMTR_A_var_pl(ij,k0,l,HNX) &
                                   + (vyt_pl(ij)+vyt_pl(ijm1)) * GMTR_A_var_pl(ij,k0,l,HNY) &
                                   + (vzt_pl(ij)+vzt_pl(ijm1)) * GMTR_A_var_pl(ij,k0,l,HNZ) &
                                   ) * (kh_pl(n,k,l)+kh_pl(ij,k,l))
          enddo

          dscl_pl(:,k,l) = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             dscl_pl(n,k,l) = dscl_pl(n,k,l) + flux_pl(v)
          enddo
          dscl_pl(n,k,l) = dscl_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_RAREA)

       enddo
       enddo

    else
       dscl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_diffusion')

    return
  end subroutine OPRT_diffusion

  !-----------------------------------------------------------------------------
  !> horizontalize 3-D vector
  subroutine OPRT_horizontalize_vec( &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl  )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_KNONE
    use mod_grd, only: &
       GRD_rscale, &
       GRD_XDIR,   &
       GRD_YDIR,   &
       GRD_ZDIR,   &
       GRD_x,      &
       GRD_x_pl
    implicit none

    real(RP), intent(inout) :: vx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: prd
    integer :: g, k, l
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_horizontalize_vec')

!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       prd = vx(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale &
           + vy(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale &
           + vz(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale

       vx(g,k,l) = vx(g,k,l) - prd * GRD_x(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
       vy(g,k,l) = vy(g,k,l) - prd * GRD_x(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
       vz(g,k,l) = vz(g,k,l) - prd * GRD_x(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          prd = vx_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale &
              + vy_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale &
              + vz_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale

          vx_pl(g,k,l) = vx_pl(g,k,l) - prd * GRD_x_pl(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
          vy_pl(g,k,l) = vy_pl(g,k,l) - prd * GRD_x_pl(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
          vz_pl(g,k,l) = vz_pl(g,k,l) - prd * GRD_x_pl(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
       enddo
       enddo
       enddo
    else
       vx_pl(:,:,:) = 0.0_RP
       vy_pl(:,:,:) = 0.0_RP
       vz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_horizontalize_vec')

    return
  end subroutine OPRT_horizontalize_vec

  !-----------------------------------------------------------------------------
  !> horizontal vorticity operator
  subroutine OPRT_vorticity( &
       scl, scl_pl, &
       vx,  vx_pl,  &
       vy,  vy_pl,  &
       vz,  vz_pl   )
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
       ADM_gmax_pl
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    implicit none

    real(RP), intent(out) :: scl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: scl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP)  :: vxt    (ADM_gall   ,TI:TJ)
    real(RP)  :: vxt_pl (ADM_gall_pl)
    real(RP)  :: vyt    (ADM_gall   ,TI:TJ)
    real(RP)  :: vyt_pl (ADM_gall_pl)
    real(RP)  :: vzt    (ADM_gall   ,TI:TJ)
    real(RP)  :: vzt_pl (ADM_gall_pl)
    real(RP)  :: flux   (ADM_gall   ,AI:AJ)
    real(RP)  :: flux_pl(ADM_gall_pl)

    integer :: nstart, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_vorticity')

    do l = 1, ADM_lall
    do k = 1, ADM_kall

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d

          vxt(n,TI) = vx(ij    ,k,l) * GMTR_T_var(n,k0,l,TI,W1) &
                    + vx(ip1j  ,k,l) * GMTR_T_var(n,k0,l,TI,W2) &
                    + vx(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TI,W3)
          vyt(n,TI) = vy(ij    ,k,l) * GMTR_T_var(n,k0,l,TI,W1) &
                    + vy(ip1j  ,k,l) * GMTR_T_var(n,k0,l,TI,W2) &
                    + vy(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TI,W3)
          vzt(n,TI) = vz(ij    ,k,l) * GMTR_T_var(n,k0,l,TI,W1) &
                    + vz(ip1j  ,k,l) * GMTR_T_var(n,k0,l,TI,W2) &
                    + vz(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TI,W3)
       enddo

       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          vxt(n,TJ) = vx(ij    ,k,l) * GMTR_T_var(n,k0,l,TJ,W1) &
                    + vx(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TJ,W2) &
                    + vx(ijp1  ,k,l) * GMTR_T_var(n,k0,l,TJ,W3)
          vyt(n,TJ) = vy(ij    ,k,l) * GMTR_T_var(n,k0,l,TJ,W1) &
                    + vy(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TJ,W2) &
                    + vy(ijp1  ,k,l) * GMTR_T_var(n,k0,l,TJ,W3)
          vzt(n,TJ) = vz(ij    ,k,l) * GMTR_T_var(n,k0,l,TJ,W1) &
                    + vz(ip1jp1,k,l) * GMTR_T_var(n,k0,l,TJ,W2) &
                    + vz(ijp1  ,k,l) * GMTR_T_var(n,k0,l,TJ,W3)
       enddo

       if ( ADM_have_sgp(l) ) then
          vxt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vxt(suf(ADM_gmin,ADM_gmin-1),TJ)
          vyt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vyt(suf(ADM_gmin,ADM_gmin-1),TJ)
          vzt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vzt(suf(ADM_gmin,ADM_gmin-1),TJ)
       endif

       nstart = suf(ADM_gmin-1,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ijm1   = n     - ADM_gall_1d

          flux(n,AI ) = 0.5_RP * ( (vxt(ijm1,TJ)+vxt(ij  ,TI)) * cinterp_HN(ij,l,AI ,1) &
                                + (vyt(ijm1,TJ)+vyt(ij  ,TI)) * cinterp_HN(ij,l,AI ,2) &
                                + (vzt(ijm1,TJ)+vzt(ij  ,TI)) * cinterp_HN(ij,l,AI ,3) )
       enddo

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d

          flux(n,AIJ) = 0.5_RP * ( (vxt(ij  ,TI)+vxt(ij  ,TJ)) * cinterp_HN(ij,l,AIJ,1) &
                                + (vyt(ij  ,TI)+vyt(ij  ,TJ)) * cinterp_HN(ij,l,AIJ,2) &
                                + (vzt(ij  ,TI)+vzt(ij  ,TJ)) * cinterp_HN(ij,l,AIJ,3) )
       enddo

       nstart = suf(ADM_gmin  ,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart,nend
          ij     = n
          im1j   = n - 1

          flux(n,AJ ) = 0.5_RP * ( (vxt(ij  ,TJ)+vxt(im1j,TI)) * cinterp_HN(ij,l,AJ ,1) &
                                + (vyt(ij  ,TJ)+vyt(im1j,TI)) * cinterp_HN(ij,l,AJ ,2) &
                                + (vzt(ij  ,TJ)+vzt(im1j,TI)) * cinterp_HN(ij,l,AJ ,3) )
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          flux(suf(ADM_gmin,ADM_gmin-1),AJ) = 0.0_RP
       endif

       do n = OPRT_nstart, OPRT_nend
          ij     = n
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d
          im1jm1 = n - 1 - ADM_gall_1d

          scl(n,k,l) = - ( flux(ij,AI ) - flux(im1j  ,AI ) &
                         + flux(ij,AIJ) - flux(im1jm1,AIJ) &
                         + flux(ij,AJ ) - flux(ijm1  ,AJ ) ) * cinterp_PRA(ij,l)
       enddo
       scl(          1:OPRT_nstart-1,k,l) = 0.0_RP
       scl(OPRT_nend+1:ADM_gall     ,k,l) = 0.0_RP
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

             vxt_pl(v) = vx_pl(n   ,k,l) * GMTR_T_var_pl(ij,k0,l,W1) &
                       + vx_pl(ij  ,k,l) * GMTR_T_var_pl(ij,k0,l,W2) &
                       + vx_pl(ijp1,k,l) * GMTR_T_var_pl(ij,k0,l,W3)

             vyt_pl(v) = vy_pl(n   ,k,l) * GMTR_T_var_pl(ij,k0,l,W1) &
                       + vy_pl(ij  ,k,l) * GMTR_T_var_pl(ij,k0,l,W2) &
                       + vy_pl(ijp1,k,l) * GMTR_T_var_pl(ij,k0,l,W3)

             vzt_pl(v) = vz_pl(n   ,k,l) * GMTR_T_var_pl(ij,k0,l,W1) &
                       + vz_pl(ij  ,k,l) * GMTR_T_var_pl(ij,k0,l,W2) &
                       + vz_pl(ijp1,k,l) * GMTR_T_var_pl(ij,k0,l,W3)
          enddo

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             flux_pl(v) = 0.5_RP * ( (vxt_pl(ijm1)+vxt_pl(ij)) * GMTR_A_var_pl(ij,k0,l,HTX) &
                                  + (vyt_pl(ijm1)+vyt_pl(ij)) * GMTR_A_var_pl(ij,k0,l,HTY) &
                                  + (vzt_pl(ijm1)+vzt_pl(ij)) * GMTR_A_var_pl(ij,k0,l,HTZ) )
          enddo

          scl_pl(:,k,l) = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + flux_pl(v)
          enddo
          scl_pl(n,k,l) = scl_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_RAREA)

       enddo
       enddo
    else
       scl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_vorticity')

    return
  end subroutine OPRT_vorticity

  !-----------------------------------------------------------------------------
  !> horizontal divergence damping
  subroutine OPRT_divdamp( &
       ddivdx, ddivdx_pl, &
       ddivdy, ddivdy_pl, &
       ddivdz, ddivdz_pl, &
       vx,     vx_pl,     &
       vy,     vy_pl,     &
       vz,     vz_pl      )
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
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    implicit none

    real(RP), intent(out) :: ddivdx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: sclt   (ADM_gall   ,TI:TJ)
    real(RP) :: sclt_pl(ADM_gall_pl)

    integer :: nstart, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: gall, gall_1d, gmin, kmin, kmax
    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_divdamp')

    gall    = ADM_gall
    gall_1d = ADM_gall_1d
    gmin    = ADM_gmin
    kmin    = ADM_kmin
    kmax    = ADM_kmax

    nstart = suf(ADM_gmin-1,ADM_gmin-1)
    nend   = suf(ADM_gmax  ,ADM_gmax  )

    do l = 1, ADM_lall
       !$omp parallel default(none),private(n,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1), &
       !$omp shared(OPRT_nstart,OPRT_nend,gall,gall_1d,gmin,kmin,kmax,nstart,nend,l,ADM_have_sgp, &
       !$omp sclt,vx,vy,vz,ddivdx,ddivdy,ddivdz,cinterp_TN,cinterp_HN,cinterp_PRA,cinterp_TRA)
       do k = kmin, kmax

          !$omp do
          do n = nstart, nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + gall_1d
             ip1jp1 = n + 1 + gall_1d

             sclt(n,TI) = ( - ( vx(ij    ,k,l) + vx(ip1j  ,k,l) ) * cinterp_TN(ij  ,l,AI ,1) &
                            - ( vx(ip1j  ,k,l) + vx(ip1jp1,k,l) ) * cinterp_TN(ip1j,l,AJ ,1) &
                            + ( vx(ip1jp1,k,l) + vx(ij    ,k,l) ) * cinterp_TN(ij  ,l,AIJ,1) &
                            - ( vy(ij    ,k,l) + vy(ip1j  ,k,l) ) * cinterp_TN(ij  ,l,AI ,2) &
                            - ( vy(ip1j  ,k,l) + vy(ip1jp1,k,l) ) * cinterp_TN(ip1j,l,AJ ,2) &
                            + ( vy(ip1jp1,k,l) + vy(ij    ,k,l) ) * cinterp_TN(ij  ,l,AIJ,2) &
                            - ( vz(ij    ,k,l) + vz(ip1j  ,k,l) ) * cinterp_TN(ij  ,l,AI ,3) &
                            - ( vz(ip1j  ,k,l) + vz(ip1jp1,k,l) ) * cinterp_TN(ip1j,l,AJ ,3) &
                            + ( vz(ip1jp1,k,l) + vz(ij    ,k,l) ) * cinterp_TN(ij  ,l,AIJ,3) &
                          ) * 0.5_RP * cinterp_TRA(ij,l,TI)

             sclt(n,TJ) = ( - ( vx(ij    ,k,l) + vx(ip1jp1,k,l) ) * cinterp_TN(ij  ,l,AIJ,1) &
                            + ( vx(ip1jp1,k,l) + vx(ijp1  ,k,l) ) * cinterp_TN(ijp1,l,AI ,1) &
                            + ( vx(ijp1  ,k,l) + vx(ij    ,k,l) ) * cinterp_TN(ij  ,l,AJ ,1) &
                            - ( vy(ij    ,k,l) + vy(ip1jp1,k,l) ) * cinterp_TN(ij  ,l,AIJ,2) &
                            + ( vy(ip1jp1,k,l) + vy(ijp1  ,k,l) ) * cinterp_TN(ijp1,l,AI ,2) &
                            + ( vy(ijp1  ,k,l) + vy(ij    ,k,l) ) * cinterp_TN(ij  ,l,AJ ,2) &
                            - ( vz(ij    ,k,l) + vz(ip1jp1,k,l) ) * cinterp_TN(ij  ,l,AIJ,3) &
                            + ( vz(ip1jp1,k,l) + vz(ijp1  ,k,l) ) * cinterp_TN(ijp1,l,AI ,3) &
                            + ( vz(ijp1  ,k,l) + vz(ij    ,k,l) ) * cinterp_TN(ij  ,l,AJ ,3) &
                          ) * 0.5_RP * cinterp_TRA(ij,l,TJ)
          enddo
          !$omp end do

          !$omp do
          do n = 1, OPRT_nstart-1
             ddivdx(n,k,l) = 0.0_RP
             ddivdy(n,k,l) = 0.0_RP
             ddivdz(n,k,l) = 0.0_RP
          enddo
          !$omp end do nowait

          !$omp do
          do n = OPRT_nend+1, gall
             ddivdx(n,k,l) = 0.0_RP
             ddivdy(n,k,l) = 0.0_RP
             ddivdz(n,k,l) = 0.0_RP
          enddo
          !$omp end do nowait

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
          !$omp end do

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

       ddivdx(:,ADM_kmin-1,l) = 0.0_RP
       ddivdx(:,ADM_kmax+1,l) = 0.0_RP
       ddivdy(:,ADM_kmin-1,l) = 0.0_RP
       ddivdy(:,ADM_kmax+1,l) = 0.0_RP
       ddivdz(:,ADM_kmin-1,l) = 0.0_RP
       ddivdz(:,ADM_kmax+1,l) = 0.0_RP
    enddo

    if ( ADM_have_pl ) then
       n = ADM_GSLF_PL

       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

             sclt_pl(v) = ( + ( vx_pl(n   ,k,l) + vx_pl(ij  ,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TNX ) &
                            + ( vy_pl(n   ,k,l) + vy_pl(ij  ,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TNY ) &
                            + ( vz_pl(n   ,k,l) + vz_pl(ij  ,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TNZ ) &
                            + ( vx_pl(ij  ,k,l) + vx_pl(ijp1,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TN2X) &
                            + ( vy_pl(ij  ,k,l) + vy_pl(ijp1,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TN2Y) &
                            + ( vz_pl(ij  ,k,l) + vz_pl(ijp1,k,l) ) * GMTR_A_var_pl(ij,  k0,l,TN2Z) &
                            - ( vx_pl(ijp1,k,l) + vx_pl(n   ,k,l) ) * GMTR_A_var_pl(ijp1,k0,l,TNX ) &
                            - ( vy_pl(ijp1,k,l) + vy_pl(n   ,k,l) ) * GMTR_A_var_pl(ijp1,k0,l,TNY ) &
                            - ( vz_pl(ijp1,k,l) + vz_pl(n   ,k,l) ) * GMTR_A_var_pl(ijp1,k0,l,TNZ ) &
                          ) * 0.5_RP * GMTR_T_var_pl(ij,k0,l,T_RAREA)
          enddo

          ddivdx_pl(:,k,l) = 0.0_RP
          ddivdy_pl(:,k,l) = 0.0_RP
          ddivdz_pl(:,k,l) = 0.0_RP

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

    call PROF_rapend('OPRT_divdamp')

    return
  end subroutine OPRT_divdamp

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_oprt
