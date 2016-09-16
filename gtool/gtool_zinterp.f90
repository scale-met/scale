!-------------------------------------------------------------------------------
!> module ZINTERP
!!
!! @par Description
!!          spacial interpolation module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-08-17 (H.Yashiro)  [new]
!!
!<
module gtool_zinterp
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use dc_log, only: &
     Log, &
#if defined(__PGI) || defined(__ES2)
     LOG_fid, &
#endif
     LOG_LMSG
  use dc_types, only: &
     SP, &
     DP
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: Zinterp_setindex
  public :: Zinterp_setcoef_HGT
  public :: Zinterp_xi2z
  public :: Zinterp_z2xi
  public :: Zinterp_setcoef_PRES
  public :: Zinterp_xi2p

  interface Zinterp_xi2z
     module procedure Zinterp_xi2zSP
     module procedure Zinterp_xi2zDP
  end interface Zinterp_xi2z

  interface Zinterp_z2xi
     module procedure Zinterp_z2xiSP
     module procedure Zinterp_z2xiDP
  end interface Zinterp_z2xi

  interface Zinterp_xi2p
     module procedure Zinterp_xi2pSP
     module procedure Zinterp_xi2pDP
  end interface Zinterp_xi2p

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
  integer,  private, allocatable :: Zinterp_xi2z_idx (:,:,:,:) !< index set   for vertical interpolation (xi->z)
  real(DP), private, allocatable :: Zinterp_xi2z_coef(:,:,:,:) !< coefficient for vertical interpolation (xi->z)
  integer,  private, allocatable :: Zinterp_z2xi_idx (:,:,:,:) !< index set   for vertical interpolation (z->xi)
  real(DP), private, allocatable :: Zinterp_z2xi_coef(:,:,:,:) !< coefficient for vertical interpolation (z->xi)

  integer,  private, allocatable :: Zinterp_xi2p_idx (:,:,:,:) !< index set   for vertical interpolation (xi->p)
  real(DP), private, allocatable :: Zinterp_xi2p_coef(:,:,:,:) !< coefficient for vertical interpolation (xi->p)

  character(len=LOG_LMSG), private :: message = ''

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine Zinterp_setindex( &
       ksize,     &
       isize,     &
       jsize,     &
       ksize_pres )
    implicit none

    integer,  intent(in) :: ksize
    integer,  intent(in) :: isize
    integer,  intent(in) :: jsize
    integer,  intent(in) :: ksize_pres
    !---------------------------------------------------------------------------

    allocate( Zinterp_xi2z_idx (ksize,isize,jsize,2) )
    allocate( Zinterp_xi2z_coef(ksize,isize,jsize,3) )
    allocate( Zinterp_z2xi_idx (ksize,isize,jsize,2) )
    allocate( Zinterp_z2xi_coef(ksize,isize,jsize,3) )

    allocate( Zinterp_xi2p_idx (ksize_pres,isize,jsize,2) )
    allocate( Zinterp_xi2p_coef(ksize_pres,isize,jsize,3) )

    return
  end subroutine Zinterp_setindex

  !-----------------------------------------------------------------------------
  subroutine Zinterp_setcoef_HGT( &
       ksize, &
       isize, &
       jsize, &
       Xi_C,  &
       Xi_F,  &
       Z_C,   &
       Z_F    )
    implicit none

    integer,  intent(in) :: ksize
    integer,  intent(in) :: isize
    integer,  intent(in) :: jsize
    real(DP), intent(in) :: Xi_C(  ksize)
    real(DP), intent(in) :: Xi_F(0:ksize)
    real(DP), intent(in) :: Z_C (  ksize,isize,jsize)
    real(DP), intent(in) :: Z_F (0:ksize,isize,jsize)

    integer :: k, i, j, kk, kp
    !---------------------------------------------------------------------------

    do j = 1, jsize
    do i = 1, isize
    do k = 1, ksize
       if ( Xi_C(k) <= Z_F(0,i,j) ) then

          Zinterp_xi2z_idx (k,i,j,1) = 1      ! dummmy
          Zinterp_xi2z_idx (k,i,j,2) = 1      ! dummmy
          Zinterp_xi2z_coef(k,i,j,1) = 0.0_DP
          Zinterp_xi2z_coef(k,i,j,2) = 0.0_DP
          Zinterp_xi2z_coef(k,i,j,3) = 1.0_DP ! set UNDEF

       elseif( Xi_C(k) <= Z_C(1,i,j) ) then

          Zinterp_xi2z_idx (k,i,j,1) = 1      ! dummmy
          Zinterp_xi2z_idx (k,i,j,2) = 1
          Zinterp_xi2z_coef(k,i,j,1) = 0.0_DP
          Zinterp_xi2z_coef(k,i,j,2) = 1.0_DP
          Zinterp_xi2z_coef(k,i,j,3) = 0.0_DP

       elseif( Xi_C(k) > Z_C(ksize,i,j) ) then

          Zinterp_xi2z_idx (k,i,j,1) = ksize
          Zinterp_xi2z_idx (k,i,j,2) = ksize   ! dummmy
          Zinterp_xi2z_coef(k,i,j,1) = 1.0_DP
          Zinterp_xi2z_coef(k,i,j,2) = 0.0_DP
          Zinterp_xi2z_coef(k,i,j,3) = 0.0_DP

       elseif( Xi_C(k) > Z_F(ksize,i,j) ) then

          Zinterp_xi2z_idx (k,i,j,1) = ksize   ! dummmy
          Zinterp_xi2z_idx (k,i,j,2) = ksize   ! dummmy
          Zinterp_xi2z_coef(k,i,j,1) = 0.0_DP
          Zinterp_xi2z_coef(k,i,j,2) = 0.0_DP
          Zinterp_xi2z_coef(k,i,j,3) = 1.0_DP ! set UNDEF

       else

          do kk = 2, ksize
             kp = kk
             if( Xi_C(k) <= Z_C(kk,i,j) ) exit
          enddo

          Zinterp_xi2z_idx (k,i,j,1) = kp - 1
          Zinterp_xi2z_idx (k,i,j,2) = kp
          Zinterp_xi2z_coef(k,i,j,1) = ( Z_C (kp,i,j) - Xi_C(k)        ) &
                                     / ( Z_C (kp,i,j) - Z_C (kp-1,i,j) )
          Zinterp_xi2z_coef(k,i,j,2) = ( Xi_C(k)      - Z_C (kp-1,i,j) ) &
                                     / ( Z_C (kp,i,j) - Z_C (kp-1,i,j) )
          Zinterp_xi2z_coef(k,i,j,3) = 0.0_DP

       endif
    enddo
    enddo
    enddo

    do j = 1, jsize
    do i = 1, isize
    do k = 1, ksize
       if ( Z_C(k,i,j) <= Xi_F(0) ) then

          Zinterp_z2xi_idx (k,i,j,1) = 1      ! dummmy
          Zinterp_z2xi_idx (k,i,j,2) = 1      ! dummmy
          Zinterp_z2xi_coef(k,i,j,1) = 0.0_DP
          Zinterp_z2xi_coef(k,i,j,2) = 0.0_DP
          Zinterp_z2xi_coef(k,i,j,3) = 1.0_DP ! set UNDEF

       elseif( Z_C(k,i,j) <= Xi_C(1) ) then

          Zinterp_z2xi_idx (k,i,j,1) = 1      ! dummmy
          Zinterp_z2xi_idx (k,i,j,2) = 1
          Zinterp_z2xi_coef(k,i,j,1) = 0.0_DP
          Zinterp_z2xi_coef(k,i,j,2) = 1.0_DP
          Zinterp_z2xi_coef(k,i,j,3) = 0.0_DP

       elseif( Z_C(k,i,j) > Xi_C(ksize) ) then

          Zinterp_z2xi_idx (k,i,j,1) = ksize
          Zinterp_z2xi_idx (k,i,j,2) = ksize  ! dummmy
          Zinterp_z2xi_coef(k,i,j,1) = 1.0_DP
          Zinterp_z2xi_coef(k,i,j,2) = 0.0_DP
          Zinterp_z2xi_coef(k,i,j,3) = 0.0_DP

       elseif( Z_C(k,i,j) > Xi_F(ksize) ) then

          Zinterp_z2xi_idx (k,i,j,1) = ksize  ! dummmy
          Zinterp_z2xi_idx (k,i,j,2) = ksize  ! dummmy
          Zinterp_z2xi_coef(k,i,j,1) = 0.0_DP
          Zinterp_z2xi_coef(k,i,j,2) = 0.0_DP
          Zinterp_z2xi_coef(k,i,j,3) = 1.0_DP ! set UNDEF

       else

          do kk = 1+1, ksize
             kp = kk
             if( Z_C(k,i,j) <= Xi_C(kk) ) exit
          enddo

          Zinterp_z2xi_idx (k,i,j,1) = kp - 1
          Zinterp_z2xi_idx (k,i,j,2) = kp
          Zinterp_z2xi_coef(k,i,j,1) = ( Xi_C(kp)    - Z_C (k,i,j) ) &
                                     / ( Xi_C(kp)    - Xi_C(kp-1)  )
          Zinterp_z2xi_coef(k,i,j,2) = ( Z_C (k,i,j) - Xi_C(kp-1)  ) &
                                     / ( Xi_C(kp)    - Xi_C(kp-1)  )
          Zinterp_z2xi_coef(k,i,j,3) = 0.0_DP

       endif
    enddo
    enddo
    enddo

    return
  end subroutine Zinterp_setcoef_HGT

  !-----------------------------------------------------------------------------
  subroutine Zinterp_xi2zSP( &
       ksize, &
       isize, &
       jsize, &
       var,   &
       var_Z  )
    use gtool_file, only: &
       RMISS
    implicit none

    integer,  intent(in)  :: ksize
    integer,  intent(in)  :: isize
    integer,  intent(in)  :: jsize
    real(SP), intent(in)  :: var  (ksize,isize,jsize)
    real(SP), intent(out) :: var_Z(ksize,isize,jsize)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = 1, jsize
    do i = 1, isize
    do k = 1, ksize
       var_Z(k,i,j) = Zinterp_xi2z_coef(k,i,j,1) * var(Zinterp_xi2z_idx(k,i,j,1),i,j) &
                    + Zinterp_xi2z_coef(k,i,j,2) * var(Zinterp_xi2z_idx(k,i,j,2),i,j) &
                    + Zinterp_xi2z_coef(k,i,j,3) * RMISS
    enddo
    enddo
    enddo

    return
  end subroutine Zinterp_xi2zSP

  !-----------------------------------------------------------------------------
  subroutine Zinterp_xi2zDP( &
       ksize, &
       isize, &
       jsize, &
       var,   &
       var_Z  )
    use gtool_file, only: &
       RMISS
    implicit none

    integer,  intent(in)  :: ksize
    integer,  intent(in)  :: isize
    integer,  intent(in)  :: jsize
    real(DP), intent(in)  :: var  (ksize,isize,jsize)
    real(DP), intent(out) :: var_Z(ksize,isize,jsize)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = 1, jsize
    do i = 1, isize
    do k = 1, ksize
       var_Z(k,i,j) = Zinterp_xi2z_coef(k,i,j,1) * var(Zinterp_xi2z_idx(k,i,j,1),i,j) &
                    + Zinterp_xi2z_coef(k,i,j,2) * var(Zinterp_xi2z_idx(k,i,j,2),i,j) &
                    + Zinterp_xi2z_coef(k,i,j,3) * RMISS
    enddo
    enddo
    enddo

    return
  end subroutine Zinterp_xi2zDP

  !-----------------------------------------------------------------------------
  subroutine Zinterp_z2xiSP( &
       ksize, &
       isize, &
       jsize, &
       var,   &
       var_Xi )
    use gtool_file, only: &
       RMISS
    implicit none

    integer,  intent(in)  :: ksize
    integer,  intent(in)  :: isize
    integer,  intent(in)  :: jsize
    real(SP), intent(in)  :: var   (ksize,isize,jsize)
    real(SP), intent(out) :: var_Xi(ksize,isize,jsize)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = 1, jsize
    do i = 1, isize
    do k = 1, ksize
       var_Xi(k,i,j) = Zinterp_z2xi_coef(k,i,j,1) * var(Zinterp_z2xi_idx(k,i,j,1),i,j) &
                     + Zinterp_z2xi_coef(k,i,j,2) * var(Zinterp_z2xi_idx(k,i,j,2),i,j) &
                     + Zinterp_z2xi_coef(k,i,j,3) * RMISS
    enddo
    enddo
    enddo

    return
  end subroutine Zinterp_z2xiSP

  !-----------------------------------------------------------------------------
  subroutine Zinterp_z2xiDP( &
       ksize, &
       isize, &
       jsize, &
       var,   &
       var_Xi )
    use gtool_file, only: &
       RMISS
    implicit none

    integer,  intent(in)  :: ksize
    integer,  intent(in)  :: isize
    integer,  intent(in)  :: jsize
    real(DP), intent(in)  :: var   (ksize,isize,jsize)
    real(DP), intent(out) :: var_Xi(ksize,isize,jsize)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = 1, jsize
    do i = 1, isize
    do k = 1, ksize
       var_Xi(k,i,j) = Zinterp_z2xi_coef(k,i,j,1) * var(Zinterp_z2xi_idx(k,i,j,1),i,j) &
                     + Zinterp_z2xi_coef(k,i,j,2) * var(Zinterp_z2xi_idx(k,i,j,2),i,j) &
                     + Zinterp_z2xi_coef(k,i,j,3) * RMISS
    enddo
    enddo
    enddo

    return
  end subroutine Zinterp_z2xiDP

  !-----------------------------------------------------------------------------
  subroutine Zinterp_setcoef_PRES( &
       ksize,      &
       isize,      &
       jsize,      &
       ksize_pres, &
       Pres_atm,   &
       Pres_sfc,   &
       Pres_out    )
    implicit none

    integer,  intent(in) :: ksize
    integer,  intent(in) :: isize
    integer,  intent(in) :: jsize
    integer,  intent(in) :: ksize_pres
    real(DP), intent(in) :: Pres_atm(ksize,isize,jsize) ! pressure in Xi coordinate [hPa]
    real(DP), intent(in) :: Pres_sfc(      isize,jsize) ! surface pressure          [hPa]
    real(DP), intent(in) :: Pres_out(ksize_pres)        ! pressure level to output  [hPa]

    real(DP) :: LnPres_atm(ksize,isize,jsize) ! (log) pressure in Xi coordinate [hPa]
    real(DP) :: LnPres_sfc(      isize,jsize) ! (log) surface pressure          [hPa]
    real(DP) :: LnPres_out(ksize_pres)        ! (log) pressure level to output  [hPa]

    integer :: k, i, j, kk, kp
    !---------------------------------------------------------------------------

    LnPres_atm(:,:,:) = dlog( Pres_atm(:,:,:) )
    LnPres_sfc(:,:)   = dlog( Pres_sfc(:,:)   )
    LnPres_out(:)     = dlog( Pres_out(:)     )

    do j = 1, jsize
    do i = 1, isize
    do k = 1, ksize_pres
       if ( LnPres_out(k) >= LnPres_sfc(i,j) ) then

          Zinterp_xi2p_idx (k,i,j,1) = 1      ! dummmy
          Zinterp_xi2p_idx (k,i,j,2) = 1      ! dummmy
          Zinterp_xi2p_coef(k,i,j,1) = 0.0_DP
          Zinterp_xi2p_coef(k,i,j,2) = 0.0_DP
          Zinterp_xi2p_coef(k,i,j,3) = 1.0_DP ! set UNDEF

       elseif( LnPres_out(k) >= LnPres_atm(1,i,j) ) then

          Zinterp_xi2p_idx (k,i,j,1) = 1      ! dummmy
          Zinterp_xi2p_idx (k,i,j,2) = 1
          Zinterp_xi2p_coef(k,i,j,1) = 0.0_DP
          Zinterp_xi2p_coef(k,i,j,2) = 1.0_DP
          Zinterp_xi2p_coef(k,i,j,3) = 0.0_DP

       elseif( LnPres_out(k) < LnPres_atm(ksize,i,j) ) then

          Zinterp_xi2p_idx (k,i,j,1) = ksize   ! dummmy
          Zinterp_xi2p_idx (k,i,j,2) = ksize   ! dummmy
          Zinterp_xi2p_coef(k,i,j,1) = 0.0_DP
          Zinterp_xi2p_coef(k,i,j,2) = 0.0_DP
          Zinterp_xi2p_coef(k,i,j,3) = 1.0_DP ! set UNDEF

       else

          do kk = 2, ksize
             kp = kk
             if( LnPres_out(k) >= LnPres_atm(kk,i,j) ) exit
          enddo

          Zinterp_xi2p_idx (k,i,j,1) = kp - 1
          Zinterp_xi2p_idx (k,i,j,2) = kp
          Zinterp_xi2p_coef(k,i,j,1) = ( LnPres_out(k)        - LnPres_atm(kp,i,j) ) &
                                     / ( LnPres_atm(kp-1,i,j) - LnPres_atm(kp,i,j) )
          Zinterp_xi2p_coef(k,i,j,2) = ( LnPres_atm(kp-1,i,j) - LnPres_out(k)      ) &
                                     / ( LnPres_atm(kp-1,i,j) - LnPres_atm(kp,i,j) )
          Zinterp_xi2p_coef(k,i,j,3) = 0.0_DP

       endif
    enddo
    enddo
    enddo

    return
  end subroutine Zinterp_setcoef_PRES

  !-----------------------------------------------------------------------------
  subroutine Zinterp_xi2pSP( &
       ksize,      &
       isize,      &
       jsize,      &
       ksize_pres, &
       var,        &
       var_P       )
    use gtool_file, only: &
       RMISS
    implicit none

    integer,  intent(in)  :: ksize
    integer,  intent(in)  :: isize
    integer,  intent(in)  :: jsize
    integer,  intent(in)  :: ksize_pres
    real(SP), intent(in)  :: var  (:,:,:)
    real(SP), intent(out) :: var_P(:,:,:)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = 1, jsize
    do i = 1, isize
    do k = 1, ksize_pres
       var_P(k,i,j) = Zinterp_xi2p_coef(k,i,j,1) * var(Zinterp_xi2p_idx(k,i,j,1),i,j) &
                    + Zinterp_xi2p_coef(k,i,j,2) * var(Zinterp_xi2p_idx(k,i,j,2),i,j) &
                    + Zinterp_xi2p_coef(k,i,j,3) * RMISS
    enddo
    enddo
    enddo

    return
  end subroutine Zinterp_xi2pSP

  !-----------------------------------------------------------------------------
  subroutine Zinterp_xi2pDP( &
       ksize,      &
       isize,      &
       jsize,      &
       ksize_pres, &
       var,        &
       var_P       )
    use gtool_file, only: &
       RMISS
    implicit none

    integer,  intent(in)  :: ksize
    integer,  intent(in)  :: isize
    integer,  intent(in)  :: jsize
    integer,  intent(in)  :: ksize_pres
    real(DP), intent(in)  :: var  (ksize     ,isize,jsize)
    real(DP), intent(out) :: var_P(ksize_pres,isize,jsize)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = 1, jsize
    do i = 1, isize
    do k = 1, ksize_pres
       var_P(k,i,j) = Zinterp_xi2p_coef(k,i,j,1) * var(Zinterp_xi2p_idx(k,i,j,1),i,j) &
                    + Zinterp_xi2p_coef(k,i,j,2) * var(Zinterp_xi2p_idx(k,i,j,2),i,j) &
                    + Zinterp_xi2p_coef(k,i,j,3) * RMISS
    enddo
    enddo
    enddo

    return
  end subroutine Zinterp_xi2pDP

end module gtool_zinterp
