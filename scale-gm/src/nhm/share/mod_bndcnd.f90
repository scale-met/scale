!-------------------------------------------------------------------------------
!> Module boundary conditions
!!
!! @par Description
!!          This module provides the subroutines for boundary conditions
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_bndcnd
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID,  &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: BNDCND_setup

  public :: BNDCND_all
  public :: BNDCND_thermo
  public :: BNDCND_rhovxvyvz
  public :: BNDCND_rhow

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
  !--- Vertical boundary condition for temperature at the top
  character(len=ADM_NSYS), private :: BND_TYPE_T_TOP    != 'TEM' : tem(kmax+1) = tem(kmax)
                                                        != 'EPL' : lagrange extrapolation

  !--- Vertical boundary condition for temperature at the ground
  character(len=ADM_NSYS), private :: BND_TYPE_T_BOTTOM != 'FIX' : tems fix
                                                        != 'TEM' : tem(kmin-1) = tem(kmin)
                                                        != 'EPL' : lagrange extrapolation

  !--- Vertical boundary condition for momentum at the top
  character(len=ADM_NSYS), private :: BND_TYPE_M_TOP    != 'RIGID' : rigid surface
                                                        != 'FREE'  : free surface

  !--- Vertical boundary condition for momentum at the ground
  character(len=ADM_NSYS), private :: BND_TYPE_M_BOTTOM != 'RIGID' : rigid surface
                                                        != 'FREE'  : free surface

  logical, private :: is_top_tem   = .true.
  logical, private :: is_top_epl   = .false.
  logical, private :: is_btm_fix   = .false.
  logical, private :: is_btm_tem   = .true.
  logical, private :: is_btm_epl   = .false.
  logical, private :: is_top_rigid = .false.
  logical, private :: is_top_free  = .true.
  logical, private :: is_btm_rigid = .true.
  logical, private :: is_btm_free  = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine BNDCND_setup
    use mod_adm, only: &
       ADM_CTL_FID, &
       ADM_proc_stop
    use mod_cnst, only: &
       CNST_TEMS0
    implicit none

    namelist / BNDCNDPARAM / &
       BND_TYPE_T_TOP,    &
       BND_TYPE_T_BOTTOM, &
       BND_TYPE_M_TOP,    &
       BND_TYPE_M_BOTTOM

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- set default
    BND_TYPE_T_TOP    = 'TEM'
    BND_TYPE_T_BOTTOM = 'TEM'
    BND_TYPE_M_TOP    = 'FREE'
    BND_TYPE_M_BOTTOM = 'RIGID'

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[bndcnd]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=BNDCNDPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** BNDCNDPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist BNDCNDPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist BNDCNDPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=BNDCNDPARAM)

    is_top_tem = .false.
    is_top_epl = .false.
    is_btm_fix = .false.
    is_btm_tem = .false.
    is_btm_epl = .false.

    if    ( BND_TYPE_T_TOP == 'TEM' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, top   ) : equal to uppermost atmosphere'
       is_top_tem = .true.
    elseif( BND_TYPE_T_TOP == 'EPL' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, top   ) : lagrange extrapolation'
       is_top_epl = .true.
    else
       write(ADM_LOG_FID,*) 'xxx Invalid BND_TYPE_T_TOP. STOP.'
       call ADM_proc_stop
    endif

    if    ( BND_TYPE_T_BOTTOM == 'FIX' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, bottom) : fixed'
       write(ADM_LOG_FID,*) '***           boundary temperature (CNST_TEMS0) : ', CNST_TEMS0
       is_btm_fix = .true.
    elseif( BND_TYPE_T_BOTTOM == 'TEM' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, bottom) : equal to lowermost atmosphere'
       is_btm_tem = .true.
    elseif( BND_TYPE_T_BOTTOM == 'EPL' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (temperature, bottom) : lagrange extrapolation'
       is_btm_epl = .true.
    else
       write(ADM_LOG_FID,*) 'xxx Invalid BND_TYPE_T_BOTTOM. STOP.'
       call ADM_proc_stop
    endif

    is_top_rigid = .false.
    is_top_free  = .false.
    is_btm_rigid = .false.
    is_btm_free  = .false.

    if    ( BND_TYPE_M_TOP == 'RIGID' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (momentum,    top   ) : rigid'
       is_top_rigid = .true.
    elseif( BND_TYPE_M_TOP == 'FREE' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (momentum,    top   ) : free'
       is_top_free  = .true.
    else
       write(ADM_LOG_FID,*) 'xxx Invalid BND_TYPE_M_TOP. STOP.'
       call ADM_proc_stop
    endif

    if    ( BND_TYPE_M_BOTTOM == 'RIGID' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (momentum,    bottom) : rigid'
       is_btm_rigid = .true.
    elseif( BND_TYPE_M_BOTTOM == 'FREE' ) then
       write(ADM_LOG_FID,*) '*** Boundary setting type (momentum,    bottom) : free'
       is_btm_free  = .true.
    else
       write(ADM_LOG_FID,*) 'xxx Invalid BND_TYPE_M_BOTTOM. STOP.'
       call ADM_proc_stop
    endif

    return
  end subroutine BNDCND_setup

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for all prognostic variables.
  subroutine BNDCND_all( &
       ijdim,      &
       rho,        &
       vx,         &
       vy,         &
       vz,         &
       w,          &
       ein,        &
       tem,        &
       pre,        &
       rhog,       &
       rhogvx,     &
       rhogvy,     &
       rhogvz,     &
       rhogw,      &
       rhoge,      &
       gsqrtgam2,  &
       phi,        &
       c2wfact,    &
       c2wfact_Gz  )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       CVdry => CNST_CV
    implicit none

    integer, intent(in)    :: ijdim           ! number of horizontal grid
    real(RP), intent(inout) :: rho(ijdim,kdim) ! density
    real(RP), intent(inout) :: vx (ijdim,kdim) ! horizontal wind (x)
    real(RP), intent(inout) :: vy (ijdim,kdim) ! horizontal wind (y)
    real(RP), intent(inout) :: vz (ijdim,kdim) ! horizontal wind (z)
    real(RP), intent(inout) :: w  (ijdim,kdim) ! vertical   wind
    real(RP), intent(inout) :: ein(ijdim,kdim) ! internal energy
    real(RP), intent(inout) :: tem(ijdim,kdim) ! temperature
    real(RP), intent(inout) :: pre(ijdim,kdim) ! pressure

    real(RP), intent(inout) :: rhog  (ijdim,kdim)
    real(RP), intent(inout) :: rhogvx(ijdim,kdim)
    real(RP), intent(inout) :: rhogvy(ijdim,kdim)
    real(RP), intent(inout) :: rhogvz(ijdim,kdim)
    real(RP), intent(inout) :: rhogw (ijdim,kdim)
    real(RP), intent(inout) :: rhoge (ijdim,kdim)

    real(RP), intent(in)    :: gsqrtgam2 (ijdim,kdim)
    real(RP), intent(in)    :: phi       (ijdim,kdim)   ! geopotential
    real(RP), intent(in)    :: c2wfact   (ijdim,kdim,2)
    real(RP), intent(in)    :: c2wfact_Gz(ijdim,kdim,6)

    integer :: ij
    !---------------------------------------------------------------------------

    !$acc  data &
    !$acc& pcopy(rho,ein,tem,pre,rhog,rhoge) &
    !$acc& pcopy(vx,vy,vz,w,rhogvx,rhogvy,rhogvz,rhogw) &
    !$acc& pcopyin(gsqrtgam2,phi,c2wfact,c2wfact_Gz)

    !--- Thermodynamical variables ( rho, ein, tem, pre, rhog, rhoge ), q = 0 at boundary
    call BNDCND_thermo( ijdim, & ! [IN]
                        tem,   & ! [INOUT]
                        rho,   & ! [INOUT]
                        pre,   & ! [INOUT]
                        phi    ) ! [IN]

    !$acc kernels pcopy(rhog,ein,rhoge) pcopyin(rho,tem,gsqrtgam2) async(0)
    do ij = 1, ijdim
       rhog (ij,kmax+1) = rho(ij,kmax+1) * gsqrtgam2(ij,kmax+1)
       rhog (ij,kmin-1) = rho(ij,kmin-1) * gsqrtgam2(ij,kmin-1)

       ein  (ij,kmax+1) = CVdry * tem(ij,kmax+1)
       ein  (ij,kmin-1) = CVdry * tem(ij,kmin-1)

       rhoge(ij,kmax+1) = rhog(ij,kmax+1) * ein(ij,kmax+1)
       rhoge(ij,kmin-1) = rhog(ij,kmin-1) * ein(ij,kmin-1)
    enddo
    !$acc end kernels

    !--- Momentum ( rhogvx, rhogvy, rhogvz, vx, vy, vz )
    call BNDCND_rhovxvyvz( ijdim,  & ! [IN]
                           rhog,   & ! [IN]
                           rhogvx, & ! [INOUT]
                           rhogvy, & ! [INOUT]
                           rhogvz  ) ! [INOUT]

    !$acc kernels pcopy(vx,vy,vz) pcopyin(rhogvx,rhogvy,rhogvz,rhog) async(0)
    do ij = 1, ijdim
       vx(ij,kmax+1) = rhogvx(ij,kmax+1) / rhog(ij,kmax+1)
       vy(ij,kmax+1) = rhogvy(ij,kmax+1) / rhog(ij,kmax+1)
       vz(ij,kmax+1) = rhogvz(ij,kmax+1) / rhog(ij,kmax+1)

       vx(ij,kmin-1) = rhogvx(ij,kmin-1) / rhog(ij,kmin-1)
       vy(ij,kmin-1) = rhogvy(ij,kmin-1) / rhog(ij,kmin-1)
       vz(ij,kmin-1) = rhogvz(ij,kmin-1) / rhog(ij,kmin-1)
    enddo
    !$acc end kernels

    !--- Momentum ( rhogw, w )
    call BNDCND_rhow( ijdim,     & ! [IN]
                      rhogvx,    & ! [IN]
                      rhogvy,    & ! [IN]
                      rhogvz,    & ! [IN]
                      rhogw,     & ! [INOUT]
                      c2wfact_Gz ) ! [IN]

    !$acc kernels pcopy(w) pcopyin(rhog,rhogw,c2wfact) async(0)
    do ij = 1, ijdim

       w(ij,kmax+1) = rhogw(ij,kmax+1) / ( c2wfact(ij,kmax+1,1) * rhog(ij,kmax+1) &
                                         + c2wfact(ij,kmax+1,2) * rhog(ij,kmax  ) )

       w(ij,kmin  ) = rhogw(ij,kmin  ) / ( c2wfact(ij,kmin  ,1) * rhog(ij,kmin  ) &
                                         + c2wfact(ij,kmin  ,2) * rhog(ij,kmin-1) )

       w(ij,kmin-1) = 0.0_RP

    enddo
    !$acc end kernels

    !$acc end data

    return
  end subroutine BNDCND_all

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for thermodynamical variables
  subroutine BNDCND_thermo( &
       ijdim, &
       tem,   &
       rho,   &
       pre,   &
       phi    )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       CNST_RAIR,  &
       CNST_EGRAV, &
       CNST_TEMS0
    implicit none

    integer, intent(in)    :: ijdim           ! number of horizontal grid
    real(RP), intent(inout) :: tem(ijdim,kdim) ! temperature
    real(RP), intent(inout) :: rho(ijdim,kdim) ! density
    real(RP), intent(inout) :: pre(ijdim,kdim) ! pressure
    real(RP), intent(in)    :: phi(ijdim,kdim) ! geopotential

    integer :: ij

    real(RP) :: z,z1,z2,z3,p1,p2,p3
    real(RP) :: lag_intpl
    lag_intpl(z,z1,p1,z2,p2,z3,p3) = ( (z-z2)*(z-z3))/((z1-z2)*(z1-z3) ) * p1 &
                                   + ( (z-z1)*(z-z3))/((z2-z1)*(z2-z3) ) * p2 &
                                   + ( (z-z1)*(z-z2))/((z3-z1)*(z3-z2) ) * p3
    !---------------------------------------------------------------------------

    !$acc data pcopy(tem,rho,pre) pcopyin(phi)

    !$acc kernels pcopy(tem) pcopyin(phi) async(0)
    do ij = 1, ijdim

       if    ( is_top_tem ) then

          tem(ij,kmax+1) = tem(ij,kmax) ! dT/dz = 0

       elseif( is_top_epl ) then

          z  = phi(ij,kmax+1) / CNST_EGRAV
          z1 = phi(ij,kmax  ) / CNST_EGRAV
          z2 = phi(ij,kmax-1) / CNST_EGRAV
          z3 = phi(ij,kmax-2) / CNST_EGRAV

          tem(ij,kmax+1) = lag_intpl( z ,                 &
                                      z1, tem(ij,kmax  ), &
                                      z2, tem(ij,kmax-1), &
                                      z3, tem(ij,kmax-2)  )

       endif

       if    ( is_btm_fix ) then
          tem(ij,kmin-1) = CNST_TEMS0
       elseif( is_btm_tem ) then
          tem(ij,kmin-1) = tem(ij,kmin) ! dT/dz = 0
       elseif( is_btm_epl ) then
          z1 = phi(ij,kmin+2) / CNST_EGRAV
          z2 = phi(ij,kmin+1) / CNST_EGRAV
          z3 = phi(ij,kmin  ) / CNST_EGRAV
          z  = phi(ij,kmin-1) / CNST_EGRAV

          tem(ij,kmin-1) = lag_intpl( z,                  &
                                      z1, tem(ij,kmin+2), &
                                      z2, tem(ij,kmin+1), &
                                      z3, tem(ij,kmin  )  )
       endif

    enddo
    !$acc end kernels

    !$acc kernels pcopy(pre,rho) pcopyin(phi,tem) async(0)
    do ij = 1, ijdim

       !--- set the boundary of pressure ( hydrostatic balance )
       pre(ij,kmax+1) = pre(ij,kmax-1) - rho(ij,kmax) * ( phi(ij,kmax+1) - phi(ij,kmax-1) )
       pre(ij,kmin-1) = pre(ij,kmin+1) - rho(ij,kmin) * ( phi(ij,kmin-1) - phi(ij,kmin+1) )

       !--- set the boundary of density ( equation of state )
       rho(ij,kmax+1) = pre(ij,kmax+1) / ( CNST_RAIR * tem(ij,kmax+1) )
       rho(ij,kmin-1) = pre(ij,kmin-1) / ( CNST_RAIR * tem(ij,kmin-1) )

    enddo
    !$acc end kernels

    !$acc end data

    return
  end subroutine BNDCND_thermo

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for horizontal momentum
  subroutine BNDCND_rhovxvyvz( &
       ijdim,  &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogvz  )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    implicit none

    integer, intent(in)    :: ijdim
    real(RP), intent(in)    :: rhog  (ijdim,kdim)
    real(RP), intent(inout) :: rhogvx(ijdim,kdim)
    real(RP), intent(inout) :: rhogvy(ijdim,kdim)
    real(RP), intent(inout) :: rhogvz(ijdim,kdim)

    integer :: ij
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(rhogvx,rhogvy,rhogvz) pcopyin(rhog) async(0)
    do ij = 1, ijdim

       if    ( is_top_rigid ) then
          rhogvx(ij,kmax+1) = -rhogvx(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
          rhogvy(ij,kmax+1) = -rhogvy(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
          rhogvz(ij,kmax+1) = -rhogvz(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
       elseif( is_top_free  ) then
          rhogvx(ij,kmax+1) =  rhogvx(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
          rhogvy(ij,kmax+1) =  rhogvy(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
          rhogvz(ij,kmax+1) =  rhogvz(ij,kmax) / rhog(ij,kmax) * rhog(ij,kmax+1)
       endif

       if    ( is_btm_rigid ) then
          rhogvx(ij,kmin-1) = -rhogvx(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
          rhogvy(ij,kmin-1) = -rhogvy(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
          rhogvz(ij,kmin-1) = -rhogvz(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
       elseif( is_btm_free  ) then
          rhogvx(ij,kmin-1) =  rhogvx(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
          rhogvy(ij,kmin-1) =  rhogvy(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
          rhogvz(ij,kmin-1) =  rhogvz(ij,kmin) / rhog(ij,kmin) * rhog(ij,kmin-1)
       endif

    enddo
    !$acc end kernels

    return
  end subroutine BNDCND_rhovxvyvz

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for vertical momentum
  subroutine BNDCND_rhow( &
       ijdim,  &
       rhogvx, &
       rhogvy, &
       rhogvz, &
       rhogw,  &
       c2wfact )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    implicit none

    integer, intent(in)    :: ijdim
    real(RP), intent(in)    :: rhogvx (ijdim,kdim)
    real(RP), intent(in)    :: rhogvy (ijdim,kdim)
    real(RP), intent(in)    :: rhogvz (ijdim,kdim)
    real(RP), intent(inout) :: rhogw  (ijdim,kdim)
    real(RP), intent(in)    :: c2wfact(ijdim,kdim,6)

    integer :: ij
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(rhogw) pcopyin(rhogvx,rhogvy,rhogvz,c2wfact) async(0)
    do ij = 1, ijdim

       if    ( is_top_rigid ) then
          rhogw(ij,kmax+1) = 0.0_RP
       elseif( is_top_free  ) then
          rhogw(ij,kmax+1) = - ( c2wfact(ij,kmax+1,1) * rhogvx(ij,kmax+1) &
                               + c2wfact(ij,kmax+1,2) * rhogvx(ij,kmax  ) &
                               + c2wfact(ij,kmax+1,3) * rhogvy(ij,kmax+1) &
                               + c2wfact(ij,kmax+1,4) * rhogvy(ij,kmax  ) &
                               + c2wfact(ij,kmax+1,5) * rhogvz(ij,kmax+1) &
                               + c2wfact(ij,kmax+1,6) * rhogvz(ij,kmax  ) )
       endif

       if    ( is_btm_rigid ) then
          rhogw(ij,kmin) = 0.0_RP
       elseif( is_btm_free  ) then
          rhogw(ij,kmin) = - ( c2wfact(ij,kmin,1) * rhogvx(ij,kmin  ) &
                             + c2wfact(ij,kmin,2) * rhogvx(ij,kmin-1) &
                             + c2wfact(ij,kmin,3) * rhogvy(ij,kmin  ) &
                             + c2wfact(ij,kmin,4) * rhogvy(ij,kmin-1) &
                             + c2wfact(ij,kmin,5) * rhogvz(ij,kmin  ) &
                             + c2wfact(ij,kmin,6) * rhogvz(ij,kmin-1) )
       endif

       rhogw(ij,kmin-1) = 0.0_RP

    enddo
    !$acc end kernels

    return
  end subroutine BNDCND_rhow

end module mod_bndcnd
