!-------------------------------------------------------------------------------
!
!+  saturation adjustment module
!
!-------------------------------------------------------------------------------
module mod_satadjust
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines for the saturation adjustment.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      11/10/24  T.Seiki, Import from NICAM
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: moist_psat_water
  public :: moist_psat_water1
  public :: moist_psat_ice
  public :: moist_psat_ice1
  public :: moist_relative_humidity
  !
  public :: moist_qsat_water
  public :: moist_qsat_water1
  public :: moist_qsat_ice
  public :: moist_qsat_ice1
  public :: moist_dqsw_dtem_rho  
  public :: moist_dqsi_dtem_rho  
  public :: moist_dqsw_dtem_dpre 
  public :: moist_dqsi_dtem_dpre 
  !
  public :: moist_dewtem 
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  integer, parameter, private :: ADM_NSYS=32
  !
  real(8), parameter, private :: QSAT_MIN  = 1.d-10
  real(8), parameter, private :: RELH_MAX  = 1.d0
  real(8), parameter, private :: DTEM_EPS0 = 1.0D-8 ! temperature convergence criterion
  real(8), parameter, private :: SMALL     = 1.0D-10
  !
  !-----------------------------------------------------------------------------
contains
  !
  subroutine moist_psat_water ( tem, psat )
    ! psat : Clasius-Clapeyron: based on CPV, CPL constant
    !
    use mod_atmos_cnst, only :        &
         CNST_RVAP,             &
         CNST_CPV,              &
         CNST_CL,               &
         CNST_LH0,              &
         CNST_LH00,             &
         CNST_PSAT0,            &
         CNST_TEM00
    !
    real(8), intent(in) :: tem(:,:)
    real(8), intent(out) :: psat(:,:)
    real(8) :: TEM_MIN = 10.d0

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_psat_water")
#endif

    
    psat(:,:) = CNST_PSAT0 &
         * ( max(tem(:,:), TEM_MIN) / CNST_TEM00 ) ** ( ( CNST_CPV - CNST_CL ) / CNST_RVAP ) &
         * exp ( CNST_LH00 / CNST_RVAP &
         * ( 1.0d0 / CNST_TEM00 - 1.0d0 / max(tem(:,:), TEM_MIN) ) )

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_psat_water")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_psat_water
  !-----------------------------------------------------------------------------
  subroutine moist_psat_ice ( tem, psat )
    ! psat : Clasius-Clapeyron: based on CPV, CPL constant
    ! for ice
    !
    use mod_atmos_cnst, only :        &
         CNST_RVAP,             &
         CNST_CPV,              &
         CNST_CI,               &
         CNST_LHS00,            &
         CNST_LHS0,             &
         CNST_PSAT0,            &
         CNST_TEM00
    !
    real(8), intent(in) :: tem(:,:)
    real(8), intent(out) :: psat(:,:)
    !
    real(8) :: TEM_MIN = 10.d0
    !

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_psat_ice")
#endif

    psat(:,:) = CNST_PSAT0 &
         * ( max(tem(:,:), TEM_MIN) / CNST_TEM00 ) ** ( ( CNST_CPV - CNST_CI ) / CNST_RVAP ) &
         * exp ( CNST_LHS00 / CNST_RVAP &
         * ( 1.0d0 / CNST_TEM00 - 1.0d0 / max(tem(:,:), TEM_MIN) ) )

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_psat_ice")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_psat_ice
  !-----------------------------------------------------------------------------
  subroutine moist_psat_ice1 ( tem, psat )
    !
    real(8), intent(in) :: tem(:)
    real(8), intent(out) :: psat(:)
    real(8) :: psat1(size(psat),1)
#ifdef _FPCOLL_
call START_COLLECTION("moist_psat_ice1")
#endif
    call moist_psat_ice ( spread(tem,2,1), psat1 )
    psat(:) = psat1(:,1)
#ifdef _FPCOLL_
call STOP_COLLECTION("moist_psat_ice1")
#endif
  end subroutine moist_psat_ice1
  ! psat : based on CCSR/NIES AGCM
  ! qsat is more accurate than that of CCSR/NIES AGCM
  subroutine moist_qsat_water ( tem, pre, qsat )
    !
    use mod_atmos_cnst, only :        &
         CNST_EPSV
    !
    real(8), intent(in) :: pre(:,:)
    real(8), intent(in) :: tem(:,:)
    real(8), intent(out) :: qsat(:,:)
    
    real(8) :: psat(size(tem,1),size(tem,2))
    

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_qsat_water")
#endif

    call moist_psat_water ( tem, psat )
    
    qsat(:,:) = &
         CNST_EPSV * psat(:,:) / ( pre(:,:) - ( 1.D0 - CNST_EPSV ) * psat(:,:) )

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_qsat_water")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_qsat_water
  !-----------------------------------------------------------------------------
  subroutine moist_qsat_ice ( tem, pre, qsat )
    !
    use mod_atmos_cnst, only :        &
         CNST_EPSV
    !
    real(8), intent(in) :: pre(:,:)
    real(8), intent(in) :: tem(:,:)
    real(8), intent(out) :: qsat(:,:)
    
    real(8) :: psat(size(tem,1),size(tem,2))
    

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_qsat_ice")
#endif

    call moist_psat_ice ( tem, psat )
    
    qsat(:,:) = &
         CNST_EPSV * psat(:,:) / ( pre(:,:) - ( 1.D0 - CNST_EPSV ) * psat(:,:) )

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_qsat_ice")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_qsat_ice
  !-----------------------------------------------------------------------------
  subroutine moist_qsat_ice1 ( tem, pre, qsat )
    implicit none
    real(8), intent(in) :: pre(:)
    real(8), intent(in) :: tem(:)
    real(8), intent(out) :: qsat(:)
    real(8) :: qsat1(size(qsat),1)
    !

#ifdef _FPCOLL_
call START_COLLECTION("moist_qsat_ice1")
#endif

    call moist_qsat_ice ( spread(tem,2,1), spread(pre,2,1), qsat1 )
    qsat(:) = qsat1(:,1)

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_qsat_ice1")
#endif

    return
  end subroutine moist_qsat_ice1
  !-----------------------------------------------------------------------------
  subroutine moist_psat_water1 ( tem, psat )
    implicit none
    real(8), intent(in) :: tem(:)
    real(8), intent(out) :: psat(:)
    real(8) :: psat1(size(psat),1)
    !

#ifdef _FPCOLL_
call START_COLLECTION("moist_psat_water1")
#endif

    call moist_psat_water ( spread(tem,2,1), psat1 )
    psat(:) = psat1(:,1)

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_psat_water1")
#endif

    return
  end subroutine moist_psat_water1
  !-----------------------------------------------------------------------------
  subroutine moist_qsat_water1 ( tem, pre, qsat )
    implicit none
    real(8), intent(in) :: pre(:)
    real(8), intent(in) :: tem(:)
    real(8), intent(out) :: qsat(:)
    real(8) :: qsat1(size(qsat),1)
    

#ifdef _FPCOLL_
call START_COLLECTION("moist_qsat_water1")
#endif

    call moist_qsat_water ( spread(tem,2,1), spread(pre,2,1), qsat1 )
    qsat(:) = qsat1(:,1)

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_qsat_water1")
#endif

    return
  end subroutine moist_qsat_water1
  !
  subroutine moist_relative_humidity(&
       rh,                           &
       rho,                          &
       tem,                          &
       qv                            &
       )
    use mod_atmos_cnst, only : &
         CNST_RVAP,      &
         CNST_TMELT
    implicit none
    real(8), intent(out) :: rh(:,:)
    real(8), intent(in) :: rho(:,:)
    real(8), intent(in) :: tem(:,:)
    real(8), intent(in) :: qv(:,:)
    !
    real(8) :: psat(size(rh,1),size(rh,2)) 
    real(8) :: rhw(size(rh,1),size(rh,2)) 
    real(8) :: rhi(size(rh,1),size(rh,2)) 
    real(8) :: delta(size(rh,1),size(rh,2)) 
    !

#ifdef _FPCOLL_
call START_COLLECTION("moist_relative_humidity")
#endif

    call moist_psat_water( &
         tem(:,:),         &
         psat(:,:)         &
         )
    rhw(:,:) = qv(:,:)/(psat(:,:)/(rho(:,:)*CNST_RVAP*tem(:,:)))
    !
    call moist_psat_ice( &
         tem(:,:),       &
         psat(:,:)       &
            )
    rhi(:,:) = qv(:,:)/(psat(:,:)/(rho(:,:)*CNST_RVAP*tem(:,:)))
    delta(:,:) = (sign(0.5D0,tem(:,:)-CNST_TMELT)+0.5D0)
    rh(:,:) = delta(:,:)*rhw(:,:) + (1.0D0-delta(:,:))*rhi(:,:)
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_relative_humidity")
#endif

    return
  end subroutine moist_relative_humidity
  !
  subroutine moist_dewtem( &
       ijdim, kdim_local,    & ! in
       wtem,                 & ! out
       tem, pre, qv, qd )
    !
    use mod_atmos_cnst, only: &
         CNST_EPSV, &
         CNST_TEM00, &
         CNST_PSAT0, &
         CNST_LH00,  &
         CNST_LHS00, &
         CNST_CPV,   &
         CNST_CL,    &
         CNST_CI,    &
         CNST_RVAP
    !    
    implicit none
    !    
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim_local
    !
    real(8), intent(out):: wtem(ijdim,kdim_local) ! dew point temperature
    !
    real(8), intent(in) :: tem(ijdim,kdim_local)
    real(8), intent(in) :: pre(ijdim,kdim_local)
    real(8), intent(in) :: qv(ijdim,kdim_local)
    real(8), intent(in) :: qd(ijdim,kdim_local)
    !
    real(8) :: prev(1:ijdim,1:kdim_local)
    real(8) :: esw(1:ijdim,1:kdim_local)
    real(8) :: lv(1:ijdim,1:kdim_local)
    !
    real(8) :: f(1:ijdim,1:kdim_local)
    real(8) :: dfdtem(1:ijdim,1:kdim_local)
    !
    logical :: flag_converge(1:ijdim,1:kdim_local)
    logical :: flag_converge_all(1:kdim_local)
    !
    real(8) :: dtem(1:ijdim)
    real(8), parameter :: PREMIN = 1.d-10
    real(8), parameter :: desw_accuracy=1.d-8
    integer, parameter :: iter_max = 20
    !
    integer :: iter
    integer :: nq
    integer :: ij, k
    !    

#ifdef _FPCOLL_
call START_COLLECTION("moist_dewtem")
#endif

    prev(:,:) = pre(:,:) * qv(:,:) &
         / ( CNST_EPSV * qd(:,:) + qv(:,:) )
    !
    wtem(:,:) = tem(:,:)*0.98d0
    flag_converge(:,:)=.false.
    do k=1, kdim_local
       flag_converge_all(k)=.false.
       do iter=1, iter_max
          if( .not. flag_converge_all(k) )then
             do ij=1, ijdim
                !
                lv(ij,k)=CNST_LH00  + (CNST_CPV - CNST_CL)*wtem(ij,k)
                esw(ij,k) = CNST_PSAT0 &
                     * ( wtem(ij,k) / CNST_TEM00 ) ** ( ( CNST_CPV - CNST_CL ) / CNST_RVAP ) &
                     * exp ( CNST_LH00 / CNST_RVAP &
                     * ( 1.0d0 / CNST_TEM00 - 1.0d0 / wtem(ij,k) ) )
                !
                f(ij,k)      = esw(ij,k) - prev(ij,k)
                dfdtem(ij,k) = lv(ij,k)*esw(ij,k)/(CNST_RVAP*wtem(ij,k)*wtem(ij,k))
                !
                if( .not. flag_converge(ij,k) )then
                   wtem(ij,k) = wtem(ij,k) - f(ij,k)/dfdtem(ij,k)
                end if
                !
                if( abs(f(ij,k)) < prev(ij,k)*desw_accuracy )then
                   flag_converge(ij,k)=.true.
                end if
                !
             end do
             !
             flag_converge_all(k)=all(flag_converge(:,k))
             !
          end if
          !
       end do
    end do
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_dewtem")
#endif

    return
  end subroutine moist_dewtem
  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water
  subroutine moist_dqsw_dtem_rho( &
       tem, rho, &
       dqsdtem  )
    use mod_atmos_cnst, only:   &
         CNST_RVAP,             &
         CNST_CPV,              &
         CNST_CL,               &
         CNST_LH0,              &
         CNST_LH00,             &
         CNST_TEM00
    implicit none
    !
    real(8), intent(in)  :: tem(:,:)
    real(8), intent(in)  :: rho(:,:)
    real(8), intent(out) :: dqsdtem(:,:)
    !
    real(8) :: psat(size(tem,1),size(tem,2)) ! saturation vapor pressure
    real(8) :: lhv(size(tem,1),size(tem,2))  ! latent heat for condensation
    !

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_dqsw_dtem_rho")
#endif

    call moist_psat_water(tem, psat)
    lhv(:,:)      = CNST_LH0 + (CNST_CPV - CNST_CL )*(tem(:,:)-CNST_TEM00)
    dqsdtem(:,:) = psat(:,:)/(rho(:,:)*CNST_RVAP*tem(:,:)*tem(:,:))&
         * (lhv(:,:)/(CNST_RVAP*tem(:,:))-1.d0)
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_dqsw_dtem_rho")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_dqsw_dtem_rho
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice
  subroutine moist_dqsi_dtem_rho( &
       tem, rho, &
       dqsdtem  )
    use mod_atmos_cnst, only:         &
         CNST_RVAP,             &
         CNST_CPV,              &
         CNST_CI,               &
         CNST_LHS0,             &
         CNST_LHS00,            &
         CNST_TEM00
    implicit none
    !
    real(8), intent(in)  :: tem(:,:)
    real(8), intent(in)  :: rho(:,:)
    real(8), intent(out) :: dqsdtem(:,:)
    !
    real(8) :: psat(size(tem,1),size(tem,2)) ! saturation vapor pressure
    real(8) :: lhv(size(tem,1),size(tem,2))  ! latent heat for condensation
    !

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_dqsi_dtem_rho")
#endif

    call moist_psat_ice(tem, psat)
    lhv(:,:)      = CNST_LHS0 + (CNST_CPV - CNST_CI )*(tem(:,:)-CNST_TEM00)
    dqsdtem(:,:) = psat(:,:)/(rho(:,:)*CNST_RVAP*tem(:,:)*tem(:,:))&
         * (lhv(:,:)/(CNST_RVAP*tem(:,:))-1.d0)
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_dqsi_dtem_rho")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_dqsi_dtem_rho
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  subroutine moist_dqsw_dtem_dpre( &
       tem, pre, &
       dqsdtem, dqsdpre )
    !
    use mod_atmos_cnst, only:   &
         CNST_EPSV,             &
         CNST_RVAP,             &
         CNST_CPV,              &
         CNST_CL,               &
         CNST_LH0,              &
         CNST_LH00,             &
         CNST_TEM00
    implicit none
    !
    real(8), intent(in)  :: tem(:,:)
    real(8), intent(in)  :: pre(:,:)
    real(8), intent(out) :: dqsdtem(:,:)
    real(8), intent(out) :: dqsdpre(:,:)
    !
    real(8) :: psat(size(tem,1),size(tem,2)) ! saturation vapor pressure
    real(8) :: lhv(size(tem,1),size(tem,2))  ! latent heat for condensation
    !
    real(8) :: den1(size(tem,1),size(tem,2)) ! denominator
    real(8) :: den2(size(tem,1),size(tem,2)) ! denominator
    !

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_dqsw_dtem_dpre")
#endif

    call moist_psat_water(tem, psat)
    !
    lhv(:,:)      = CNST_LH0 + (CNST_CPV - CNST_CL )*(tem(:,:)-CNST_TEM00)
    den1(:,:)     = (pre(:,:)-(1.d0-CNST_EPSV)*psat(:,:))*(pre(:,:)-(1.d0-CNST_EPSV)*psat(:,:)) 
    den2(:,:)     = den1(:,:)*CNST_RVAP*tem(:,:)*tem(:,:)
    dqsdpre(:,:)  = - CNST_EPSV*                  psat(:,:)/den1(:,:)
    dqsdtem(:,:)  =   CNST_EPSV*pre(:,:)*lhv(:,:)*psat(:,:)/den2(:,:)
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_dqsw_dtem_dpre")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_dqsw_dtem_dpre
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  subroutine moist_dqsi_dtem_dpre( &
       tem, pre, &
       dqsdtem, dqsdpre )
    !
    use mod_atmos_cnst, only:   &
         CNST_EPSV,             &
         CNST_RVAP,             &
         CNST_CPV,              &
         CNST_CI,               &
         CNST_LHS0,             &
         CNST_LHS00,            &
         CNST_TEM00
    implicit none
    !
    real(8), intent(in)  :: tem(:,:)
    real(8), intent(in)  :: pre(:,:)
    real(8), intent(out) :: dqsdtem(:,:)
    real(8), intent(out) :: dqsdpre(:,:)
    !
    real(8) :: psat(size(tem,1),size(tem,2)) ! saturation vapor pressure
    real(8) :: lhv(size(tem,1),size(tem,2))  ! latent heat for condensation
    !
    real(8) :: den1(size(tem,1),size(tem,2)) ! denominator
    real(8) :: den2(size(tem,1),size(tem,2)) ! denominator
    !

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_dqsi_dtem_dpre")
#endif

    call moist_psat_ice(tem, psat)
    !
    lhv(:,:)      = CNST_LHS0 + (CNST_CPV - CNST_CI )*(tem(:,:)-CNST_TEM00)
    den1(:,:)     = (pre(:,:)-(1.d0-CNST_EPSV)*psat(:,:))*(pre(:,:)-(1.d0-CNST_EPSV)*psat(:,:)) 
    den2(:,:)     = den1(:,:)*CNST_RVAP*tem(:,:)*tem(:,:)
    dqsdpre(:,:)  = - CNST_EPSV*                  psat(:,:)/den1(:,:)
    dqsdtem(:,:)  =   CNST_EPSV*pre(:,:)*lhv(:,:)*psat(:,:)/den2(:,:)
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_dqsi_dtem_dpre")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_dqsi_dtem_dpre

end module mod_satadjust
