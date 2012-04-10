!-------------------------------------------------------------------------------
!
!+  Thermodynamics variables module
!
!-------------------------------------------------------------------------------
module mod_thrmdyn
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines in which the thermodyanics variables
  !       are calculated.
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.00      2011/10/24 T.Seiki:  Import from NICAM
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  use mod_atmos_cnst,  only :&
       CNST_RAIR,      &
       CNST_RVAP,      &
       CNST_CV,        &
       CNST_CP,        &
       CNST_CPV,       &
       CNST_CL,        &
       CNST_CI,        &
       CNST_PRE00,     &
       CNST_KAPPA,     &
       CNST_EPSV,      &
       CNST_TEM00,     &
       CNST_PSAT0,     &
       CNST_LHF0,      &
       CNST_LH0,       & 
       CNST_LHF00,     &
       CNST_LH00,      &
       NQW_STR,NQW_END,   &
       I_QV,              &
       I_QC, I_QR,        &
       I_QI, I_QS, I_QG,  &
       CVW, CPW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: thrmdyn_rho
  public :: thrmdyn_ein
  public :: thrmdyn_eth
  public :: thrmdyn_tem
  public :: thrmdyn_pre
  public :: thrmdyn_th
  public :: thrmdyn_tempreth
  public :: thrmdyn_tempre
  public :: thrmdyn_tempre2
  public :: thrmdyn_cv
  public :: thrmdyn_cp
  public :: thrmdyn_qd
  public :: thrmdyn_ent
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_rho( &
       rho,               & !--- OUT : density     
       pre,               & !--- IN  : pressure    
       tem,               & !--- IN  : temperature 
       qd,                & !--- IN  : dry concentration 
       qv )                 !--- IN  : water concentration 
    !------ 
    !------ Density calculation in all regions
    !------    1. calculation region of rho
    !------                   : (:,:,:)
    !------ 
    implicit none
    !
    real(8), intent(out) :: rho(:,:)
    real(8), intent(in)  :: pre(:,:)
    real(8), intent(in)  :: tem(:,:)
    real(8), intent(in)  :: qd(:,:)
    real(8), intent(in)  :: qv(:,:)
    !

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_rho")
#endif

    rho(:,:) = pre(:,:) / tem(:,:) &
         / ( qd(:,:)*CNST_RAIR+qv(:,:)*CNST_RVAP )

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_rho")
#endif

    return
  end subroutine thrmdyn_rho
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_ein( &
       ein,               &  !--- OUT : internal energy
       tem,               &  !--- IN  : temperature
       qd,                &  !--- IN  : dry concentration 
       q )                   !--- IN  : water concentration 
    !------ 
    !------ Internal energy calculation in all regions
    !------    1. calculation region of ein
    !------                   : (:,:)
    !------ 
    !------ CAUTION : ein = CV*T*qd + CVV*T*qv + CPL*T*qc
    !------           ein_moist = CV*T*qd + (CVV*T+LH00)*qv + CPL*T*qc
    !              
    implicit none
    !
    real(8), intent(out) :: ein(:,:)
    real(8), intent(in)  :: tem(:,:)
    real(8), intent(in)  :: qd(:,:)
    real(8), intent(in)  :: q(:,:,:)
    !
    real(8) :: cv(size(tem,1),size(tem,2))
    !

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_ein")
#endif

    call thrmdyn_cv( &
         cv,         & !--- out
         q,          & !--- in
         qd )          !--- in
    ein(:,:) = cv(:,:) * tem(:,:)

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_ein")
#endif

    return
  end subroutine thrmdyn_ein
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_eth( &
       eth,               &  !--- OUT : enthalpy
       ein,               &  !--- IN  : internal energy
       pre,               &  !--- IN  : pressure
       rho )                 !--- IN  : density
    implicit none
    !
    real(8), intent(out) :: eth(:,:)
    real(8), intent(in)  :: ein(:,:)
    real(8), intent(in)  :: pre(:,:)
    real(8), intent(in)  :: rho(:,:)
    !

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_eth")
#endif

    eth(:,:) = ein(:,:) + pre(:,:) / rho(:,:)

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_eth")
#endif

    return
  end subroutine thrmdyn_eth
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tem( &
       tem,               &  !--- OUT  : temperature       
       ein,               &  !--- IN : internal energy
       qd,                &  !--- IN  : dry concentration 
       q )                   !--- IN  : water concentration 
    !
    !------ CAUTION : ein = CV*T*qd + CVV*T*qv + CPL*T*qc
    !------           ein_moist = CV*T*qd + (CVV*T+LH00)*qv + CPL*T*qc
    implicit none
    !
    real(8), intent(out) :: tem(:,:)
    real(8), intent(in)  :: ein(:,:)
    real(8), intent(in)  :: qd(:,:)
    real(8), intent(in)  :: q(:,:,:)
    !
    real(8)  :: cv(size(tem,1),size(tem,2))
    !

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_tem")
#endif

    call thrmdyn_cv( &
         cv,         & !--- out
         q,          & !--- in
         qd )          !--- in
    tem(:,:) = ein(:,:) / cv(:,:)

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_tem")
#endif

    return
  end subroutine thrmdyn_tem
  ! xxxxx [Add] T.Seiki
  subroutine thrmdyn_tempre2( &
       tem,               &  !--- OUT  : temperature       
       pre,               &  !--- OUT  : pressure
       rho,               &  !--- IN  : 
       th,                &  !--- IN  : 
       qd,                &  !--- IN  : dry concentration 
       qv )                  !--- IN  : vapor concentration 
    !
    implicit none
    !
    real(8), intent(out) :: tem(:,:)
    real(8), intent(out) :: pre(:,:)
    real(8), intent(in)  :: rho(:,:)
    real(8), intent(in)  :: th(:,:)
    real(8), intent(in)  :: qd(:,:)
    real(8), intent(in)  :: qv(:,:)
    !
    real(8) :: rpres0, wkappa

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_tempre2")
#endif

    !
    ! pre = rho*tem*(qd*rair+qv*rvap)
    ! tem = th*(pre/pre0)**kappa
    !
    rpres0   = 1.d0/CNST_PRE00
    wkappa   = 1.d0/(1.d0-CNST_KAPPA)
    !
    tem(:,:) = ( th(:,:) &
         * (rho(:,:)*(qd(:,:)*CNST_RAIR+qv(:,:)*CNST_RVAP)*rpres0 )**CNST_KAPPA )**wkappa
    pre(:,:) = rho(:,:)*tem(:,:)*(qd(:,:)*CNST_RAIR+qv(:,:)*CNST_RVAP)
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_tempre2")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_tempre2
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_pre( &
       pre,               &  !--- OUT : pressure
       tem,               &  !--- IN  : temperature       
       rho,               &  !--- IN  : density
       qd,                &  !--- IN  : dry concentration 
       q )                   !--- IN  : water concentration 
    implicit none
    !
    real(8), intent(out) :: pre(:,:)
    real(8), intent(in)  :: tem(:,:)
    real(8), intent(in)  :: rho(:,:)
    real(8), intent(in)  :: qd(:,:)
    real(8), intent(in)  :: q(:,:,:)
    !

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_pre")
#endif

    pre(:,:) = rho(:,:) * tem(:,:) &
         * ( qd(:,:)*CNST_RAIR+q(:,:,I_QV)*CNST_RVAP )

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_pre")
#endif

    return
  end subroutine thrmdyn_pre
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_th( &
       th,               &  !--- OUT : potential temperature
       tem,              &  !--- IN  : temperature       
       pre )                !--- IN  : pressure
    !
    implicit none
    !
    real(8), intent(out) :: th(:,:)
    real(8), intent(in)  :: tem(:,:)
    real(8), intent(in)  :: pre(:,:)
    !
    real(8) :: p0k
    !

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_th")
#endif

    p0k=CNST_PRE00**CNST_KAPPA
    !
    th(:,:) =tem(:,:)*(abs(pre(:,:))**(-CNST_KAPPA))*p0k

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_th")
#endif

    return
  end subroutine thrmdyn_th
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tempreth( &
       tem,                    &  !--- OUT : temperature              
       pre,                    &  !--- OUT : pressure
       th,                     &  !--- OUT : potential temperature
       ein,                    &  !--- IN  : internal energy
       rho,                    &  !--- IN  : density
       qd,                     &  !--- IN  : dry concentration 
       q )                        !--- IN  : water concentration 
    implicit none
    !
    real(8), intent(out) :: tem(:,:)
    real(8), intent(out) :: pre(:,:)
    real(8), intent(out) :: th(:,:)
    real(8), intent(in)  :: ein(:,:)
    real(8), intent(in)  :: rho(:,:)
    real(8), intent(in)  :: qd(:,:)
    real(8), intent(in)  :: q(:,:,:)
    !
    real(8) :: cv(size(tem,1),size(tem,2))
    !
    real(8) :: p0k
    !

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_tempreth")
#endif

    p0k=CNST_PRE00**CNST_KAPPA
    !
    call thrmdyn_cv( &
         cv,         & !--- out
         q,          & !--- in
         qd )          !--- in
    !
    tem(:,:) = ein(:,:) / cv(:,:)
    pre(:,:) = rho(:,:) * tem(:,:) &
         * ( qd(:,:)*CNST_RAIR+q(:,:,I_QV)*CNST_RVAP )
    th(:,:) =tem(:,:)*(abs(pre(:,:))**(-CNST_KAPPA))*p0k

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_tempreth")
#endif

    return
  end subroutine thrmdyn_tempreth
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tempre( &
       tem,                  &  !--- OUT : temperature              
       pre,                  &  !--- OUT : pressure
       ein,                  &  !--- IN  : internal energy
       rho,                  &  !--- IN  : density
       qd,                   &  !--- IN  : dry concentration 
       q )                      !--- IN  : water concentration 
    !
    implicit none
    !
    real(8), intent(out) :: tem(:,:)
    real(8), intent(out) :: pre(:,:)
    real(8), intent(in)  :: ein(:,:)
    real(8), intent(in)  :: rho(:,:)
    real(8), intent(in)  :: qd(:,:)
    real(8), intent(in)  :: q(:,:,:)
    !
    real(8) :: cv(size(tem,1),size(tem,2))
    !

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_tempre")
#endif

    call thrmdyn_cv( &
         cv,         & !--- out
         q,          & !--- in
         qd )          !--- in
    !
    tem(:,:) = ein(:,:) / cv(:,:)
    pre(:,:) = rho(:,:) * tem(:,:) &
         * ( qd(:,:)*CNST_RAIR+q(:,:,I_QV)*CNST_RVAP )

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_tempre")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_tempre
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_cv( &
       cva,              & !--- OUT : specific heat
       q,                & !--- IN  : mass concentration
       qd)                 !--- IN  : dry mass concentration

    implicit none
    !
    real(8), intent(out) :: cva(:,:)
    real(8), intent(in) :: q(:,:,:)
    real(8), intent(in) :: qd(:,:)
    !
    integer :: nq
    !

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_cv")
#endif

    cva(:,:) = qd(:,:) * CNST_CV
    do nq = NQW_STR,NQW_END
       cva(:,:) = cva(:,:) + q(:,:,nq) * CVW(nq)
    end do
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_cv")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_cv
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_cp( &
       cpa,              & !--- OUT : specific heat
       q,                & !--- IN  : mass concentration
       qd )                !--- IN  : dry mass concentration
    implicit none
    !
    real(8), intent(out) :: cpa(:,:)
    real(8), intent(in) :: q(:,:,:)
    real(8), intent(in) :: qd(:,:)
    !
    integer :: nq
    !

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_cp")
#endif

    cpa(:,:) = qd(:,:) * CNST_CP
    do nq = NQW_STR,NQW_END
       cpa(:,:) = cpa(:,:) + q(:,:,nq) * CPW(nq)
    end do
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_cp")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_cp
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_qd( &
       qd,               & !--- OUT  : dry mass concentration
       q )                 !--- IN  : mass concentration
    implicit none
    !
    real(8), intent(out) :: qd(:,:)
    real(8), intent(in) :: q(:,:,:)
    !
    integer :: nq
    !

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_qd")
#endif

    qd(:,:) = 1.0D0
    do nq = NQW_STR,NQW_END
       qd(:,:) = qd(:,:) - q(:,:,nq)
    end do

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_qd")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_qd
  !-----------------------------------------------------------------------------
  subroutine thrmdyn_ent( &
       ent,               & !--- OUT :
       tem,               & !--- IN
       pre,               & !--- IN
       q,                 & !--- IN
       qd                 & !--- IN
       )
    !
    implicit none
    !
    real(8), intent(out):: ent(:,:)
    real(8), intent(in) :: tem(:,:)
    real(8), intent(in) :: pre(:,:)
    real(8), intent(in) :: q(:,:,:)
    real(8), intent(in) :: qd(:,:)
    !
    real(8) :: pred(size(tem,1),size(tem,2))
    real(8) :: prev(size(tem,1),size(tem,2))
    real(8) :: cpa(size(tem,1),size(tem,2))
    !
    real(8), parameter :: PREMIN = 1.d-10
    !
    integer :: nq

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_ent")
#endif

    !
    ! entropy
    ! ent = qd * sd + qv * sv + qc * sc + qr * sr ...
    pred(:,:) = pre(:,:) * CNST_EPSV * qd(:,:) &
         / ( CNST_EPSV * qd(:,:) + q(:,:,I_QV) )
    prev(:,:) = pre(:,:) * q(:,:,I_QV) &
         / ( CNST_EPSV * qd(:,:) + q(:,:,I_QV) )
    pred(:,:) = max( pred(:,:), PREMIN )
    prev(:,:) = max( prev(:,:), PREMIN )
    !
    !  T.Mitsui, strict definition of entropy is given by 
    ! (8),(11),(13) Satoh (2003), Mon.Wea.Rev.
    !
    ent(:,:) = qd(:,:) * CNST_CP   * log ( tem(:,:)  / CNST_TEM00 )   &
         -     qd(:,:) * CNST_RAIR * log ( pred(:,:) / CNST_PRE00 )   &
         + q(:,:,I_QV) * CNST_CPV  * log ( tem(:,:)  / CNST_TEM00 )   &
         - q(:,:,I_QV) * CNST_RVAP * log ( prev(:,:) / CNST_PSAT0 )   &
         + q(:,:,I_QV) * CNST_LH0 / CNST_TEM00 
    do nq = NQW_STR, NQW_END
       if      ( (nq==I_QC) .or. (nq==I_QR) )then
          ent(:,:) = ent(:,:) + q(:,:,nq) * CNST_CL * log( tem(:,:) / CNST_TEM00 )
       else if ( (nq==I_QI) .or. (nq==I_QS) .or. (nq==I_QG)  ) then
          ent(:,:) = ent(:,:) + q(:,:,nq) * CNST_CI * log( tem(:,:) / CNST_TEM00 ) &
               - q(:,:,nq) * CNST_LHF0 / CNST_TEM00
       end if
    end do
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_ent")
#endif

    return
  end subroutine thrmdyn_ent
  !
end module mod_thrmdyn
!-------------------------------------------------------------------------------
