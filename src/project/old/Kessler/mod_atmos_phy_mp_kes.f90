!-------------------------------------------------------------------------------
!
!+  NICAM Double moment Water 6 scheme
!
!-------------------------------------------------------------------------------
module mod_atmos_phy_mp
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines for the Kessler parametrization.
  !
  !       
  !++ Current Corresponding Author : T.Seiki
  ! 
  !++ History: NDW6
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !        0.00   12/01/14 Y.Miyamoto
  !               
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_setup
  public :: ATMOS_PHY_MP
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=32), save :: WLABEL(11)
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Wrapper for KESSLER'

    WLABEL( 1) = "VAPOR"
    WLABEL( 2) = "CLOUD"
    WLABEL( 3) = "RAIN"
    WLABEL( 4) = "ICE"
    WLABEL( 5) = "SNOW"
    WLABEL( 6) = "GRAUPEL"
    WLABEL( 7) = "CLOUD_NUM"
    WLABEL( 8) = "RAIN_NUM"
    WLABEL( 9) = "ICE_NUM"
    WLABEL(10) = "SNOW_NUM"
    WLABEL(11) = "GRAUPEL_NUM"

    return
  end subroutine ATMOS_PHY_MP_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP
    use mod_const, only: &
       Rdry   => CONST_Rdry,   &
       CPdry  => CONST_CPdry,  &
       CVdry  => CONST_CVdry,  &
       RovCP  => CONST_RovCP,  &
       CPovCV => CONST_CPovCV, &
       Pstd   => CONST_PRE00,  &
       Rvap   => CONST_Rvap,   &
       CPvap  => CONST_CPvap,  &
       CVvap  => CONST_CVvap,  &
       Lv0    => CONST_LH0,    &
       cl     => CONST_CL,     &
       t0     => CONST_TEM00
    use mod_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP, &
       TIME_NOWSTEP
    use mod_grid, only : &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       GRID_CDZ
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait, &
       COMM_total
    use mod_atmos_vars, only: &
       var => atmos_var, &
       A_NAME,      &
       VA  => A_VA, &
       QA  => A_QA, &
       I_DENS,      &
       I_MOMX,      &
       I_MOMY,      &
       I_MOMZ,      &
       I_RHOT,      &
       I_QV,        &
       I_QC,        &
       I_QR
    implicit none

    ! local
    real(8) :: pott  (KA,IA,JA)      ! potential temperature [K]
    real(8) :: dens_l(KA,IA,JA)      ! density [kg/m3]
    real(8) :: momz_l(KA,IA,JA)      ! momentum (z) [kg/m3 * m/s]
    real(8) :: pott_l(KA,IA,JA)      ! potential temperature [K]
    real(8) :: qtrc_l(KA,IA,JA,3)    ! tracer mixing ratio   [kg/kg],[1/m3]

    integer :: k, i, j, iv, cc, re, dz(KA)
    real(8), parameter :: tt1=1.7269d1, tt2=3.585d1
    real(8) :: dens_s(KA)
    real(8) :: dens_sf, ptot
    real(8) :: dt
    real(8) :: qvs
    real(8) :: tem
    real(8) :: Lv, cpm, cvm, rm
    real(8) :: tvel1, tvel2, pp, dtqv, es, et
    real(8) :: th2, th_old, qc_old, ql_old
    real(8) :: ar, ac, er, fr, c3, aa

    !---------------------------------------------------------------------------

    dz(:) = GRID_CDZ(:)
    dt    = TIME_DTSEC_ATMOS_PHY_MP

    pott  (:,:,:)   = var(:,:,:,I_RHOT) / var(:,:,:,I_DENS) 
    dens_l(:,:,:)   = var(:,:,:,I_DENS)
    momz_l(:,:,:)   = var(:,:,:,I_MOMZ)
    pott_l(:,:,:)   = pott(:,:,:)
    qtrc_l(:,:,:,1) = var(:,:,:,5+I_QV)
    qtrc_l(:,:,:,2) = var(:,:,:,5+I_QC)
    qtrc_l(:,:,:,3) = var(:,:,:,5+I_QR)

    do k = 1, KA
       dens_s(k) = sum(dens_l(k,:,:))/dble(IA*JA)
    end do
    dens_sf = 1.16711355984244d0

    do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             cc  = 0
             aa  = 0.D0
             ar  = 0.D0
             ac  = 0.D0
             er  = 0.D0
             Lv  = Lv0
             cpm = CPdry + CPvap * qtrc_l(k,i,j,1)
             cvm = CVdry + CVvap * qtrc_l(k,i,j,1)
             rm  = Rdry  +  Rvap * qtrc_l(k,i,j,1)
             th2 = 0.D0
             dtqv= 0.D0
             pp  = Pstd * ( dens_l(k,i,j) * pott_l(k,i,j) * Rdry / Pstd )**CPovCV / 1.D2
             tem = pp * 1.D2 / ( dens_l(k,i,j) * Rdry )
             qvs = 6.22D-1 * 6.1078D0/pp * dexp ( tt1 * ( tem - t0 ) / ( tem - tt2 ) )
             ptot= tem / pott(k,i,j)

             if ( qtrc_l(k,i,j,2) < 0.D0 ) write(IO_FID_LOG,*) "minus qc0",k,i,j,qtrc_l(k,i,j,2)
             if ( qtrc_l(k,i,j,3) < 0.D0 ) write(IO_FID_LOG,*) "minus qr0",k,i,j,qtrc_l(k,i,j,3)

!             if ( qtrc_l(k,i,j,1) > qvs .or. qtrc_l(k,i,j,2) > 0.D0 ) then
             if ( qtrc_l(k,i,j,1) > qvs ) then
 91             continue
                cpm = CPdry + CPvap * qtrc_l(k,i,j,1)
                cvm = CVdry + CVvap * qtrc_l(k,i,j,1)
                rm  = Rdry  +  Rvap * qtrc_l(k,i,j,1)
                pp  = Pstd * ( dens_l(k,i,j) * pott_l(k,i,j) * Rdry / Pstd )**CPovCV / 1.D2
                tem = pp * 1.D2 / ( dens_l(k,i,j) * Rdry )
                th2 = CVdry * Lv / ( cvm * CPdry * ptot ) &
                    - pott_l(k,i,j) * Rvap / cvm * ( 1.D0 - ( Rdry * cpm ) / ( CPdry * rm ) )
! --- Tetens' formula ---
                qvs = 6.22d-1 * 6.1078D0/pp * dexp ( tt1 * ( tem - t0 ) / ( tem - tt2 ) )
! --- if thed air is supersaturated ---
                dtqv = ( qtrc_l(k,i,j,1) - qvs ) / &
                       ( 1.D0 + tt1 * ( t0 - tt2 ) * qvs * Lv / CPdry / ( tem - tt2 ) )
!                dtqv = dmin1( dtqv , qtrc_l(k,i,j,2) )
!                dtqv = dmin1( dtqv , -1.D0*qtrc_l(k,i,j,2) , 0.d0 )
!                if ( qtrc_l(k,i,j,1) + dtqv < 1.0D-20 ) dtqv = 1.D0-20 - qtrc_l(k,i,j,1)
!                goto 92
!                if (cc > 1) write(IO_FID_LOG,*) "demand 1",qtrc(k,i,j,1),qtrc(k,i,j,2),qtrc(k,i,j,3),dtqv
                th_old          = pott_l(k,i,j  )
                qc_old          = qtrc_l(k,i,j,2)
                pott_l(k,i,j  ) = pott_l(k,i,j  ) + th2 * dtqv
                qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) - dtqv
                qtrc_l(k,i,j,2) = qtrc_l(k,i,j,2) + dtqv
                cc = cc + 1
!                if (cc > 2) write(IO_FID_LOG,*) "demand 2",k,i,j,cc,pott_l(k,i,j  ),qtrc_l(k,i,j,1),qtrc_l(k,i,j,2),pp,tem,th2,dtqv
                if ( qtrc_l(k,i,j,2) < 0.D0 ) then
                   write(IO_FID_LOG,*) "minus qc1",k,i,j,qtrc_l(k,i,j,2),qc_old,dtqv,t0,qvs
                   pott_l(k,i,j  ) = pott_l(k,i,j  ) + th2 * qtrc_l(k,i,j,2)
                   qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) - qtrc_l(k,i,j,2)
                   qtrc_l(k,i,j,2) = 0.D0
                   goto 92
                endif
                if ( dabs( pott_l(k,i,j) - th_old ) <= 1.D-5 ) then
!                    write(IO_FID_LOG,*) "converged th1",k,i,j,dabs( pott_l(k,i,j) - th_old ),th_old,pott_l(k,i,j),dtqv,th2,qtrc_l(k,i,j,2),qtrc_l(k,i,j,3),pp
                   goto 92
                else
                   if ( cc > 300 ) then
                      write(IO_FID_LOG,*) "not converged th1", k,i,j,dabs( pott_l(k,i,j) - th_old ),th_old,pott_l(k,i,j), pott(k,i,j)
                      goto 92
                   end if
                   goto 91
                end if
             end if
 92          continue
             qc_old          = qtrc_l(k,i,j,2)
             ql_old          = qtrc_l(k,i,j,3)
             if ( qtrc_l(k,i,j,2) < 0.D0 ) write(IO_FID_LOG,*) "minus qc0s",k,i,j,qtrc_l(k,i,j,2),tem,Lv,cpm,cvm,rm,tem,qvs
             if ( qtrc_l(k,i,j,3) < 0.D0 ) write(IO_FID_LOG,*) "minus qr0s",k,i,j,qtrc_l(k,i,j,3),dtqv
!             if (cc > 1) write(IO_FID_LOG,*) "returdation 0",k,i,j,pott_l(k,i,j  ),qtrc_l(k,i,j,1),qtrc_l(k,i,j,2),pp,tem,th2,dtqv
! --- Terminal velocity ---
! --- Durran and Klemp (1983) ---
             pp  = Pstd * ( dens_l(k,i,j) * pott_l(k,i,j) * Rdry / Pstd )**CPovCV / 1.D2
             tem = pp * 1.D2 / ( dens_l(k,i,j) * Rdry )
!             if ( k<15 .and. qtrc_l(k,i,j,2) > 0 ) write(IO_FID_LOG,*) "testa 01",k,i,j,pp,tem,qvs,qtrc_l(k,i,j,1),qtrc_l(k,i,j,2),qtrc_l(k,i,j,3)
!             if ( k<15 .and. qtrc_l(k,i,j,3) > 0 ) write(IO_FID_LOG,*) "testb 02",pp,tem,qvs,qtrc_l(k,i,j,1),qtrc_l(k,i,j,2),qtrc_l(k,i,j,3)
             if ( qtrc_l(k,i,j,3) > 0.D0 ) then
! --- ARPS ---
                if ( k >= KS .and. k < KE ) then
                   tvel1 = 3.634D1 * ( 1.0D-3 * dens_s(k+1) * dabs(qtrc_l(k+1,i,j,3)) )**(1.346D-1) &
                         * dsqrt( dens_sf/dens_s(k+1) )
                   tvel2 = 3.634D1 * ( 1.0D-3 * dens_s(k  ) * dabs(qtrc_l(k  ,i,j,3)) )**(1.346D-1) &
                         * dsqrt( dens_sf/dens_s(k  ) )
                   fr = ( dens_l(k+1,i,j)*tvel1*dabs(qtrc_l(k+1,i,j,3)) &
                        - dens_l(k  ,i,j)*tvel2*dabs(qtrc_l(k  ,i,j,3)) ) / dz(k) 
!                write(IO_FID_LOG,*) "rain 00",k,i,j,tvel1,tvel2,fr,dens_s(k+1), qtrc_l(k+1,i,j,3),( 1.0D-3 * dens_s(k+1) * dabs(qtrc_l(k+1,i,j,3)) )**(1.346D-1),dsqrt( dens_s(KS)/dens_s(k+1) )
                else if ( k == KS ) then
                   tvel1 = 3.634D1 * ( 1.0D-3 * dens_s(k+1) * dabs(qtrc_l(k+1,i,j,3)) )**(1.346D-1) &
                         * dsqrt( dens_sf/dens_s(k+1) )
                   tvel2 = 3.634D1 * ( 1.0D-3 * dens_s(k  ) * dabs(qtrc_l(k  ,i,j,3)) )**(1.346D-1) &
                         * dsqrt( dens_sf/dens_s(k  ) )
                   fr = ( dens_l(k+1,i,j)*tvel1*dabs(qtrc_l(k+1,i,j,3)) &
                        - dens_l(k  ,i,j)*tvel2*dabs(qtrc_l(k  ,i,j,3)) ) / dz(k)
!                write(IO_FID_LOG,*) "rain 01",k,i,j,tvel1,tvel2,fr,dens_s(k+1), qtrc_l(k+1,i,j,3),( 1.0D-3 * dens_s(k+1) * dabs(qtrc_l(k+1,i,j,3)) )**(1.346D-1),dsqrt( dens_s(KS)/dens_s(k+1) )
                else if ( k == KE ) then
                   tvel1 = 3.634D1 * ( 1.0D-3 * dens_s(k  ) * dabs(qtrc_l(k  ,i,j,3)) )**(1.346D-1) &
                         * dsqrt( dens_sf/dens_s(k  ) )
                   tvel2 = 3.634D1 * ( 1.0D-3 * dens_s(k  ) * dabs(qtrc_l(k  ,i,j,3)) )**(1.346D-1) &
                         * dsqrt( dens_sf/dens_s(k  ) )
!                   fr = ( dens_l(k  ,i,j)*tvel1*dabs(qtrc_l(k  ,i,j,3)) &
                   fr = ( 0.d0 &
                        - dens_l(k  ,i,j)*tvel2*dabs(qtrc_l(k  ,i,j,3)) ) / dz(k)
!                write(IO_FID_LOG,*) "rain 02",k,i,j,tvel1,tvel2,fr,dens_s(k),qtrc_l(k,i,j,3),( 1.0D-3 * dens_s(k) * dabs(qtrc_l(k,i,j,3)) )**(1.346D-1),dsqrt( dens_s(KS)/dens_s(k) )
                endif
             else
                fr = 0.D0
             endif
             if ( qtrc_l(k,i,j,3) > 0.D0 ) then
                qvs = 6.22D-1 * 6.1078D0/pp * dexp ( tt1 * ( tem - t0 ) / ( tem - tt2 ) )
!                write(IO_FID_LOG,*) "rain 1",k,i,j,pott_l(k,i,j  ),qtrc_l(k,i,j,1),qtrc_l(k,i,j,2),qtrc_l(k,i,j,3)
                if ( qtrc_l(k,i,j,1) < qvs ) then
                   c3 = 1.6D0 + 3.03922D1 *( ( dens_l(k,i,j) * qtrc_l(k,i,j,3) )**(2.046D-1) )
! --- evaporation rate (Ogura and Takahashi, 1971) ---
                   er = ( 1.D0 - qtrc_l(k,i,j,1) / qvs ) * c3 * ( ( dens_l(k,i,j) * qtrc_l(k,i,j,3) )**(5.25D-1) ) &
                      / ( dens_l(k,i,j) * ( 2.03D4 + 9.584D6 / ( pp*1.D2 * qvs ) ) )
!                   er = dmin1( er , qtrc_l(k,i,j,3)/dt )
!                   if ( qtrc_l(k,i,j,1) + er * dt > qvs ) er = ( qvs - qtrc_l(k,i,j,1) ) / dt
!                   write(IO_FID_LOG,*) "rain 2",c3,er
                endif
             endif
! --- Kessler (1969) parameterization --- 
! --- Autoconversion ---
!       if (qc(i,k,4).ge.0.5d-3) then 
!        ar=1.0d-3*(qc(i,k,4)-0.5d-3) 
! --- ARPS ---
             if ( qtrc_l(k,i,j,2) >= 1.0d-3 ) then
                ar = 1.0D-3 * ( qtrc_l(k,i,j,2) - 1.0D-3 )
!                if (k < 7) write(IO_FID_LOG,*) "cloud 1",qtrc_l(k,i,j,2),qtrc_l(k,i,j,3),ar
             endif
! --- Accretion (Collection) ---
             if ( qtrc_l(k,i,j,2) > 0.D0 .and. qtrc_l(k,i,j,3) > 0.D0 ) then
                ac = 2.2D0 * qtrc_l(k,i,j,2) * ( qtrc_l(k,i,j,3)**(8.75D-1) )
!                if (k < 7) write(IO_FID_LOG,*) "cloud 2",qtrc_l(k,i,j,2),qtrc_l(k,i,j,3),ac
             endif
             if ( ar + ac > dabs(qtrc_l(k,i,j,2)) / dt ) then
                aa = ar + ac
!                ar = qtrc_l(k,i,j,2) * ar / aa
!                ac = qtrc_l(k,i,j,2) * ac / aa
                write(IO_FID_LOG,*) "cloud 3",k,i,j,qtrc_l(k,i,j,2),aa,ar,ac,dt
             endif
             th2 = CVdry * Lv / ( cvm * CPdry * ptot) &
                 - pott_l(k,i,j) * Rvap / cvm * ( 1.D0 - ( Rdry * cpm ) / ( CPdry * rm ) )
!             pott_l(k,i,j  ) = pott_l(k,i,j  ) - th2 * ( dtqv + dt * er )
!             qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) + dtqv + dt * er
!             qtrc_l(k,i,j,2) = qtrc_l(k,i,j,2) - dtqv - dt * ( ar + ac )
             qc_old          = qtrc_l(k,i,j,2)
             ql_old          = qtrc_l(k,i,j,3)
             pott_l(k,i,j  ) = pott_l(k,i,j  ) - th2 * ( dt * er )
             qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) + dt * er
             qtrc_l(k,i,j,2) = qtrc_l(k,i,j,2) - dt * ( ar + ac )
             qtrc_l(k,i,j,3) = qtrc_l(k,i,j,3) + dt * ( ar + ac - er ) + dt * ( fr / dens_l(k,i,j) )
             qtrc_l(k,i,j,1) =  max ( qtrc_l(k,i,j,1) , 0.D0 )
             if ( qtrc_l(k,i,j,2) < 0.D0 ) then
                write(IO_FID_LOG,*) "minus qc2",k,i,j,qtrc_l(k,i,j,2), th2, qc_old,ar,ac,er,dt
!                pott_l(k,i,j  ) = pott_l(k,i,j  ) - th2 * qc_old
!                qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) + qc_old
                pott_l(k,i,j  ) = pott_l(k,i,j  ) + th2 * qtrc_l(k,i,j,2)
                qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) - qtrc_l(k,i,j,2)
                qtrc_l(k,i,j,2) = 0.D0
             endif
             if ( qtrc_l(k,i,j,3) < 0.D0 ) then
                write(IO_FID_LOG,*) "minus qr2",k,i,j,qtrc_l(k,i,j,3), th2, ql_old,ar,ac,er,fr
!                pott_l(k,i,j  ) = pott_l(k,i,j  ) - th2 * ql_old
!                qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) + ql_old
                pott_l(k,i,j  ) = pott_l(k,i,j  ) + th2 * qtrc_l(k,i,j,3)
                qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) - qtrc_l(k,i,j,3)
                qtrc_l(k,i,j,3) = 0.D0
             endif
             dtqv = 0.d0
             cc   = 0.d0
!             if ( k<15 .and. qtrc_l(k,i,j,2) > 0 ) write(IO_FID_LOG,*) "cloudb 1",k,i,j,qtrc_l(k,i,j,2),qtrc_l(k,i,j,3),ar,ac,er,dt,fr
!             if ( k<15 .and. qtrc_l(k,i,j,3) > 0 ) write(IO_FID_LOG,*) "cloudi 1",k,i,j,qtrc_l(k,i,j,2),qtrc_l(k,i,j,3),ar,ac,er,dt,fr
!             if (k < 7 .and. qtrc_l(k,i,j,2) /= 0) write(IO_FID_LOG,*) "returdation 1",k,i,j,pott_l(k,i,j  ),qtrc_l(k,i,j,1),qtrc_l(k,i,j,2),qtrc_l(k,i,j,3),dtqv,th2,ar,ac,aa
!             goto 94
             pp  = Pstd * ( dens_l(k,i,j) * pott_l(k,i,j) * Rdry / Pstd )**CPovCV / 1.D2
             tem = pp * 1.D2 / ( dens_l(k,i,j) * Rdry )
             qvs = 6.22D-1 * 6.1078D0/pp * dexp ( tt1 * ( tem - t0 ) / ( tem - tt2 ) )
!             if ( qtrc_l(k,i,j,1) > qvs .or. qtrc_l(k,i,j,2) > 0.D0 ) then
             if ( qtrc_l(k,i,j,1) > qvs ) then
!                aa  = th2 * tem/pott(k,i,j) * 1.767D1 * 2.435D2
 93             continue
                cpm = CPdry + CPvap * qtrc_l(k,i,j,1)
                cvm = CVdry + CVvap * qtrc_l(k,i,j,1)
                rm  = Rdry  +  Rvap * qtrc_l(k,i,j,1)
                pp  = Pstd * ( dens_l(k,i,j) * pott_l(k,i,j) * Rdry / Pstd )**CPovCV / 1.D2
                tem = pp * 1.D2 / ( dens_l(k,i,j) * Rdry )
! --- Tetens' formula ---
                qvs = 6.22D-1 * 6.1078D0/pp * dexp ( tt1 * ( tem - t0 ) / ( tem - tt2 ) )
                th2 = CVdry * Lv / ( cvm * CPdry * tem/pott(k,i,j) ) &
                    - pott_l(k,i,j) * Rvap / cvm * ( 1.D0 - ( Rdry * cpm ) / ( CPdry * rm ) )
!
! --- if thed air is supersaturated ---
                dtqv = ( qtrc_l(k,i,j,1) - qvs ) / &
                       ( 1.D0 + tt1 * ( t0 - tt2 ) * qvs * Lv / CPdry / ( tem - tt2 ) )
!                dtqv = dmin1( dtqv , qtrc_l(k,i,j,2) )
!                dtqv = ( qvs - qtrc_l(k,i,j,1) ) / ( 1.D0 + aa * qvs / ( ( tem - 2.965D1 )**2 ) )
                th_old          = pott_l(k,i,j  )
                qc_old          = qtrc_l(k,i,j,2)
                pott_l(k,i,j  ) = pott_l(k,i,j  ) + th2 * dtqv
                qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) - dtqv
                qtrc_l(k,i,j,2) = qtrc_l(k,i,j,2) + dtqv
                cc = cc + 1
                if ( qtrc_l(k,i,j,2) < 0.D0 ) then
                   write(IO_FID_LOG,*) "minus qc3",k,i,j,qtrc_l(k,i,j,2), th2, qc_old, dtqv
                   pott_l(k,i,j  ) = pott_l(k,i,j  ) + th2 * qtrc_l(k,i,j,2)
                   qtrc_l(k,i,j,1) = qtrc_l(k,i,j,1) - qtrc_l(k,i,j,2)
                   qtrc_l(k,i,j,2) = 0.D0
                   goto 94
                endif
                if ( dabs( pott_l(k,i,j) - th_old ) <= 1.D-5 ) then
!                   write(IO_FID_LOG,*) "converged th2",k,i,j,dabs( pott_l(k,i,j) - th_old ),th_old,pott_l(k,i,j)
                   goto 94
                else
                   if ( cc > 300 ) then
                      write(IO_FID_LOG,*) "not converged th2", k,i,j,dabs( pott_l(k,i,j) - th_old ),th_old,pott_l(k,i,j),pott(k,i,j)
                      goto 94
                   end if
                   goto 93
                endif
             endif
 94          continue
          enddo
       enddo
    enddo

    var(:,:,:,I_RHOT) = pott_l(:,:,:)*dens_l(:,:,:)
    var(:,:,:,I_DENS) = dens_l(:,:,:)
    var(:,:,:,I_MOMZ) = momz_l(:,:,:)
    var(:,:,:,5+I_QV) = qtrc_l(:,:,:,1)
    var(:,:,:,5+I_QC) = qtrc_l(:,:,:,2)
    var(:,:,:,5+I_QR) = qtrc_l(:,:,:,3)


    ! fill IHALO & JHALO
    do iv = 1, VA
       call COMM_vars8( var(:,:,:,iv), iv )
       call COMM_wait ( var(:,:,:,iv), iv )
    enddo

    ! check total mass
    call COMM_total( var(:,:,:,:), A_NAME(:) )

    return
  end subroutine ATMOS_PHY_MP

end module mod_atmos_phy_mp 
