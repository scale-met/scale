!-------------------------------------------------------------------------------
!> MACRO for Thermodynamics
!!
!! @par Description
!!          Thermodynamics macros
!!
!! @author Theam Scale
!!
!! @par History
!! @li      2012-12-22 (S.Nishizawa)  [new]
!!
!<
!-------------------------------------------------------------------------------
#define CALC_QDRY(qdry, q, k, i, j, iqw) \
    qdry = 1.0_RP; do iqw=QQS,QQE; qdry = qdry-q(k,i,j,iqw); enddo

#define CALC_CP(cp, qdry, q, k, i, j, iqw, CPdry, AQ_CP) \
    cp = qdry*CPdry; do iqw=QQS,QQE; cp = cp+q(k,i,j,iqw)*AQ_CP(iqw); enddo

#define CALC_CV(cv, qdry, q, k, i, j, iqw, CVdry, AQ_CV) \
    cv = qdry*CVdry; do iqw=QQS,QQE; cv = cv+q(k,i,j,iqw)*AQ_CV(iqw); enddo

#define CALC_R(r, qv, qdry, Rdry, Rvap) \
    r = Rdry * qdry + Rvap * qv

#define CALC_PRE(pre, dens, pott, r, cp, P00) \
  pre = P00 * ( dens*r*pott/P00 )**(cp/(cp-r))
