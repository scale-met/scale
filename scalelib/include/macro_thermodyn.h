!-------------------------------------------------------------------------------
!> MACRO for Thermodynamics
!!
!! @par Description
!!          Thermodynamics macros
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-22 (S.Nishizawa)  [new]
!!
!<
!-------------------------------------------------------------------------------
#ifdef DRY
#define CALC_QDRY(qdry, q, mass, k, i, j, iqw)	\
    qdry = 1.0_RP
#else
#define CALC_QDRY(qdry, q, mass, k, i, j, iqw)				\
    qdry = 1.0_RP; do iqw=1,QA; qdry = qdry-q(k,i,j,iqw)*mass(iqw); enddo
#endif

#ifdef DRY
#define CALC_R(r, qdry, q, k, i, j, iqw, Rdry, Rq) \
    r = Rdry
#else
#define CALC_R(r, qdry, q, k, i, j, iqw, Rdry, Rq) \
    r = qdry*Rdry; do iqw=1,QA; r = r+q(k,i,j,iqw)*Rq(iqw); enddo
#endif

#ifdef DRY
#define CALC_CV(cv, qdry, q, k, i, j, iqw, CVdry, CVq) \
    cv = CVdry
#else
#define CALC_CV(cv, qdry, q, k, i, j, iqw, CVdry, CVq) \
    cv = qdry*CVdry; do iqw=1,QA; cv = cv+q(k,i,j,iqw)*CVq(iqw); enddo
#endif

#ifdef DRY
#define CALC_CP(cp, qdry, q, k, i, j, iqw, CPdry, CPq) \
    cp = CPdry
#else
#define CALC_CP(cp, qdry, q, k, i, j, iqw, CPdry, CPq) \
    cp = qdry*CPdry; do iqw=1,QA; cp = cp+q(k,i,j,iqw)*CPq(iqw); enddo
#endif


#define CALC_PRE(pre, dens, pott, r, cp, P00) \
  pre = P00 * ( dens*r*pott/P00 )**(cp/(cp-r))
