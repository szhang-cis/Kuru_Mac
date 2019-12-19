      SUBROUTINE LAHE01T(POSGPT,SHAPDT,DVOLDT,PROPST,EHISTT,ELDIST,
     .                   SHAPET,WSTIFT,VELCMT,WSTI1T,WSTI2T,
     .                    frrrc, frrr1, fdvol,
     .                   WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,
     .                   LNODST,DISIMT,HACHET,WHADEL,CENTRO,ADVEMT,
     .                   TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE PHASE CHANGE MATRIX (LATENT HEAT 
C     CONTRIBUTION INTO THE JACOBIAN MATRIX FOR ELEMENT NO. 5 ) IN A
C     UNIPHASE UNIDIMENSIONAL SUBELEMENT.
C
C     (TWO-NODED ELEMENTS)
C    
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION POSGPT(NDIMET,*), SHAPDT(NNODLT,*),
     .          DVOLDT(*),        WSTIFT(*)
      DIMENSION PROPST(*),
     .          EHISTT(NHISTT,*), ELDIST(NDOFCT,*),
     .          SHAPET(NNODLT,*), VELCMT(*)
      DIMENSION ELCODT(NDIMET,*)
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .          XJACMT(NDIMET,*), CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*)
      DIMENSION LNODST(*),        DISIMT(*)
      DIMENSION WSTI1T(*),        WSTI2T(*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          DVOLIT(*)
C
      CALL SHAP01T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,SHAPDT,
     .             DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,EHISTT,ELDIST,
     .              frrrc, frrr1, fdvol,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT)
C
      CALL MAUF05T(DVOLDT,PROPST,SHAPDT,WSTIFT,EHISTT,VELCMT,WSTI1T,
     .             WSTI2T,ELDIST,DISIMT,WHADEL,TEINIT,FPCHLT,DVOLIT)
C
      RETURN
      END
