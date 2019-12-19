      SUBROUTINE PHCH00T(POSGPT,SHAPDT,DVOLDT,PROPST,EHISTT,ELDIST,
     .                   SHAPET,ELELMT,VELCMT,ELEL1T,ELEL2T,INPCCT,
     .                    frrrc, frrr1, fdvol,
     .                   WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,
     .                   LNODST,DISIMT,ELELTT,HACHET,WHADEL,CENTRO,
     .                   ADVEMT,TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE PHASE-CHANGE RESIDUAL "THERMAL FORCES"
C     IN A UNIPHASE UNIDIMENSIONAL SUBELEMENT FOR THE CASE:
C
C     NPHCH=1
C
C     (ISOTHERMAL PHASE-CHANGE)
C    
C     (UNIDIMENSIONAL TWO-NODED ELEMENTS)
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
     .          DVOLDT(*)
      DIMENSION PROPST(*),        VELCMT(*),
     .          ELELMT(*),        EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*), DISIMT(*),
     .          ELELTT(*),        SHAPET(NNODLT,*),
     .          ELEL1T(*),        ELEL2T(*)
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .          XJACMT(NDIMET,*), CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*)
      DIMENSION ELCODT(NDIMET,*), LNODST(*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          DVOLIT(*)
C
      CALL SHAP00T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,SHAPDT,
     .             DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,EHISTT,ELDIST,
     .              frrrc, frrr1, fdvol,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT)
C
      CALL FRUF05T(PROPST,VELCMT,ELELMT,DVOLDT,SHAPDT,EHISTT,ELDIST,
     .             DISIMT,ELELTT,ELEL1T,ELEL2T,WHADEL,INPCCT,TEINIT,
     .             FPCHLT,DVOLIT)
C
      RETURN
      END
