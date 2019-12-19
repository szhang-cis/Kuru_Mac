      SUBROUTINE LAHE05T(POSGPT,SHAPDT,DVOLDT,PROPST,EHISTT,ELDIST,
     .                   SHAPET,WSTIFT,VELCMT,WSTI1T,WSTI2T,
     .                    frrrc, frrr1, frrr2, frr12,
     .                    fsssc, fsss1, fsss2,
     .                    fss12, fdvoc, fdvo1, fdvo2,
     .                   WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,
     .                   LNODST,DISIMT,HACHET,WHADEL,CENTRO,ADVEMT,
     .                   TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE PHASE CHANGE MATRIX (LATENT HEAT 
C     CONTRIBUTION INTO THE JACOBIAN MATRIX FOR ELEMENT NO. 5 ) IN A
C     UNIPHASE BIDIMENSIONAL SUBELEMENT.
C
C     (FOUR-NODED ELEMENTS)
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
      DIMENSION DVOLDT(*),        POSGPT(NDIMET,*),
     .          SHAPDT(NNODLT,*)
      DIMENSION WSTIFT(*)
C
      DIMENSION PROPST(*),        EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*), SHAPET(NNODLT,*)
      DIMENSION ELCODT(NDIMET,*)
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .          XJACMT(NDIMET,*), CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*)
      DIMENSION LNODST(*),        VELCMT(*),
     .          DISIMT(*)
      DIMENSION WSTI1T(*),        WSTI2T(*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          DVOLIT(*)
C
      CALL SHAP05T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,SHAPDT,
     .             DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,EHISTT,ELDIST,
     .              frrrc, frrr1, frrr2, frr12,
     .              fsssc, fsss1, fsss2, fss12,
     .              fdvoc, fdvo1, fdvo2,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT)
C
      CALL MAUF05T(DVOLDT,PROPST,SHAPDT,WSTIFT,EHISTT,VELCMT,WSTI1T,
     .             WSTI2T,ELDIST,DISIMT,WHADEL,TEINIT,FPCHLT,DVOLIT)
C
      RETURN
      END
