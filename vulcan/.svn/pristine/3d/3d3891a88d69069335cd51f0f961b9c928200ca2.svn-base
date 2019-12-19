      SUBROUTINE LAHC01T(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .                   PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .                   DMATXT,SIGMAT,XJACMT,
     .                   SHAPDT,DVOLDT,
     .                   POSGPT,WEIGPT,
     .                   DERIDT,
     .                   CARDDT,GPCDDT,
     .                   VELCMT,WSTI1T,
     .                   WSTI2T,
     .                   LNODST,ELCODT,
     .                    frrrc, frrr1, fdvol,HACHET,WHADEL,CENTRO,
     .                   ADVEMT,TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE ADVECTIVE PHASE-CHANGE MATRIX 
C     (CONVECTIVE LATENT HEAT CONTRIBUTION INTO THE JACOBIAN MATRIX FOR
C     ELEMENT NO. 5 ) IN A UNIPHASE UNIDIMENSIONAL SUBELEMENT
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
      DIMENSION CARTDT(NDIMET,NNODLT,*), DVOLUT(*),
     .          EHISTT(NHISTT,*),        ELDIST(NDOFCT,*),
     .          EPMTXT(*),               GPCODT(NDIMET,*),
     .          PROPST(*),               SHAPET(NNODLT,*),
     .          STRSGT(NSTR1T,*),        ESTIFT(*),
     .          HSTIFT(NEVABT,NNODET)
      DIMENSION BMATXT(NSTR1T,*),        DMATXT(NSTR1T,*),
     .          SIGMAT(*),               XJACMT(NDIMET,*)
C
      DIMENSION DVOLDT(*),        POSGPT(NDIMET,*),
     .          SHAPDT(NNODLT,*)
C
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .                            CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*), VELCMT(*)
      DIMENSION WSTI1T(*),        WSTI2T(*)
      DIMENSION LNODST(*)
      DIMENSION ELCODT(NDIMET,*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*)
      DIMENSION CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*)
C
      CALL SHAP01T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,SHAPDT,
     .             DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,EHISTT,ELDIST,
     .              frrrc, frrr1, fdvol,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT)
C
      CALL STUC05T(CARDDT,DVOLDT,EHISTT,ELDIST,EPMTXT,GPCDDT,
     .             PROPST,SHAPDT,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .             DMATXT,SIGMAT,XJACMT,
     .             VELCMT,WSTI1T,WSTI2T,WHADEL,ADVEMT,TEINIT,FPCHLT)
C
      RETURN
      END
