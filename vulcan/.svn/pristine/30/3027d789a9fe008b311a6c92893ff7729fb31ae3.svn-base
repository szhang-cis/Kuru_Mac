      SUBROUTINE STMC05T(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .                   PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .                   DMATXT,SIGMAT,XJACMT,SHAPDT,DVOLDT,
     .                   POSGPT,WEIGPT,DERIDT,CARDDT,GPCDDT,
     .                   VELCMT,WSTI1T,
     .                   WSTI2T,ELCODT,HACHET,WHADEL,CENTRO,ADVEMT,
     .                   TEINIT,FPCHLT,LNODST)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE JACOBIAN MATRIX DUE TO ADVECTIVE
C     PHASE-CHANGE EFFECTS FOR MULTIPHASE ELEMENTS
C     ( ELEMENT TYPE NO. 5 )
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
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*)
C
      GOTO(1,2,3), NDIMET
C
    1 CALL STMC1DT(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .             PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .             DMATXT,SIGMAT,XJACMT,
     .             SHAPDT,DVOLDT,
     .             POSGPT,WEIGPT,
     .             DERIDT,
     .             CARDDT,GPCDDT,
     .             VELCMT,WSTI1T,
     .             WSTI2T,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .             LNODST,ELCODT,FPCHLT)
      RETURN
C
    2 CALL STMC2DT(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .             PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .             DMATXT,SIGMAT,XJACMT,
     .             SHAPDT,DVOLDT,
     .             POSGPT,WEIGPT,
     .             DERIDT,
     .             CARDDT,GPCDDT,
     .             VELCMT,WSTI1T,
     .             WSTI2T,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .             LNODST,ELCODT,FPCHLT)
      RETURN
C
    3 CALL STMC3DT(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .             PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .             DMATXT,SIGMAT,XJACMT,
     .             SHAPDT,DVOLDT,
     .             POSGPT,WEIGPT,
     .             DERIDT,
     .             CARDDT,GPCDDT,
     .             VELCMT,WSTI1T,
     .             WSTI2T,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .             LNODST,ELCODT,FPCHLT)
      RETURN
C
      END
