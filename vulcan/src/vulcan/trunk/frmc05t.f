      SUBROUTINE FRMC05T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                   GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                   STRSGT,STRA0T,STRS0T,
     .                   BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                   SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                   SHAPDT,DVOLDT,POSGPT,WEIGPT,DERIDT,
     .                   CARDDT,GPCDDT,
     .                   ELEL1T,ELEL2T,
     .                   HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE INTERNAL RESISTING HEATS DUE TO
C     PHASE-CHANGE CONVECTION EFFECTS
C     ( FOR ELEMENT NO. 5 )
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
     .          EHISTT(NHISTT,*),        ELCODT(NDIMET,*),
     .          ELDIST(NDOFCT,*),
     .          EPMTXT(*),               GPCODT(NDIMET,*), 
     .          LNODST(*),               PROPST(*),
     .          RMAT1T(NDIMET,*),        SHAPET(NNODLT,*),
     .          STRANT(NSTR1T,*),        STRSGT(NSTR1T,*),
     .          STRA0T(NSTR1T,*),        STRS0T(NSTR1T,*)
      DIMENSION BMATXT(NSTR1T,*),        BMSIGT(*),
     .          DESIGT(*),               DMATXT(NSTR1T,*),
     .          DSTRAT(*),               PRESGT(*),
     .          SGTOTT(*),               SIGMAT(*),
     .          TSTRAT(*),               XJACMT(NDIMET,*),
     .          ELELTT(*),               VELCMT(*)
C
      DIMENSION DVOLDT(*),        SHAPDT(NNODLT,*),
     .          POSGPT(NDIMET,*)
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .                            CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*)
      DIMENSION ELEL1T(*),        ELEL2T(*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*),
     .          ADVEMT(NDIMET,*), TEINIT(NDOFCT,*),
     .          FPCHLT(NFPCH,*)
C
      GOTO(1,2,3), NDIMET
C
    1 CALL FRMC1DT(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .             GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .             STRSGT,STRA0T,STRS0T,
     .             BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .             SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .             SHAPDT,DVOLDT,
     .             POSGPT,WEIGPT,
     .             DERIDT,
     .             CARDDT,GPCDDT,
     .             ELEL1T,ELEL2T,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT)
      RETURN
C
    2 CALL FRMC2DT(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .             GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .             STRSGT,STRA0T,STRS0T,
     .             BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .             SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .             SHAPDT,DVOLDT,
     .             POSGPT,WEIGPT,
     .             DERIDT,
     .             CARDDT,GPCDDT,
     .             ELEL1T,ELEL2T,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT)
      RETURN
C
    3 CALL FRMC3DT(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .             GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .             STRSGT,STRA0T,STRS0T,
     .             BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .             SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .             SHAPDT,DVOLDT,
     .             POSGPT,WEIGPT,
     .             DERIDT,
     .             CARDDT,GPCDDT,
     .             ELEL1T,ELEL2T,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT)
      RETURN
C
      END
