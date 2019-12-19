      SUBROUTINE PHCX12T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                   GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                   STRSGT,STRA0T,STRS0T,
     .                   BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                   SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                   SHAPDT,DVOLDT,
     .                   POSGPT,WEIGPT,
     .                   DERIDT,
     .                   CARDDT,GPCDDT,
     .                   ELEL1T,ELEL2T,
     .                    frrrc, frrr1, frrr2, frrr3,
     .                    fsssc, fsss1, fsss2, fsss3, 
     .                    ftttc, fttt1, fttt2, fttt3, 
     .                    fdvol,
     .                   HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE ADVECTIVE PHASE-CHANGE RESIDUAL
C     "THERMAL FORCES" IN A UNIPHASE TRIDIMENSIONAL SUBELEMENT FOR
C     THE CASE:
C
C     (NON-ISOTHERMAL PHASE-CHANGE)
C
C     (TETRAHEDRA)
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
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*)
      DIMENSION CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*)
C
      CALL SHAP12T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,SHAPDT,
     .             DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,EHISTT,ELDIST,
     .             DISIMT,
     .              frrrc, frrr1, frrr2, frrr3,
     .              fsssc, fsss1, fsss2, fsss3,
     .              ftttc, fttt1, fttt2, fttt3,
     .              fdvol,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT)
C
      CALL FRUX05T(CARDDT,DVOLDT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .             GPCDDT,LNODST,PROPST,RMAT1T,SHAPDT,STRANT,
     .             STRSGT,STRA0T,STRS0T,
     .             BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .             SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .             WHADEL,ADVEMT,TEINIT,FPCHLT)
C
      RETURN
      END
