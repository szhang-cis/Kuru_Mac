      SUBROUTINE PHCC08T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                   GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                   STRSGT,STRA0T,STRS0T,
     .                   BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                   SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                   SHAPDT,DVOLDT,
     .                   POSGPT,WEIGPT,
     .                   DERIDT,
     .                   CARDDT,GPCDDT,
     .                   ELEL1T,ELEL2T,
     .                  frrrc, frrr1, frrr2, frrr3, frr12, frr13, frr23,
     .                  fr123,
     .                  fsssc, fsss1, fsss2, fsss3, fss12, fss13, fss23,
     .                  fs123,
     .                  ftttc, fttt1, fttt2, fttt3, ftt12, ftt13, ftt23,
     .                  ft123,
     .                  fdvoc, fdvo1, fdvo2, fdvo3, fdv12, fdv13, fdv23,
     .                  fdv11, fdv22, fdv33, fd123, fd112, fd113, fd122,
     .                  fd133, fd223, fd233, f1123, f1223, f1233,
     .                  HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE ADVECTIVE PHASE-CHANGE RESIDUAL
C     "THERMAL FORCES" IN A UNIPHASE/MULTIPHASE TRIDIMENSIONAL
C     SUBELEMENT FOR THE CASE:
C
C     NON-ISOTHERMAL PHASE-CHANGE
C    
C     (NON-ISOTHERMAL PHASE-CHANGE)
C    
C     (EIGHT-NODED ELEMENTS)
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
      CALL SHAP08T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,SHAPDT,
     .             DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,EHISTT,ELDIST,
     .             DISIMT,
     .              frrrc, frrr1, frrr2, frrr3, frr12, frr13, frr23,
     .              fr123,
     .              fsssc, fsss1, fsss2, fsss3, fss12, fss13, fss23,
     .              fs123,
     .              ftttc, fttt1, fttt2, fttt3, ftt12, ftt13, ftt23,
     .              ft123,
     .              fdvoc, fdvo1, fdvo2, fdvo3, fdv12, fdv13, fdv23,
     .              fdv11, fdv22, fdv33, fd123, fd112, fd113, fd122,
     .              fd133, fd223, fd233, f1123, f1223, f1233,
     .             HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT)
C
      CALL FRUC05T(CARDDT,DVOLDT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .             GPCDDT,LNODST,PROPST,RMAT1T,SHAPDT,STRANT,
     .             STRSGT,STRA0T,STRS0T,
     .             BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .             SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .             WHADEL,ADVEMT,TEINIT,FPCHLT)
C
      RETURN
      END
