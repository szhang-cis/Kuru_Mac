      SUBROUTINE PHCC01T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                   GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                   STRSGT,STRA0T,STRS0T,
     .                   BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                   SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                   SHAPDT,DVOLDT,
     .                   POSGPT,WEIGPT,
     .                   DERIDT,
     .                   CARDDT,GPCDDT,
     .                   ELEL1T,ELEL2T,
     .                    frrrc, frrr1, fdvol,HACHET,WHADEL,CENTRO,
     .                   ADVEMT,TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE PHASE-CHANGE RESIDUAL "THERMAL FORCES"
C     IN A MULTIPHASE UNIDIMENSIONAL SUBELEMENT FOR THE CASE
C
C     (NON-ISOTHERMAL PHASE-CHANGE)
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
      CALL SHAP01T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,SHAPDT,
     .             DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,EHISTT,ELDIST,
     .              frrrc, frrr1, fdvol,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT)
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
