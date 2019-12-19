      SUBROUTINE PHCH11T(PROPST,VELCMT,ELELMT,SHAPDT,DVOLDT,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,
     .                    frrrc, frrr1, frrr2, frrr3,
     .                    fsssc, fsss1, fsss2, fsss3,
     .                    ftttc, fttt1, fttt2, fttt3,
     .                    fdvol,
     .                   WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,
     .                   LNODST,POSGPT,ELEL1T,ELEL2T,INPCCT,HACHET,
     .                   WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE PHASE-CHANGE RESIDUAL "THERMAL FORCES"
C     IN A UNIPHASE TRIDIMENSIONAL SUBELEMENT FOR THE CASE:
C
C     (ISOTHERMAL PHASE-CHANGE)  not implemented yet !!!!
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
      DIMENSION SHAPDT(NNODLT,*),        DVOLDT(*)
      DIMENSION PROPST(*),               VELCMT(*),
     .          ELELMT(*),               EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*),        DISIMT(*), 
     .          ELELTT(*)
      DIMENSION ELCODT(NDIMET,*)
      DIMENSION WEIGPT(*),               POSGPT(NDIMET,*),
     .          DERIDT(NDIMET,NNODLT,*), XJACMT(NDIMET,*),
     .          CARDDT(NDIMET,NNODLT,*), GPCDDT(NDIMET,*)
      DIMENSION LNODST(*)
      DIMENSION ELEL1T(*),               ELEL2T(*)
      DIMENSION HACHET(NNODLT),          WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*),        ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*),        FPCHLT(NFPCH,*),
     .          DVOLIT(*)
C
      CALL SHAP11T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,SHAPDT,
     .             DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,EHISTT,ELDIST,
     .             DISIMT,
     .              frrrc, frrr1, frrr2, frrr3,
     .              fsssc, fsss1, fsss2, fsss3,
     .              ftttc, fttt1, fttt2, fttt3,
     .              fdvol,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT)
C
      CALL FRUF05T(PROPST,VELCMT,ELELMT,DVOLDT,SHAPDT,EHISTT,ELDIST,
     .             DISIMT,ELELTT,ELEL1T,ELEL2T,WHADEL,INPCCT,TEINIT,
     .             FPCHLT,DVOLIT)
C
      RETURN
      END
