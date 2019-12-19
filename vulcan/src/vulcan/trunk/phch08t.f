      SUBROUTINE PHCH08T(PROPST,VELCMT,ELELMT,SHAPDT,DVOLDT,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,WEIGPT,DERIDT,XJACMT,
     .                   CARDDT,GPCDDT,ELCODT,LNODST,POSGPT,
     .                    frrrc, frrr1, frrr2, frrr3, frr12, frr13,
     .                    frr23, fr123,
     .                    fsssc, fsss1, fsss2, fsss3, fss12, fss13,
     .                    fss23, fs123,
     .                    ftttc, fttt1, fttt2, fttt3, ftt12, ftt13,
     .                    ftt23, ft123,
     .                    fdvoc, fdvo1, fdvo2, fdvo3, fdv12, fdv13,
     .                    fdv23, fdv11, fdv22, fdv33, fd123, fd112,
     .                    fd113, fd122, fd133, fd223, fd233, f1123,
     .                    f1223, f1233,
     .                   ELEL1T,ELEL2T,INPCCT,HACHET,WHADEL,CENTRO,
     .                   ADVEMT,TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE PHASE-CHANGE RESIDUAL "THERMAL FORCES"
C     IN A UNIPHASE/MULTIPHASE TRIDIMENSIONAL SUBELEMENT FOR THE CASE:
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
      DIMENSION SHAPDT(NNODLT,*),        DVOLDT(*)
C
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
      CALL FRUF05T(PROPST,VELCMT,ELELMT,DVOLDT,SHAPDT,EHISTT,ELDIST,
     .             DISIMT,ELELTT,ELEL1T,ELEL2T,WHADEL,INPCCT,TEINIT,
     .             FPCHLT,DVOLIT)
C
      RETURN
      END
