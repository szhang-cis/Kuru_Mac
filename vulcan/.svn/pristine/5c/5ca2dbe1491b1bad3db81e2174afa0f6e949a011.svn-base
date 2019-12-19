      SUBROUTINE LAHE08T(POSGPT,SHAPDT,DVOLDT,PROPST,EHISTT,ELDIST,
     .                   SHAPET,WSTIFT,VELCMT,WSTI1T,WSTI2T,
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
     .                   WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,
     .                   LNODST,DISIMT,HACHET,WHADEL,CENTRO,ADVEMT,
     .                   TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE PHASE CHANGE MATRIX (LATENT HEAT 
C     CONTRIBUTION INTO THE JACOBIAN MATRIX FOR ELEMENT NO. 5 ) IN A
C     UNIPHASE TRIDIMENSIONAL SUBELEMENT.
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
      CALL MAUF05T(DVOLDT,PROPST,SHAPDT,WSTIFT,EHISTT,VELCMT,WSTI1T,
     .             WSTI2T,ELDIST,DISIMT,WHADEL,TEINIT,FPCHLT,DVOLIT)
C
      RETURN
      END
