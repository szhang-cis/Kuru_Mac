      SUBROUTINE LAHC08T(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .                   PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .                   DMATXT,SIGMAT,XJACMT,SHAPDT,DVOLDT,
     .                   POSGPT,WEIGPT,DERIDT,CARDDT,GPCDDT,
     .                   VELCMT,WSTI1T,WSTI2T,LNODST,ELCODT,
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
     .                   HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE ADVECTIVE PHASE-CHANGE MATRIX 
C     (ADVECTIVE LATENT HEAT CONTRIBUTION INTO THE JACOBIAN MATRIX FOR
C     ELEMENT NO. 5 ) IN A UNIPHASE TRIDIMENSIONAL SUBELEMENT
C
C     (8-NODED ELEMENTS)
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
      CALL STUC05T(CARDDT,DVOLDT,EHISTT,ELDIST,EPMTXT,GPCDDT,
     .             PROPST,SHAPDT,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .             DMATXT,SIGMAT,XJACMT,
     .             VELCMT,WSTI1T,WSTI2T,WHADEL,ADVEMT,TEINIT,FPCHLT)
C
      RETURN
      END
