      SUBROUTINE SHAP08T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,
     .                   SHAPDT,DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,
     .                   EHISTT,ELDIST,DISIMT,
     .                    frrrc, frrr1, frrr2, frrr3, frr12, frr13,
     .                    frr23, fr123,
     .                    fsssc, fsss1, fsss2, fsss3, fss12, fss13,
     .                    fss23, fs123,
     .                    ftttc, fttt1, fttt2, fttt3, ftt12, ftt13,
     .                    ftt23, ft123,
     .                    fdvoc, fdvo1, fdvo2, fdvo3, fdv12, fdv13,
     .                    fdv23, fdv11, fdv22, fdv33, fd123, fd112,
     .                    fd113, fd122, fd133, fd223, fd233, f1123,
     .                    f1223, f1233,HACHET,WHADEL,CENTRO,ADVEMT,
     .                   TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE SHAPE FUNCTIONS
C     IN A BIPHASE/MULTUPHASE TRIDIMENSIONAL ELEMENT FOR THE CASE:
C
C     NON-ISOTHERMAL PHASE-CHANGE
C    
C     (EIGHT-NODED BRICK ELEMENTS)
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
      COMMON/THICKNESS/THICKT
C
      DIMENSION CARDDT(NDIMET,NNODLT,*), DVOLDT(*),
     .          ELCODT(NDIMET,*),        GPCDDT(NDIMET,*),
     .          LNODST(*),               PROPST(*),
     .          SHAPDT(NNODLT,*)
      DIMENSION DERIDT(NDIMET,NNODLT,*), POSGPT(NDIMET,*),
     .          WEIGPT(*),               XJACMT(NDIMET,*)
C
      DIMENSION VELCMT(*),               EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*),        DISIMT(*)
      DIMENSION HACHET(NNODLT),          WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*),        ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*),        FPCHLT(NFPCH,*)
C
      TWOPIT=6.283185307179586D0
C
C**** IDENTIFY SAMPLING POINTS AND WEIGHTS
C
      CALL RULEPW(NDIMET,NGAULT,NRULET,POSGPT,WEIGPT)
C
C**** LOOP ON INTEGRATION POINTS
C
      DO 100 IGAUST=1,NGAULT
C
C**** COMPUTE SHAPE FUNCTIONS AND DERIVATIVES
C
      EXISP=POSGPT(1,IGAUST)
      ETASP=0.0D0
      IF(NDIMET.GE.2) ETASP=POSGPT(2,IGAUST)
      EZETA=0.0D0
      IF(NDIMET.EQ.3) EZETA=POSGPT(3,IGAUST)
C
      RDISCT=1.0D0*frrrc+frrr1*exisp+frrr2*etasp+frrr3*ezeta+
     .       frr12*exisp*etasp+frr13*exisp*ezeta+frr23*etasp*ezeta+
     .       fr123*exisp*etasp*ezeta
      SDISCT=1.0D0*fsssc+fsss1*exisp+fsss2*etasp+fsss3*ezeta+
     .       fss12*exisp*etasp+fss13*exisp*ezeta+fss23*etasp*ezeta+
     .       fs123*exisp*etasp*ezeta
      TDISCT=1.0D0*ftttc+fttt1*exisp+fttt2*etasp+fttt3*ezeta+
     .       ftt12*exisp*etasp+ftt13*exisp*ezeta+ftt23*etasp*ezeta+
     .       ft123*exisp*etasp*ezeta
C
      CALL SHAFUN(DERIDT(1,1,IGAUST),RDISCT,SDISCT,TDISCT,NDIMET,NNODLT,
     .            NQUTRT,     0,SHAPDT(1,IGAUST))
C
C**** CARTESIAN DERIVATIVES 
C
      CALL JACOBST(CARDDT(1,1,IGAUST), DERIDT(1,1,IGAUST), DETJMT,
     .             ELCODT,             GPCDDT(1,IGAUST),   IELEMT,
     .             NDIMET,NNODLT,      SHAPDT(1,IGAUST),   XJACMT,
     .             LUREST,LUPRIT)
C
C**** MAKE A CONTROL ON THE CORRECTNESS OF THE SHAPE FUNCTION AND THEIR 
C     DERIV.
C
      ICHEKT=0
      IF(ICHEKT.EQ.1)CALL TESTM2(CARDDT(1,1,IGAUST),DERIDT(1,1,IGAUST),
     .                           SHAPDT(1,IGAUST),IGAUST,NDIMET,NNODLT)
C
C**** INTEGRATION WEIGHTS
C
      fdnou=fdvoc+fdvo1*exisp+fdvo2*etasp+fdvo3*ezeta+
     .      fdv12*exisp*etasp+fdv13*exisp*ezeta+fdv23*etasp*ezeta+
     .      fdv11*exisp*exisp+fdv22*etasp*etasp+fdv33*ezeta*ezeta+
     .      fd123*exisp*etasp*ezeta+fd112*exisp*exisp*etasp+
     .      fd113*exisp*exisp*ezeta+fd122*exisp*etasp*etasp+
     .      fd133*exisp*ezeta*ezeta+fd223*etasp*etasp*ezeta+
     .      fd233*etasp*ezeta*ezeta+
     .      f1123*exisp*exisp*etasp*ezeta+
     .      f1223*exisp*etasp*etasp*ezeta+
     .      f1233*exisp*etasp*ezeta*ezeta
C
      if(fdnou.le.0.0D0)
     . CALL RUNENDT('SHAP08T: DVOLD=0.0 !!!!!           ')
C
                      DVOLDT(IGAUST)=WEIGPT(IGAUST)*DETJMT*fdnou
      IF(NTYPET.EQ.3) DVOLDT(IGAUST)=DVOLDT(IGAUST)*TWOPIT*
     .                               GPCDDT(1,IGAUST)
                      DVOLDT(IGAUST)=DVOLDT(IGAUST)*THICKT
C
C**** DEFINES THE WEIGHT FUNCTION WHEN IGALET=0
C
      IF(IGALET.EQ.0) THEN
       DO INODLT=1,NNODLT
        WHADEL(INODLT,IGAUST)=SHAPDT(INODLT,IGAUST)
       ENDDO
       GO TO 100
      ENDIF
C
C**** DEFINES THE UPWIND FUNCTION
C
      NSUPE=NSUBD(1,1)
      CALL HACHEC(ELCODT,PROPST,ELDIST,VELCMT,HACHET,CENTRO,
     .            ADVEMT,TEINIT,SHAPDT(1,IGAUST),NSUPE,FPCHLT,
     .            EHISTT(1,IGAUST),CARDDT(1,1,IGAUST))
C
C**** COMPUTES PERTURBATION FUNCTIONS AND DERIVATIVES
C
      CALL WHAFUN(DERIDT(1,1,IGAUST),RDISCT,SDISCT,TDISCT,NDIMET,NNODLT,
     .            WHADEL(1,IGAUST),HACHET,XJACMT,IPERT)
C
C**** COMPUTES WEIGHT FUNCTION (SHAPE FUNCTION + PERTURBATION FUNCTION)
C
      DO INODLT=1,NNODLT
       WHADEL(INODLT,IGAUST)=SHAPDT(INODLT,IGAUST)+
     .                       WHADEL(INODLT,IGAUST)
      ENDDO
C
  100 CONTINUE
C
      RETURN
      END
