      SUBROUTINE SHAP09T(CARDDT,DVOLDT,ELCODT,GPCDDT,LNODST,PROPST,
     .                   SHAPDT,DERIDT,POSGPT,WEIGPT,XJACMT,VELCMT,
     .                   EHISTT,ELDIST,DISIMT,
     .                    frrrc, frrr1, frrr2, frr12,
     .                    fsssc, fsss1, fsss2, fss12,
     .                    fdvoc, fdvo1, fdvo2,HACHET,WHADEL,CENTRO,
     .                   ADVEMT,TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE SHAPE FUNCTIONS
C     IN A BIPHASE BIDIMENSIONAL ELEMENT FOR THE CASES:
C
C     NPHCH=7
C     NPHCH=8
C     NPHCH=9
C
C     (QUADRILATERAL ELEMENTS)
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
C**** IDENTIFY SAMPLING POINTS
C
      posgpt(1,1)=0.0D0
      posgpt(2,1)=0.0D0
      posgpt(1,2)=0.0D0
      posgpt(2,2)=0.0D0
      posgpt(1,3)=0.0D0
      posgpt(2,3)=0.0D0
c
      weigpt(1)=2.0D0
      weigpt(2)=0.0D0
      weigpt(3)=0.0D0
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
      RDISCT=1.0D0*frrrc+exisp*frrr1+etasp*frrr2+frr12*exisp*etasp
      SDISCT=1.0D0*fsssc+exisp*fsss1+etasp*fsss2+fss12*exisp*etasp
      TDISCT=0.0D0
C
      CALL SHAFUN(DERIDT(1,1,IGAUST),RDISCT,SDISCT,TDISCT,NDIMET,NNODLT,
     .            NQUTRT,     0,SHAPDT(1,IGAUST))
C
C**** CARTESIAN DERIVATIVES 
C
      CALL JACOBST(CARDDT(1,1,IGAUST),DERIDT(1,1,IGAUST),DETJMT,
     .             ELCODT,            GPCDDT(1,IGAUST),   IELEMT,
     .             NDIMET,NNODLT,      SHAPDT(1,IGAUST),  XJACMT,
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
                      DVOLDT(IGAUST)=WEIGPT(IGAUST)*DETJMT*(fdvoc+
     .                               exisp*fdvo1+etasp*fdvo2)
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
