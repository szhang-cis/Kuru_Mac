      SUBROUTINE FR1D01T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,
     .                   POSGPT,ELCODT,WEIGPT,DERIDT,XJACMT,CARDDT,
     .                   GPCDDT,LNODST,ELEL1T,ELEL2T,INPCCT,TEMPLT,
     .                   HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT,
     .                   DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES" DUE TO 
C     PHASE-CHANGE WHEN IT HAS TWO PHASES ONLY FOR 1D ELEMENTS
C
C     FOR THE CASE: NPHCHT=1
C    
C***********************************************************************
C***********************************************************************
C     CASE
C
C NDIMET   NNODET   NISOTT             CASE               NPHCHT
C***********************************************************************
C
C   1        2        1       *---------R----------*        1
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
      DIMENSION PROPST(*),        VELCMT(*),
     .          ELELMT(*),        DVOLUT(*),
     .          SHAPET(NNODLT,*), EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*), DISIMT(*),
     .          ELELTT(*)
      DIMENSION DVOLDT(*),        SHAPDT(NNODLT,*),
     .          POSGPT(NDIMET,*)
      DIMENSION GPCODT(NDIMET,*), ELCODT(NDIMET,*)
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .          XJACMT(NDIMET,*), CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*)
      DIMENSION LNODST(*)
      DIMENSION ELEL1T(*),        ELEL2T(*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          DVOLIT(*)
C
C**** CALCULATE THE INTERFACE COORDINATES AND THE SHAPE FUNCTIONS FOR 
C     BIPHASE ELEMENTS
C
      IF(INPCCT.EQ.1) THEN
       ELDI1=ELDIST(1,1)-VELCMT(1)*DTIMET
       ELDI2=ELDIST(1,2)-VELCMT(2)*DTIMET
       CALL COIN01T(ELDI2,ELDI1,TEMPLT,RCOOR,ISOLI)
      ELSE
       ELDI1=ELDIST(1,1)
       ELDI2=ELDIST(1,2)
       CALL COIN01T(ELDI2,ELDI1,TEMPLT,RCOOR,ISOLI)
      ENDIF
C
      r1=-1.0
      r2=rcoor
      frrrc=0.5*(r1+r2)
      frrr1=0.5*(r2-r1)
      fdvol=frrr1
      CALL PHCH00T(POSGPT,SHAPDT,DVOLDT,PROPST,EHISTT,ELDIST,SHAPET,
     .             ELELMT,VELCMT,ELEL1T,ELEL2T,INPCCT,
     .              frrrc, frrr1, fdvol,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,LNODST,
     .             DISIMT,ELELTT,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .             FPCHLT,DVOLIT)
C
      r1=rcoor
      r2=1.0
      frrrc=0.5*(r1+r2)
      frrr1=0.5*(r2-r1)
      fdvol=frrr1
      CALL PHCH00T(POSGPT,SHAPDT,DVOLDT,PROPST,EHISTT,ELDIST,SHAPET,
     .             ELELMT,VELCMT,ELEL1T,ELEL2T,INPCCT,
     .              frrrc, frrr1, fdvol,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,LNODST,
     .             DISIMT,ELELTT,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .             FPCHLT,DVOLIT)
      RETURN
      END
