      SUBROUTINE FRTF05T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,
     .                   NPHCHT,POSGPT,ELCODT,WEIGPT,DERIDT,XJACMT,
     .                   CARDDT,GPCDDT,ELEL1T,ELEL2T,HACHET,
     .                   WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT,DVOLIT,
     .                   LNODST,INPCCT,TEMPLT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES" DUE TO 
C     PHASE-CHANGE FOR BIPHASED ELEMENTS
C     ( ELEMENT  NO 5 )
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
      DIMENSION LNODST(*),        ELEL1T(*),
     .          ELEL2T(*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          DVOLIT(*)
C
      GOTO(1,2,3), NDIMET
C
    1 CALL FRTF1DT(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,NPHCHT,POSGPT,
     .             ELCODT,WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,
     .             ELEL1T,ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,
     .             ADVEMT,TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    2 CALL FRTF2DT(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,NPHCHT,POSGPT,
     .             ELCODT,WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,
     .             ELEL1T,ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,
     .             ADVEMT,TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    3 CALL RUNENDT('FRTF05T:NDIME=3                    ')
      RETURN
C
      END
