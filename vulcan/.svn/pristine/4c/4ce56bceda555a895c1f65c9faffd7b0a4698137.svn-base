      SUBROUTINE FRTF1DT(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,
     .                   NPHCHT,POSGPT,ELCODT,WEIGPT,DERIDT,XJACMT,
     .                   CARDDT,GPCDDT,LNODST,ELEL1T,ELEL2T,INPCCT,
     .                   TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .                   FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES" DUE TO 
C     PHASE-CHANGE WHEN IT HAS TWO PHASES ONLY FOR 1D ELEMENTS
C     (2 & 3 NODED)
C    
C***********************************************************************
C     EHIST(1) = Density
C     EHIST(2) = Temperature Derivative of Density
C     EHIST(3) = Specific Heat coefficient)
C     EHIST(4) = Temperature Derivative of the Specific Heat coefficient
C***********************************************************************
C***********************************************************************
C    CASE
C
C NDIME   NNODE   NISOT              CASE              NPHCH
C***********************************************************************
C
C   1       2       1       *---------R----------*       1
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
C     BIPHASE ELEMENTS (2 & 3 NODED)
C
      GOTO(1,2), NPHCHT
C
    1 CONTINUE
      CALL FR1D01T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    2 CONTINUE
cc
cc not implemented yet (for 3-noded elements)
cc
cc    2 CALL FR1D02T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
cc       .           DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
cc       .           WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
cc       .           ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
cc     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
      END
