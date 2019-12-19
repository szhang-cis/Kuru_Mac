      SUBROUTINE FRTF2DT(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,
     .                   NPHCHT,POSGPT,ELCODT,WEIGPT,DERIDT,XJACMT,
     .                   CARDDT,GPCDDT,LNODST,ELEL1T,ELEL2T,INPCCT,
     .                   TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .                   FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES" DUE TO 
C     PHASE-CHANGE WHEN IT HAS TWO PHASES ONLY FOR 2D ELEMENTS
C    
C***********************************************************************
C     EHIST(1) = Density
C     EHIST(2) = Temperature Derivative of Density
C     EHIST(3) = Specific Heat coefficient)
C     EHIST(4) = Temperature Derivative of the Specific Heat coefficient
C***********************************************************************
C***********************************************************************
C    CASES
C
C NDIME   NNODE   NISOT                    CASE              NPHCH
C***********************************************************************
C
C                                 4                    3
C   2       4       1             *--------------------*       1
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 S                    !
C                                 !                    !
C                                 !                    !
C                                 *-------R------------*
C                                 1                    2
C
C                                 4                    3
C   2       4       1             *--------------------*       2
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 !                    S
C                                 !                    !
C                                 !                    !
C                                 *-------R------------*
C                                 1                    2
C
C                                 4                    3
C   2       4       1             *-------R------------*       3
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 !                    S
C                                 !                    !
C                                 !                    !
C                                 *--------------------*
C                                 1                    2
C
C                                 4                    3
C   2       4       1             *-------R------------*       4
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 S                    !
C                                 !                    !
C                                 !                    !
C                                 *--------------------*
C                                 1                    2
C
C                                 4                    3
C   2       4       1             *-------R------------*       5
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 *-------R------------*
C                                 1                    2
C
C                                 4                    3
C   2       4       1             *--------------------*       6
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 S                    S
C                                 !                    !
C                                 !                    !
C                                 *--------------------*
C                                 1                    2
C
C                                 4                    3
C   2       4       1             *-------R------------*       10
C                                 !                    !
C                                 !                    !
C                                 !                    !
C                                 S                    S
C                                 !                    !
C                                 !                    !
C                                 *-------R------------*
C                                 1                    2
C
C   2       3       1            3*                            7
C                                 !  . 
C                                 !     . 
C                                 !        X 
C                                 !           . 
C                                 !              . 
C                                 !                 . 
C                                 *--------R-----------*
C                                 1                    2
C
C   2       3       1            3*                            8
C                                 !  . 
C                                 !     . 
C                                 S        . 
C                                 !           . 
C                                 !              . 
C                                 !                 . 
C                                 *--------R-----------*
C                                 1                    2
C
C   2       3       1            3*                            9
C                                 !  . 
C                                 !     . 
C                                 S        X 
C                                 !           . 
C                                 !              . 
C                                 !                 . 
C                                 *--------------------*
C                                 1                    2
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
      GOTO(1,2,3,4,5,6,7,8,9,10), NPHCHT
C
    1 CALL FR2D01T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    2 CALL FR2D02T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    3 CALL FR2D03T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    4 CALL FR2D04T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    5 CALL FR2D05T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    6 CALL FR2D06T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    7 CALL FR2D07T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    8 CALL FR2D08T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
    9 CALL FR2D09T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
   10 CALL FR2D10T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,ELDIST,
     .             DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,POSGPT,ELCODT,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,LNODST,ELEL1T,
     .             ELEL2T,INPCCT,TEMPLT,HACHET,WHADEL,CENTRO,ADVEMT,
     .             TEINIT,FPCHLT,DVOLIT)
      RETURN
C
      END
