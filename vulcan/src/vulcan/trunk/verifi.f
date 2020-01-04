
      SUBROUTINE VERIFI(DESMN,DESJ2,DESJ3,DIFFE,EHIST,NCRIT,PREMN,
     .                  PREJ2,PREJ3,PROPS,REDUC,WCOJ3,YCOJ3,ZCOJ2)
C***********************************************************************
C
C****THIS ROUTINE CHECKS IF THE TRIAL REDUCTION SATISFIES THE
C    BOUNDARY  SURFACE.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION EHIST(*), PROPS(*)
      IF(NCRIT.EQ.3)GOTO 10
C
C***LINEAR PLASTIC MODEL
C
      CALL VERIF1(DESMN,DESJ2,DESJ3,DIFFE,PREMN,PREJ2,PREJ3,
     .            PROPS,REDUC,WCOJ3,YCOJ3,ZCOJ2)
C
      RETURN
C
C***CAM-CLAY
C
   10 CALL VERIF2(DESMN,DESJ2,DESJ3,DIFFE,EHIST,NCRIT,PREMN,
     .            PREJ2,PREJ3,PROPS,REDUC,WCOJ3,YCOJ3,ZCOJ2)
C
      RETURN
      END