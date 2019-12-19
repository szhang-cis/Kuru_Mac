      SUBROUTINE VERIPCT(ELEL2,ELEL1,NINDI,IVEPC)
C***********************************************************************
C
C**** THIS ROUTINE CONTROLS THAT ELEL2 >= ELEL1 (IVEPC=0)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ELEL1(*), ELEL2(*)
C
      IVEPC=0
      DO IINDI=1,NINDI
       WS1=ELEL1(IINDI)
       WS2=ELEL2(IINDI)
       IF(WS1.GT.WS2) IVEPC=IVEPC+1
      ENDDO
C
      RETURN
      END
