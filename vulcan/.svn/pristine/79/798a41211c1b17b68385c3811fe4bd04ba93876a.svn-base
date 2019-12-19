      SUBROUTINE TESTM1(ESTIF,IELEM,NEVAB)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE SUM OF ALL TERMS IN A MATRIX
C     STORED AS AN ARRAY
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ESTIF(*)
C
      SUMAI=0.
      DO 10 IRAWS=1,NEVAB
      KOUNT=IRAWS
        DO 20 JRAWS=1,IRAWS
          SUMAI=SUMAI+ESTIF(KOUNT)
          KOUNT=KOUNT+(NEVAB-JRAWS)
   20   CONTINUE
C
      KOUNT=KOUNT-(NEVAB-(JRAWS-1))
        DO 30 KRAWS=JRAWS,NEVAB
          KOUNT=KOUNT+1
          SUMAI=SUMAI+ESTIF(KOUNT)
   30   CONTINUE
   10 CONTINUE
c     PRINT *, IELEM,SUMAI
      RETURN
      END
