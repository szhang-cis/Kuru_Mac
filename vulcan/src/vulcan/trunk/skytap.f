      SUBROUTINE SKYTAP(KRESL,KSYMM,GSTDI,GSTLO,GSTUP,IDSK2,NLAST,
     .                  NEQNS)
C***********************************************************************
C
C*** THIS ROUTINE STORES A NEW FACTORISED GLOBAL MATRIX OR READS BACK
C    FROM TAPE AN OLD ONE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION GSTDI(NEQNS), GSTLO(NLAST), GSTUP(NLAST)
C
      IF(KRESL.EQ.1)THEN
        IF(KSYMM.EQ.0) WRITE(IDSK2) GSTDI,GSTLO,GSTUP
        IF(KSYMM.EQ.1) WRITE(IDSK2) GSTDI,GSTUP
      ELSE
        IF(KSYMM.EQ.0) READ (IDSK2) GSTDI,GSTLO,GSTUP
        IF(KSYMM.EQ.1) READ (IDSK2) GSTDI,GSTUP
      ENDIF
C
      RETURN
      END
