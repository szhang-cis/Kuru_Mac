      SUBROUTINE SKYCEK(GSTUP,NTERM,SAVAL)
C***********************************************************************
C
C**** THIS ROUTINE TESTS FOR RANK
C
C.... INPUTS
C
C       GSTUP(NTERM) -  COLUMN TO OF UNREDUCED ELEMENTS IN ARRAY
C       NTERM        -  NUMBER OF ELEMENTS IN COLUMN
C
C.... OUTPUTS
C
C       SAVAL        -  SUM OF ABSOLUTE VALUES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION GSTUP(*)
C
      SAVAL=0.
      DO 100 ITERM=1,NTERM
  100 SAVAL=SAVAL+DABS(GSTUP(ITERM))
C
      RETURN
      END
