      SUBROUTINE WRITDK(NUNIT,MATRI,LENR1,NRCUR)
C***********************************************************************
C
C**** THIS ROUTINE WRITES TO DISK 'NUNIT' A MATRIX WHICH LENGTH
C     IS FIXED AS LENR1
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION MATRI(LENR1)
C
      WRITE(NUNIT,REC=NRCUR) MATRI
C
      RETURN
      END
