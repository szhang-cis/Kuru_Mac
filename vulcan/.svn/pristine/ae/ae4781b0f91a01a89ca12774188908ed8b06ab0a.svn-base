      DOUBLE PRECISION FUNCTION VECDOT(AVECT,BVECT,NTERM)
C***********************************************************************
C
C**** THIS FUNCTION PERFORMS THE VECTOR DOT PRODUCT
C
C.... INPUT PARAMETERS
C       AVECT(NTERM) -FIRST  VECTOR
C       BVECT(NTERM) -SECOND VECTOR
C       NTERM        -NUMBER OF TERMS IN THE VECTORS
C
C.... OUPUT PARAMETER
C       VECDOT       -DOT PRODUCT OF AVECT AND BVECT
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AVECT(*), BVECT(*)
C
       VECDOT=0.
       DO 100 ITERM=1,NTERM
  100  VECDOT=VECDOT+AVECT(ITERM)*BVECT(ITERM)
C
      RETURN
      END
