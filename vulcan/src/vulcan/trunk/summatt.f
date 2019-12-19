      SUBROUTINE SUMMATT(ELELM,ELEL2,ELEL1,NINDI)
C***********************************************************************
C
C****THIS ROUTINE SUM ELEL1 AND ELEL2
C    
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ELELM(*),  ELEL2(*),  ELEL1(*)
C
      DO IINDI=1,NINDI
       ELELM(IINDI)=ELELM(IINDI)+ELEL2(IINDI)-ELEL1(IINDI)
      ENDDO
C
      RETURN
      END
