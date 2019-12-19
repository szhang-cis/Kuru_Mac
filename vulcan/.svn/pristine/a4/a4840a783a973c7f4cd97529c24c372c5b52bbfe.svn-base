      SUBROUTINE TOTSIG(ELDIS,KPORE,NDOFC,NNODE,NSTR1,PROPS,SIGMA,
     .                  SHAPE)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE TOTAL STRESS AT GAUSSIAN POINT
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ELDIS(NDOFC,*), SHAPE(*), SIGMA(*)
      DIMENSION DELTA(6),       PROPS(*)
C
ctm      DATA DELTA/1.,1.,0.,1.,0.,0./
C
      DELTA(1)=1.
      DELTA(2)=1.
      DELTA(3)=0.
      DELTA(4)=1.
      DELTA(5)=0.
      DELTA(6)=0.
C
      COMPW=PROPS(19)                      ! ?????
      IF(COMPW.EQ.0.0) COMPW=1.0D00
C
      PORES=0.0D00
      DO 20 INODE=1,NNODE
 20   PORES=PORES+ELDIS(NDOFC,INODE)*SHAPE(INODE)
      PORES=COMPW*PORES
C
C**** CALULATE THE TOTAL STRESS COMPONENTS
C
      DO 30 ISTR1=1,NSTR1
   30 SIGMA(ISTR1)=SIGMA(ISTR1)-DELTA(ISTR1)*PORES
C
      RETURN
      END
