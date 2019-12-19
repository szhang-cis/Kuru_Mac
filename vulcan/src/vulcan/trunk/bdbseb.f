      SUBROUTINE BDBSEB(BMATX,CARTD,GPCOD,NDIME,NDOFN,NEVAB,NNODE,NSTR1,
     .                  NTYPE,SHAPE)
C**********************************************************************
C
C****THIS ROUTINE CALCULATE B-MATRIX
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION CARTD(NDIME,*), BMATX(NSTR1,*), SHAPE(*)
C
      DO ISTR1=1,NSTR1
       DO IEVAB=1,NEVAB
        BMATX(ISTR1,IEVAB)=0.0
       ENDDO
      ENDDO
C
      GOTO (1,1,3,4,5),NTYPE
C
C***PLANE STRESS OR STRAIN
C
    1 CONTINUE
C
      I1=0
      DO 10 I=1,NNODE
      I1=(I-1)*NDOFN+1
      I2=I1+1
C
      BMATX(1,I1)=CARTD(1,I)
      BMATX(2,I2)=CARTD(2,I)
      BMATX(3,I1)=CARTD(2,I)
      BMATX(3,I2)=CARTD(1,I)
C
   10 CONTINUE
C
      RETURN
C
C***AXIAL-SYMMETRIC CASE
C
    3 CONTINUE
C
      I1=0
      DO 20 I=1,NNODE
      I1=(I-1)*NDOFN+1
      I2=I1+1
      CONSI=SHAPE(I)/GPCOD
C
      BMATX(1,I1)=CARTD(1,I)
      BMATX(2,I2)=CARTD(2,I)
      BMATX(3,I1)=CARTD(2,I)
      BMATX(3,I2)=CARTD(1,I)
      BMATX(4,I1)=CONSI
C
   20 CONTINUE
C
      RETURN
C
C***3-D CASE
C
    4 CONTINUE
C
      I1=0
      DO 40 I=1,NNODE
      I1=(I-1)*NDOFN+1
      I2=I1+1
      I3=I1+2
C
      BMATX(1,I1)=CARTD(1,I)
      BMATX(2,I2)=CARTD(2,I)
      BMATX(3,I1)=CARTD(2,I)
      BMATX(3,I2)=CARTD(1,I)
      BMATX(4,I3)=CARTD(3,I)
      BMATX(5,I1)=CARTD(3,I)
      BMATX(5,I3)=CARTD(1,I)
      BMATX(6,I2)=CARTD(3,I)
      BMATX(6,I3)=CARTD(2,I)
C
   40 CONTINUE
C
      RETURN
C
C***1-D CASE
C
    5 CONTINUE
C
      I1=0
      DO 50 I=1,NNODE
      I1=(I-1)*NDOFN+1
C
      BMATX(1,I1)=CARTD(1,I)
C
   50 CONTINUE
C
      RETURN
      END
