      SUBROUTINE RMATR2(NDIME,NSTRS,RMAT1,RMAT2)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE TENSOR TRASFORMATION MATRIX 
C     RMAT2(NSTRS,NSTRS) FROM THE GLOBAL CARTESIAN REFERENCE 
C     SYSTEM TO A LOCAL CARTESIAN REFERENCE SYSTEM OF THE CONSTITUTIVE 
C     TENSOR.
C  
C                      ,        T
C                 DMATX  = RMAT2 * DMATX *RMAT2 ;
C
C          2        2                   2
C        L1       M1        L1M1      N1        L1N1        M1N1
C          2        2                   2
C        L2       M2        L2M2      N2        L2N2        M2N2
C
C    2*L1L2   2*M1M2  (L1M2+L2M1) 2*N1N2  (L1N2+L2N1)  (M1N2+M2N1)
C          2        2                  2
C        L3       M3        L3M3     N3        L3N3        M3N3
C
C    2*L1L3   2*M1M3  (L1M3+L3M1) 2*N1N3  (L1N3+L3N1)  (M1N3+M3N1)
C
C    2*L2L3   2*M2M3  (L2M3+L3M2) 2*N2N3  (L2N3+L3N2)  (M2N3+M3N2)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION RMAT1(NDIME,*), RMAT2(NSTRS,*)
C
C**** CASE OF NSTRS=3
C
      AL1=RMAT1(1,1)
      AM1=RMAT1(1,2)
      AL2=RMAT1(2,1)
      AM2=RMAT1(2,2)
C
      RMAT2(1,1)=   AL1*AL1
      RMAT2(1,2)=   AM1*AM1
      RMAT2(1,3)=   AL1*AM1
C
      RMAT2(2,1)=   AL2*AL2
      RMAT2(2,2)=   AM2*AM2
      RMAT2(2,3)=   AL2*AM2
C
      RMAT2(3,1)=2.*AL1*AL2
      RMAT2(3,2)=2.*AM1*AM2
      RMAT2(3,3)=AL1*AM2+AL2*AM1
C
C**** CASE OF NSTRS=4.
C
      IF(NSTRS.EQ.4)THEN
       RMAT2(4,1)=0.
       RMAT2(4,2)=0.
       RMAT2(4,3)=0.
       RMAT2(4,4)=1.
       RMAT2(1,4)=0.
       RMAT2(2,4)=0.
       RMAT2(3,4)=0.
      END IF
C
C**** 3D-CASE
C
      IF(NSTRS.EQ.6)THEN
C
      AN1=RMAT1(1,3)
      AN2=RMAT1(2,3)
      AL3=RMAT1(3,1)
      AM3=RMAT1(3,2)
      AN3=RMAT1(3,3)
C
      RMAT2(1,4)=   AN1*AN1
      RMAT2(1,5)=   AL1*AN1 
      RMAT2(1,6)=   AM1*AN1
C
      RMAT2(2,4)=   AN2*AN2
      RMAT2(2,5)=   AL2*AN2 
      RMAT2(2,6)=   AM2*AN2
C
      RMAT2(3,4)=2.*AN1*AN2
      RMAT2(3,5)=AL1*AN2+AL2*AN1 
      RMAT2(3,6)=AM1*AN2+AM2*AN1
C
      RMAT2(4,1)=   AL3*AL3
      RMAT2(4,2)=   AM3*AM3
      RMAT2(4,3)=   AL3*AM3
      RMAT2(4,4)=   AN3*AN3
      RMAT2(4,5)=   AL3*AN3 
      RMAT2(4,6)=   AM3*AN3
C
      RMAT2(5,1)=2.*AL1*AL3
      RMAT2(5,2)=2.*AM1*AM3
      RMAT2(5,3)=AL1*AM3+AL3*AM1
      RMAT2(5,4)=2.*AN1*AN3
      RMAT2(5,5)=AL1*AN3+AL3*AN1 
      RMAT2(5,6)=AM1*AN3+AM3*AN1
C
      RMAT2(6,1)=2.*AL2*AL3
      RMAT2(6,2)=2.*AM2*AM3
      RMAT2(6,3)=AL2*AM3+AL3*AM2
      RMAT2(6,4)=2.*AN2*AN3
      RMAT2(6,5)=AL2*AN3+AL3*AN2 
      RMAT2(6,6)=AM2*AN3+AM3*AN2
C
      END IF
C
      RETURN
      END
