      SUBROUTINE RMATR1(NDIME,TABPO,RMAT1)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE TENSOR TRASFORMATION MATRIX 
C     RMAT1(NDIME,NDIME) 
C          FROM THE GLOBAL CARTESIAN REFERENCE SYSTEM 
C          TO   THE  LOCAL CARTESIAN REFERENCE SYSTEM
C  
C                                             ! L1  M1  N1 !
C     VECTOR  = RMAT1 * VECTOR        RMAT1 = ! L2  M2  N2 !
C           L                 G               ! L3  M3  N3 !
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION TABPO(*), RMAT1(NDIME,*)
C
C***2-D CASE
C
      IF(NDIME.EQ.2) THEN
        THETA=TABPO(1)
        RMAT1(1,1)= DCOS(THETA)
        RMAT1(1,2)= DSIN(THETA)
        RMAT1(2,1)=-RMAT1(1,2)
        RMAT1(2,2)= RMAT1(1,1)
C
        RETURN
      ENDIF
C
C***3-D CASE
C
      IF(NDIME.EQ.3) THEN
        A1=TABPO(1)
        A2=TABPO(2)
        A3=TABPO(3)
C
        RMAT1(1,1)= DCOS(A3)*DCOS(A2)-DSIN(A3)*DSIN(A2)*DCOS(A1)
        RMAT1(1,2)= DCOS(A3)*DSIN(A2)+DSIN(A3)*DCOS(A2)*DCOS(A1)
        RMAT1(1,3)= DSIN(A1)*DSIN(A3)
        RMAT1(2,1)=-DSIN(A3)*DCOS(A2)-DCOS(A3)*DSIN(A2)*DCOS(A1)
        RMAT1(2,2)=-DSIN(A3)*DSIN(A2)+DCOS(A3)*DCOS(A2)*DCOS(A1)
        RMAT1(2,3)= DSIN(A1)*DCOS(A3)
        RMAT1(3,1)= DSIN(A1)*DSIN(A2)
        RMAT1(3,2)=-DSIN(A1)*DCOS(A2)
        RMAT1(3,3)= DCOS(A1)
C
        RETURN
C
      ENDIF
C
      END
