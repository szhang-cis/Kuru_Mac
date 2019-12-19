      SUBROUTINE LENJNT(AMAXL,ELCOD,RMAT1,NDIME,NNODE)
C*******************************************************************
C
C****THIS ROUTINE COMPUTES THE CHARACTERISTIC JOINT LENGTH
C
C*******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RMAT1(NDIME,*), ELCOD(NDIME,*)
      DIMENSION V(3)
C
C***FIND LONGEST ARRAY IN ELEMENT
C
      AMAXL=0.0
C
      DO 100 INODE=      1,NNODE
      DO 100 JNODE=INODE+1,NNODE
C
C..............CONNECTING VECTOR
C
CC$DIR SCALAR
        DO I=1,NDIME
          V(I)=ELCOD(I,INODE)-ELCOD(I,JNODE)
        END DO
C
C..............PROYECT V ON JOINT DIRECTION
C
      ALENG=0.
CC$DIR SCALAR
      DO JDIME=1,NDIME
        ALENG=ALENG+RMAT1(1,JDIME)*V(JDIME)  
      END DO
      ALENG=ABS(ALENG)
C
C..............LOOK FOR MAXIMUM
C
      IF(ALENG.GT.AMAXL) AMAXL=ALENG
C
 100  CONTINUE
C
      RETURN
      END
