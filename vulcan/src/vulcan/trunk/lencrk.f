      SUBROUTINE LENCRK(ALENG,BETAM,ELCOD,DVOLU,NCRAK,NDIME,NNODE,
     .                  NTYPE,THICK,IELEM)
C*******************************************************************
C
C****THIS ROUTINE COMPUTES THE CHARACTERISTIC LENGTHS
C
C*******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ALENG(NDIME), BETAM(NDIME,NDIME), ELCOD(NDIME,NNODE)
      DIMENSION V(3),     VM(3),      ALM(3)
      DATA TWOPI/6.283185307179586/
C
CC$DIR SCALAR
      DO I=1,NDIME
        ALM(I)=0.0
      ENDDO
C
C***FIND LONGEST ARRAY IN ELEMENT
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
C..............PROYECT V ON CRACK DIRECTIONS
C
CC$DIR SCALAR
      DO ICRAK=1,NDIME
        VM(ICRAK)=0.
CC$DIR SCALAR
        DO JDIME=1,NDIME
          VM(ICRAK)=VM(ICRAK)+BETAM(ICRAK,JDIME)*V(JDIME)  
        END DO
      END DO
C
C..............LOOK FOR MAXIMUM LENGHTS
C
CC$DIR SCALAR
      DO IDIME=1,NDIME
        IF(ABS(VM(IDIME)).GT.ALM(IDIME)) ALM(IDIME)=ABS(VM(IDIME))
      ENDDO
C
 100  CONTINUE
C
C***COMPUTE VOLUME OF THE ELEMENT
C
      TVOLU=1.0
CC$DIR SCALAR
      DO I=1,NDIME
        TVOLU=TVOLU*ALM(I)
      ENDDO
C
C***CORRECT FOR PLANE STRESS
C
      IF(NTYPE.NE.4) TVOLU=TVOLU*THICK
C
C***CORRECT FOR AXISYMMETRIC
C
      IF(NTYPE.EQ.3) THEN
        RADUS=0.0D00
        DO INODE=1,NNODE
          RADUS=RADUS+ELCOD(1,NNODE)
        ENDDO
        RADUS=RADUS/FLOAT(NNODE)
        TVOLU=TVOLU*TWOPI*RADUS
      ENDIF
C
C***SCALE TO POINT VOLUME
C
      EXP=1.0D+00/DFLOAT(NDIME)
      FAC=(DVOLU/TVOLU)**EXP
C
CC$DIR SCALAR
      DO I=1,NCRAK
        ALENG(I)=ALM(I)*FAC
      END DO
C
C***CORRECT FOR THIRD CRACK IN 2D
C
      IF(NDIME.EQ.2.AND.NCRAK.EQ.3) THEN
        ALENG(3)=1.0D00
        IF(NTYPE.NE.4) ALENG(3)=THICK
        IF(NTYPE.EQ.3) ALENG(3)=TWOPI*RADUS*THICK
      ENDIF
C
      RETURN
      END
