      SUBROUTINE RULIRO(NGAUS,POSGP,WEIGP)
C***********************************************************************
C
C****THIS ROUTINE SETS UP THE IRONS' RULES INTEGRATION CONSTANTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION POSGP(3,*), WEIGP(*)
C
C***DEFINE SAMPLING WEIGHTS
C
      IF(NGAUS.EQ.6) THEN
        A1=8.0/6.0
        DO 10 IGAUS=1,NGAUS
   10   WEIGP(IGAUS)=A1
      ENDIF
      IF(NGAUS.EQ.14) THEN
        A1=0.886426592797784
        A2=0.335180055401662
        DO 20 IGAUS=1,6
   20   WEIGP(IGAUS)=A1
        DO 30 IGAUS=7,NGAUS
   30   WEIGP(IGAUS)=A2
      ENDIF
      IF(NGAUS.EQ.15) THEN
        A1=1.56444444444444444
        A2=0.35555555555555556
        A3=0.53777777777777778
        WEIGP(1)=A1
        DO 40 IGAUS=2,7
   40   WEIGP(IGAUS)=A2
        DO 50 IGAUS=8,NGAUS
   50   WEIGP(IGAUS)=A3
      ENDIF
C
C***DEFINE SAMPLING POSITIONS
C
      DO 55 IGAUS=1,NGAUS
      DO 55 IDIME=1,3
   55 POSGP(IDIME,IGAUS)=0.0
C
      IF(NGAUS.EQ.6) THEN
        B=1.0
        IGAUS=1
        DO 60 IDIME=1,3
        POSGP(IDIME,IGAUS  )= B
        POSGP(IDIME,IGAUS+1)=-B
        IGAUS=IGAUS+2
   60 CONTINUE
      ENDIF
      IF(NGAUS.EQ.14.OR.NGAUS.EQ.15) THEN
        IF(NGAUS.EQ.14) THEN
          B=0.795822426
          C=0.758786911
          IGAUS=1
        ENDIF
        IF(NGAUS.EQ.15) THEN
          B=1.0
          C=(1.0/(9.0*A3))**0.25
          IGAUS=2
        ENDIF
        DO 70 IDIME=1,3
        POSGP(IDIME,IGAUS  )= B
        POSGP(IDIME,IGAUS+1)=-B
        IGAUS=IGAUS+2
   70   CONTINUE
        IF(NGAUS.EQ.14) IGAUS=7
        IF(NGAUS.EQ.15) IGAUS=8
        C1=C
        DO 80 J=1,2
        C2=C
        DO 90 K=1,2
        C3=C
        DO 95 L=1,2
        POSGP(1,IGAUS)=C1
        POSGP(2,IGAUS)=C2
        POSGP(3,IGAUS)=C3
        C3=-C
   95   IGAUS=IGAUS+1
        C2=-C
   90   CONTINUE
        C1=-C
   80   CONTINUE
      ENDIF
      RETURN
      END
