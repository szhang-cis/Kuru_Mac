      SUBROUTINE RESTVM(BUFFE,IFLAG,ARRAY,LENGH,NRCUR)
C***********************************************************************
C
C**** THIS ROUTINE WRITES (IFLAG=1) OR READS (IFLAG=2)
C     ONE ARRAY ('ARRAY' WHICH DIMENSION IS 'LENGH') TO OR FROM
C     'NRCUR'   CURRENT POSITION POINTER OF 'BUFFE' ARRAY
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ARRAY(LENGH), BUFFE(*)
C
      MRCUR=NRCUR-1
C
      IF(IFLAG.EQ.1) THEN
        DO 10 ILENG=1,LENGH
          BUFFE(MRCUR+ILENG)=ARRAY(ILENG)
   10   CONTINUE
      ELSE
        DO 20 ILENG=1,LENGH
          ARRAY(ILENG)=BUFFE(MRCUR+ILENG)
   20   CONTINUE
       ENDIF
C
      RETURN
      END
