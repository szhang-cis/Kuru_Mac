      SUBROUTINE SMOESQ(ELVAL,NDIME,NNODL)
C***********************************************************************
C
C****THIS ROUTINE INTERPOLATES VALUES FOR THE CORNER NODES
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ELVAL(NNODL)
C
      IF(NDIME.EQ.1) THEN
        CALL RUNEND('ERROR: SMOESQ > NOT IMPLEMENTED')
        IF(NNODL.EQ.3) ELVAL(3)=0.5*(ELVAL(1)+ELVAL(2))
      ELSE IF(NDIME.EQ.2) THEN
        IF(NNODL.EQ.6) THEN
          CALL RUNEND('ERROR: SMOESQ > NOT IMPLEMENTED')
          ELVAL(4)=0.5*(ELVAL(1)+ELVAL(2))
          ELVAL(5)=0.5*(ELVAL(2)+ELVAL(3))
          ELVAL(6)=0.5*(ELVAL(3)+ELVAL(1))
        ELSE IF(NNODL.EQ.8) THEN
          ELVAL(1)=ELVAL(1)+0.5*(ELVAL(8)+ELVAL(5))
          ELVAL(2)=ELVAL(2)+0.5*(ELVAL(5)+ELVAL(6))
          ELVAL(3)=ELVAL(3)+0.5*(ELVAL(6)+ELVAL(7))
          ELVAL(4)=ELVAL(4)+0.5*(ELVAL(7)+ELVAL(8))
        ELSE IF(NNODL.EQ.9) THEN
          CALL RUNEND('ERROR: SMOESQ > NOT IMPLEMENTED')
          ELVAL(5)=0.5*(ELVAL(1)+ELVAL(2))
          ELVAL(6)=0.5*(ELVAL(2)+ELVAL(3))
          ELVAL(7)=0.5*(ELVAL(3)+ELVAL(4))
          ELVAL(8)=0.5*(ELVAL(4)+ELVAL(1))
          ELVAL(9)=0.25*(ELVAL(1)+ELVAL(2)+ELVAL(3)+ELVAL(4))
        END IF
      ELSE IF(NDIME.EQ.3) THEN
        IF(NNODL.EQ.10) THEN
          CALL RUNEND('ERROR: SMOESQ > NOT IMPLEMENTED')
          ELVAL( 5)=0.5*(ELVAL(1)+ELVAL(2))
          ELVAL( 6)=0.5*(ELVAL(2)+ELVAL(3))
          ELVAL( 7)=0.5*(ELVAL(3)+ELVAL(1))
          ELVAL( 8)=0.5*(ELVAL(4)+ELVAL(1))
          ELVAL( 9)=0.5*(ELVAL(4)+ELVAL(2))
          ELVAL(10)=0.5*(ELVAL(4)+ELVAL(3))
      ELSE IF(NNODL.EQ.20) THEN
          ELVAL(1)=ELVAL(1)+0.5*(ELVAL(12)+ELVAL( 9)+ELVAL(13))
          ELVAL(2)=ELVAL(2)+0.5*(ELVAL( 9)+ELVAL(10)+ELVAL(14))
          ELVAL(3)=ELVAL(3)+0.5*(ELVAL(10)+ELVAL(11)+ELVAL(15))
          ELVAL(4)=ELVAL(4)+0.5*(ELVAL(11)+ELVAL(12)+ELVAL(16))
          ELVAL(5)=ELVAL(5)+0.5*(ELVAL(20)+ELVAL(17)+ELVAL(13))
          ELVAL(6)=ELVAL(6)+0.5*(ELVAL(17)+ELVAL(18)+ELVAL(14))
          ELVAL(7)=ELVAL(7)+0.5*(ELVAL(18)+ELVAL(19)+ELVAL(15))
          ELVAL(8)=ELVAL(8)+0.5*(ELVAL(19)+ELVAL(20)+ELVAL(16))
        ELSE IF(NNODL.EQ.27) THEN
          CALL RUNEND('ERROR: SMOESQ > NOT IMPLEMENTED')
          ELVAL( 9)=0.5*(ELVAL(1)+ELVAL(2))
          ELVAL(10)=0.5*(ELVAL(2)+ELVAL(3))
          ELVAL(11)=0.5*(ELVAL(3)+ELVAL(4))
          ELVAL(12)=0.5*(ELVAL(4)+ELVAL(1))
          ELVAL(13)=0.5*(ELVAL(1)+ELVAL(5))
          ELVAL(14)=0.5*(ELVAL(2)+ELVAL(6))
          ELVAL(15)=0.5*(ELVAL(3)+ELVAL(7))
          ELVAL(16)=0.5*(ELVAL(4)+ELVAL(8))
          ELVAL(17)=0.5*(ELVAL(5)+ELVAL(6))
          ELVAL(18)=0.5*(ELVAL(6)+ELVAL(7))
          ELVAL(19)=0.5*(ELVAL(7)+ELVAL(8))
          ELVAL(20)=0.5*(ELVAL(8)+ELVAL(5))
          ELVAL(21)=0.25*(ELVAL(1)+ELVAL(2)+ELVAL(3)+ELVAL(4))
          ELVAL(22)=0.25*(ELVAL(1)+ELVAL(2)+ELVAL(6)+ELVAL(5))
          ELVAL(23)=0.25*(ELVAL(2)+ELVAL(3)+ELVAL(7)+ELVAL(6))
          ELVAL(24)=0.25*(ELVAL(3)+ELVAL(4)+ELVAL(8)+ELVAL(7))
          ELVAL(25)=0.25*(ELVAL(4)+ELVAL(1)+ELVAL(5)+ELVAL(8))
          ELVAL(26)=0.25*(ELVAL(5)+ELVAL(6)+ELVAL(7)+ELVAL(8))
          ELVAL(27)=0.125*(ELVAL(1)+ELVAL(2)+ELVAL(3)+ELVAL(4)+
     .                     ELVAL(5)+ELVAL(6)+ELVAL(7)+ELVAL(8))
        END IF
      END IF
C
      RETURN
      END
