      SUBROUTINE SMOMID(ELVAL,NDIME,NNODL,NQUTR)
C***********************************************************************
C
C**** THIS ROUTINE INTERPOLATES VALUES FOR THE MID SIDE NODES
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ELVAL(NNODL)
C
      IF(NDIME.EQ.1) THEN
        IF(NNODL.EQ.3) ELVAL(3)=0.5*(ELVAL(1)+ELVAL(2))
      ELSE IF(NDIME.EQ.2) THEN
        IF(NNODL.EQ.4) THEN
          IF(NQUTR.EQ.2) ELVAL(4)=1.0/3.0*(ELVAL(1)+ELVAL(2)+ELVAL(3))
        ENDIF
        IF(NNODL.EQ.6) THEN
          ELVAL(4)=0.5*(ELVAL(1)+ELVAL(2))
          ELVAL(5)=0.5*(ELVAL(2)+ELVAL(3))
          ELVAL(6)=0.5*(ELVAL(3)+ELVAL(1))
        ELSE IF(NNODL.EQ.8) THEN
          ELVAL(5)=0.5*(ELVAL(1)+ELVAL(2))
          ELVAL(6)=0.5*(ELVAL(2)+ELVAL(3))
          ELVAL(7)=0.5*(ELVAL(3)+ELVAL(4))
          ELVAL(8)=0.5*(ELVAL(4)+ELVAL(1))
        ELSE IF(NNODL.EQ.9) THEN
          ELVAL(5)=0.5*(ELVAL(1)+ELVAL(2))
          ELVAL(6)=0.5*(ELVAL(2)+ELVAL(3))
          ELVAL(7)=0.5*(ELVAL(3)+ELVAL(4))
          ELVAL(8)=0.5*(ELVAL(4)+ELVAL(1))
          ELVAL(9)=0.25*(ELVAL(1)+ELVAL(2)+ELVAL(3)+ELVAL(4))
        END IF
      ELSE IF(NDIME.EQ.3) THEN
        IF(NNODL.EQ.5) THEN
          ELVAL(5)=0.25*(ELVAL(1)+ELVAL(2)+ELVAL(3)+ELVAL(4))
        ELSE IF(NNODL.EQ.10) THEN
          ELVAL( 5)=0.5*(ELVAL(1)+ELVAL(2))
          ELVAL( 6)=0.5*(ELVAL(2)+ELVAL(3))
          ELVAL( 7)=0.5*(ELVAL(3)+ELVAL(1))
          ELVAL( 8)=0.5*(ELVAL(4)+ELVAL(1))
          ELVAL( 9)=0.5*(ELVAL(4)+ELVAL(2))
          ELVAL(10)=0.5*(ELVAL(4)+ELVAL(3))
        ELSE IF(NNODL.EQ.20) THEN
          ELVAL( 9)=0.5*(ELVAL(1)+ELVAL(2))
          ELVAL(10)=0.5*(ELVAL(2)+ELVAL(3))
          ELVAL(11)=0.5*(ELVAL(3)+ELVAL(4))
          ELVAL(12)=0.5*(ELVAL(4)+ELVAL(1))
          ELVAL(13)=0.5*(ELVAL(5)+ELVAL(6))
          ELVAL(14)=0.5*(ELVAL(6)+ELVAL(7))
          ELVAL(15)=0.5*(ELVAL(7)+ELVAL(8))
          ELVAL(16)=0.5*(ELVAL(8)+ELVAL(5))
          ELVAL(17)=0.5*(ELVAL(1)+ELVAL(5))
          ELVAL(18)=0.5*(ELVAL(2)+ELVAL(6))
          ELVAL(19)=0.5*(ELVAL(3)+ELVAL(7))
          ELVAL(20)=0.5*(ELVAL(4)+ELVAL(8))
        ELSE IF(NNODL.EQ.27) THEN                      ! to be revised !
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
