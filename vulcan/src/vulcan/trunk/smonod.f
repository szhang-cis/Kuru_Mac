      SUBROUTINE SMONOD(NDIME,NNODL,NNODS,NQUTR)
C***********************************************************************
C
C**** THIS ROUTINE DETERMINES THE NUMBER OF LOW ORDER NODES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      NNODS=NNODL
      IF(NDIME.EQ.1) THEN
        IF(NNODL.EQ.3) NNODS=2
      ELSE IF(NDIME.EQ.2) THEN
        IF(NNODL.EQ.4) THEN
          IF(NQUTR.EQ.2) NNODS=3
        ENDIF
        IF(NNODL.EQ.6) NNODS=3
        IF(NNODL.EQ.8) NNODS=4
        IF(NNODL.EQ.9) NNODS=4
      ELSE IF(NDIME.EQ.3) THEN
        IF(NNODL.EQ. 5) NNODS=4
        IF(NNODL.EQ.10) NNODS=4
        IF(NNODL.EQ.20) NNODS=8
        IF(NNODL.EQ.27) NNODS=8
      END IF
C
      RETURN
      END
