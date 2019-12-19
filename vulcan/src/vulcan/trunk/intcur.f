      SUBROUTINE INTCUR(XJACM,DETJM,NDIME,NDIML)
C***********************************************************************
C
C***THIS ROUTINE CALCULATES THE "DETERMINAT" FOR CURVILINEAR INTEGRALS
C
C     XJACM  -  GIVEN NDIME*NDIML JACOBIAN MATRIX FOR BOUNDARY ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XJACM(NDIME,NDIML)
C
      IF(NDIME.EQ.2.AND.NDIML.EQ.1) THEN
        DETJM=DSQRT(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1))
      ELSE IF(NDIME.EQ.3.AND.NDIML.EQ.1) THEN
        DETJM=DSQRT(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)+
     .                                    XJACM(3,1)*XJACM(3,1))
      ELSE ! NDIME=3 & NDIML=2
        DETJM=DSQRT( (XJACM(1,1)*XJACM(2,2)-XJACM(2,1)*XJACM(1,2))*
     .               (XJACM(1,1)*XJACM(2,2)-XJACM(2,1)*XJACM(1,2))+
     .               (XJACM(1,1)*XJACM(3,2)-XJACM(3,1)*XJACM(1,2))*
     .               (XJACM(1,1)*XJACM(3,2)-XJACM(3,1)*XJACM(1,2))+
     .               (XJACM(2,1)*XJACM(3,2)-XJACM(3,1)*XJACM(2,2))*
     .               (XJACM(2,1)*XJACM(3,2)-XJACM(3,1)*XJACM(2,2)) )
      ENDIF
      RETURN
      END
