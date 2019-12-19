      SUBROUTINE CENTROI(ELCODT,NDIMET,NNODLT,CENTRO)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE COORDINATES OF THE ELEMENTAL CENTROID
C     MEASURED WITH REGARD TO THE NODES
C     ( ELEMENT NO. 5 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ELCODT(NDIMET,*), CENTRO(NDIMET,*)
C
      GO TO (1,2,3),NDIMET
C
    1 CONTINUE
      IF(NNODLT.EQ.2) THEN
       CENTRO(1,1)=(ELCODT(1,2)-ELCODT(1,1))/2.0
       CENTRO(1,2)=-CENTRO(1,1)
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.3) THEN
       call runendt('centroi: 1D 3-noded not implemented yet')
       RETURN
      ENDIF
C
    2 CONTINUE
      IF(NNODLT.EQ.3) THEN
       DO IDIMET=1,NDIMET
        CENTRO(IDIMET,1)=(-2.0*ELCODT(IDIMET,1)+ELCODT(IDIMET,2)+
     .                   ELCODT(IDIMET,3))/3.0
        CENTRO(IDIMET,2)=(ELCODT(IDIMET,1)-2.0*ELCODT(IDIMET,2)+
     .                   ELCODT(IDIMET,3))/3.0
        CENTRO(IDIMET,3)=(ELCODT(IDIMET,1)+ELCODT(IDIMET,2)-
     .                   2.0*ELCODT(IDIMET,3))/3.0
       ENDDO
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.4) THEN
       DO IDIMET=1,NDIMET
        CENTRO(IDIMET,1)=(-3.0*ELCODT(IDIMET,1)+ELCODT(IDIMET,2)+
     .                   ELCODT(IDIMET,3)+ELCODT(IDIMET,4))/4.0
        CENTRO(IDIMET,2)=(ELCODT(IDIMET,1)-3.0*ELCODT(IDIMET,2)+
     .                   ELCODT(IDIMET,3)+ELCODT(IDIMET,4))/4.0
        CENTRO(IDIMET,3)=(ELCODT(IDIMET,1)+ELCODT(IDIMET,2)-
     .                   3.0*ELCODT(IDIMET,3)+ELCODT(IDIMET,4))/4.0
        CENTRO(IDIMET,4)=(ELCODT(IDIMET,1)+ELCODT(IDIMET,2)+
     .                   ELCODT(IDIMET,3)-3.0*ELCODT(IDIMET,4))/4.0
       ENDDO
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.6) THEN
       call runendt('centroi: 2D 6-noded not implemented yet')
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.8) THEN
       call runendt('centroi: 2D 8-noded not implemented yet')
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.9) THEN
       call runendt('centroi: 2D 9-noded not implemented yet')
       RETURN
      ENDIF
C
    3 CONTINUE
      call runendt('centroi: 3D not implemented yet')
C
      RETURN
C
      END
