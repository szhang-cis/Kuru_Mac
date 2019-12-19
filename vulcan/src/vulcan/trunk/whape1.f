      SUBROUTINE WHAPE1(DERIV,S,NDIME,NNODE,SHAPE,HACHE,
     .                  XJACM,IPERT,XJACI,CARTT,SHAPT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES PERTURBATIONS FUNCTIONS AND THEIR
C     DERIVATIVES FOR LINEAR AND QUADRATIC ISOPARAMETRIC 1-D ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DERIV(NDIME,*), SHAPE(*), HACHE(*), XJACM(NDIME,*)
C
      DIMENSION XJACI(NDIME,*), CARTT(NDIME,*), SHAPT(NDIME,*)
C
C**** 2 NODED ELEMENT
C
      IF(NNODE.EQ.2) THEN
       GO TO (1,2,3) IPERT
C
    1  SHAPE(1)=(1.-S*S)*HACHE(1)
       SHAPE(2)=(1.-S*S)*HACHE(2)
C
       DERIV(1,1)=-2.0*S*HACHE(1)
       DERIV(1,2)=-2.0*S*HACHE(2)
       RETURN
C
    2  SHAPE(1)=(1.-S*S)*HACHE(1)
       SHAPE(2)=(1.-S*S)*HACHE(2)
C
       DERIV(1,1)=0.0
       DERIV(1,2)=0.0
       RETURN
C
    3  CONTINUE
C
C**** CALCULATE THE DETERMINANT AND INVERS OF JACOBIAN MATRIX
C
       CALL INVMTX(XJACM,XJACI,DETJM,NDIME)
C
       SHAPT(1,1)=-0.5
       SHAPT(1,2)= 0.5
C
C**** CALCULATE THE CARTESIAN DERIVATIVES CORRESPONDING TO "SHAPT"
C
C      Note: instead of PROMA1 & PROMA3, it is possible to use:
C
C       DO IDIME=1,NDIME
C        DO INODE=1,NNODE
C         CARTT(IDIME,INODE)=0.0
C         DO JDIME=1,NDIME
C          CARTT(IDIME,INODE)=CARTT(IDIME,INODE)+
C     .                       XJACI(IDIME,JDIME)*SHAPT(JDIME,INODE)
C         ENDDO
C        ENDDO
C       ENDDO
C
C       DO INODE=1,NNODE
C        SHAPE(INODE)=0.0
C        DO IDIME=1,NDIME
C         SHAPE(INODE)=SHAPE(INODE)+CARTT(IDIME,INODE)*HACHE(IDIME)
C        ENDDO
C       ENDDO
C
       CALL PROMA1(CARTT,XJACI,SHAPT,NDIME,NNODE,NDIME)
       CALL PROMA3(SHAPE,CARTT,HACHE,NNODE,NDIME,    1)
C
       DERIV(1,1)=0.0
       DERIV(1,2)=0.0
       RETURN
      ENDIF
C
C**** 3 NODED ELEMENT
C
      call runendt('whape1: 3 noded not implemented')
      IF(NNODE.EQ.3) THEN
        SHAPE(1)=S*(-1.+S)*0.5
        SHAPE(2)=S*( 1.+S)*0.5
        SHAPE(3)=1.-S*S
C
        DERIV(1,1)= -0.5+S
        DERIV(1,2)=  0.5+S
        DERIV(1,3)= -2.0*S
C
        RETURN
      ENDIF
C
      END
