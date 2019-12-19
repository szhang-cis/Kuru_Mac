      SUBROUTINE BDBSEM(C,D,E,GPCOD,KSYMM,NDIME,NDOFN,NEVAB,
     .                  NNODE,NSTRS,NTYPE,SHAPE)
C***********************************************************************
C                                 T
C**** THIS ROUTINE PERFORMS  E = B D B  WHEN B IS THE BMATX FOR SMALL
C     DISPLACEMENT FORMULATION
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION C(NDIME,*), D(NSTRS,*), E(NEVAB,*), SHAPE(*)
C
      GOTO (1,1,3,4,5),NTYPE
C
C**** PLANE STRESS OR STRAIN
C
    1 DO I=1,NNODE
       I1=(I-1)*NDOFN+1
       I2=I1+1
       IN=I                ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) IN=1 ! UNSYMMETRIC CASE
C$DIR NO_RECURRENCE
       DO J=IN,NNODE
        J1=(J-1)*NDOFN+1
        J2=J1+1
        E(I1,J1)=C(1,I)*(D(1,1)*C(1,J)+D(1,3)*C(2,J))
     .          +C(2,I)*(D(3,1)*C(1,J)+D(3,3)*C(2,J))
        E(I1,J2)=C(1,I)*(D(1,2)*C(2,J)+D(1,3)*C(1,J))
     .          +C(2,I)*(D(3,2)*C(2,J)+D(3,3)*C(1,J))
        E(I2,J1)=C(2,I)*(D(2,1)*C(1,J)+D(2,3)*C(2,J))
     .          +C(1,I)*(D(3,1)*C(1,J)+D(3,3)*C(2,J))
        E(I2,J2)=C(2,I)*(D(2,2)*C(2,J)+D(2,3)*C(1,J))
     .          +C(1,I)*(D(3,2)*C(2,J)+D(3,3)*C(1,J))
       END DO
      END DO
C
      RETURN
C
C**** AXIAL-SYMMETRIC CASE
C
    3 DO I=1,NNODE
       I1=(I-1)*NDOFN+1
       I2=I1+1
       CONSI=SHAPE(I)/GPCOD
       IN=I                ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) IN=1 ! UNSYMMETRIC CASE
C$DIR NO_RECURRENCE
       DO J=IN,NNODE
        J1=(J-1)*NDOFN+1
        J2=J1+1
        CONSJ=SHAPE(J)/GPCOD
        E(I1,J1)=C(1,I)*(D(1,1)*C(1,J)+D(1,3)*C(2,J)+D(1,4)*CONSJ)
     .          +C(2,I)*(D(3,1)*C(1,J)+D(3,3)*C(2,J)+D(3,4)*CONSJ)
     .          +CONSI*(D(4,1)*C(1,J)+D(4,3)*C(2,J)+D(4,4)*CONSJ)
        E(I1,J2)=C(1,I)*(D(1,2)*C(2,J)+D(1,3)*C(1,J))
     .          +C(2,I)*(D(3,2)*C(2,J)+D(3,3)*C(1,J))
     .          +CONSI*(D(4,2)*C(2,J)+D(4,3)*C(1,J))
        E(I2,J1)=C(2,I)*(D(2,1)*C(1,J)+D(2,3)*C(2,J)+D(2,4)*CONSJ)
     .          +C(1,I)*(D(3,1)*C(1,J)+D(3,3)*C(2,J)+D(3,4)*CONSJ)
        E(I2,J2)=C(2,I)*(D(2,2)*C(2,J)+D(2,3)*C(1,J))
     .          +C(1,I)*(D(3,2)*C(2,J)+D(3,3)*C(1,J))
       END DO
      END DO
C
      RETURN
C
C**** 3-D CASE
C
    4 DO I=1,NNODE
       I1=(I-1)*NDOFN+1
       I2=I1+1
       I3=I1+2
       IN=I                ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) IN=1 ! UNSYMMETRIC CASE
C$DIR NO_RECURRENCE
       DO J=IN,NNODE
        J1=(J-1)*NDOFN+1
        J2=J1+1
        J3=J1+2
        E(I1,J1)=C(1,I)*(D(1,1)*C(1,J)+D(1,3)*C(2,J)+D(1,5)*C(3,J))
     .          +C(2,I)*(D(3,1)*C(1,J)+D(3,3)*C(2,J)+D(3,5)*C(3,J))
     .          +C(3,I)*(D(5,1)*C(1,J)+D(5,3)*C(2,J)+D(5,5)*C(3,J))
        E(I1,J2)=C(1,I)*(D(1,2)*C(2,J)+D(1,3)*C(1,J)+D(1,6)*C(3,J))
     .          +C(2,I)*(D(3,2)*C(2,J)+D(3,3)*C(1,J)+D(3,6)*C(3,J))
     .          +C(3,I)*(D(5,2)*C(2,J)+D(5,3)*C(1,J)+D(5,6)*C(3,J))
        E(I1,J3)=C(1,I)*(D(1,4)*C(3,J)+D(1,5)*C(1,J)+D(1,6)*C(2,J))
     .          +C(2,I)*(D(3,4)*C(3,J)+D(3,5)*C(1,J)+D(3,6)*C(2,J))
     .          +C(3,I)*(D(5,4)*C(3,J)+D(5,5)*C(1,J)+D(5,6)*C(2,J))
        E(I2,J1)=C(2,I)*(D(2,1)*C(1,J)+D(2,3)*C(2,J)+D(2,5)*C(3,J))
     .          +C(1,I)*(D(3,1)*C(1,J)+D(3,3)*C(2,J)+D(3,5)*C(3,J))
     .          +C(3,I)*(D(6,1)*C(1,J)+D(6,3)*C(2,J)+D(6,5)*C(3,J))
        E(I2,J2)=C(2,I)*(D(2,2)*C(2,J)+D(2,3)*C(1,J)+D(2,6)*C(3,J))
     .          +C(1,I)*(D(3,2)*C(2,J)+D(3,3)*C(1,J)+D(3,6)*C(3,J))
     .          +C(3,I)*(D(6,2)*C(2,J)+D(6,3)*C(1,J)+D(6,6)*C(3,J))
        E(I2,J3)=C(2,I)*(D(2,4)*C(3,J)+D(2,5)*C(1,J)+D(2,6)*C(2,J))
     .          +C(1,I)*(D(3,4)*C(3,J)+D(3,5)*C(1,J)+D(3,6)*C(2,J))
     .          +C(3,I)*(D(6,4)*C(3,J)+D(6,5)*C(1,J)+D(6,6)*C(2,J))
        E(I3,J1)=C(3,I)*(D(4,1)*C(1,J)+D(4,3)*C(2,J)+D(4,5)*C(3,J))
     .          +C(1,I)*(D(5,1)*C(1,J)+D(5,3)*C(2,J)+D(5,5)*C(3,J))
     .          +C(2,I)*(D(6,1)*C(1,J)+D(6,3)*C(2,J)+D(6,5)*C(3,J))
        E(I3,J2)=C(3,I)*(D(4,2)*C(2,J)+D(4,3)*C(1,J)+D(4,6)*C(3,J))
     .          +C(1,I)*(D(5,2)*C(2,J)+D(5,3)*C(1,J)+D(5,6)*C(3,J))
     .          +C(2,I)*(D(6,2)*C(2,J)+D(6,3)*C(1,J)+D(6,6)*C(3,J))
        E(I3,J3)=C(3,I)*(D(4,4)*C(3,J)+D(4,5)*C(1,J)+D(4,6)*C(2,J))
     .          +C(1,I)*(D(5,4)*C(3,J)+D(5,5)*C(1,J)+D(5,6)*C(2,J))
     .          +C(2,I)*(D(6,4)*C(3,J)+D(6,5)*C(1,J)+D(6,6)*C(2,J))
       END DO
      END DO
C
      RETURN
C
C**** 1-D CASE
C
    5 DO I=1,NNODE
       I1=(I-1)*NDOFN+1
       IN=I                ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) IN=1 ! UNSYMMETRIC CASE
C$DIR NO_RECURRENCE
       DO J=IN,NNODE
        J1=(J-1)*NDOFN+1
        E(I1,J1)=C(1,I)*D(1,1)*C(1,J)
       END DO
      END DO
C
      RETURN
C
      END
