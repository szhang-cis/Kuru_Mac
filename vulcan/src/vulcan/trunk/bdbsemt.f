      SUBROUTINE BDBSEMT(C,D,E,GPCOD,KSYMM,NDIME,NDOFN,NEVAB,
     .                   NNODE,NSTRS,NTYPE,SHAPE,W)
C***********************************************************************
C                                 T
C**** THIS ROUTINE PERFORMS  E = W  D B  WHEN W & B ARE THE CARTESIAN
C     DERIVATIVE OF THE WEIGHTING & SHAPE FUNCTIONS, RESPECTIVELY
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION C(NDIME,*), D(NSTRS,*), E(NEVAB,*), SHAPE(*)
      DIMENSION W(NDIME,*)
C
      DO INODE=1,NNODE
       INDEX=INODE                ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) INDEX=1     ! UNSYMMETRIC CASE
       DO JNODE=INDEX,NNODE
        E(INODE,JNODE)=0.0D0
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          E(INODE,JNODE)=E(INODE,JNODE)+
     .                   W(IDIME,INODE)*D(IDIME,JDIME)*C(JDIME,JNODE)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C
      RETURN
      END
