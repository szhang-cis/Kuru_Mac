      SUBROUTINE BDBSEMC(C,D,E,GPCOD,KSYMM,NDIME,NDOFN,NEVAB,
     .                   NNODE,NSTRS,NTYPE,SHAPE)
C**********************************************************************
C                                 T
C**** THIS ROUTINE PERFORMS  E = N  D B  WHEN B IS THE BMATX IS THE
C     CARTESIAN DERIVATIVE OF THE SHAPE FUNCTIONS
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION C(NDIME,*), D(NSTRS,*), E(NEVAB,*), SHAPE(*)
C
C**** INITIALISES E
C
      DO IEVAB=1,NEVAB            ! nevab=nnode
       DO JEVAB=1,NEVAB
        E(IEVAB,JEVAB)=0.0
       ENDDO
      ENDDO
C
      DO INODE=1,NNODE
       INDEX=INODE                ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) INDEX=1     ! UNSYMMETRIC CASE
       DO JNODE=INDEX,NNODE
        DO IDIME=1,NDIME
         E(INODE,JNODE)=E(INODE,JNODE)+
     .                  SHAPE(INODE)*D(IDIME,1)*C(IDIME,JNODE)
        ENDDO
       ENDDO
      ENDDO
C
      RETURN
      END
