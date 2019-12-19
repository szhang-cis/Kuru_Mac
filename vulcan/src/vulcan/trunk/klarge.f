      SUBROUTINE KLARGE(CARTD,DVOLU,ESTIF,GPCOD,KSYMM,NDIME,NDOFN,
     .                  NEVAB,NNODE,NTYPE,SHAPE,SIGMA)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE ADDITIONAL STIFFNESS TERMS DUE TO
C     THE FINITE LAGRANGE STRAIN FORMULATION
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION CARTD(NDIME,*), ESTIF(*), SHAPE(*), SIGMA(*)
      DIMENSION GPCOD(*)
C
      GOTO (1,1,3,4,5),NTYPE
C
C**** PLANE STRAIN & STRESS
C
    1 DO INODE=1,NNODE
       KEVAB=(INODE-1)*NDOFN
       DO JNODE=INODE,NNODE
        JEVAB=(JNODE-1)*NDOFN
        ACOEF=SIGMA(1)*CARTD(1,JNODE)+SIGMA(3)*CARTD(2,JNODE)
        BCOEF=SIGMA(3)*CARTD(1,JNODE)+SIGMA(2)*CARTD(2,JNODE)
        CONST=(ACOEF*CARTD(1,INODE)+BCOEF*CARTD(2,INODE))*DVOLU
        DO IDOFN=1,NDOFN
         IEVAB=KEVAB+IDOFN
         JEVAB=JEVAB+1
         KLOCA=(2*NEVAB-IEVAB)*(IEVAB-1)/2+JEVAB
         IF(KSYMM.EQ.0) KLOCA=(IEVAB-1)*NEVAB+JEVAB  ! unsymmetric
         ESTIF(KLOCA)=ESTIF(KLOCA)+CONST
        ENDDO
       ENDDO
      ENDDO
C
      GO TO 100
C
C**** AXI-SYMMETRIC
C
    3 GPCO2=GPCOD(1)*GPCOD(1)
      DO INODE=1,NNODE
       KEVAB=(INODE-1)*NDOFN
       DO JNODE=INODE,NNODE
        JEVAB=(JNODE-1)*NDOFN
        ACOEF=SIGMA(1)*CARTD(1,JNODE)+SIGMA(3)*CARTD(2,JNODE)
        BCOEF=SIGMA(3)*CARTD(1,JNODE)+SIGMA(2)*CARTD(2,JNODE)
        CONST=(ACOEF*CARTD(1,INODE)+BCOEF*CARTD(2,INODE))*DVOLU
        CCOEF=SIGMA(4)*SHAPE(INODE)*SHAPE(JNODE)/GPCO2
        DO IDOFN=1,NDOFN
         IEVAB=KEVAB+IDOFN
         JEVAB=JEVAB+1
         KLOCA=(2*NEVAB-IEVAB)*(IEVAB-1)/2+JEVAB
         IF(KSYMM.EQ.0) KLOCA=(IEVAB-1)*NEVAB+JEVAB
         DCOEF=0.
         IF(IDOFN.EQ.1) DCOEF=CCOEF*DVOLU
         ESTIF(KLOCA)=ESTIF(KLOCA)+CONST+DCOEF 
        ENDDO
       ENDDO
      ENDDO
C
      GO TO 100
C
C**** 3-D CASE
C
    4 DO INODE=1,NNODE
       KEVAB=(INODE-1)*NDOFN
       DO JNODE=INODE,NNODE
        JEVAB=(JNODE-1)*NDOFN
        ACOEF=SIGMA(1)*CARTD(1,JNODE)+SIGMA(3)*CARTD(2,JNODE)+
     .        SIGMA(5)*CARTD(3,JNODE)
        BCOEF=SIGMA(3)*CARTD(1,JNODE)+SIGMA(2)*CARTD(2,JNODE)+
     .        SIGMA(6)*CARTD(3,JNODE)
        CCOEF=SIGMA(5)*CARTD(1,JNODE)+SIGMA(6)*CARTD(2,JNODE)+
     .        SIGMA(4)*CARTD(3,JNODE)
        CONST=(ACOEF*CARTD(1,INODE)+
     .         BCOEF*CARTD(2,INODE)+CCOEF*CARTD(3,INODE))*DVOLU
        DO IDOFN=1,NDOFN
         IEVAB=KEVAB+IDOFN
         JEVAB=JEVAB+1
         KLOCA=(2*NEVAB-IEVAB)*(IEVAB-1)/2+JEVAB
         IF(KSYMM.EQ.0) KLOCA=(IEVAB-1)*NEVAB+JEVAB
         ESTIF(KLOCA)=ESTIF(KLOCA)+CONST
        ENDDO
       ENDDO
      ENDDO
C
      GO TO 100
C
C**** 1-D CASE
C
    5 DO INODE=1,NNODE
       KEVAB=(INODE-1)*NDOFN
       DO JNODE=INODE,NNODE
        JEVAB=(JNODE-1)*NDOFN
        CONST=SIGMA(1)*CARTD(1,JNODE)*CARTD(1,INODE)*DVOLU
        DO IDOFN=1,NDOFN
         IEVAB=KEVAB+IDOFN
         JEVAB=JEVAB+1
         KLOCA=(2*NEVAB-IEVAB)*(IEVAB-1)/2+JEVAB
         IF(KSYMM.EQ.0) KLOCA=(IEVAB-1)*NEVAB+JEVAB
         ESTIF(KLOCA)=ESTIF(KLOCA)+CONST
        ENDDO
       ENDDO
      ENDDO
C
C**** COMPLETE MATRIX FOR RECTANGULAR CASE
C
  100 IF(KSYMM.EQ.0) THEN
       DO IEVAB=1,NEVAB
        DO JEVAB=IEVAB,NEVAB
         KLOCS=(IEVAB-1)*NEVAB+JEVAB
         KLOCI=(JEVAB-1)*NEVAB+IEVAB
         ESTIF(KLOCI)=ESTIF(KLOCS)
        END DO
       END DO
      END IF
C
      RETURN
      END
