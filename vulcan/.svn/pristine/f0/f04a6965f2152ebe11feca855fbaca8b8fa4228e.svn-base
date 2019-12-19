      SUBROUTINE EQHEATT(BMATX,BMSIG,CARTD,DVOLU,ELDIS,GPCOD,
     .                   LARGE,NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTR1,
     .                   NTYPE,SHAPE,SGTOT,XJACM,SIGMA,elelt,
     .                   WARTD)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE EQUIVALENT NODAL "HEATS FORCES" DUE TO
C     THE EQUILIBRATED HEATS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION BMATX(NSTR1,*), BMSIG(*), CARTD(NDIME,*),
     .          ELDIS(NDOFC,*), SHAPE(*), SGTOT(*),
     .          XJACM(NDIME,*), SIGMA(*), elelt(*)
      DIMENSION WARTD(NDIME,*)
C
      DO 10 IDIME=1,NDIME
      DO 10 IEVAB=1,NNODE
      BMSIG(IEVAB)=BMSIG(IEVAB)+WARTD(IDIME,IEVAB)*sgtot(IDIME)*DVOLU
      elelt(ievab)=elelt(ievab)+wartd(idime,ievab)*sgtot(idime)*dvolu
   10 CONTINUE
C
      RETURN
      END
