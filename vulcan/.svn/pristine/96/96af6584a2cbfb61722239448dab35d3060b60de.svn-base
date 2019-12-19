      SUBROUTINE EQLOAD(BMATX,BMSIG,CARTD,DVOLU,ELDIS,GPCOD,
     .                  LARGE,NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .                  NTYPE,SHAPE,SGTOT,
     .                  XJACM)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE EQUIVALENT NODAL FORCES DUE TO
C     THE EQUILIBRATED STRESS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION BMATX(NSTRE,*), BMSIG(*), CARTD(NDIME,*), 
     .          ELDIS(NDOFC,*), SHAPE(*), GPCOD(*),       SGTOT(*),
     .          XJACM(NDIME,*)
C
      NEVAL=NNODE*NDOFN      ! local NEVAB
C
      IF(LARGE.EQ.0) THEN
       CALL BMASIG(BMSIG,CARTD,DVOLU,GPCOD,NDIME,NDOFN,NNODE,
     .             NTYPE,SHAPE,SGTOT)
      ELSE
       CALL BLARGE(BMATX,CARTD,ELDIS,GPCOD,LARGE,NDIME,NDOFC,
     .             NNODE,NSTRE,NTYPE,SHAPE,XJACM)
C
       DO 10 ISTRE=1,NSTRE
       DO 10 IEVAB=1,NEVAL
       BMSIG(IEVAB)=BMSIG(IEVAB)+BMATX(ISTRE,IEVAB)*SGTOT(ISTRE)*DVOLU
   10  CONTINUE
      ENDIF
C
      RETURN
      END
