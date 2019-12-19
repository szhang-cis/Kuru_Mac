      SUBROUTINE KMATRI(BMATX,CARTD,DMATX,DVOLU,ESTIF,GPCOD,LARGE,
     .                  KSYMM,NDIME,NDOFN,NEVAB,NKOVA,NNODE,NSTRE,
     .                  NSTRS,NTYPE,SHAPE,NSTR1,NDOFC,ELDIS,EMATX,
     .                  XJACM)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE MATERIAL STIFFNESS MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION BMATX(NSTRE,*), CARTD(NDIME,*), DMATX(NSTRS,*),
     .          ESTIF(*),       GPCOD(*),       SHAPE(*)
      DIMENSION ELDIS(NDOFC,*), XJACM(NDIME,*), EMATX(NEVAB,*)
C
      NEVAL=NNODE*NDOFN      ! local NEVAB
C
C                   T
C**** PERFORME THE B   D  B PRODUCT
C
      IF(LARGE.EQ.0) THEN
       CALL BDBSEM(CARTD,DMATX,EMATX,GPCOD,KSYMM,NDIME,NDOFN,
     .             NEVAB,NNODE,NSTRS,NTYPE,SHAPE)
      ELSE
       CALL BLARGE(BMATX,CARTD,ELDIS,GPCOD,
     .             LARGE,NDIME,NDOFC,NNODE,NSTRE,NTYPE,SHAPE,
     .             XJACM)
C
       CALL BDBCO1(BMATX,DMATX,EMATX,KSYMM,NEVAB,NSTRE,NSTRS,NEVAL,
     .             NSTRE)
      ENDIF
C
C**** ADD THE CONTRIBUTION TO THE STIFFNESS MATRIX
C
      CALL SQTOTA(ESTIF,EMATX,DVOLU,NEVAB,KSYMM)
C
      RETURN
      END
