      SUBROUTINE KMATRIT(BMATX,CARTD,DMATX,DVOLU,ESTIF,GPCOD,LARGE,
     .                   KSYMM,NDIME,NDOFN,NEVAB,NKOVA,NNODE,NSTRE,
     .                   NSTRS,NTYPE,SHAPE,WARTD,EMATX)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE MATERIAL CONDUCTIVITY MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION BMATX(NSTRE,*), CARTD(NDIME,*), DMATX(NSTRS,*),
     .          ESTIF(*),       SHAPE(*)
c     DIMENSION EMATX(20,20)
      DIMENSION WARTD(NDIME,*), EMATX(NEVAB,*)
C
C                    T
C**** PERFORMS THE W  D  B PRODUCT (W=W_X & B=CARTD)
C
      IF(LARGE.EQ.0.OR.LARGE.EQ.3) THEN
       CALL BDBSEMT(CARTD,DMATX,EMATX,GPCOD,KSYMM,NDIME,NDOFN,
     .              NEVAB,NNODE,NSTRS,NTYPE,SHAPE,WARTD)
      ELSE
c      CALL RUNENDT('ERROR IN KMATRIT: LARGE=1,2 NOT IMPLEMENTED')
c      CALL BDBCO1T(BMATX,DMATX,EMATX,KSYMM,NEVAB,NSTRE,NSTRS)
      ENDIF
C
C**** ADD THE CONTRIBUTION TO THE CONDUCTIVITY MATRIX
C
      CALL SQTOTA(ESTIF,EMATX,DVOLU,NEVAB,KSYMM)
C
      RETURN
      END
