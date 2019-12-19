      SUBROUTINE HMATRI(BMATX,CARTD,DVOLU,GPCOD,HSTIF,IFLAG,NDIME,
     .                  NEVAB,NNODE,NSTRE,NTYPE,PROPS,SHAPE)
C**********************************************************************
C
C****THIS ROUTINE EVALUATES THE ELEMENT SOIL-FLUID COUPLED MATRIX
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BMATX(NEVAB,*), CARTD(NDIME,*), DELTA(4), 
     .          HSTIF(NEVAB,*), PROPS(*),       SHAPE(*)
      DATA DELTA/1.,1.,0.,1./
C
      CONST=PROPS(19)
      IF(IFLAG.EQ.0)THEN
C
C***SMALL DEFORMATIONS
C
        IEVAB=0
        DO 10 INODE=1,NNODE
        DO 10 IDIME=1,NDIME
        IEVAB=IEVAB+1
C
        FMULT=CARTD(IDIME,INODE)
        IF(NTYPE.EQ.3.AND.IDIME.EQ.1)FMULT=FMULT+SHAPE(INODE)/GPCOD
C
        DO 10 JNODE=1,NNODE
        HSTIF(IEVAB,JNODE)=HSTIF(IEVAB,JNODE)
     .                    +CONST*FMULT*SHAPE(JNODE)*DVOLU
   10   CONTINUE
C
      ELSE
C
C***LARGE DEFORMATIONS ############???????????????????????????????
C
        NTIME=NSTRE
        IF(NTYPE.EQ.4)NTIME=4
C
        DO 30 IEVAB=1,NEVAB
C
        FMULT=0.
        DO 40 ITIME=1,NTIME
   40   FMULT=FMULT+BMATX(ITIME,IEVAB)*DELTA(ITIME) ! >>>>>>>> MAL
C
        DO 30 JNODE=1,NNODE
        HSTIF(IEVAB,JNODE)=HSTIF(IEVAB,JNODE)
     .                    +CONST*FMULT*SHAPE(JNODE)*DVOLU
   30   CONTINUE
C       
      ENDIF
      RETURN
      END
