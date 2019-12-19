      SUBROUTINE BINCOM(BMATX,CARTD,ELDIS,GPCOD,LARGE,NDIME,NDOFC,
     .                  NNODE,NSTRE,NTYPE,SHAPE,XJACM,DMATX,NSTRS)
C********************************************************************
C
C**** THIS ROUTINE EVALUATES THE B-MATRIX & D-MATRIX FOR THE
C     INCOMPRESSIBILITY CONDITION
C
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION BMATX(NSTRE,*), CARTD(NDIME,*), ELDIS(NDOFC,*),
     .          GPCOD(*),       SHAPE(*),       XJACM(NDIME,*),
     .          DMATX(NSTRS,*)
C
C**** DEFINES DMATX
C
      DMATX(1,1)=1.0
C
C**** 1-D CASE
C
      IF(NTYPE.EQ.5) THEN
       DO INODE=1,NNODE
        BMATX(1,INODE)=CARTD(1,INODE)
       ENDDO
       RETURN
      ENDIF
C
C**** BUILDS UP THE BMATX
C
      LGASH=0
      DO 30 INODE=1,NNODE
      MGASH=LGASH+1
      NGASH=MGASH+1
      LGASH=NGASH
      BMATX(1,MGASH)=CARTD(1,INODE)
      BMATX(1,NGASH)=CARTD(2,INODE)
C
C**** COMPLETES FOR THE AXI-SYMMETRIC CASE (nothing in this case)
C
      IF(NTYPE.NE.3)GOTO 40
C

C
C**** COMPLETES FOR THE 3-D CASE
C
   40 IF(NTYPE.NE.4)GOTO 30
C
      LGASH=NGASH+1
      BMATX(1,LGASH)=CARTD(3,INODE)
C
   30 CONTINUE
C
      RETURN
      END
