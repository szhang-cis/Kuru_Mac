      SUBROUTINE BBMATWL(CARTD,ELDIS,GPCOD,
     .                   NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .                   NTYPE,SHAPE,NSTR1,NSTRS,BMATX,CMEAN,
     .                   XJACM,XJA3M,DETJM,
     .                   XJACN,XJA3N,DETJN)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE B-SHEAR STRAIN-DISPLACEMENT MATRIX
C     FOR LARGE STRAINS (LARGE=1)
C
C     Notes:
C
C     This routine is based on blarge.f
C     CMEAN in this routine has the dimensions of BSBAR
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION CARTD(NDIME,*), 
     .          ELDIS(NDOFC,*), SHAPE(*),
     .          XJACM(NDIME,*), XJACN(NDIME,*)
      DIMENSION BMATX(NSTRE,*), CMEAN(NSTR1,*)
C
      B13=1.0D0/3.0D0
      B23=2.0D0/3.0D0
C
C**** 1-D CASE
C
      IF(NTYPE.EQ.5) THEN
       DO INODE=1,NNODE
        BMATX(1,INODE)=CARTD(1,INODE)*XJACM(1,1)
       ENDDO
      ENDIF
C
C**** CALCULATE CONSTANT TERM FOR AXI-SYMMETRIC CASE IN FINITE
C     DEFORMATIONS
C
      FMULT=1.0D0
      IF(NTYPE.NE.3) GOTO 10
      FMULT=0.0D0
      DO 20 JNODE=1,NNODE
   20 FMULT=FMULT+ELDIS(1,JNODE)*SHAPE(JNODE)
      FMULT=FMULT/GPCOD                         ! =XJA3M
   10 CONTINUE
C
C**** BUILD UP THE BMATX
C
      LGASH=0
      DO 30 INODE=1,NNODE
      MGASH=LGASH+1
      NGASH=MGASH+1
      LGASH=NGASH
      BMATX(1,MGASH)=B23*CARTD(1,INODE)*XJACM(1,1)+
     .               B13*CMEAN(1,INODE)*XJACN(1,1)+
     .               B13*CMEAN(2,INODE)*XJACN(1,2)-
     .               B13*CARTD(2,INODE)*XJACM(1,2)
      BMATX(1,NGASH)=B23*CARTD(1,INODE)*XJACM(2,1)+
     .               B13*CMEAN(1,INODE)*XJACN(2,1)+
     .               B13*CMEAN(2,INODE)*XJACN(2,2)-
     .               B13*CARTD(2,INODE)*XJACM(2,2)
      BMATX(2,MGASH)=B23*CARTD(2,INODE)*XJACM(1,2)+
     .               B13*CMEAN(2,INODE)*XJACN(1,2)+
     .               B13*CMEAN(1,INODE)*XJACN(1,1)-
     .               B13*CARTD(1,INODE)*XJACM(1,1)
      BMATX(2,NGASH)=B232*CARTD(2,INODE)*XJACM(2,2)+
     .               B13*CMEAN(2,INODE)*XJACN(2,2)+
     .               B13*CMEAN(1,INODE)*XJACN(2,1)-
     .               B13*CARTD(1,INODE)*XJACM(2,1)
      BMATX(3,MGASH)=CARTD(2,INODE)*XJACM(1,1)+
     .               CARTD(1,INODE)*XJACM(1,2)
      BMATX(3,NGASH)=CARTD(2,INODE)*XJACM(2,1)+
     .               CARTD(1,INODE)*XJACM(2,2)
C
C**** COMPLETE FOR THE AXI-SYMMETRIC CASE
C
      IF(NTYPE.NE.3) GOTO 40
C
      BMATX(1,MGASH)=BMATX(1,MGASH)+B13*CMEAN(3,INODE)*XJA3N-
     .                              B13*SHAPE(INODE)*FMULT/GPCOD
      BMATX(2,MGASH)=BMATX(2,MGASH)+B13*CMEAN(3,INODE)*XJA3N-
     .                              B13*SHAPE(INODE)*FMULT/GPCOD
C
      BMATX(4,MGASH)=B23*SHAPE(INODE)*FMULT/GPCOD+
     .               B13*CMEAN(3,INODE)*XJA3N+
     .               B13*CMEAN(1,INODE)*XJACN(1,1)-
     .               B13*CARTD(1,INODE)*XJACM(1,1)+
     .               B13*CMEAN(2,INODE)*XJACN(1,2)-
     .               B13*CARTD(2,INODE)*XJACM(1,2)
      BMATX(4,NGASH)=B13*CMEAN(1,INODE)*XJACN(2,1)-
     .               B13*CARTD(1,INODE)*XJACM(2,1)+
     .               B13*CMEAN(2,INODE)*XJACN(2,2)-
     .               B13*CARTD(2,INODE)*XJACM(2,2)
C
C**** COMPLETE FOR THE 3-D CASE         ! to be revised !!!!
C
   40 IF(NTYPE.NE.4) GOTO 30
C
      LGASH=NGASH+1
      BMATX(1,LGASH)=CARTD(1,INODE)*XJACM(3,1)
      BMATX(2,LGASH)=CARTD(2,INODE)*XJACM(3,2)
      BMATX(3,LGASH)=CMEAN(2,INODE)*XJACN(3,1)+
     .               CMEAN(1,INODE)*XJACN(3,2)
      BMATX(4,MGASH)=CARTD(3,INODE)*XJACM(1,3)
      BMATX(4,NGASH)=CARTD(3,INODE)*XJACM(2,3)
      BMATX(4,LGASH)=CARTD(3,INODE)*XJACM(3,3)
      BMATX(5,MGASH)=CMEAN(1,INODE)*XJACN(1,3)+
     .               CMEAN(3,INODE)*XJACN(1,1)
      BMATX(5,NGASH)=CMEAN(1,INODE)*XJACN(2,3)+
     .               CMEAN(3,INODE)*XJACN(2,1)
      BMATX(5,LGASH)=CMEAN(1,INODE)*XJACN(3,3)+
     .               CMEAN(3,INODE)*XJACN(3,1)
      BMATX(6,MGASH)=CMEAN(2,INODE)*XJACN(1,3)+
     .               CMEAN(3,INODE)*XJACN(1,2)
      BMATX(6,NGASH)=CMEAN(2,INODE)*XJACN(2,3)+
     .               CMEAN(3,INODE)*XJACN(2,2)
      BMATX(6,LGASH)=CMEAN(2,INODE)*XJACN(3,3)+
     .               CMEAN(3,INODE)*XJACN(3,2)
C
   30 CONTINUE
C
      RETURN
      END
