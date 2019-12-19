      SUBROUTINE BBMATXL(CARTD,ELDIS,GPCOD,
     .                   NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .                   NTYPE,SHAPE,NSTR1,NSTRS,BMATX,CMEAN,
     .                   XJACM,XJA3M,DETJM,DEFRA,
     .                   XJACN,XJA3N,DETJN,
     .                   KPARB)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE B-BAR STRAIN-DISPLACEMENT MATRIX
C     FOR LARGE STRAINS (LARGE=1)
C
C     Notes:
C
C     This routine is based on blarge.f
C     CMEAN in this routine has the dimensions of BSBAR
C     For KPARB=19-21 (only useful for NTYPE=3), B-BAR_F33 is computed
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION CARTD(NDIME,*), 
     .          ELDIS(NDOFC,*), SHAPE(*),
     .          XJACM(NDIME,*), XJACN(NDIME,*)
      DIMENSION BMATX(NSTRE,*), CMEAN(NSTR1,*)
      DIMENSION RCGTT(6)
C
C**** PRELIMINARY COMPUTATIONS
C
      DEFR2=DEFRA*DEFRA
      B13=1.0D0/3.0D0
      B23=2.0D0/3.0D0
C
      RCGTT(1)=XJACM(1,1)*XJACM(1,1)/DEFR2   ! right Cauchy-Green tensor
      DXNU1=1.0D0
      DXMU1=1.0D0
      IF(NTYPE.NE.5) THEN
       RCGTT(1)=(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1))/DEFR2
       RCGTT(2)=(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2))/DEFR2
       RCGTT(3)=(XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2))/DEFR2
       RCGTT(4)=1.0D0
       DXNU1= XJACN(2,2)
       DXNU2=-XJACN(2,1)
       DXMU1= XJACM(2,2)/DEFRA
       DXMU2=-XJACM(2,1)/DEFRA
       DXNV1=-XJACN(1,2)
       DXNV2= XJACN(1,1)
       DXMV1=-XJACM(1,2)/DEFRA
       DXMV2= XJACM(1,1)/DEFRA
       IF(NTYPE.EQ.3) THEN
        RCGTT(4)=XJA3M*XJA3M/DEFR2
        DXNU1=DXNU1*XJA3N
        DXNU2=DXNU2*XJA3N
        DXMU1=DXMU1*XJA3M/DEFRA
        DXMU2=DXMU2*XJA3M/DEFRA
        DXNV1=DXNV1*XJA3N
        DXNV2=DXNV2*XJA3N
        DXMV1=DXMV1*XJA3M/DEFRA
        DXMV2=DXMV2*XJA3M/DEFRA
        DE2DN=DETJN/XJA3N
        DE2DM=DETJM/(XJA3M/DEFRA)
       ENDIF
       IF(NTYPE.EQ.4) THEN
        RCGTT(1)=(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)+
     .            XJACM(3,1)*XJACM(3,1))/DEFR2
        RCGTT(2)=(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)+
     .            XJACM(3,2)*XJACM(3,2))/DEFR2
        RCGTT(3)=(XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2)+
     .            XJACM(3,1)*XJACM(3,2))/DEFR2
        RCGTT(4)=(XJACM(1,3)*XJACM(1,3)+XJACM(2,3)*XJACM(2,3)+
     .            XJACM(3,3)*XJACM(3,3))/DEFR2
        RCGTT(5)=(XJACM(1,1)*XJACM(1,3)+XJACM(2,1)*XJACM(2,3)+
     .            XJACM(3,1)*XJACM(3,3))/DEFR2
        RCGTT(6)=(XJACM(1,2)*XJACM(1,3)+XJACM(2,2)*XJACM(2,3)+
     .            XJACM(3,2)*XJACM(3,3))/DEFR2
        DXNU1= (XJACN(2,2)*XJACN(3,3)-XJACN(3,2)*XJACN(2,3))
        DXNU2=-(XJACN(2,1)*XJACN(3,3)-XJACN(3,1)*XJACN(2,3))
        DXNU3= (XJACN(2,1)*XJACN(3,2)-XJACN(3,1)*XJACN(2,2))
        DXMU1= (XJACM(2,2)*XJACM(3,3)-XJACM(3,2)*XJACM(2,3))/DEFR2
        DXMU2=-(XJACM(2,1)*XJACM(3,3)-XJACM(3,1)*XJACM(2,3))/DEFR2
        DXMU3= (XJACM(2,1)*XJACM(3,2)-XJACM(3,1)*XJACM(2,2))/DEFR2
        DXNV1=-(XJACN(1,2)*XJACN(3,3)-XJACN(3,2)*XJACN(1,3))
        DXNV2= (XJACN(1,1)*XJACN(3,3)-XJACN(1,3)*XJACN(3,1))
        DXNV3=-(XJACN(1,1)*XJACN(3,2)-XJACN(3,1)*XJACN(1,2))
        DXMV1=-(XJACM(1,2)*XJACM(3,3)-XJACM(3,2)*XJACM(1,3))/DEFR2
        DXMV2= (XJACM(1,1)*XJACM(3,3)-XJACM(1,3)*XJACM(3,1))/DEFR2
        DXMV3=-(XJACM(1,1)*XJACM(3,2)-XJACM(3,1)*XJACM(1,2))/DEFR2
        DXNW1= (XJACN(1,2)*XJACN(2,3)-XJACN(2,2)*XJACN(1,3))
        DXNW2=-(XJACN(1,1)*XJACN(2,3)-XJACN(2,1)*XJACN(1,3))
        DXNW3= (XJACN(1,1)*XJACN(2,2)-XJACN(2,1)*XJACN(1,2))
        DXMW1= (XJACM(1,2)*XJACM(2,3)-XJACM(2,2)*XJACM(1,3))/DEFR2
        DXMW2=-(XJACM(1,1)*XJACM(2,3)-XJACM(2,1)*XJACM(1,3))/DEFR2
        DXMW3= (XJACM(1,1)*XJACM(2,2)-XJACM(2,1)*XJACM(1,2))/DEFR2
       ENDIF           ! ntype.eq.4
      ENDIF            ! ntype.ne.5
C
C**** 1-D CASE
C
      IF(NTYPE.EQ.5) THEN
       DO INODE=1,NNODE
        AU=(DXNU1*CMEAN(1,INODE))/DETJN-
     .     (DXMU1*CARTD(1,INODE))/DETJM
        BMATX(1,INODE)=(CARTD(1,INODE)*XJACM(1,1)/DEFRA+
     .                                           B13*RCGTT(1)*AU)*DEFRA2
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
      FMULT=FMULT/GPCOD                              ! = XJA3M
   10 CONTINUE
C
C**** BUILD UP THE BMATX
C
      LGASH=0
      DO 30 INODE=1,NNODE
      MGASH=LGASH+1
      NGASH=MGASH+1
      LGASH=NGASH
C
      AU=(DXNU1*CMEAN(1,INODE)+DXNU2*CMEAN(2,INODE))/DETJN-
     .   (DXMU1*CARTD(1,INODE)+DXMU2*CARTD(2,INODE))/DETJM
      AV=(DXNV1*CMEAN(1,INODE)+DXNV2*CMEAN(2,INODE))/DETJN-
     .   (DXMV1*CARTD(1,INODE)+DXMV2*CARTD(2,INODE))/DETJM
C
      BMATX(1,MGASH)=(CARTD(1,INODE)*XJACM(1,1)/DEFRA+
     .                                            B13*RCGTT(1)*AU)*DEFR2
      BMATX(1,NGASH)=(CARTD(1,INODE)*XJACM(2,1)/DEFRA+
     .                                            B13*RCGTT(1)*AV)*DEFR2
      BMATX(2,MGASH)=(CARTD(2,INODE)*XJACM(1,2)/DEFRA+
     .                                            B13*RCGTT(2)*AU)*DEFR2
      BMATX(2,NGASH)=(CARTD(2,INODE)*XJACM(2,2)/DEFRA+
     .                                            B13*RCGTT(2)*AV)*DEFR2
      BMATX(3,MGASH)=(CARTD(2,INODE)*XJACM(1,1)/DEFRA+
     .                CARTD(1,INODE)*XJACM(1,2)/DEFRA+
     .                                            B23*RCGTT(3)*AU)*DEFR2
      BMATX(3,NGASH)=(CARTD(2,INODE)*XJACM(2,1)/DEFRA+
     .                CARTD(1,INODE)*XJACM(2,2)/DEFRA+
     .                                            B23*RCGTT(3)*AV)*DEFR2
C
C**** COMPLETE FOR THE AXI-SYMMETRIC CASE
C
      IF(NTYPE.NE.3) GOTO 40
C
      AU33=DE2DN*CMEAN(3,INODE)/DETJN-
     .    (DE2DM*SHAPE(INODE)/GPCOD)/DETJM
C
      ACO33=1.0D0
c     IF(KPARB.EQ.19.OR.KPARB.EQ.20.OR.KPARB.EQ.21) ACO33=0.0D0
      IF(KPARB.EQ.19.OR.KPARB.EQ.20) ACO33=0.0D0    ! approxim.
C
      BMATX(1,MGASH)=BMATX(1,MGASH)+B13*RCGTT(1)*AU33*DEFR2*ACO33
      BMATX(2,MGASH)=BMATX(2,MGASH)+B13*RCGTT(2)*AU33*DEFR2*ACO33
      BMATX(3,MGASH)=BMATX(3,MGASH)+B23*RCGTT(3)*AU33*DEFR2*ACO33
C
      BMATX(4,MGASH)=(SHAPE(INODE)*FMULT/GPCOD+B13*RCGTT(4)*AU+
     .                                         B13*RCGTT(4)*AU33)*ACO33+
     .  (1.0D0-ACO33)*CMEAN(3,INODE)*XJA3N
      BMATX(4,NGASH)=B13*RCGTT(4)*AV*DEFR2*ACO33
C
C**** COMPLETE FOR THE 3-D CASE
C
   40 IF(NTYPE.NE.4) GOTO 30
C
      LGASH=NGASH+1
C
      AUX=(DXNU3*CMEAN(3,INODE))/DETJN-
     .    (DXMU3*CARTD(3,INODE))/DETJM
      AVX=(DXNV3*CMEAN(3,INODE))/DETJN-
     .    (DXMV3*CARTD(3,INODE))/DETJM
      AU=AU+AUX
      AV=AV+AVX
      AW=(DXNW1*CMEAN(1,INODE)+DXNW2*CMEAN(2,INODE)+
     .    DXNW3*CMEAN(3,INODE))/DETJN-
     .   (DXMW1*CARTD(1,INODE)+DXMW2*CARTD(2,INODE)+
     .    DXMW3*CARTD(3,INODE))/DETJM
C
      BMATX(1,MGASH)=BMATX(1,MGASH)+B13*RCGTT(1)*AUX*DEFR2
      BMATX(1,NGASH)=BMATX(1,NGASH)+B13*RCGTT(1)*AVX*DEFR2
      BMATX(2,MGASH)=BMATX(2,MGASH)+B13*RCGTT(2)*AUX*DEFR2
      BMATX(2,NGASH)=BMATX(2,NGASH)+B13*RCGTT(2)*AVX*DEFR2
      BMATX(3,MGASH)=BMATX(3,MGASH)+B23*RCGTT(3)*AUX*DEFR2
      BMATX(3,NGASH)=BMATX(3,NGASH)+B23*RCGTT(3)*AVX*DEFR2
C
      BMATX(1,LGASH)=(CARTD(1,INODE)*XJACM(3,1)/DEFRA+
     .                                            B13*RCGTT(1)*AW)*DEFR2
      BMATX(2,LGASH)=(CARTD(2,INODE)*XJACM(3,2)/DEFRA+
     .                                            B13*RCGTT(2)*AW)*DEFR2
      BMATX(3,LGASH)=(CARTD(2,INODE)*XJACM(3,1)/DEFRA+
     .                CARTD(1,INODE)*XJACM(3,2)/DEFRA+
     .                                            B23*RCGTT(3)*AW)*DEFR2
      BMATX(4,MGASH)=(CARTD(3,INODE)*XJACM(1,3)/DEFRA+
     .                                            B13*RCGTT(4)*AU)*DEFR2
      BMATX(4,NGASH)=(CARTD(3,INODE)*XJACM(2,3)/DEFRA+
     .                                            B13*RCGTT(4)*AV)*DEFR2
      BMATX(4,LGASH)=(CARTD(3,INODE)*XJACM(3,3)/DEFRA+
     .                                            B13*RCGTT(4)*AW)*DEFR2
      BMATX(5,MGASH)=(CARTD(1,INODE)*XJACM(1,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(1,1)/DEFRA+
     .                                            B23*RCGTT(5)*AU)*DEFR2
      BMATX(5,NGASH)=(CARTD(1,INODE)*XJACM(2,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(2,1)/DEFRA+
     .                                            B23*RCGTT(5)*AV)*DEFR2
      BMATX(5,LGASH)=(CARTD(1,INODE)*XJACM(3,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(3,1)/DEFRA+
     .                                            B23*RCGTT(5)*AW)*DEFR2
      BMATX(6,MGASH)=(CARTD(2,INODE)*XJACM(1,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(1,2)/DEFRA+
     .                                            B23*RCGTT(6)*AU)*DEFR2
      BMATX(6,NGASH)=(CARTD(2,INODE)*XJACM(2,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(2,2)/DEFRA+
     .                                            B23*RCGTT(6)*AV)*DEFR2
      BMATX(6,LGASH)=(CARTD(2,INODE)*XJACM(3,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(3,2)/DEFRA+
     .                                            B23*RCGTT(6)*AW)*DEFR2
C
   30 CONTINUE
C
      RETURN
      END
