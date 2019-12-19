      SUBROUTINE OUTCAU(CARTD,ELDIS,STRSG,STRAN,XJACM,GPCOD,SHAPE,XJACI,
     .                  BSBAR)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE CAUCHY STRESS FROM THE SECOND PIOLA
C     KIRCHHOFF STRESS
C     ( ELEMENT NO. 1 )
C
C
C     Note:
C
C     FPS33 is the zz component of the deformation gradient for plane
C     stress conditions. It is useful, only for large strains, to
C     compute the Cauchy stress tensor.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      DIMENSION CARTD(NDIME,NNODS,*), ELDIS(NDIME,*)
      DIMENSION XJACM(NDIME,*),       XJACI(NDIME,*)
      DIMENSION STRSG(NSTR1,*),       STRAN(NSTR1,*),
     .          SIGMA(3,3),           AUXIL(3,3)
      DIMENSION GPCOD(NDIME,*),       SHAPE(NNODS,*)
      DIMENSION BSBAR(NSTR1,NEVAB,*)
      DIMENSION XJACN(3,3),           XJANI(3,3)
C
      KPARB=INT(PPARB)
C
C**** LOOP ON INTEGRATION POINTS
C
      DO 100 IGAUS=1,NGAUL
C
C**** CALCULATE THE JACOBIAN MATRIX
C
C     Note: NNODS instead of NNODL has been considered
C
      IF(KPARB.LE.0) THEN
       CALL LINEAR(CARTD(1,1,IGAUS), DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .             GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODS,NSTR1,
     .             PROPS,            SHAPE(1,IGAUS),DUMMY,
     .             DUMMY,DUMMY,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,DETJB,
     .             ELCOD,
     .                 1)
      ELSE
       CALL LINEAB(BSBAR(1,1,IGAUS), DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .             GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODS,NSTR1,
     .             PROPS,            SHAPE(1,IGAUS),DUMMY,
     .             DUMMY,DUMMY,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .             ELCOD,CARTD(1,1,IGAUS),NEVAB,
     .             XJACN,XJANI,XJA3N,XJ3NI,DETJN,DETJB,
     .                 1)
      ENDIF
C
C**** BUILD STRESS TENSOR
C
      SIGMA(1,1)=STRSG(1,IGAUS)
      IF(NTYPE.EQ.4)THEN
       SIGMA(1,2)=STRSG(3,IGAUS)
       SIGMA(1,3)=STRSG(5,IGAUS)
       SIGMA(2,2)=STRSG(2,IGAUS)
       SIGMA(2,3)=STRSG(6,IGAUS)
       SIGMA(3,3)=STRSG(4,IGAUS)
      ELSE IF(NTYPE.NE.5) THEN
       SIGMA(1,2)=STRSG(3,IGAUS)
       SIGMA(2,2)=STRSG(2,IGAUS)
       IF(NTYPE.NE.1) SIGMA(3,3)=STRSG(4,IGAUS)
      ENDIF
      DO 60  I=1,NDIME
      DO 60  J=I,NDIME
   60 SIGMA(J,I)=SIGMA(I,J)
C                                           t
C**** COMPUTE CAUCHY STRESS: S = 1/J ( F:S:F )
C
      DO 70 I=1,NDIME
      DO 70 J=1,NDIME
      AUXIL(I,J)=0.D0
      DO 70 K=1,NDIME
      DO 70 L=1,NDIME
   70 AUXIL(I,J)=AUXIL(I,J)+XJACM(I,K)*SIGMA(K,L)*XJACM(J,L)
      DO 80 I=1,NDIME
      DO 80 J=1,NDIME
   80 SIGMA(I,J)=AUXIL(I,J)  
      IF(NTYPE.EQ.3) SIGMA(3,3)=XJA3M*SIGMA(3,3)*XJA3M
C
      IF(NTYPE.EQ.1) THEN
       FPS33=DSQRT(2.0D0*STRAN(4,IGAUS)+1.0D0)
       DETJB=DETJB*FPS33
      ENDIF
C
      STRSG(1,IGAUS)=SIGMA(1,1)/DETJB
      IF(NTYPE.EQ.4)THEN
       STRSG(3,IGAUS)=SIGMA(1,2)/DETJB
       STRSG(5,IGAUS)=SIGMA(1,3)/DETJB
       STRSG(2,IGAUS)=SIGMA(2,2)/DETJB
       STRSG(6,IGAUS)=SIGMA(2,3)/DETJB
       STRSG(4,IGAUS)=SIGMA(3,3)/DETJB
      ELSE IF(NTYPE.NE.5) THEN
       STRSG(3,IGAUS)=SIGMA(1,2)/DETJB
       STRSG(2,IGAUS)=SIGMA(2,2)/DETJB
       IF(NTYPE.NE.1) STRSG(4,IGAUS)=SIGMA(3,3)/DETJB
      ENDIF
C
  100 CONTINUE
C
      RETURN
      END
