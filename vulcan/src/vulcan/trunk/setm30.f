      SUBROUTINE SETM30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,
     .                  PROPS,RMAT1,SHAPE,STIFH,THICK,
     .                  DERIV,POSGP,WEIGP,XJACM,RMAT2,CMEAN,SHAPR,
     .                  WSTIR,VNORL,BSBAR,CARTS,SHAPS,GPCOS,ITASI,ITASK)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME NEEDED CONSTANT MATRICES FOR
C     FUTURE USE ( ELEMENT NO. 30 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      COMMON/JACOBSA/IEROR,KEROR
C
      DIMENSION CARTD(NDIME,NNODS,*), DVOLU(*),
     .          ELCOD(NDIME,*),       EPMTX(*),
     .          GPCOD(NDIME,*),       LNODS(*),
     .          PROPS(*),             RMAT1(NDIME,*), 
     .          SHAPE(NNODS,*),       EMASS(*),
     .          STIFH(*)
      DIMENSION DERIV(NDIME,NNODS,*), 
     .          POSGP(NDIME,*),       WEIGP(*),
     .          XJACM(NDIME,*),       RMAT2(NSTRS,*)
      DIMENSION CMEAN(3,*),           SHAPR(NNODS,*),
     .          WSTIR(NNODS,*),       VNORL(NDIME,*),
     .          BSBAR(NSTR1,NEVAB,*)
      DIMENSION CARTS(NDIME,NNODS,*), SHAPS(NNODS,*), GPCOS(NDIME,*)
C
      TWOPI=6.283185307179586D0
C
C**** INITIALIZATION
C
      DO IDIME=1,NDIME
       DO IGAUS=1,NGAUS
        GPCOD(IDIME,IGAUS)=0.0D0
       ENDDO
      ENDDO
      DO IGAUS=1,NGAUS
       WEIGP(IGAUS)=0.0D0
      ENDDO
C
C**** EVALUATE QUANTITIES ASSOCIATED WITH THE GAUSSIAN POINTS
C
c     IFLAG=INT(PROPS(1))
C
C**** ROTATION MATRIX 
C
c     IF(IFLAG.EQ.1)
c    . CALL ROTATE(IELEM,NDIME,NSTRS,PROPS,RMAT1,RMAT2)
C
C**** ELASTIC CONSTANT MATRIX
C
ctm   CALL ELASTI(NSTRS,NTYPE,EPMTX,PROPS,RMAT2,ncrit)
C
C**** IDENTIFY SAMPLING POINTS AND WEIGHTS
C
      CALL RULEPW(NDIME,NGAUL,NRULE,POSGP,WEIGP)
C
C**** LOOP ON INTEGRATION POINTS
C
      DO 100 IGAUS=1,NGAUL
C
C**** COMPUTE SHAPE FUNCTIONS AND DERIVATIVES
C
      EXISP=POSGP(1,IGAUS)
      ETASP=0.0D0
      IF(NDIME.GE.2) ETASP=POSGP(2,IGAUS)
      EZETA=0.0D0
      IF(NDIME.EQ.3) EZETA=POSGP(3,IGAUS)
      CALL SHAFUN(DERIV(1,1,IGAUS),EXISP,ETASP,EZETA,NDIME,NNODS,
     .            NQUTR,    0,SHAPE(1,IGAUS))
C
C**** CARTESIAN DERIVATIVES 
C
      CALL JACOBS(CARTD(1,1,IGAUS), DERIV(1,1,IGAUS), DETJM,
     .            ELCOD,            GPCOD(1,IGAUS),   IELEM,
     .            NDIME,NNODS,      SHAPE(1,IGAUS),   XJACM,
     .            LURES,LUPRI)
      IF(IEROR.NE.0) GO TO 101
C
C**** MAKE A CONTROL ON THE CORRECTNESS OF THE SHAPE FUNCTION AND
C     THEIR DERIV.
C
      ICHEK=0
      IF(ICHEK.EQ.1)CALL TESTM2(CARTD(1,1,IGAUS),DERIV(1,1,IGAUS),
     .                          SHAPE(1,  IGAUS),IGAUS,NDIME,NNODS)
C
C**** INTEGRATION WEIGHTS
C
                       DVOLU(IGAUS)=WEIGP(IGAUS)*DETJM
      IF(NTYPE.EQ.3)   DVOLU(IGAUS)=DVOLU(IGAUS)*TWOPI*GPCOD(1,IGAUS)
      IF(NTYPE.NE.4)   DVOLU(IGAUS)=DVOLU(IGAUS)*THICK
C
  100 CONTINUE
C
C**** CHECKS GAUSSIAN VOLUME FOR AXISYMMETRIC PROBLEMS FOR OUTPUT OPER.
C
      IF(NTYPE.EQ.3) THEN
       CALL GAUCEK(DVOLU,   30,ITASK,NGAUL,NNODL,ICEKE)
       IF(ICEKE.NE.0) CALL RUNEND('ERROR: WRONG NGAUS IN gaucek.f-30')
      ENDIF
C
      IF(ITASI.EQ.1) GO TO 101
C
C**** EVALUATION OF B-BAR-SHEAR
C
C     Notes:
C
C     A constant B-bar-shear matrix (BSBAR) can only be evaluated for
C     small strains. For large strains, this matrix is not constant and,
C     therefore, can not be computed in this routine.
C
C     For small strains, CMEAN is computed within a NGAUS loop (except 
C     for KPARB=4, 5 and 6). For large strains, it is necessary to
C     calculate CMEAN for every NGAUS and, therefore, this array is
C     stored in BSBAR. 
C
C     BSBAR is not evaluated for output operations (ITASK=11-14) when
C     LARGE=0 (for LARGE > 0 it is necessary to compute BSBAR)
C
C     For LARGE > 0, KPARB=1,3,4,6,7,9 compute B-bar considering the
C     following decomposition: F-bar=Fdev*Fvol-bar
C     For LARGE > 0, KPARB=2,3,5,6,8,9 compute B-(bar+shear) considering
C     constant Green-Lagrange shear strains
C     For LARGE > 0, KPARB=10-12 is another option to compute B-bar
C     For LARGE > 0, KPARB=13-15 uses B-e_shear which is based on the
C     shear Almansi strains
C     For LARGE > 0, KPARB=16-18 uses B-(e_shear-bar) which is based on
C     the shear Almansi strains and F-bar=Fdev*Fvol-bar
C     For LARGE > 0, KPARB=19-21 uses B-bar_F33 which is identical to
C     B-bar considering F_33=F_33-bar
C     For LARGE > 0, KPARB=22-24 uses J-bar
C
      ITASP=0
      IF((ITASK.GE.11.AND.ITASK.LE.14).AND.LARGE.EQ.0) ITASP=1
C
      IF(ITASP.EQ.0) THEN
C
C**** COMPUTES PRELIMINARY ARRAYS
C
       KPARB=INT(PPARB)
       IF(KPARB.NE.0) THEN
        IF(KPARB.LT.0) KPARB=-KPARB
C
C**** COMPUTE "REDUCED" SHAPE FUNCTIONS
C
        IF((KPARB.GE.1.AND.KPARB.LE.3).OR.
     .     (KPARB.GE.7.AND.KPARB.LE.9).OR.
     .     (KPARB.EQ.10.OR.KPARB.EQ.12).OR.
     .     (KPARB.EQ.13.OR.KPARB.EQ.15).OR.
     .     (KPARB.EQ.16.OR.KPARB.EQ.18).OR.
     .     (KPARB.EQ.19.OR.KPARB.EQ.21).OR.
     .     (KPARB.EQ.22.OR.KPARB.EQ.24)) THEN
         DO IGAUS=1,NGAUL
          EXISP=POSGP(1,IGAUS)
          ETASP=0.0D0
          IF(NDIME.GE.2) ETASP=POSGP(2,IGAUS)
          EZETA=0.0D0
          IF(NDIME.EQ.3) EZETA=POSGP(3,IGAUS)
          CALL SHAFUR(EXISP,ETASP,EZETA,NDIME,NNODS,NGAUR,
     .                NQUTR,SHAPR(1,IGAUS))
         ENDDO
        ENDIF                  ! (kparb.ge.1.and.kparb.le.3)...
C
C**** COMPUTE "REDUCED" MASS MATRIX
C
        IF((KPARB.GE.1.AND.KPARB.LE.3).OR.
     .      KPARB.EQ.10.OR.
     .      KPARB.EQ.13.OR.
     .      KPARB.EQ.16.OR.
     .      KPARB.EQ.19.OR.
     .      KPARB.EQ.22) THEN
         CALL WMATRR(WSTIR,SHAPR,DVOLU,NGAUR,NGAUL,NNODS)
        ENDIF                  ! (kparb.ge.1.and.kparb.le.3)...
C
C**** COMPUTE "REDUCED" SHAPE FUNCTIONS AND DERIVATIVES
C
        IF((KPARB.GE.7.AND.KPARB.LE.9).OR.
     .      KPARB.EQ.12.OR.
     .      KPARB.EQ.15.OR.
     .      KPARB.EQ.18.OR.
     .      KPARB.EQ.21.OR.
     .      KPARB.EQ.24) THEN
         DO IDIME=1,NDIME
          DO IGAUS=1,NGAUR
           GPCOS(IDIME,IGAUS)=0.0D0
          ENDDO
         ENDDO
C
         NRULX=NRULE
         IF(NDIME.EQ.3.AND.NNODL.EQ.5) NRULX=6
C
         CALL RULEPW(NDIME,NGAUR,NRULX,POSGP,WEIGP)
C
         DO IGAUS=1,NGAUR
          EXISP=POSGP(1,IGAUS)
          ETASP=0.0D0
          IF(NDIME.GE.2) ETASP=POSGP(2,IGAUS)
          EZETA=0.0D0
          IF(NDIME.EQ.3) EZETA=POSGP(3,IGAUS)
          CALL SHAFUN(DERIV(1,1,IGAUS),EXISP,ETASP,EZETA,NDIME,NNODS,
     .                NQUTR,    0,SHAPS(1,IGAUS))
C
          CALL JACOBS(CARTS(1,1,IGAUS), DERIV(1,1,IGAUS), DETJM,
     .                ELCOD,            GPCOS(1,IGAUS),   IELEM,
     .                NDIME,NNODS,      SHAPS(1,IGAUS),   XJACM,
     .                LURES,LUPRI)
         ENDDO
        ENDIF                  ! (kparb.ge.7.and.kparb.le.9)...
C
        IF(LARGE.EQ.0) THEN
C
C**** B-(BAR-SHEAR) MATRIX
C
         GO TO (1,2,3,4,5,6,7,8,9), KPARB
C
C**** B-BAR (SMOOTHING TECHNIQUE)
C
    1    DO IGAUS=1,NGAUL
          CALL BBMEAR(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,SHAPR,WSTIR,NGAUR,
     .                IGAUS)
          CALL BBMATX(BSBAR(1,1,IGAUS),CMEAN,CARTD(1,1,IGAUS),
     .                GPCOD(1,IGAUS),SHAPE(1,IGAUS))
         ENDDO
         GO TO 14
C
C**** B-SHEAR (SMOOTHING TECHNIQUE)
C
    2    DO IGAUS=1,NGAUL
          CALL BBMEAR(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,SHAPR,WSTIR,NGAUR,
     .                IGAUS)
          CALL BBMATY(BSBAR(1,1,IGAUS),CMEAN,CARTD(1,1,IGAUS),
     .                GPCOD(1,IGAUS),SHAPE(1,IGAUS))
         ENDDO
         GO TO 14
C
C**** B-BAR+SHEAR (SMOOTHING TECHNIQUE)
C
    3    DO IGAUS=1,NGAUL
          CALL BBMEAR(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,SHAPR,WSTIR,NGAUR,
     .                IGAUS)
          CALL BBMATZ(BSBAR(1,1,IGAUS),CMEAN,CARTD(1,1,IGAUS),
     .                GPCOD(1,IGAUS),SHAPE(1,IGAUS))
         ENDDO
         GO TO 14
C
C**** B-BAR (AVERAGE TECHNIQUE)
C
    4    CALL BBMEAN(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,    3)
         DO IGAUS=1,NGAUL
          CALL BBMATX(BSBAR(1,1,IGAUS),CMEAN,CARTD(1,1,IGAUS),
     .                GPCOD(1,IGAUS),SHAPE(1,IGAUS))
         ENDDO
         GO TO 14
C
C**** B-SHEAR (AVERAGE TECHNIQUE)
C
    5    CALL BBMEAN(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,    3)
         DO IGAUS=1,NGAUL
          CALL BBMATY(BSBAR(1,1,IGAUS),CMEAN,CARTD(1,1,IGAUS),
     .                GPCOD(1,IGAUS),SHAPE(1,IGAUS))
         ENDDO
         GO TO 14
C
C**** B-BAR+SHEAR (AVERAGE TECHNIQUE)
C
    6    CALL BBMEAN(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,    3)
         DO IGAUS=1,NGAUL
          CALL BBMATZ(BSBAR(1,1,IGAUS),CMEAN,CARTD(1,1,IGAUS),
     .                GPCOD(1,IGAUS),SHAPE(1,IGAUS))
         ENDDO
         GO TO 14
C
C**** B-BAR (SELECTIVE INTEGRATION)
C
    7    DO IGAUS=1,NGAUL
          CALL BBMEAS(CARTS,CMEAN,GPCOS,SHAPS,DVOLU,SHAPR,WSTIR,NGAUR,
     .                IGAUS)
          CALL BBMATX(BSBAR(1,1,IGAUS),CMEAN,CARTD(1,1,IGAUS),
     .                GPCOD(1,IGAUS),SHAPE(1,IGAUS))
         ENDDO
         GO TO 14
C
C**** B-SHEAR (SELECTIVE INTEGRATION)
C
    8    DO IGAUS=1,NGAUL
          CALL BBMEAS(CARTS,CMEAN,GPCOS,SHAPS,DVOLU,SHAPR,WSTIR,NGAUR,
     .                IGAUS)
          CALL BBMATY(BSBAR(1,1,IGAUS),CMEAN,CARTD(1,1,IGAUS),
     .                GPCOD(1,IGAUS),SHAPE(1,IGAUS))
         ENDDO
         GO TO 14
C
C**** B-BAR+SHEAR (SELECTIVE INTEGRATION)
C
    9    DO IGAUS=1,NGAUL
          CALL BBMEAS(CARTS,CMEAN,GPCOS,SHAPS,DVOLU,SHAPR,WSTIR,NGAUR,
     .                IGAUS)
          CALL BBMATZ(BSBAR(1,1,IGAUS),CMEAN,CARTD(1,1,IGAUS),
     .                GPCOD(1,IGAUS),SHAPE(1,IGAUS))
         ENDDO
         GO TO 14
C
   14    CONTINUE
C
        ELSE                 ! large strains
C
C**** CARTESIAN DERIVATIVE OF THE SHAPE FUNCTIONS
C 
         GO TO (21,21,21,22,22,22,23,23,23,
     .          21,22,23,21,22,23,21,22,23,
     .          21,22,23,21,22,23), KPARB
C
C**** CMEAN (SMOOTHING TECHNIQUE)
C
   21    DO IGAUS=1,NGAUL
          CALL BBMEAR(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,SHAPR,WSTIR,NGAUR,
     .                IGAUS)
          DO IDIME=1,3
           DO INODE=1,NNODL
            BSBAR(IDIME,INODE,IGAUS)=CMEAN(IDIME,INODE)
           ENDDO
          ENDDO
         ENDDO
         GO TO 34
C
C**** CMEAN (AVERAGE TECHNIQUE)
C
   22    DO IGAUS=1,NGAUL
          CALL BBMEAN(CARTD,BSBAR(1,1,IGAUS),GPCOD,SHAPE,DVOLU,NSTR1)
         ENDDO
         GO TO 34
C
C**** CMEAN (SELECTIVE INTEGRATION)
C
   23    DO IGAUS=1,NGAUL
          CALL BBMEAS(CARTS,CMEAN,GPCOS,SHAPS,DVOLU,SHAPR,WSTIR,NGAUR,
     .                IGAUS)
          DO IDIME=1,3
           DO INODE=1,NNODL
            BSBAR(IDIME,INODE,IGAUS)=CMEAN(IDIME,INODE)
           ENDDO
          ENDDO
         ENDDO
         GO TO 34
C
   34    CONTINUE
C
        ENDIF              ! large.eq.0
C
       ENDIF               ! kparb.ne.0
      ENDIF                ! itasp.eq.0
C
C**** COMPUTE OUTWARD UNIT NORMAL (FOR CONTACT PROBLEM)
C
C     Note: only necessary when using node to node contact element
C           (ITYPE=4). For the side to side contact element (ITYPE=32),
C           the normal vector is computed in setm32.f
C
      IF(NTYPE.EQ.3) THEN
       GPCOA=0.0D0                                      ! average radius
       DO IGAUL=1,NGAUL
        GPCOA=GPCOA+GPCOD(1,IGAUL)
       ENDDO
       GPCOA=GPCOA/NGAUL
      ENDIF
C
      DO IDIME=1,NDIME
       DO INODL=1,NNODS
        VNORL(IDIME,INODL)=0.0D+00
        DO IGAUS=1,NGAUL
         DVOLI=DVOLU(IGAUS)
         IF(NTYPE.EQ.3) THEN
          IF(GPCOD(1,IGAUS).GT.(1.0E-10*GPCOA))       ! see GPCOA in *.f
     .                         DVOLI=DVOLU(IGAUS)/(TWOPI*GPCOD(1,IGAUS))
         ENDIF
         VNORL(IDIME,INODL)=VNORL(IDIME,INODL)+
     .                      CARTD(IDIME,INODL,IGAUS)*DVOLI
        ENDDO
       ENDDO
      ENDDO
C
C**** COMPUTE MASS MATRIX TO PERFORM SMOOTHING OF GAUSSIAN TENSORS
C
      IF(KSGAU.NE.0.OR.IGALE.EQ.1) THEN
       CALL SMOMTX(ELCOD,DVOLU,SHAPE,EMASS)
      ENDIF
C
C**** COMPUTE STABILIZATION MATRIX FOR HOURGLASS CONTROL
C
      IF(NHOUR.NE.0) THEN
       CALL HOURMX(CARTD,EPMTX,DVOLU,ELCOD,GPCOD,PROPS,STIFH,THICK)
      ENDIF
C
  101 CONTINUE
      RETURN
      END
