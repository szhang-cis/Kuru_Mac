      SUBROUTINE STIF30(CARTD,DVOLU,EHIST,ELCOD,ELDIS,EPMTX,GPCOD,
     .                  PROPS,SHAPE,STRSG,ESTIF,HSTIF,EMASS,
     .                  BMATX,DMATX,SIGMA,XJACM,TENOD,BSBAR,DTENO,
     .                  XJACN,AUXS1)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS & COUPLING MATRICES
C         ( ELEMENT TYPE NO. 30 )
C
C     Notes:
C
C     STRATEGY PARAMETERS COMING FROM A COMMON OF TUNING PARAMETERS
C
C     KPART: Indicator for Temperature Interpolation)
C             = 0 (T=cte)  = 1 (T=linear)
C
C     KPARB: Indicator for B-cmean (see setm30.f)
C            =  0 (without B-cmean)
C            =  1 (with B-bar; smoothing technique)
C            =  2 (with B-shear; smoothing technique)
C            =  3 (with B-(bar+shear); smoothing technique)
C            =  4 (with B-bar; average technique)
C            =  5 (with B-shear; average technique)
C            =  6 (with B-(bar+shear); average technique)
C            =  7 (with B-bar; selective integration)
C            =  8 (with B-shear; selective integration)
C            =  9 (with B-(bar+shear); selective integration)
C
C            Only for LARGE > 0
C            = 10 (with another B-bar; smoothing technique)
C            = 11 (with another B-bar; average technique)
C            = 12 (with another B-bar; selective integration)
C            = 13 (with B-e_shear; smoothing technique)
C            = 14 (with B-e_shear; average technique)
C            = 15 (with B-e_shear; selective integration)
C            = 16 (with B-(e_shear+bar); smoothing technique)
C            = 17 (with B-(e_shear+bar); average technique)
C            = 18 (with B-(e_shear+bar); selective integration)
C            = 19 (with B-bar & F_33-bar; smoothing technique)
C            = 20 (with B-bar & F_33-bar; average technique)
C            = 21 (with B-bar & F_33-bar; selective integration)
C            = 22 (with J-bar; smoothing technique)
C            = 23 (with J-bar; average technique)
C            = 24 (with J-bar; selective integration)
C
C     KPARB < 0: same as KPARB > 0 but the stiffness is computed with 
C                the standard B (corresponding to KPARB=0) and, 
C                therefore, only the residual vector is calculated with
C                the B-mean chosen (corresponding to INT(KPARB))
C
C     KPARI: Indicator for Stress Print (see chek01.f)
C             = 0 (not smoothed)  = 1 (smoothed)
C
C     ESTAB: Artificial Shear Modulus (used in liquid phase only)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
      COMMON/JACOBSA/IEROR,KEROR
C
      DIMENSION CARTD(NDIME,NNODL,*), DVOLU(*),          EHIST(NHIST,*),
     .          ELCOD(NDIME,*),       ELDIS(NDOFC,*),    EPMTX(*),
     .          GPCOD(NDIME,*),       PROPS(*),          SHAPE(NNODL,*),
     .          STRSG(NSTR1,*),       ESTIF(NKOVA),
     .          HSTIF(NEVAB,*),
     .          EMASS(*)
      DIMENSION BMATX(NSTRE,*),       DMATX(NSTRS,*),     SIGMA(*),
     .          XJACM(NDIME,*)
      DIMENSION BSBAR(NSTR1,NEVAB,*), TENOD(*), DTENO(*)
      DIMENSION XJACN(NDIME,*),       AUXS1(NEVAB,*)
      DIMENSION XJACI(3,3),           XJANI(3,3)
C
      IF(IEROR.NE.0) RETURN
C
C**** INCOMPRESSIBILITY CONDITION
C
      IF(KPROB.EQ.5) THEN
C
C**** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
C
       DO IGAUS=1,NGAUL
        CALL BINCOM(BMATX,CARTD(1,1,IGAUS), ELDIS,GPCOD(1,IGAUS),
     .              LARGE,NDIME,NDOFC,NNODL,NSTRE,NTYPE,SHAPE(1,IGAUS),
     .              XJACM,DMATX,NSTRS)
        CALL KMATRB(BMATX,       CARTD(1,1,IGAUS), DMATX,
     .              DVOLU(IGAUS),ESTIF,            GPCOD(1,IGAUS),
     .              LARGE,KSYMM, NDIME,NDOFN,NEVAB,NKOVA,NNODL,
     .              NSTRE,NSTRS,NTYPE,SHAPE(1,IGAUS),
     .              NSTRE,NDOFC,ELDIS,AUXS1,    ! NSTRE instead of NSTR1
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .              BMATX,
     .              XJACN,XJANI,XJA3N,XJ3NI,DETJN)
       ENDDO
       RETURN
      ENDIF
C
C**** LOOK FOR MODEL DEFINITION (not used now)
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IPEP1=100
      IF(ISOTR.EQ.1) IPEP1=200
      IPEP2=INT(PROPS(2))   ! 10=elastic; 20=plastic; 30=viscoplastic;
C                           ! 40=hyperelastic
      IPEP3=INT(PROPS(3))   ! 1=no temp.-dependent; 2=temp.-depend.
      IPEPE=IPEP1+IPEP2+IPEP3
C
      KPART=INT(PPART)
      KPARB=INT(PPARB)
      KPARI=INT(PPARI)
C
C**** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
C
      DO 10 IGAUS=1,NGAUL
C
C**** SET UP THE ELASTIC MATRIX (only for isotropic materials !!!)
C
      IF(NMEMOM.EQ.1) THEN
       KOUNT=0
       DO ISTRS=1,NSTRS
        DO JSTRS=ISTRS,NSTRS
         KOUNT=KOUNT+1
         DMATX(ISTRS,JSTRS)=EPMTX(KOUNT)
         DMATX(JSTRS,ISTRS)=DMATX(ISTRS,JSTRS)
        ENDDO
       ENDDO
      ENDIF
C
C**** MODIFY D-MATRIX ACCORDING TO SELECTED NONLINEAR CONSTITUTIVE MODEL
C
C     GENERAL FORMULATION FOR ELASTO-PLASTICITY-DAMAGE
C     (only for symmetric constitutive tensors)
C
      IF(NALGO.GT.1.AND.NCRIT.NE.0) THEN
       IF(NCRIT.GE.31.AND.NCRIT.LE.42) THEN
        CALL PLBMTX(DMATX,EHIST(IPLAS(2),IGAUS),NSTRS,NSTR1,KSYMM)
       ENDIF     ! ncrit.ge.31.and.ncrit.le.41
      ENDIF      ! nalgo.gt.1.and.ncrit.ne.0
C
C**** COMPUTE THE ADDITIONAL TERMS FOR FINITE DISPLACEMENT
C
      IF(LARGE.NE.0) THEN
C
C**** CALCULATE THE JACOBIAN MATRIX
C
       IF(KPARB.LE.0) THEN
        CALL LINEAR(CARTD(1,1,IGAUS), DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .              GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODL,NSTR1,
     .              PROPS,            SHAPE(1,IGAUS),DUMMY,
     .              DUMMY,DUMMY,
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,DETJB,
     .              ELCOD,
     .                  1)
       ELSE
        CALL LINEAB(BSBAR(1,1,IGAUS), DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .              GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODL,NSTR1,
     .              PROPS,            SHAPE(1,IGAUS),DUMMY,
     .              DUMMY,DUMMY,
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .              ELCOD,CARTD(1,1,IGAUS),NEVAB,
     .              XJACN,XJANI,XJA3N,XJ3NI,DETJN,DETJB,
     .                  1)
       ENDIF
C
       IF(KSTIG.EQ.1) THEN
        DO ISTR1=1,NSTR1
         SIGMA(ISTR1)=STRSG(ISTR1,IGAUS)
        ENDDO
C
        IF(KPARB.LE.0) THEN
         CALL KLARGE(CARTD(1,1,IGAUS), DVOLU(IGAUS),   ESTIF,
     .               GPCOD(1,IGAUS),   KSYMM, NDIME,NDOFN,
     .               NEVAB,NNODL,NTYPE,SHAPE(1,IGAUS), SIGMA)
        ELSE
         CALL KLARGB(CARTD(1,1,IGAUS), DVOLU(IGAUS),   ESTIF,
     .               GPCOD(1,IGAUS),   KSYMM, NDIME,NDOFN,
     .               NEVAB,NNODL,NTYPE,SHAPE(1,IGAUS), SIGMA,
     .             BSBAR(1,1,IGAUS),NSTR1,BMATX,NSTRE,NSTRS,ELDIS,NDOFC,
     .               XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .               XJACN,XJANI,XJA3N,XJ3NI,DETJN,
     .               KPARB)
        ENDIF
       ENDIF            ! kstig.eq.1
      ENDIF             ! large.ne.0
C
C**** CALCULATE THE MATERIAL STIFFNESS MATRIX
C
      IF(KPARB.LE.0) THEN
       CALL KMATRI(BMATX       ,CARTD(1,1,IGAUS), DMATX,
     .             DVOLU(IGAUS),ESTIF,            GPCOD(1,IGAUS),
     .             LARGE,KSYMM, NDIME,NDOFN,NEVAB,NKOVA,NNODL,
     .             NSTRE,NSTRS,NTYPE,SHAPE(1,IGAUS),
     .             NSTR1,NDOFC,ELDIS,AUXS1,
     .             XJACM)
      ELSE
       CALL KMATRB(BMATX,       CARTD(1,1,IGAUS), DMATX,
     .             DVOLU(IGAUS),ESTIF,            GPCOD(1,IGAUS),
     .             LARGE,KSYMM, NDIME,NDOFN,NEVAB,NKOVA,NNODL,
     .             NSTRE,NSTRS,NTYPE,SHAPE(1,IGAUS),
     .             NSTR1,NDOFC,ELDIS,AUXS1,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .             BSBAR(1,1,IGAUS),
     .             XJACN,XJANI,XJA3N,XJ3NI,DETJN)
      ENDIF
C
   10 CONTINUE
C
C**** ADD HOURGLASS CONTROL
C
      IF(NHOUR.NE.0) THEN
       IF((NDIME.EQ.2).AND.((NNODL.NE.4).OR.(NGAUL.NE.1))) RETURN
       IF((NDIME.EQ.3).AND.((NNODL.NE.8).OR.(NGAUL.NE.1))) RETURN
       FMULT=1.0
       IF(KELAS.NE.1) THEN
        IF(NCRIT.GE.20.AND.NCRIT.LE.22) FMULT=1.0-EHIST(11,1)
       ENDIF  
       DO IKOVA=1,NKOVA
        ESTIF(IKOVA)=ESTIF(IKOVA)+FMULT*EMASS(IKOVA)
       ENDDO
      ENDIF
C
      RETURN
      END
