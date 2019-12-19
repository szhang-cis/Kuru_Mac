      SUBROUTINE FRIN30(CARTD,DVOLU,EHIST,ELCOD,ELDIS,EMASS,EPMTX,
     .                  GPCOD,LNODS,PROPS,RMAT1,SHAPE,STRAN,STRSG,
     .                  STRA0,STRS0,THICK,VANIS,
     .                  BMSIG,BMATX,DESIG,DMATX,DSTRA,PRESG,SGTOT,
     .                  SIGMA,TSTRA,XJACM,TENOD,DTENO,
     .                  PWOEL,PREAL,TGAPL,TENOI,BSBAR,XJACN,FPCHL)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING FORCES
C     ( FOR ELEMENT NO. 30 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
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
      DIMENSION CARTD(NDIME,NNODL,*), DVOLU(*),
     .          EHIST(NHIST,*),       ELCOD(NDIME,*),
     .          ELDIS(NDOFC,*),       EMASS(*),
     .          EPMTX(*),             GPCOD(NDIME,*), 
     .          LNODS(*),             PROPS(*),
     .          RMAT1(NDIME,*),       SHAPE(NNODL,*),
     .          STRAN(NSTR1,*),       STRSG(NSTR1,*),
     .          STRA0(NSTR1,*),       STRS0(NSTR1,*), VANIS(NANIV,*)
      DIMENSION BMATX(NSTRE,*),       BMSIG(*),
     .          DESIG(*),             DMATX(NSTRS,*), 
     .          DSTRA(*),             XJACM(NDIME,*), 
     .          PRESG(*),             SGTOT(*), 
     .          SIGMA(*),             TSTRA(*),
     .          TENOD(*),             DTENO(*)
      DIMENSION PWOEL(*),             PREAL(*),
     .          TGAPL(*),             TENOI(*)
      DIMENSION BSBAR(NSTR1,NEVAB,*),
     .          XJACN(NDIME,*),
     .          FPCHL(NFPCH,*)
      DIMENSION XJACI(3,3),           XJANI(3,3)
      DIMENSION STRS1(6)
C
      IF(IEROR.NE.0) RETURN
C
C**** INCOMPRESSIBILITY CONDITION
C
      IF(KPROB.EQ.5) THEN
C
C**** LOOP OVER INTEGRATION POINTS
C
       DO IGAUS=1,NGAUL
        CALL BINCOM(BMATX,CARTD(1,1,IGAUS), ELDIS,GPCOD(1,IGAUS),
     .              LARGE,NDIME,NDOFC,NNODL,NSTRE,NTYPE,SHAPE(1,IGAUS),
     .              XJACM,DMATX,NSTRS)
        CALL LINCOM(BMATX, DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .              GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODL,NSTR1,
     .              PROPS,            SHAPE(1,IGAUS),STRAN(1,IGAUS),
     .              STRA0(1,IGAUS),   TSTRA,            XJACM,SGTOT)
        CALL EQLOAB(BMATX,BMSIG,CARTD(1,1,IGAUS),
     .              DVOLU(IGAUS),   ELDIS,  GPCOD(1,IGAUS),   LARGE,
     .              NDIME,NDOFC,NDOFN,NEVAB,NNODL,NSTRE,      NTYPE,
     .              SHAPE(1,IGAUS), SGTOT,
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .              NSTR1,NSTRS,BMATX,
     .              XJACN,XJANI,XJA3N,XJ3NI,DETJN)
       ENDDO
       RETURN
      ENDIF
C
C**** INITIALISES MECHANICAL COUPLING TERM
C
      IF(ITERME.GT.0) THEN              ! bidirectional coupled problems
       DO INODE=1,NNODL
        PWOEL(INODE)=0.0D0
        PREAL(INODE)=0.0D0
        TGAPL(INODE)=0.0D0
       END DO
      ENDIF
C
C**** LOOK FOR MODEL DEFINITION
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IPEP1=100
      IF(ISOTR.EQ.1) IPEP1=200
      IPEP2=INT(PROPS(2))   ! 10=elastic; 20=plastic; 30=viscoplastic;
C                           ! 40=hyperelastic
      IPEP3=INT(PROPS(3))   ! 1=no temp.-dependent; 2=temp.-depend.
      IPEP4=INT(PROPS(4))   ! standard; non-standard
      IPEPE=IPEP1+IPEP2+IPEP3
C
      KPART=INT(PPART)
      KPARB=INT(PPARB)
      KPARI=INT(PPARI)
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO 100 IGAUS=1,NGAUL
C
C**** SET UP THE ELASTIC MATRIX
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
C**** EVALUATES INCREMENTAL STRAINS & INCREMENTAL ELASTIC STRESSES
C
      IF(KPARB.EQ.0) THEN
       CALL LINEAR(CARTD(1,1,IGAUS), DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .             GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODL,NSTR1,
     .             PROPS,            SHAPE(1,IGAUS),STRAN(1,IGAUS),   
     .             STRA0(1,IGAUS),   TSTRA,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,DETJB,
     .             ELCOD,
     .                 0)
      ELSE
       CALL LINEAB(BSBAR(1,1,IGAUS), DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .             GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODL,NSTR1,
     .             PROPS,            SHAPE(1,IGAUS),STRAN(1,IGAUS),   
     .             STRA0(1,IGAUS),   TSTRA,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .             ELCOD,CARTD(1,1,IGAUS),NEVAB,
     .             XJACN,XJANI,XJA3N,XJ3NI,DETJN,DETJB,
     .                 0)
      ENDIF
C
C**** EVALUATE EFFECTIVE STRESSES FOR LINEAR ELASTICITY
C
      IF(NCRIT.EQ.0) THEN
       DO ISTR1=1,NSTR1
        STRSG(ISTR1,IGAUS)=STRSG(ISTR1,IGAUS)+DESIG(ISTR1)
        IF(NMEMOM.EQ.0) THEN
         SGTOT(ISTR1)=STRSG(ISTR1,IGAUS)
        ELSE
         SGTOT(ISTR1)=STRSG(ISTR1,IGAUS)-STRS0(ISTR1,IGAUS)
        ENDIF
       ENDDO
      ELSE
C
C**** (A) EVALUATES TOTAL 'TRIAL' STRESSES (EFFECTIVE) FOR OTHER
C         CONSTITUTIVE EQUATIONS
C
       DO ISTR1=1,NSTR1
        SGTOT(ISTR1)=0.0D0
        PRESG(ISTR1)=STRSG(ISTR1,IGAUS)
        SIGMA(ISTR1)=PRESG(ISTR1)+DESIG(ISTR1)
        STRS1(ISTR1)=0.0D0
       ENDDO
C
C**** (B) EVALUATES TOTAL 'TRUE' STRESSES (EFFECTIVE)
C
       IF(NCRIT.GE.31.AND.NCRIT.LE.42) THEN
        TEMPG=0.0D0
        DTEMG=0.0D0
        TEINI=0.0D0
#ifndef restricted
        IF(IPEP3.EQ.2)            ! temperature-dependent models
     .   CALL GATECA(TEMPG,DTEMG,TEINI,
     .               TENOD,DTENO,TENOI,
     .               SHAPE(1,IGAUS),KPART)
#endif
        DBOLU=DVOLU(IGAUS)
        IF(INITV.EQ.1) THEN       ! non-standard initial conditions
         IF(NPRE2.GT.0) THEN
          DO ISTR1=1,NSTR1
           STRS1(ISTR1)=STRS0(ISTR1,IGAUS)
          ENDDO
         ENDIF
        ENDIF
C
        IF(IPEP4.EQ.0)            ! standard models
     .   CALL CEPL36(SIGMA,DESIG,EHIST(IPLAS(1),IGAUS),PRESG,
     .               EHIST(IPLAS(3),IGAUS),DMATX,PROPS,SGTOT,
     .               DVOLU,EHIST(IPLAS(2),IGAUS),GPCOD(1,IGAUS),
     .               STRAN(1,IGAUS),DSTRA,
     .               XJACM,XJACI,XJA3M,XJA3I,DETJB,
     .               TEMPG,DTEMG,VANIS,
     .               PWOEL,SHAPE(1,IGAUS),DBOLU,THICK,TEINI,
     .               FPCHL,EHIST(IPLAS(6),IGAUS),
     .               EHIST(IPLAS(4),IGAUS),EHIST(IPLAS(5),IGAUS),
     .               EHIST(IPLAS(7),IGAUS),      ! dim. to be improved
     .               STRS1)
C
#ifndef restricted
        IF(IPEP4.EQ.1)            ! non-standard model: SG cast iron
     .   CALL CEPL51(SIGMA,DESIG,EHIST(IPLAS(1),IGAUS),PRESG,
     .               EHIST(IPLAS(3),IGAUS),DMATX,PROPS,SGTOT,
     .               DVOLU,EHIST(IPLAS(2),IGAUS),GPCOD(1,IGAUS),
     .               STRAN(1,IGAUS),DSTRA,
     .               XJACM,XJACI,XJA3M,XJA3I,DETJB,
     .               TEMPG,DTEMG,
     .               PWOEL,SHAPE(1,IGAUS),DBOLU,THICK,TEINI,
     .               FPCHL,EHIST(IPLAS(6),IGAUS),
     .               EHIST(IPLAS(4),IGAUS),EHIST(IPLAS(5),IGAUS),
     .               EHIST(IPLAS(7),IGAUS),
     .               STRS1)
#endif
C
        IF(IPEP4.EQ.2)            ! non-standard model: green sand
     .   call runend('non-standard model 2 not implemented yet')
C
#ifndef restricted
        IF(IPEP4.EQ.3)            ! non-standard model: dual-phase steel
     .   CALL CEPL53(SIGMA,DESIG,EHIST(IPLAS(1),IGAUS),PRESG,
     .               EHIST(IPLAS(3),IGAUS),DMATX,PROPS,SGTOT,
     .               DVOLU,EHIST(IPLAS(2),IGAUS),GPCOD(1,IGAUS),
     .               STRAN(1,IGAUS),DSTRA,
     .               XJACM,XJACI,XJA3M,XJA3I,DETJB,
     .               TEMPG,DTEMG,VANIS,
     .               PWOEL,SHAPE(1,IGAUS),DBOLU,THICK,TEINI,
     .               FPCHL,EHIST(IPLAS(6),IGAUS),
     .               EHIST(IPLAS(4),IGAUS),EHIST(IPLAS(5),IGAUS),
     .               EHIST(IPLAS(7),IGAUS),      ! dim. to be improved
     .               STRS1,EHIST(IPLAS(8),IGAUS))
#endif
C
       ENDIF       ! ncrit.ge.31.or.ncrit.le.41
C
C**** (C) STORES FINAL STRESSES IN THE PROPER ARRAY
C
       DO ISTR1=1,NSTR1
        STRSG(ISTR1,IGAUS)=SGTOT(ISTR1)
        IF(NMEMOM.EQ.1) THEN
         SGTOT(ISTR1)=STRSG(ISTR1,IGAUS)-STRS0(ISTR1,IGAUS)
        ENDIF
       ENDDO
C
      ENDIF        ! ncrit.eq.0
C
C**** INTEGRATE THE STRESSES INTO THE INTERNAL NODAL FORCES
C
      IF(KPARB.EQ.0) THEN
       CALL EQLOAD(BMATX,      BMSIG,CARTD(1,1,IGAUS), 
     .             DVOLU(IGAUS),   ELDIS,  GPCOD(1,IGAUS),   LARGE,
     .             NDIME,NDOFC,NDOFN,NEVAB,NNODL,NSTRE,      NTYPE,
     .             SHAPE(1,IGAUS), SGTOT,
     .             XJACM)
      ELSE
       CALL EQLOAB(BSBAR(1,1,IGAUS),BMSIG,CARTD(1,1,IGAUS), 
     .             DVOLU(IGAUS),   ELDIS,  GPCOD(1,IGAUS),   LARGE,
     .             NDIME,NDOFC,NDOFN,NEVAB,NNODL,NSTRE,      NTYPE,
     .             SHAPE(1,IGAUS), SGTOT,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .             NSTR1,NSTRS,BMATX,
     .             XJACN,XJANI,XJA3N,XJ3NI,DETJN)
      ENDIF
  100 CONTINUE
C
C**** CORRECTS COUPLING TERM AT CORNER NODES FOR INTERC=1
C
#ifndef restricted
      IF(IPEP3.EQ.2) THEN                 ! temperature-dependent models
       IF(ITERME.GT.0) THEN               ! bidirectional coupled
        IF(ITERMP.GT.0) THEN              ! coupling term
         IF(INTERC.EQ.1) THEN
          CALL SMOESQ(PWOEL,NDIME,NNODL)
         ENDIF
        ENDIF
       ENDIF
      ENDIF
#endif
C
C**** ADD HOUGLASS CONTROL
C
      IF(NHOUR.NE.0) THEN
       IF((NDIME.EQ.2).AND.((NNODL.NE.4).OR.(NGAUL.NE.1))) RETURN
       IF((NDIME.EQ.3).AND.((NNODL.NE.8).OR.(NGAUL.NE.1))) RETURN
       FMULT=1.0
       IF(KELAS.NE.1) THEN
        IF(NCRIT.GE.20.AND.NCRIT.LE.22) FMULT=1.0-EHIST(11,1)
       ENDIF 
       CALL HOURFC(BMSIG,ELDIS,EMASS,FMULT,NEVAB) 
      ENDIF
C
      RETURN
      END
