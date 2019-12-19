      SUBROUTINE FRIN01(CARTD,DVOLU,EHIST,ELCOD,ELDIS,EMASS,EPMTX,
     .                  GPCOD,LNODS,PROPS,RMAT1,SHAPE,STRAN,STRSG,
     .                  STRA0,STRS0,THICK,
     .                  BMSIG,BMATX,DESIG,DMATX,DSTRA,PRESG,SGTOT,
     .                  SIGMA,TSTRA,XJACM,TENOD,DTENO,PWOEL,PREAL,
     .                  TGAPL)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING FORCES OF THE SOLID PHASE
C               ( FOR ELEMENT NO. 1 )
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
      DIMENSION CARTD(NDIME,NNODE,*), DVOLU(*),
     .          EHIST(NHIST,*),       ELCOD(NDIME,*),
     .          ELDIS(NDOFC,*),       EMASS(*),
     .          EPMTX(*),             GPCOD(NDIME,*), 
     .          LNODS(*),             PROPS(*),
     .          RMAT1(NDIME,*),       SHAPE(NNODE,*),
     .          STRAN(NSTR1,*),       STRSG(NSTR1,*),
     .          STRA0(NSTR1,*),       STRS0(NSTR1,*)
      DIMENSION BMATX(NSTRE,*),       BMSIG(*),
     .          DESIG(*),             DMATX(NSTRS,*),
     .          DSTRA(*),             XJACM(NDIME,*), 
     .          PRESG(*),             SGTOT(*), 
     .          SIGMA(*),             TSTRA(*),
     .          TENOD(*),             DTENO(*)
      DIMENSION PWOEL(*),             PREAL(*),
     .          TGAPL(*)
C
C**** INITIALISE MECHANICAL COUPLING TERM
C
      DO INODE=1,NNODE
       PWOEL(INODE)=0.0
       PREAL(INODE)=0.0
       TGAPL(INODE)=0.0
      END DO
C
C***UPDATE GEOMETRY FOR LARGE DISPLACEMENTS
C
      IF(LARGE.EQ.1) THEN
        DO IDIME=1,NDIME
          DO INODE=1,NNODL
             ELDIS(IDIME,INODE)=ELCOD(IDIME,INODE)+ELDIS(IDIME,INODE)
          ENDDO
        ENDDO
      ENDIF
C
C***LOOP OVER INTEGRATION POINTS
C
      DO 100 IGAUS=1,NGAUL
C
C***EVALUATE INCREMENTAL STRAINS & INCREMENTAL ELASTIC STRESSES
C
C
C***SET UP THE ELASTIC MATRIX
C
      KOUNT=0
      DO 15 ISTRS=1,    NSTRS
      DO 15 JSTRS=ISTRS,NSTRS
      KOUNT=KOUNT+1
      DMATX(ISTRS,JSTRS)=EPMTX(KOUNT)
   15 DMATX(JSTRS,ISTRS)=DMATX(ISTRS,JSTRS)
C
      CALL LINEAR(CARTD(1,1,IGAUS), DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .            GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODL,NSTR1,
     .            PROPS,            SHAPE(1,IGAUS),STRAN(1,IGAUS),   
     .            STRA0(1,IGAUS),   TSTRA,            XJACM)
C
C***EVALUATE EFFECTIVE STRESSES
C
C
C     LINEAR ELASTICITY
C
      IF(NCRIT.EQ.0) THEN                 
        DO 40 ISTR1=1,NSTR1
        STRSG(ISTR1,IGAUS)=STRSG(ISTR1,IGAUS)+DESIG(ISTR1)
   40   SGTOT(ISTR1)=STRSG(ISTR1,IGAUS)-STRS0(ISTR1,IGAUS)
      ELSE
C
C     OTHER CONSTITUTIVE EQUATIONS
C
C      (a)  evaluate total 'trial' stresses ( effective )
C
        DO 50 ISTR1=1,NSTR1
        SGTOT(ISTR1)=0.0D00
        PRESG(ISTR1)=STRSG(ISTR1,IGAUS)
   50   SIGMA(ISTR1)=PRESG(ISTR1)+DESIG(ISTR1)
C
C      (b)  evaluate total 'true' stresses  ( effective )
C
c      IF(NCRIT.GT.0.AND.NCRIT.LE.4)
c     .  CALL CEPL14(DEPSV,DESIG,DMATX,EHIST(1,IGAUS),IGAUS,KUNLD,
c     .              NBACK,NSTR1,PRESG,PROPS,         SGTOT,SIGMA)
c      IF(NCRIT.EQ.5)
c     .  CALL CECRA5(DVOLU(IGAUS),EHIST(1,IGAUS),ELCOD,IGAUS,NDIME,
c     .              NNODE,NSTR1 ,PROPS,         SGTOT,SIGMA,TSTRA,
c     .              THICK)
c      IF(NCRIT.EQ.10)
c     .  CALL CEJN10(DMATX,ELCOD,EHIST(1,IGAUS),EHIST(2,IGAUS),
c     .              IGAUS,NDIME,NNODE,NSTR1,PRESG,PROPS,RMAT1,
c     .              SGTOT,SIGMA,TSTRA,EHIST(3,IGAUS),ehist(4,igaus),
c     .              tenod(3),tenod(4),dstra)
c
c*  tenod(3) and tenod(4) define the element state
c
c      IF(NCRIT.GE.20.AND.NCRIT.LE.30)
c     .   CALL CEDM20(DESIG,DMATX,DVOLU,EHIST(1,IGAUS),IGAUS,
c     .               NDIME,NNODE,NSTR1,PRESG,PROPS,SGTOT,
c     .               SIGMA,TSTRA,ELCOD,GPCOD(1,1),THICK)
      IF(NCRIT.GE.31.AND.NCRIT.LE.38)THEN
c       III59=INT(PROPS(59))
C
c      CALL IDENPR(PROPS,
c     .            NYOUN,VYOUN,NPOIS,VPOIS,NALPH,VALPH,NCCER,VCCER,
c     .      IVIFL,NVISC,VVISC,NEXPO,VEXPO,
c     .      ISOTT,NCCOE,VCCOE,NCCOB,VCCOB,NCCOQ,VCCOQ,
c     .            TEREF,TEMRO,TEMPL)
C
c       IF(III59.EQ.3) CALL GATECA(TEMPG,DTEMG,SHAPE(1,IGAUS),TENOD,
c     .                                  DTENO,TEMPL,KPART)
C
       DBOLU=DVOLU(IGAUS)
c         CALL CEPL35(SIGMA,DESIG,EHIST(1,IGAUS),PRESG,EHIST(7,IGAUS),
c     .               DMATX,PROPS,SGTOT,DVOLU,EHIST(18,IGAUS),
c     .               GPCOD(1,1),STRAN(1,IGAUS),DSTRA,XJACM,
c     .               TEMPG,DTEMG,EHIST(39,IGAUS),PWOEL,SHAPE(1,IGAUS),
c     .               DBOLU,
c     .               EHIST(45,IGAUS),EHIST(46,IGAUS),EHIST(47,IGAUS),
c     .               EHIST(48,IGAUS),EHIST(49,IGAUS),EHIST(50,IGAUS),
c     .               THICK)
      ENDIF
C
C      (c)  store final stresses in the proper array
C
        DO 60 ISTR1=1,NSTR1
        STRSG(ISTR1,IGAUS)=SGTOT(ISTR1)
   60   SGTOT(ISTR1)=STRSG(ISTR1,IGAUS)-STRS0(ISTR1,IGAUS)
C
      ENDIF
C
C***TRANSFORMA EL ESTADO TENSIONAL DE:   PIOLA ---> CAUCHY
C                  (solo en grandes deformaciones)
C      IF(LARGE.EQ.1) CALL CAUCHY(CARTD,ELDIS,STRSG,XJACM)
C
C***CALCULATE THE TOTAL STRESS AT GAUSSIAN POINT
C
ctm       IF(KPORE.EQ.2)
ctm     .    CALL TOTSIG(ELDIS,KPORE,NDOFC,NNODL,NSTR1,PROPS,SGTOT,
ctm     .                SHAPE(1,IGAUS))
C
C***INTEGRATE THE STRESSES INTO THE INTERNAL NODAL FORCES
C
      CALL EQLOAD(BMATX,          BMSIG,  CARTD(1,1,IGAUS), 
     .            DVOLU(IGAUS),   ELDIS,  GPCOD(1,IGAUS),   LARGE,
     .            NDIME,NDOFC,NDOFN,NEVAB,NNODL,NSTRE,      NTYPE,
     .            SHAPE(1,IGAUS), SGTOT,  XJACM)
  100 CONTINUE
C
C***ADD HOUGLASS CONTROL
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
