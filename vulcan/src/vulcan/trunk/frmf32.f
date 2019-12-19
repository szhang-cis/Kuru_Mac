      SUBROUTINE FRMF32(PROPS,LNODS,ELDIS,BMSIG,DVOLU,SHAPE,EHIST,TENOD,
     .                  ELCOD,PWOEL,PREAL,TGAPL,VNORL,VTANL,DISIL,DISPL,
     .                  POSGP,WEIGP,XJACM,ELCO2,DVOLD,GPCDD,SHAPD,DERID,
     .                  STRSG,THICK)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE CONTACT FORCES FOR ELEMENT 32 WHEN
C     CONTACT ELEMENT SUBDIVISION IS CONSIDERED
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
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
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM,PSUBC
C
      DIMENSION PROPS(*),               ELDIS(NDOFC,*),
     .          BMSIG(*),               DVOLU(*),
     .          SHAPE(NNODN,*),         EHIST(NHIST,*)
      DIMENSION TENOD(*),               ELCOD(NDIME,*),
     .          PWOEL(*),               PREAL(*),
     .          TGAPL(*),               LNODS(*)
      DIMENSION VNORL(NDIME,*),
     .          VTANL(NDIME,NDIME-1,*),
     .          DISIL(NDOFC,*),         DISPL(NDOFC,*)
      DIMENSION POSGP(NDIME,*),         WEIGP(*),
     .          XJACM(NDIME,*),         ELCO2(NDIME,*)
      DIMENSION DERID(NDIML,NNODN,*),   DVOLD(*),
     .          GPCDD(NDIME,*),         SHAPD(NNODN,*)
      DIMENSION STRSG(NSTR1,*)
C
      NNOBO=NNODL/2
      IF(NOCOL.EQ.1) NNOBO=NNODL                   ! non-coincident mesh
      NEVBO=NNOBO*NDOFN
C
C**** PROPERTIES FOR THE CONTACT PROBLEM
C
      RIGIN=PROPS(2)
      TEMPL=PROPS(5)
      INOTE=INT(PROPS(6))
      ICOMO=INT(PROPS(7))
      IITEN=INT(RITEN)      ! 0: tun. not considered; 1: tun. considered
      NSUBC=INT(PSUBC)
C
      IF(ICONC.EQ.0) THEN
       IF(IITEN.EQ.0) THEN
        INP17=INT(PROPS(17))
        IF(INP17.EQ.1) THEN
         NSUBC=INT(PROPS(18))
        ENDIF
       ENDIF
      ENDIF
cc    IF(NSUBC.EQ.1) RETURN
C
C**** INITIALISES MECHANICAL COUPLING TERM (only for friction problems)
C
      IF(ITERME.GT.0) THEN   ! bidirectional coupling
       DO INODL=1,NNODL
        PWOEL(INODL)=0.0D0
       END DO
      ENDIF
C
C**** PERFORMS THE CONTACT ELEMENT SUBDIVISION
C
      IF(NDIML.EQ.1) THEN
       IF(NNOBO.EQ.2) THEN
        nr=nsubc
        dr=2.00D0/nr
        rir=-1.00D0
        rirm=-1.00D0
        do ir=1,nr
         rir=rirm
         rirm=rir+dr
         r1=rir
         r2=rirm
         frrrc=0.50D0*(r1+r2)
         frrr1=0.50D0*(r2-r1)
         fdvol=frrr1
C
         CALL SHAP01(DVOLD,ELCO2,GPCDD,LNODS,PROPS,SHAPD,THICK,
     .               DERID,POSGP,WEIGP,XJACM,VNORL,STRSG,VTANL,
     .               frrrc,frrr1,fdvol)
C
         CALL FRIN32(PROPS,LNODS,ELDIS,BMSIG,DVOLD,SHAPD,EHIST,TENOD,
     .               ELCOD,PWOEL,PREAL,TGAPL,VNORL,VTANL,DISIL,DISPL,
     .                   1)
        enddo
       ENDIF                 ! nnobo.eq.2
      ENDIF                  ! ndiml.eq.1
C
      IF(NDIML.EQ.2) THEN
       call runend('error: contact element subdiv. not impl. in 3D')
      ENDIF                  ! ndiml.eq.2
C
      IF(NOCOL.EQ.0) THEN
       DO IEVAB=1,NEVBO
        BMSIG(IEVAB+NEVBO)=-BMSIG(IEVAB)
       END DO
      ENDIF                         ! nocoi.eq.0
C
C**** COMPUTES NODAL NORMAL GAP & PRESSURE TO CONVECTION-RADIATION
C     COEFFICIENT OR JUST TO PRINT (see outgap.f) ONLY FOR
C     COINCIDENT MESHES (NOCOL=0)
C     (ASSUMPTION: VNORL IS THE SAME FOR EVERY GAUSS POINT)
C
      IF(NOCOL.EQ.0) THEN
       DO INODL=1,NNOBO
        DGAPN=0.0D0
        DO IDOFN=1,NDOFN
         ELDI1=ELDIS(IDOFN,INODL)
         ELDI2=ELDIS(IDOFN,INODL+NNOBO)
         DGAPN=DGAPN+(ELDI1-ELDI2)*VNORL(IDOFN,1)
        ENDDO
C
        IF(TGAUS.LT.TEMPL) THEN
         IF(DGAPN.GE.0.0D0) THEN    ! contact
          PRESN=RIGIN*DGAPN         ! simplification for ICOMO=0
          DGAPN=0.0D0
         ELSE                       ! a normal gap is produced
          PRESN=0.0D0
          DGAPN=-DGAPN
         ENDIF
C
         TGAPL(INODL)=DGAPN
         TGAPL(INODL+NNOBO)=DGAPN
         PREAL(INODL)=PRESN
         PREAL(INODL+NNOBO)=PRESN
        ELSE                        ! always in contact
         PRESN=RIGIN*DGAPN          ! simplification for ICOMO=0
         DGAPN=0.0D0
         IF(PRESN.LT.0.0D0) PRESN=-PRESN
C
         TGAPL(INODL)=0.0D0
         TGAPL(INODL+NNOBO)=0.0D0
         PREAL(INODL)=PRESN
         PREAL(INODL+NNOBO)=PRESN
        ENDIF
       END DO
      ENDIF                         ! nocoi.eq.0
C
      RETURN
      END
