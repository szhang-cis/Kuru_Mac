      SUBROUTINE MASS32(DVOLU,PROPS,LNODS,SHAPE,WSTIF,EHIST,ELDIS,TENOD,
     .                  DISPL,ELCOD,VNORL,DMATX,AUXS1,AUXS2,BMATX,VTANL)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS MATRIX FOR THE GAP ELEMENT
C     ( ELEMENT NO. 32 )
C
C
C     Note: SHAPE, ELDIS, TENOD & DISPL must be dimensioned
C           with NNODN due to non-coincident contact mesh (NOCOL=1)
C
C***********************************************************************
C
C**** PROPERTIES FOR THE CONTACT PROBLEM
C
C     PROPS(2)=RIGIN            : NORNAL STIFFNESS
C     PROPS(5)=TEMPL            : LIQUIDUS TEMPERATURE
C     PROPS(6)=INOTE            : DENOTES WHICH NODE BELONGS TO PIECE
C
C
C     REGULARIZATION OF CONTACT JACOBIAN MATRIX: 3 FORMS
C
C     Note: if two or more forms are input, the one with higher
C           numbering will be considered
C
C     1) JACOBIAN_REGULARIZATION_FACTORS (DEFINED FOR EACH CONTACT
C        MATERIAL FOR THE WHOLE ANALYSIS)
C
C        FIXED PARAMETER:
C        IITEF=10
C
C        PROPS(10)=INP10: INDEX FOR "JACOBIAN_REGULARIZATION_FACTORS"
C        READS: TRUPL, TRUPM, NNN20
C        PROPS(11)=TRUPL: FIRST JACOB. REGULAR. FACTOR
C        PROPS(12)=TRUPM: SECOND JACOB. REGULAR. FACTOR
C        PROPS(13)=NNN20: DO LOOP INDEX
C
C        PROPS(14)=INP14: INDEX FOR "MIN_&_MAX_ADMISSIBLE_CONTACT_GAPS"
C        READS: TOLGA, TOLGAM
C        PROPS(15)=TOLGA: MINIMUM ADMISSIBLE CONTACT GAP
C        PROPS(16)=TOLGAM: MAXIMUM ADMISSIBLE CONTACT GAP
C
C        DEFAULTS VALUES (defined in inpset.f):
C        IITEF=10
C        TRUPL=1.0E-06
C        TRUPM=0.0
C        TOLGA=-1.0E-08
C        TOLGAM=1.0E+20
C        DEFAULT VALUE (defined in prop04.f):
C        NNN20=20
C
C
C     2) TUNING_PARAMETERS (DEFINED FOR EACH CONTACT SET AND THEY CAN 
C        BE INPUT AT MATERIAL PROPERTIES AND/OR INTERVAL DATA LEVEL)
C
C        FIXED PARAMETER:
C        NNN20=20
C
C        IITEN: INDEX FOR TUNING PARAMETERS
C        READS: IITEF, TRUPL, TRUPM, TOLGA, TOLGAM
C               (COMING FROM A COMMON OF TUNING PARAMETERS)
C        IITEF: ITERATION FROM WHICH TRUPL ACTS
C        TRUPL: FIRST JACOB. REGULAR. FACTOR
C        TRUPM: SECOND JACOB. REGULAR. FACTO
C        TOLGA: MINIMUM CONTACT ADMISIBLE GAP
C        TOLGAM: MAXIMUM CONTACT ADMISSIBLE GAP
C
C        DEFAULTS VALUES (defined in inpset.f):
C        IITEF=10
C        TRUPL=1.0E-06
C        TRUPM=0.0
C        TOLGA=-1.0E-08
C        TOLGAM=1.0E+20
C        DEFAULT VALUE (defined in mass32.f):
C        NNN20=20
C
C
C     3) CONTACT_CONTROL_FACTORS (THE SAME FOR ALL CONTACT SETS AND
C        THEY HAVE TO BE DEFINED FOR EVERY TIME STEP)
C
C        FIXED PARAMETERS:
C        IITEF=10
C        NNN20=20
C        TOLGA=-1.0E-08
C        TOLGAM=1.0E+20
C
C        ICONC: CONTACT CONTROL INDEX
C        READS: TRUPL, TRUPM
C        TRUPL: FIRST JACOB. REGULAR. FACTOR
C        TRUPM: SECOND JACOB. REGULAR. FACTOR
C
C        DEFAULTS VALUES (defined in inpset.f):
C        IITEF=10
C        TRUPL=1.0E-06
C        TRUPM=0.0
C        TOLGA=-1.0E-08
C        TOLGAM=1.0E+20
C        DEFAULT VALUE (defined in mass32.f):
C        NNN20=20
C
C***********************************************************************
C
C**** PROPERTIES FOR THE LINKING PROBLEM
C
C     PROPS(2-4)=RIGIN          : STIFFNESSES
C     PROPS(5)=TEMPL            : LIQUIDUS TEMPERATURE
C     PROPS(6)=INOTE            : DENOTES WHICH NODE BELONGS TO PIECE
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
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM
C
      DIMENSION DVOLU(*),       PROPS(*),
     .          SHAPE(NNODN,*), WSTIF(*),
     .          EHIST(NHIST,*), ELDIS(NDOFC,*),
     .          TENOD(*),       DISPL(NDOFN,*),
     .          ELCOD(NDIME,*), VNORL(NDIME,*),
     .          LNODS(*)
C
      DIMENSION AUXS1(NEVAB,*), AUXS2(NEVAB,*),
     .          DMATX(NDIME,*), BMATX(NDIME,*), ! NDIME instead of NSTR1
     .          VTANL(NDIME,NDIME-1,*)
C
C**** DETERMINES BETWEEN CONTACT OR LINKING
C
      ICOLI=INT(PROPS(1))
      IF(ICOLI.EQ.1) THEN
C
       NNOBO=NNODL/2
       IF(NOCOL.EQ.1) NNOBO=NNODL                  ! non-coincident mesh
       NEVBO=NNOBO*NDOFN
C
       RIGIN=PROPS(2)
       TEMPL=PROPS(5)
       INOTE=INT(PROPS(6))
       ICOMO=INT(PROPS(7))
       IITEN=INT(RITEN)     ! 0: tun. not considered; 1: tun. considered
       IITEF=INT(RITEF)
C
C**** DEFINES JACOBIAN REGULARIZATION FACTORS & MINIMUM & MAXIMUM
C     ADDMISSIBLE CONTACT GAPS
C
       NNN20=20
       IF(ICONC.EQ.0) THEN
        IF(IITEN.EQ.0) THEN
         INP10=INT(PROPS(10))
         IF(INP10.EQ.1) THEN
          IITEF=10
          TRUPL=PROPS(11)
          TRUPM=PROPS(12)
          NNN20=INT(PROPS(13))
         ENDIF
         INP14=INT(PROPS(14))
         IF(INP14.EQ.1) THEN
          TOLGA=PROPS(15)
          TOLGAM=PROPS(16)
         ENDIF
        ENDIF
       ENDIF
C
       IFRIC=INT(PROPS(51))
       IF(IFRIC.EQ.1) THEN
        IFRIM=INT(PROPS(52))
        IF(IFRIM.EQ.1) THEN
         IFRCO=INT(PROPS(53))
         IF(IFRCO.EQ.1) THEN
          RIGIT=PROPS(54)
          FRIMU=PROPS(55)
         ENDIF
        ENDIF
        IF(IFRIM.EQ.2) THEN
         RIGIT=PROPS(54)
        ENDIF
       ENDIF
C
C**** LOOP ON INTEGRATION POINTS
C
       DO 100 IGAUS=1,NGAUL
C
       IF(NOCOL.EQ.1) THEN                              ! checks contact
        IAUXY=(IGAUS-1)*NNODL
        IF(LNODS(NNOBO).EQ.LNODS(NNODL+NNOBO+IAUXY)) GO TO 100
       ENDIF
C
C**** COMPUTES NORMAL GAP (CURRENT)
C
C     DGAPN  > 0: CONTACT PROBLEM
C     DGAPN <= 0: NO CONTACT PROBLEM
C
       DGAPN=0.0D0
       DO IDOFN=1,NDOFN
        DO INODL=1,NNOBO
         ELDI1=ELDIS(IDOFN,INODL)
         IF(NOCOL.EQ.0) THEN
          ELDI2=ELDIS(IDOFN,INODL+NNOBO)
          DGAPN=DGAPN+SHAPE(INODL,IGAUS)*(ELDI1-ELDI2)*
     .                                                VNORL(IDOFN,IGAUS)
         ELSE
          ELDI2=ELDIS(IDOFN,INODL+NNOBO+IAUXY)
          DGAPN=DGAPN+(SHAPE(INODL,IGAUS)*ELDI1-
     .                 SHAPE(INODL+NNOBO,IGAUS)*ELDI2)*
     .                                                VNORL(IDOFN,IGAUS)
         ENDIF
        END DO
       END DO
C
C**** COMPUTES NORMAL GAP (LAST CONVERGED)
C
       DGAPA=0.0D0
       DO IDOFN=1,NDOFN
        DO INODL=1,NNOBO
         ELDI1=ELDIS(IDOFN,INODL)-DISPL(IDOFN,INODL)
         IF(NOCOL.EQ.0) THEN
          ELDI2=ELDIS(IDOFN,INODL+NNOBO)-DISPL(IDOFN,INODL+NNOBO)
          DGAPA=DGAPA+SHAPE(INODL,IGAUS)*(ELDI1-ELDI2)*
     .                                                VNORL(IDOFN,IGAUS)
         ELSE
          ELDI2=ELDIS(IDOFN,INODL+NNOBO+IAUXY)-
     .                                    DISPL(IDOFN,INODL+NNOBO+IAUXY)
          DGAPA=DGAPA+(SHAPE(INODL,IGAUS)*ELDI1-
     .                 SHAPE(INODL+NNOBO,IGAUS)*ELDI2)*
     .                                                VNORL(IDOFN,IGAUS)
         ENDIF
        END DO
       END DO
C
C**** COMPUTES TEMPERATURE
C
       TGAUS=0.0D0
       IF(ITERME.GE.0) THEN                 ! coupled problems
        IF(INOTE.EQ.1) THEN
         DO INODL=1,NNOBO
          TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL)
         END DO
        ELSE
         DO INODL=1,NNOBO
          IF(NOCOL.EQ.0) THEN
           TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL+NNOBO)
          ELSE
           CALL RUNEND('ERROR: INOTE=2 NOT IMPLEMENTED YET (mass32.f)')
          ENDIF
         END DO
        ENDIF
       ENDIF
C
C**** COMPUTES CONTACT RIGIDITY
C
       IF(ICOMO.EQ.1) THEN
        RIGIN=2.0D0*RIGIN*DABS(DGAPN)
       ENDIF
       IF(ICOMO.EQ.2) THEN
        DGABA=PROPS(8)
        IF(DGAPN.GT.0.0D0.AND.DGAPN.LT.DGABA) THEN
         RIGIN=RIGIN/DGABA*DGAPN
        ENDIF
       ENDIF
       IF(ICOMO.EQ.3) THEN
        RIGI0=PROPS(9)
        RIGIN=2.0D0*RIGIN*DABS(DGAPN)+RIGI0
       ENDIF
       IF(ICOMO.EQ.4) THEN
        DGABA=PROPS(8)
        RIGI0=PROPS(9)
        IF(DGAPN.GT.0.0D0.AND.DGAPN.LT.DGABA) THEN
         RIGIN=(RIGIN-RIGI0)/DGABA*DGAPN+RIGI0
        ENDIF
       ENDIF
       IF(ICOMO.EQ.5) THEN                     ! mode I debounding model
        IF(DGAPN.GT.0.0D0) THEN                ! compression
         RIGIN=PROPS(2)                        ! K_n_compression
        ELSE
         RIGIN=PROPS(17)*PROPS(17)/(2.0D0*PROPS(9)*PROPS(8))   ! K_n
         DAMAN=EHIST(NDIME+NDIME*NDIME,IGAUS)                  ! d_n
         RIGIN=(1.0D0-DAMAN)*RIGIN
        ENDIF
       ENDIF
C
C**** "STABILIZED" CONTACT RIGIDITY
C
       RIGIX=RIGIN
       IF(TGAUS.LT.TEMPL) THEN
        IF(DGAPN.LE.0.0D0) THEN
         DO III20=1,NNN20
          IF(IITER.GT.(III20*IITEF)) RIGIX=RIGIX*TRUPL
         ENDDO
         IF(DGAPA.LE.TOLGA) RIGIX=RIGIX*TRUPM
        ENDIF      ! dgapn.le.0.0
       ENDIF       ! tgaus.lt.templ
C
C**** INITIALISES DMATX
C
       DO IDIME=1,NDIME
        DO JDIME=1,NDIME
         DMATX(IDIME,JDIME)=0.0D+00
        ENDDO
       ENDDO
C
C**** COMPUTES CONTACT CONSTITUTIVE TENSOR
C
       KDIME=NDIME-1
       DO IDIME=1,NDIME
        PROYI=VNORL(IDIME,IGAUS)
        DO JDIME=1,NDIME
         PROYJ=VNORL(JDIME,IGAUS)
         DMATX(IDIME,JDIME)=RIGIX*PROYI*PROYJ
         IF(IFRIC.EQ.1) THEN            ! "stabilized" friction rigidity
          KDIME=KDIME+1
          DMATX(IDIME,JDIME)=DMATX(IDIME,JDIME)+
     .                                    RIGIX/RIGIN*EHIST(KDIME,IGAUS)
         ENDIF
        ENDDO
       ENDDO
C
       IF(NOCOL.EQ.0) THEN
        CALL WBOUND(NDIME,NEVBO,NNOBO,
     .              SHAPE(1,IGAUS),WSTIF,DMATX,BMATX,
     .              KSYMM,AUXS1)
C
        CALL WINTERT(DVOLU(IGAUS),NEVAB,NEVBO,WSTIF,KSYMM,
     .               AUXS1,AUXS2)
       ELSE
        CALL WBOUNC(DVOLU(IGAUS),NDIME,NDOFN,NEVBO,NNOBO,PROPS,
     .              SHAPE(1,IGAUS),WSTIF,DMATX,BMATX,ELDIS,NDOFC,KSYMM,
     .              AUXS1,NEVAB,AUXS2,IAUXY,LNODS)
       ENDIF
C
  100  CONTINUE
C
      ELSE              ! icoli=2 (linking)
C
       NNOBO=NNODL/2
       NEVBO=NNOBO*NDOFN
C
C**** LOOP ON INTEGRATION POINTS
C
       DO 101 IGAUS=1,NGAUL
C
C**** PROPS(1+IDOFN)=PENALTIES PARAMETERS IN X-Y-Z DIRECTIONS
C     PROPS(5)=LIQUIDUS TEMPERATURE
C
       TEMPL=PROPS(5)
       INOTE=INT(PROPS(6))
C
C**** COMPUTES TEMPERATURE
C
       TGAUS=0.0D0
       IF(ITERME.GE.0) THEN                 ! coupled problems
        IF(INOTE.EQ.1) THEN
         DO INODL=1,NNOBO
          TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL)
         END DO
        ELSE
         DO INODL=1,NNOBO
          TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL+NNOBO)
         END DO
        ENDIF
       ENDIF
C
C**** INITIALISES DMATX
C
       DO IDIME=1,NDIME
        DO JDIME=1,NDIME
         DMATX(IDIME,JDIME)=0.0D+00
        ENDDO
       ENDDO
C
C**** COMPUTES LINKING CONSTITUTIVE TENSOR
C
       IF(TGAUS.GT.TEMPL) THEN
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          IF(IDIME.EQ.JDIME) DMATX(IDIME,JDIME)=PROPS(1+IDIME)
         ENDDO
        ENDDO
       ENDIF                       ! tgaus.gt.templ
C
       CALL WBOUND(NDIME,NEVBO,NNOBO,
     .             SHAPE(1,IGAUS),WSTIF,DMATX,BMATX,
     .             KSYMM,AUXS1)
C
       CALL WINTERT(DVOLU(IGAUS),NEVAB,NEVBO,WSTIF,KSYMM,
     .              AUXS1,AUXS2)
C
  101  CONTINUE
C
      ENDIF                      ! icoli.eq.1
C
      RETURN
      END
