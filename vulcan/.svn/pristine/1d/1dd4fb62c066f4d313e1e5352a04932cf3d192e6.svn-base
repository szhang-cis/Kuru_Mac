      SUBROUTINE PLVARI(PROPS,ALENG,SIGMA,DEEPP,DEEPL,CAPAP,
     .                  VECTG,NSTRS,ANGFI,ANGSI,
     .                  DFTEQ,DFAFI,YIELD,CAPAD,NSTR1,YIELP,
     .                  COEPD,VADEG,DMATP,COUTH,
     .                  DMATX,DMATE,RETEN,DM,CCERO,NCRIP,EFFPL,EFFP1,DP,
     .                  CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,DISTM,
     .                  BACKA,BACMA,BACKS,BACMS,STRAP,STRAL,ISOTT,IKINE,
     .                  POROS,CPOEF,CPOE1,CPOE2,IPORO,DESIG,DESIL,
     .                  A2COE,A3COE,PREYA,DTIME,DAMAG,EFFP2,
     .                  LARGE,IFREN,NDIME,XJACM,XJA3M,XJACI,XJA3I,IREKI,
     .                  IRECR,RECRY,SRECR,ERECR,ASREC,ALREC,TRMTX,TLAMD,
     .                  NAUXI)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES:
C
C     THE INTERNAL VARIABLES:
C
C               CAPAP: plastic hardening func. in solidification prob.
C               EFFPL: effective plastic strain
C               ANGFI: internal friction
C               ANGSI: dilatancy
C
C
C     Notes:
C
C     For plastic & viscoplastic models:
C     PLASTIC HARDENING COEFFICIENT: CCOEF (h_Cp)
C
C     For plastic models:
C     PLASTIC HARDENING MODULUS: HARDS (H_Cp)     (see plhard.f)
C
C     When performing the product DEEPL:DEEPL, the fact that factor "2"
C     is included in VECTG must be taken into account.
C
C
C
C     MODELS FOR THE COMPUTATION OF THE PLASTIC HARDENING FUNCTION
C
C     ISOTT defines the isotropic hardening model considered:
C
C     MODEL 1 (work hardening; i.e., thesis):
C             CCOEF is the hardening coefficient (h_Cp in thesis)
C             Evolution equation: \dot Cp = \dot lambda h_Cp \sigma : R
C
C     MODEL 2 (linear strain hardening):
C             CCOEF is the hardening coefficient (h_Cp in thesis)
C             Evolution equation: \dot Cp = h_Cp dot \barEp where
C             \barEp=effec. pl. str.
C
C     MODEL 3 (V2 version):
C             CCOEB & CCOEQ: b & Q parameters of V2 viscoplastic model
C             Closed form: Cp = Q (1-exp(-b \barEp)) where
C             \barEp=effec. pl. str.
C
C     MODEL 4 (linear version of MODEL 3)
C             CCOEF is the hardening coefficient (secant)
C             Closed form: Cp = h_Cp \barEp
C
C     MODEL 5 (linear strain hardening):
C             CCOEF is the hardening coefficient (h_Cp in thesis)
C             Evolution equation: Cp = h_Cp \barEp where
C             \barEp=effec. pl. str. computed without factor 2/3 & sqrt
C             (idem MODEL=4 with other definition for \barEp)
C             This model is not consistent with the Von Mises yield
C             function
C
C     MODEL 6 (non-linear strain hardening; non-linear version of
C             MODEL 2):
C             CCOEB & CCOEQ: n & A parameters
C             Closed form: Cp = A \barEp^n where \barEp=effec. pl. str.
C
C     MODEL 7 (non-linear strain hardening; very similar to MODEL=6):
C             CCOEB & CCOEQ: n & A parameters
C             Closed form: Cp = C = A (\barEp_o+\barEp)^n where
C             \barEp=effec. pl. str. and A*\barEp_o^n=Cth
C
C     MODEL 8 (model 3 + model 4):
C             CCOEF, CCOEB & CCOEQ: h_Cp, b & Q parameters
C             Closed form: Cp = h_Cp \barEp + Q (1-exp(-b \barEp)) where
C             \barEp=effec. pl. str.
C
C     MODEL 9 Closed form: C = f(\barEp,T), where f is an explicit
C             function of the effective plastic strain and the
C             temperature
C
C     MODEL 10 Johnson & Cook hardening model (strain, strain rate &
C              temperature-dependent)
C              CCOEF, CCOEB & CCOEQ: H, n & A parameters
C              Closed form: Cp = C = (Cth + A \barEp^n)*
C                                                  (1 + H ln \dot\barEp)
C              where \barEp=effec. pl. str.
C                    \dot\barEp=effec. pl. str. rate
C              The temperature-dependency is included in Cth & A
C              H & n are normally temperature-independent
C
C     MODEL 11 Idem MODEL 10 with the strain hardening term of MODEL 7
C
C     MODEL 12 Idem MODEL 7 with the consideration of a critical
C              effective plastic strain in order to model the effect
C              of plastification without hardening which is, e.g.,
C              a typical phenomenon of low carbon steels (Luders band
C              formation)
C              CCOEF, CCOEB & CCOEQ: \barEp_c, n & A parameters
C              Closed form: Cp = C = A (\barEp_o+<\barEp-\barEp_c>)^n
C              where \barEp=effec. pl. str., \barEp_c=critical effec.
C              pl. str. and A*\barEp_o^n=Cth
C
C     MODEL 13 Idem MODEL 6 with the consideration of a critical 
C              effective plastic strain in order to model the effect
C              of plastification without hardening (Idem MODEL 12)
C              CCOEF, CCOEB & CCOEQ: \barEp_c, n & A parameters
C              Closed form: Cp = A (<\barEp-\barEp_c>)^n
C              where \barEp=effec. pl. str., \barEp_c=critical effec.
C
C     These models can be applied to plastic or viscoplastic analysis.
C     They can include damage effects through the factor (1-d^p).
C
C
C     MODELS FOR THE COMPUTATION OF THE BACK STRESS TENSOR
C
C     IKINE defines the kinematic hardening model considered:
C
C     MODEL 3 (Armstrong-Frederick model):
C             Two coefficients: CKOEB & CKOEQ
C             L_e(x)=CKOEB d^p - CKOEQ \dot\barEp x
C             L_e: Lie derivative using the strain-like push-forward &
C                  & pull-back transformations
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION PROPS(*),    SIGMA(*),       DEEPP(*),   VECTG(6),
     .          HCAPA(6),    GE(2),          PRSIG(3)
      DIMENSION SA(3),SB(3), SC(3)
      DIMENSION DMATP(6,*),  DMATX(NAUXI,*), DMATE(6,*), DM(6,6)
      DIMENSION BACKA(*),    BACKS(*),       STRAP(*),   DP(6)
      DIMENSION UNOMA(6),    DESIG(*),       DEEPL(*),   STRAL(*),
     .          DESIL(*),    XJACM(NDIME,*), XJACI(NDIME,*),
     .          BACMA(*),    BACMS(*)
      DIMENSION DBACS(6),    DBACM(6),       RECRY(*)
      DIMENSION TRMTX(NDIME,*), STRAX(6), DEEPX(6)
C
      IF(LARGE.EQ.0) THEN
       DO ISTR1=1,NSTR1
        STRAL(ISTR1)=STRAP(ISTR1)
        DEEPL(ISTR1)=DEEPP(ISTR1)
       ENDDO
      ELSE
C
C**** TRANSFORMS TOTAL & INCREMENTAL PLASTIC DEFORMATION
C
       DETJX=1.0D0
       IF(IFREN.EQ.1.OR.IFREN.EQ.4) THEN                     ! no change
        CALL PUFOBA(STRAP,STRAL,XJACI,XJA3I,DETJX,    3)
        CALL PUFOBA(DEEPP,DEEPL,XJACI,XJA3I,DETJX,    3)
       ENDIF
       IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.IFREN.EQ.5.OR.
     .    IFREN.EQ.6.OR.IFREN.EQ.7.OR.IFREN.EQ.8.OR.
     .    IFREN.EQ.9.OR.IFREN.EQ.10) THEN
        CALL PUFOBA(STRAP,STRAL,XJACI,XJA3I,DETJX,    2)     ! Almansi
        CALL PUFOBA(DEEPP,DEEPL,XJACI,XJA3I,DETJX,    2)
       ENDIF
      ENDIF        ! large.eq.0
C
C**** ROTATES TOTAL & INCREMENTAL PLASTIC DEFORMATION FROM GLOBAL TO
C     LOCAL (MATERIAL) SYSTEM OF COORDINATES (THIS IS NEEDED FOR
C     ORTHOTROPIC MATERIALS)
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IF(ISOTR.NE.0) THEN
       DO ISTR1=1,NSTR1
        STRAX(ISTR1)=STRAL(ISTR1)
        DEEPX(ISTR1)=DEEPL(ISTR1)
       ENDDO
       IF(IFREN.EQ.2.OR.IFREN.EQ.3) THEN
        XJA3I=1.0D0
        DETJX=1.0D0
        CALL PUFOBA(STRAX,STRAL,TRMTX,XJA3I,DETJX,    1)
        CALL PUFOBA(DEEPX,DEEPL,TRMTX,XJA3I,DETJX,    1)
       ELSE
        CALL RUNEND('ERROR: IFREN NE 2 OR 3 NOT IMPLEMENTED')
       ENDIF
      ENDIF
C
      DO ISTR1=1,NSTR1
       UNOMA(ISTR1)=0.0D0
      ENDDO
      UNOMA(1)=1.0D0
      UNOMA(2)=1.0D0
      UNOMA(4)=1.0D0
C
      IREKI=0
C
C**** FLOW POTENTIAL CHOICE
C
      GO TO (31,32,33,34,35,36,37,38,39,40,41,42) (NCRIP-30)
C
C**** TRESCA
C
   31 CONTINUE
      CALL RUNEND('ERROR IN PLVARI       ')
      RETURN
C
C**** VON MISES
C
   32 CONTINUE
C
C**** COMPUTATION OF TOTAL EFFECTIVE PLASTIC DEFORMATION
C
      EFFP1=0.0D+00
      DO ISTR1=1,NSTR1
       IF(ISTR1.EQ.3.OR.ISTR1.EQ.5.OR.ISTR1.EQ.6) THEN
        EFFP1=EFFP1+0.5D0*DEEPL(ISTR1)*DEEPL(ISTR1)
       ELSE
        EFFP1=EFFP1+DEEPL(ISTR1)*DEEPL(ISTR1)
       ENDIF
      ENDDO
      EFFP1=DSQRT(2.0D0/3.0D0*EFFP1)
      EFFPL=EFFPL+EFFP1
C
C**** COMPUTES THE PLASTIC HARDENING FUNCTION
C
      IF(ISOTT.EQ.1) THEN
       DCAPA=0.0D0
       DO ISTRS=1,NSTRS
        DCAPA=DCAPA+CCOEF*SIGMA(ISTRS)*DEEPL(ISTRS)
       ENDDO
       CAPAP=CAPAP+DCAPA*(1.0D0-DAMAG)
c      COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.2) THEN
       DCAPA=CCOEF*EFFP1
       CAPAP=CAPAP+DCAPA*(1.0D0-DAMAG)
c      COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.3) THEN
       CAPAP=CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))*(1.0D0-DAMAG)
c      COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.4) THEN
       CAPAP=CCOEF*EFFPL*(1.0D0-DAMAG)
       COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.5) THEN
       EFFPL=EFFPL-EFFP1      ! changes the standard definition of EFFPL
       EFFP1=EFFP1*EFFP1*3.0D0/2.0D0
       EFFPL=EFFPL+EFFP1
       CAPAP=CCOEF*EFFPL*(1.0D0-DAMAG)
c      COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.6) THEN
       CAPAP=CCOEQ*(EFFPL**CCOEB)*(1.0D0-DAMAG)
       COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.7) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPAP=CCOEQ*((EFFPL0+EFFPL)**CCOEB)*(1.0D0-DAMAG)
c      COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.8) THEN
       CAPAP=CCOEF*EFFPL+CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))*(1.0D0-DAMAG)
c      COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.9) THEN
       CALL CCOEFET(CCOEFE,DCCOEFE,
     .              TEMPG,EFFPL,PROPS)
       CAPAP=CCOEFE*(1.0D0-DAMAG)
c      COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.10) THEN
       CAPA1=(CCERO+CCOEQ*(EFFPL**CCOEB))
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
c      COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.11) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPA1=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
c      COUTH=-CAPAP*EFFP1     ! isotropic hardening dissipation
      ENDIF
      IF(ISOTT.EQ.12) THEN
       CAPAP=CCERO
       IF(EFFPL.GE.CCOEF) THEN
        EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
        CAPAP=CCOEQ*((EFFPL0+EFFPL-CCOEF)**CCOEB)
        CAPAP=CAPAP*(1.0D0-DAMAG)
c       COUTH=-CAPAP*EFFP1    ! isotropic hardening dissipation
       ENDIF
      ENDIF
      IF(ISOTT.EQ.13) THEN
       IF(EFFPL.GT.CCOEF) THEN
        CAPAP=CCOEQ*((EFFPL-CCOEF)**CCOEB)
        CAPAP=CAPAP*(1.0D0-DAMAG)
c       COUTH=-CAPAP*EFFP1    ! isotropic hardening dissipation
       ENDIF
      ENDIF
C
C**** COMPUTES THE BACK STRESS TENSOR
C
      IF(IKINE.EQ.1) THEN
       IREKI=1
       DO ISTR1=1,NSTR1
        BACKS(ISTR1)=CKOEF*STRAL(ISTR1)
c       COUTH=-CAPAP*EFFP1    ! kinematic hardening dissipation
       ENDDO
      ENDIF
      IF(IKINE.EQ.2) THEN
       CALL RUNEND('ERROR: IKINE=2 NOT IMPLEMENTED YET - plvari.f')
      ENDIF
      IF(IKINE.EQ.3) THEN
       IREKI=1
       IFOBA=3
       IF(LARGE.NE.0) IFOBA=2
       DO ISTR1=1,NSTR1
        DBACS(ISTR1)=CKOEB*DEEPL(ISTR1)-CKOEQ*EFFP1*BACKS(ISTR1)
       ENDDO
       DETJX=1.0D0
       CALL PUFOBA(DBACS,DBACM,XJACM,XJA3M,DETJX,IFOBA)   ! pull-back
       CALL PUFOBA(BACKS,BACMS,XJACM,XJA3M,DETJX,IFOBA)   ! pull-back
       DO ISTR1=1,NSTR1
        BACMS(ISTR1)=BACMS(ISTR1)+DBACM(ISTR1)
       ENDDO
       CALL PUFOBA(BACMS,BACKS,XJACI,XJA3I,DETJX,IFOBA)   ! push-forward
c      COUTH=-CAPAP*EFFP1     ! kinematic hardening dissipation
      ENDIF
      IF(IKINE.EQ.4) THEN
       IREKI=1
       IFOBA=3
       IF(LARGE.NE.0) IFOBA=1
       DO ISTR1=1,NSTR1
        DBACS(ISTR1)=CKOEB*DEEPL(ISTR1)-CKOEQ*EFFP1*BACKS(ISTR1)
       ENDDO
       DETJX=1.0D0
       CALL PUFOBA(DBACS,DBACM,XJACI,XJA3I,DETJX,IFOBA)   ! pull-back
       CALL PUFOBA(BACKS,BACMS,XJACI,XJA3I,DETJX,IFOBA)   ! pull-back
       DO ISTR1=1,NSTR1
        BACMS(ISTR1)=BACMS(ISTR1)+DBACM(ISTR1)
       ENDDO
       CALL PUFOBA(BACMS,BACKS,XJACM,XJA3M,DETJX,IFOBA)   ! push-forward
c      COUTH=-CAPAP*EFFP1     ! kinematic hardening dissipation
      ENDIF
      IF(IKINE.EQ.5) THEN
       CALL RUNEND('ERROR: IKINE=5 NOT IMPLEMENTED YET - plvari.f')
      ENDIF
C
C**** COMPUTES DEFORMATION RESISTANCE & GRAIN SIZE (RECRYSTALLIZATION)
C
      IF(IRECR.EQ.1) THEN
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
C
       EASTE=ERECR*1000.0D0      ! initial guess for L* (from mum to mm)
       IF(EFFPX.GT.0.0D0) THEN             ! to avoid zero**RECRY(6)
        PP= 1.0D0/RECRY(18)
        PM=-1.0D0/RECRY(18)
        AA=EFFPX/RECRY(1)*(EASTE/RECRY(8))**RECRY(2)*RECRY(21)   ! A
        AB=(RECRY(19)/(RECRY(5)*RECRY(11)*AA**RECRY(12)))**PP
        AC1=AA**RECRY(6)+DSQRT(AA**(2.0D0*RECRY(6))+1.0D0)
        AC2=DLOG(AC1)
        AC=AC2**PM
        AR=AB*AC-EASTE
        I=0
        DO WHILE ((DABS(AR).GT.ERECR/10000.0D0).AND.(I.LT.100))
         I=I+1
         DADL=RECRY(2)*AA/EASTE
         DBDL=-RECRY(12)/RECRY(18)*AB/AA*DADL
         DCDL=-AC/(RECRY(18)*AC2*AC1)*RECRY(6)/AA*DADL*(AA**RECRY(6)+
     .        AA**(2.0D0*RECRY(6))/DSQRT(AA**(2.0D0*RECRY(6))+1.0D0))
         AJ=-DBDL*AC-AB*DCDL+1.0D0
         DEASTE=AR/AJ
         EASTE=EASTE+DEASTE
         AA=EFFPX/RECRY(1)*(EASTE/RECRY(8))**RECRY(2)*RECRY(21)
         AB=(RECRY(19)/(RECRY(5)*RECRY(11)*AA**RECRY(12)))**PP
         AC1=AA**RECRY(6)+DSQRT(AA**(2.0D0*RECRY(6))+1.0D0)
         AC2=DLOG(AC1)
         AC=AC2**PM
         AR=AB*AC-EASTE
        ENDDO
       ENDIF
       EASTE=EASTE/1000.0D0                             ! from mum to mm
C
       SASTE0=SRECR                           ! rough estimation
       IF(EFFPX.GT.0.0D0)                     ! to avoid zero**RECRY(12)
     .  SASTE0=RECRY(11)*(EFFPX/RECRY(1)*RECRY(21))**RECRY(12)   ! So*
       SASTE=SASTE0*(EASTE/RECRY(8))**(RECRY(2)*RECRY(12))       ! S*
C
       EMAC1=RECRY(14)*SASTE-RECRY(7)
       IF(EMAC1.LT.0.0D0) EMAC1=0.0D0
       EPSIC=2.0D0/DSQRT(3.0D0)*RECRY(17)/DISTM*EMAC1            ! ep_C
C
       EMAC2=(ERECR/RECRY(8))**(RECRY(2)*RECRY(12))-
     .       (EASTE/RECRY(8))**(RECRY(2)*RECRY(12))
       IF(EMAC2.LT.0.0D0) EMAC2=0.0D0
       EPSIR=2.0D0/DSQRT(3.0D0)*RECRY(16)/DISTM*SASTE0*EMAC2     ! ep_R
C
       SF=SRECR                                                  ! Sf
       IF((EFFPL-EPSIC).LE.EPSIR) THEN
        AS1=1.0D0-RECRY(10)
        AS2=1.0D0/AS1
        EMAC3=RECRY(9)/SASTE0*(RECRY(8)/EASTE)**(RECRY(2)*RECRY(12))*
     .     (RECRY(10)-1.0D0)*(EPSIC+EPSIR)+(1.0D0-
     .     RECRY(7)/SASTE0*(RECRY(8)/EASTE)**(RECRY(2)*RECRY(12)))**AS1
        IF(EMAC3.GT.0.0D0) THEN                     ! to avoid zero**AS2
         EMAC3=EMAC3**AS2
        ELSE
         EMAC3=0.0D0
        ENDIF
        SF=SASTE0*(EASTE/RECRY(8))**(RECRY(2)*RECRY(12))*(1.0D0-EMAC3)
       ENDIF
C
       IF(EPSIR.GT.0.0D0) THEN
        EMAC4=(EFFPL-EPSIC)/EPSIR
        IF(EMAC4.LT.0.0D0) EMAC4=0.0D0
        XR=1.0D0-DEXP(-RECRY(15)*RECRY(8)/EASTE*EMAC4)           ! XR
       ELSE
        XR=1.0D0
        IF(EFFPL.LT.EPSIC) XR=0.0D0
       ENDIF
C
       AS3=1.0D0-(RECRY(8)/EASTE)**(RECRY(2)*RECRY(12))*SRECR/SASTE0
       AS4=1.0D0
       IF(AS3.LT.0.0D0) THEN
        AS3=-AS3
        AS4=-1.0D0
       ENDIF
       ASREC=RECRY(9)*AS3**RECRY(10)*AS4-RECRY(13)*XR*(SRECR-SF)
       DSRECR=ASREC*EFFPX                                        ! dot S
C
       XRC=1.0D0-DEXP(-RECRY(15)*RECRY(8)/EASTE)                 ! XRC
C
       EMAC5=ERECR-EASTE
       IF(EMAC5.LT.0.0D0) EMAC5=0.0D0
       EMAC6=XR-XRC
       IF(EMAC6.LT.0.0D0) EMAC6=0.0D0
       ALREC=0.0D0
       DERECR=0.0D0                                              ! dot L
       IF(EFFPL.GT.EPSIC) THEN
        ALREC=-RECRY(13)*XR*EMAC5
        DERECR=ALREC*EFFPX+
     .         RECRY(20)*(1.0D0-DEXP(-EMAC6))/RECRY(21)*RECRY(8)/ERECR
       ENDIF
C
       SRECR=SRECR+DSRECR*DTIME                            ! integration
       ERECR=ERECR+DERECR*DTIME
      ENDIF                   ! irecr.eq.1
      RETURN
C
C**** MOHR-COULOMB & VERSION  WITH "TENSION CUT-OFF"
C
   33 CONTINUE
      CALL RUNEND('ERROR IN PLVARI       ')
      RETURN
C
C**** DRUCKER-PRAGER
C
   34 CONTINUE
      CALL RUNEND('ERROR IN PLVARI       ')
      RETURN
C
C**** J. LUBLINER'S THEORY
C
   35 CONTINUE
      CALL RUNEND('ERROR IN PLVARI       ')
      RETURN
C
C**** ABOUAF'S MODEL
C
   36 CONTINUE
      CALL RUNEND('ERROR IN PLVARI       ')
      RETURN
C
C**** WEBER BROWN'S MODEL
C
   37 CONTINUE
      CALL RUNEND('ERROR IN PLVARI       ')
      RETURN
C
C**** SG CAST IRON
C
   38 CONTINUE
C
C**** COMPUTATION OF TOTAL EFFECTIVE PLASTIC DEFORMATION
C
      EFFP1=0.0D+00
      DO ISTR1=1,NSTR1
       IF(ISTR1.EQ.3.OR.ISTR1.EQ.5.OR.ISTR1.EQ.6) THEN
        EFFP1=EFFP1+0.5D0*DEEPL(ISTR1)*DEEPL(ISTR1)
       ELSE
        EFFP1=EFFP1+DEEPL(ISTR1)*DEEPL(ISTR1)
       ENDIF
      ENDDO
      EFFP1=DSQRT(2.0D0/3.0D0*EFFP1)
      EFFPL=EFFPL+EFFP1
C
C**** COMPUTES THE PLASTIC HARDENING FUNCTION
C
      IF(ISOTT.EQ.1) THEN
       DCAPA=0.0D0
       DO ISTRS=1,NSTRS
        DCAPA=DCAPA+CCOEF*SIGMA(ISTRS)*DEEPL(ISTRS)
       ENDDO
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.2) THEN
       DCAPA=CCOEF*EFFP1
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.3) THEN
       CAPAP=CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.4) THEN
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.5) THEN
       EFFPL=EFFPL-EFFP1
       EFFP1=EFFP1*EFFP1*3.0D0/2.0D0
       EFFPL=EFFPL+EFFP1
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.6) THEN
       CAPAP=CCOEQ*(EFFPL**CCOEB)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPAP=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       CAPAP=CCOEF*EFFPL+CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.9) THEN
       CALL CCOEFET(CCOEFE,DCCOEFE,
     .              TEMPG,EFFPL,PROPS)
       CAPAP=CCOEFE
      ENDIF
      IF(ISOTT.EQ.10) THEN
       CAPA1=(CCERO+CCOEQ*(EFFPL**CCOEB))
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.11) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPA1=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.12) THEN
       CAPAP=CCERO
       IF(EFFPL.GT.CCOEF) THEN
        EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
        CAPAP=CCOEQ*((EFFPL0+EFFPL-CCOEF)**CCOEB)
       ENDIF
      ENDIF
C
      IF(IKINE.EQ.1) THEN
       DO ISTR1=1,NSTR1
        BACKA(ISTR1)=CKOEF*STRAL(ISTR1)
        BACKS(ISTR1)=BACKA(ISTR1)
       ENDDO
      ENDIF
      RETURN
C
C**** GREEN SAND
C
   39 CONTINUE
C
C**** COMPUTATION OF TOTAL EFFECTIVE PLASTIC DEFORMATION
C
      EFFP1=0.0D+00
      DO ISTR1=1,NSTR1
       IF(ISTR1.EQ.3.OR.ISTR1.EQ.5.OR.ISTR1.EQ.6) THEN
        EFFP1=EFFP1+0.5D0*DEEPL(ISTR1)*DEEPL(ISTR1)
       ELSE
        EFFP1=EFFP1+DEEPL(ISTR1)*DEEPL(ISTR1)
       ENDIF
      ENDDO
      EFFP1=DSQRT(2.0D0/3.0D0*EFFP1)
      EFFPL=EFFPL+EFFP1
C
C**** COMPUTES THE PLASTIC HARDENING FUNCTION
C
      IF(ISOTT.EQ.1) THEN
       DCAPA=0.0D0
       DO ISTRS=1,NSTRS
        DCAPA=DCAPA+CCOEF*SIGMA(ISTRS)*DEEPL(ISTRS)
       ENDDO
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.2) THEN
       DCAPA=CCOEF*EFFP1
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.3) THEN
       CAPAP=CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.4) THEN
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.5) THEN
       EFFPL=EFFPL-EFFP1
       EFFP1=EFFP1*EFFP1*3.0D0/2.0D0
       EFFPL=EFFPL+EFFP1
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.6) THEN
       CAPAP=CCOEQ*(EFFPL**CCOEB)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPAP=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       CAPAP=CCOEF*EFFPL+CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.9) THEN
       CALL CCOEFET(CCOEFE,DCCOEFE,
     .              TEMPG,EFFPL,PROPS)
       CAPAP=CCOEFE
      ENDIF
      IF(ISOTT.EQ.10) THEN
       CAPA1=(CCERO+CCOEQ*(EFFPL**CCOEB))
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.11) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPA1=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.12) THEN
       CAPAP=CCERO
       IF(EFFPL.GT.CCOEF) THEN
        EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
        CAPAP=CCOEQ*((EFFPL0+EFFPL-CCOEF)**CCOEB)
       ENDIF
      ENDIF
C
      IF(IKINE.EQ.1) THEN
       DO ISTR1=1,NSTR1
        BACKA(ISTR1)=CKOEF*STRAL(ISTR1)
        BACKS(ISTR1)=BACKA(ISTR1)
       ENDDO
      ENDIF
      RETURN
C
C**** GREEN SAND
C
   40 CONTINUE
C
C**** COMPUTATION OF TOTAL EFFECTIVE PLASTIC DEFORMATION
C
      EFFP1=0.0D+00
      DO ISTR1=1,NSTR1
       IF(ISTR1.EQ.3.OR.ISTR1.EQ.5.OR.ISTR1.EQ.6) THEN
        EFFP1=EFFP1+0.5D0*DEEPL(ISTR1)*DEEPL(ISTR1)
       ELSE
        EFFP1=EFFP1+DEEPL(ISTR1)*DEEPL(ISTR1)
       ENDIF
      ENDDO
      EFFP1=DSQRT(2.0D0/3.0D0*EFFP1)
      EFFPL=EFFPL+EFFP1
C
C**** COMPUTES THE PLASTIC HARDENING FUNCTION
C
      IF(ISOTT.EQ.1) THEN
       DCAPA=0.0D0
       DO ISTRS=1,NSTRS
        DCAPA=DCAPA+CCOEF*SIGMA(ISTRS)*DEEPL(ISTRS)
       ENDDO
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.2) THEN
       DCAPA=CCOEF*EFFP1
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.3) THEN
       CAPAP=CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.4) THEN
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.5) THEN
       EFFPL=EFFPL-EFFP1
       EFFP1=EFFP1*EFFP1*3.0D0/2.0D0
       EFFPL=EFFPL+EFFP1
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.6) THEN
       CAPAP=CCOEQ*(EFFPL**CCOEB)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPAP=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       CAPAP=CCOEF*EFFPL+CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.9) THEN
       CALL CCOEFET(CCOEFE,DCCOEFE,
     .              TEMPG,EFFPL,PROPS)
       CAPAP=CCOEFE
      ENDIF
      IF(ISOTT.EQ.10) THEN
       CAPA1=(CCERO+CCOEQ*(EFFPL**CCOEB))
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.11) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPA1=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.12) THEN
       CAPAP=CCERO
       IF(EFFPL.GT.CCOEF) THEN
        EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
        CAPAP=CCOEQ*((EFFPL0+EFFPL-CCOEF)**CCOEB)
       ENDIF
      ENDIF
C
      IF(IKINE.EQ.1) THEN
       DO ISTR1=1,NSTR1
        BACKA(ISTR1)=CKOEF*STRAL(ISTR1)
        BACKS(ISTR1)=BACKA(ISTR1)
       ENDDO
      ENDIF
      RETURN
C
C**** GURSON
C
   41 CONTINUE
C
      IF(LARGE.EQ.0) THEN
       DO ISTR1=1,NSTR1
        DESIL(ISTR1)=DESIG(ISTR1)
       ENDDO
      ELSE
C
C**** TRANSFORMS INCREMENTAL STRESSES
C
       IF(IFREN.EQ.1.OR.IFREN.EQ.4) THEN                     ! no change
        DO ISTR1=1,NSTR1
         DESIL(ISTR1)=DESIG(ISTR1)
        ENDDO
       ENDIF
       IF(IFREN.EQ.2.OR.IFREN.EQ.5.OR.IFREN.EQ.7.OR.         ! Kirchhoff
     .    IFREN.EQ.8.OR.IFREN.EQ.9.OR.IFREN.EQ.10) THEN
        DETJX=1.0D0
        CALL PUFOBA(DESIG,DESIL,XJACM,XJA3M,DETJX,    1)
       ENDIF
       IF(IFREN.EQ.3.OR.IFREN.EQ.6) THEN                     ! Cauchy
        DETJX=DETJM
        CALL PUFOBA(DESIG,DESIL,XJACM,XJA3M,DETJX,    1)
       ENDIF
      ENDIF        ! large.eq.0
C
C**** COMPUTATION OF TOTAL EFFECTIVE PLASTIC DEFORMATION
C
      PREYS=A3COE*CCERO+A2COE*CAPAP          ! as in pltoha.f
      IFVER=INT(PROPS(35))
      IF(IFVER.EQ.4) THEN
       IF(CAPAP.EQ.0.0D0) PREYS=CCERO
      ENDIF
C
      EFFP1=0.0D+00
      DO ISTR1=1,NSTR1
c      IF(ISTR1.EQ.3.OR.ISTR1.EQ.5.OR.ISTR1.EQ.6) THEN    ??
c       EFFP1=EFFP1+0.5D0*SIGMA(ISTR1)*DEEPL(ISTR1)
c      ELSE
        EFFP1=EFFP1+SIGMA(ISTR1)*DEEPL(ISTR1)
c      ENDIF
      ENDDO
c     EFFP1=EFFP1/((1.0D0-POROS)*YIELD)
      EFFP1=EFFP1/((1.0D0-POROS)*PREYS)
      EFFPL=EFFPL+EFFP1
C
C**** COMPUTES THE PLASTIC HARDENING FUNCTION
C
      IF(ISOTT.EQ.1) THEN
       DCAPA=0.0D0
       DO ISTRS=1,NSTRS
        DCAPA=DCAPA+CCOEF*SIGMA(ISTRS)*DEEPL(ISTRS)
       ENDDO
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.2) THEN
       DCAPA=CCOEF*EFFP1
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.3) THEN
       CAPAP=CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.4) THEN
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.5) THEN
       EFFPL=EFFPL-EFFP1
       EFFP1=EFFP1*EFFP1*3.0D0/2.0D0
       EFFPL=EFFPL+EFFP1
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.6) THEN
       CAPAP=CCOEQ*(EFFPL**CCOEB)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPAP=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       CAPAP=CCOEF*EFFPL+CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.9) THEN
       CALL CCOEFET(CCOEFE,DCCOEFE,
     .              TEMPG,EFFPL,PROPS)
       CAPAP=CCOEFE
      ENDIF
      IF(ISOTT.EQ.10) THEN
       CAPA1=(CCERO+CCOEQ*(EFFPL**CCOEB))
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.11) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPA1=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.12) THEN
       CAPAP=CCERO
       IF(EFFPL.GT.CCOEF) THEN
        EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
        CAPAP=CCOEQ*((EFFPL0+EFFPL-CCOEF)**CCOEB)
       ENDIF
      ENDIF
C
      IF(IKINE.EQ.1) THEN
       DO ISTR1=1,NSTR1
        BACKA(ISTR1)=CKOEF*STRAL(ISTR1)
        BACKS(ISTR1)=BACKA(ISTR1)
       ENDDO
      ENDIF
C
C**** POROSITY EVOLUTION
C
C     MODEL 1 Gurson model: only nucleation
C
C     MODEL 2 Gurson model: only growth
C
C     MODEL 3 Gurson model: nucleation & growth
C
C     MODEL 4 Gurson-Tvergaard model: nucleation & growth
C
C
C     Notes:
C
C     PREYA: last converged total hardening function (it is necessary
C            to be stored when nucleation is considered. Here, it is
C            always stored)
C
      IF(IPORO.EQ.1.OR.IPORO.EQ.3) THEN                     ! nucleation
       DPORO=0.0D0
       HEAVI=0.0D0
       DO ISTR1=1,NSTR1
        HEAVI=HEAVI+UNOMA(ISTR1)*DESIL(ISTR1)
       ENDDO
       HEAVI=HEAVI/3.0D0+PREYS-PREYA
       IF(HEAVI.GT.0.0D0) DPORO=CPOEF/PREYS*HEAVI
       POROS=POROS+DPORO
      ENDIF
      IF(IPORO.EQ.4) THEN                                   ! nucleation
       DPORO=CPOEF*EFFP1
       POROS=POROS+DPORO
      ENDIF
      IF(IPORO.EQ.2.OR.IPORO.EQ.3.OR.IPORO.EQ.4) THEN       ! growth
       DPORO=0.0D0
       POAUX=1.0D0-POROS
       DO ISTR1=1,NSTR1
        DPORO=DPORO+POAUX*UNOMA(ISTR1)*DEEPL(ISTR1)
       ENDDO
       POROS=POROS+DPORO
      ENDIF
      IF(POROS.GT.0.99) THEN                       ! approximation
       CALL RUNMEN('POROS = 1 IN SOME POINTS  ')
       POROS=0.99D0
      ENDIF
      RETURN
C
C**** HILL 48
C
   42 CONTINUE
C
C**** COMPUTATION OF TOTAL EFFECTIVE PLASTIC DEFORMATION
C
      EFFP1=TLAMD*YIELP/YIELD
      EFFPL=EFFPL+EFFP1
C
C**** COMPUTES THE PLASTIC HARDENING FUNCTION
C
      IF(ISOTT.EQ.1) THEN
       DCAPA=0.0D0
       DO ISTRS=1,NSTRS
        DCAPA=DCAPA+CCOEF*SIGMA(ISTRS)*DEEPL(ISTRS)
       ENDDO
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.2) THEN
       DCAPA=CCOEF*EFFP1
       CAPAP=CAPAP+DCAPA
      ENDIF
      IF(ISOTT.EQ.3) THEN
       CAPAP=CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.4) THEN
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.5) THEN
       EFFPL=EFFPL-EFFP1      ! changes the standard definition of EFFPL
       EFFP1=EFFP1*EFFP1*3.0D0/2.0D0
       EFFPL=EFFPL+EFFP1
       CAPAP=CCOEF*EFFPL
      ENDIF
      IF(ISOTT.EQ.6) THEN
       CAPAP=CCOEQ*(EFFPL**CCOEB)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPAP=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       CAPAP=CCOEF*EFFPL+CCOEQ*(1.0D0-DEXP(-CCOEB*EFFPL))
      ENDIF
      IF(ISOTT.EQ.9) THEN
       CALL CCOEFET(CCOEFE,DCCOEFE,
     .              TEMPG,EFFPL,PROPS)
       CAPAP=CCOEFE
      ENDIF
      IF(ISOTT.EQ.10) THEN
       CAPA1=(CCERO+CCOEQ*(EFFPL**CCOEB))
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.11) THEN
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       CAPA1=CCOEQ*((EFFPL0+EFFPL)**CCOEB)
       EFFPX=EFFP1/DTIME            ! effective plastic deformation rate
       IF(EFFPX.GT.EFFP2) EFFP2=EFFPX       ! maximum eff. pl. def. rate
       IF(EFFP2.GT.0.0D0) CAPA2=1.0D0+CCOEF*DLOG(EFFP2)
       IF(CAPA2.GT.1.0D0) CAPA1=CAPA1*CAPA2
       CAPAP=CAPA1*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.12) THEN
       CAPAP=CCERO
       IF(EFFPL.GT.CCOEF) THEN
        EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
        CAPAP=CCOEQ*((EFFPL0+EFFPL-CCOEF)**CCOEB)
       ENDIF
      ENDIF
C
      IF(IKINE.EQ.1) THEN
       DO ISTR1=1,NSTR1
        BACKA(ISTR1)=CKOEF*STRAL(ISTR1)
        BACKS(ISTR1)=BACKA(ISTR1)
       ENDDO
      ENDIF
      RETURN
C
      END
