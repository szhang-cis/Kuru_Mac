      SUBROUTINE DAVARI(PROPS,ALENG,SIGMA,DEEPP,CAPAP,
     .                  VECTG,NSTRS,ANGFI,ANGSI,
     .                  DFTEQ,DFAFI,YIELD,CAPAD,NSTR1,
     .                  COEPD,VADEG,DMATP,COUTH,
     .                  DMATX,DMATE,RETEN,DM,CCERO,NCRIP,EFFPL,EFFP1,DP,
     .                  CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .                  BACKA,BACKS,STRAP,ISOTT,IKINE,NTYPE,
     .                  DTIME,TLAMD,DLAMD,DAMAG,CDAMA,CFRAC,CCOB9,
     .                  CCOB10,CCOB11,CCOB12,IDAMG,PMEAN,NAUXI)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES:
C
C     THE INTERNAL VARIABLES:
C
C               DAMAG: damage
C
C     IDAMG: damage indicator
C
C     =0     no damage
C
C
C     Damage models governed by plasticity or viscoplasticity
C
C     =1     damage variable is an internal var. with evolution law
C     =2     damage variable is an internal var. without evolution law
C     =3     damage variable is an internal variable of green sand model
C     =4     Lemaitre model
C     =5     Lemaitre model with no evolution for negative pressure
C
C
C     Damage models governed by a damage criterion
C
C     =21    it only evaluates a damage criterion
C     =22    not implemented
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/PROPSC/POISC1,POISC2,POISC3,STRRA,STRR4,DAMAP1,DAMAP2
C
      DIMENSION PROPS(*),   SIGMA(*),       DEEPP(6),   VECTG(6),
     .          HCAPA(6),   GE(2),          PRSIG(3)
      DIMENSION SA(3),SB(3),SC(3)
      DIMENSION DMATP(6,6), DMATX(NAUXI,*), DMATE(6,6), DM(6,6)
      DIMENSION BACKA(*),   BACKS(*),       STRAP(*),   DP(6),
     .          STRSP(3)
C
      IF(IDAMG.EQ.0) RETURN
C
      IF(IDAMG.LE.20) THEN   ! damage models governed by pl. or viscopl.
C
C**** FLOW POTENTIAL CHOICE
C
       GO TO (31,32,33,34,35,36,37,38,39,40,41) (NCRIP-30)
C
C**** TRESCA
C
   31  CONTINUE
       CALL RUNEND('ERROR IN DAVARI       ')
       RETURN
C 
C**** VON MISES
C
   32  CONTINUE
       IF(IDAMG.EQ.1) THEN
c       DAMAG=DAMAG+DLAMD*CDAMA*DTIME        ! ojo con viscoplas.=>DTIME
c       DAMAG=DAMAG+TLAMD*CDAMA
        DAMAG=DAMAG+CDAMA*EFFP1
c       COUTH=COUTH+XXX                   ! damage dissipation
       ENDIF
       IF(IDAMG.EQ.2.OR.IDAMG.EQ.3) THEN
        call runend('ERROR: idamg=2,3 not implemented in davari.f')
       ENDIF
       IF(IDAMG.EQ.4.OR.(IDAMG.EQ.5.AND.PMEAN.GT.0.0D0)) THEN
        IF(DAMAG.LT.CCOB9) THEN
         IF(EFFPL.GT.CDAMA) THEN
          CAPAP=CAPAP/(1.0D0-DAMAG)       ! undamaged hardening function
          DAMAG=DAMAG+CCOB9/(CFRAC-CDAMA)*
     .              (DAMAP1+DAMAP2*(PMEAN/YIELD)*(PMEAN/YIELD))*
     .               EFFP1
c    .              (EFFPL**(2.0D00*CCOEB))*EFFP1 ! to be revised
          IF(DAMAG.GT.CCOB10) THEN
           DAMAG=CCOB10
           CALL RUNMEN('WARNING: DAMAGE VARIABLE SET TO FINAL VALUE')
          ENDIF
          CAPAP=(1.0D0-DAMAG)*CAPAP       ! damaged hardening function
         ENDIF
        ELSE
         IF(DAMAG.LT.CCOB10) THEN
          CAPAP=CAPAP/(1.0D0-DAMAG)       ! undamaged hardening function
          CCOB9X=CCOB10*(1.0D0-DAMAG)**CCOB11           ! proposed slope
          DAMAG=DAMAG+CCOB9X/(CFRAC-CDAMA)*
     .              (DAMAP1+DAMAP2*(PMEAN/YIELD)*(PMEAN/YIELD))*
     .               EFFP1
c    .              (EFFPL**(2.0D00*CCOEB))*EFFP1 ! to be revised
          IF(DAMAG.GT.CCOB10) THEN
           DAMAG=CCOB10
           CALL RUNMEN('WARNING: DAMAGE VARIABLE SET TO FINAL VALUE')
          ENDIF
          CAPAP=(1.0D0-DAMAG)*CAPAP       ! damaged hardening function
         ENDIF
        ENDIF
c       COUTH=COUTH+XXX                   ! damage dissipation
       ENDIF
       IF(IDAMG.EQ.6) THEN
        IF(DAMAG.LT.CCOB12) THEN
         CAPAP=CAPAP/(1.0D0-DAMAG)        ! undamaged hardening function
         IF(NTYPE.EQ.4) THEN
          CALL PRINSI(SIGMA,STRSP)
         ELSE
          CALL PRIDIR(NTYPE,SIGMA,STRSP)
         ENDIF
         STRSM=STRSP(1)
         IF(STRSP(2).GT.STRSM) STRSM=STRSP(2)
         IF(STRSP(3).GT.STRSM) STRSM=STRSP(3)
         IF(STRSM.LT.0.0D0) STRSM=0.0D0
         PMEAX=PMEAN
         IF(PMEAX.LT.0.0D0) PMEAX=0.0D0
         YVARI=CCOB10*PMEAX+CCOB11*STRSM+(1.0D0-CCOB10-CCOB11)*YIELD
         DAMAG=DAMAG+(YVARI/CDAMA)**CFRAC/(1.0D0-DAMAG)**CCOB9*EFFP1
         IF(DAMAG.GT.CCOB12) THEN
          DAMAG=CCOB12
          CALL RUNMEN('WARNING: DAMAGE VARIABLE SET TO FINAL VALUE')
         ENDIF
         CAPAP=(1.0D0-DAMAG)*CAPAP        ! damaged hardening function
        ENDIF
c       COUTH=COUTH+XXX                   ! damage dissipation
       ENDIF
       IF(IDAMG.EQ.11) THEN               ! index-based models
        DAMAG=DAMAG+YIELD*EFFP1
       ENDIF
       RETURN
C
C**** MOHR-COULOMB & VERSION  WITH "TENSION CUT-OFF"
C
   33  CONTINUE
       CALL RUNEND('ERROR IN DAVARI       ')
       RETURN
C
C**** DRUCKER-PRAGER
C
   34  CONTINUE
       CALL RUNEND('ERROR IN DAVARI       ')
       RETURN
C
C**** J. LUBLINER'S THEORY
C
   35  CONTINUE
       CALL RUNEND('ERROR IN DAVARI       ')
       RETURN
C
C**** ABOUAF'S MODEL
C
   36  CONTINUE
       CALL RUNEND('ERROR IN DAVARI       ')
       RETURN
C
C**** WEBER BROWN'S MODEL
C
   37  CONTINUE
       CALL RUNEND('ERROR IN DAVARI       ')
       RETURN
C
C**** SG CAST IRON
C
   38  CONTINUE         ! no damage for sg cast iron
       RETURN
C
C**** GREEN SAND
C
   39  CONTINUE
       IF(IDAMG.EQ.1) THEN
        DAMAG=DAMAG+DLAMD*CDAMA*DTIME
       ENDIF
       IF(DAMAG.GT.1.0) THEN
        DAMAG=1.0D0
        CALL RUNMEN('WARNING: DAMAG SET TO 1.0')
       ENDIF
       RETURN
C
C**** GREEN SAND
C
   40  CONTINUE
       IF(IDAMG.EQ.1) THEN
        DAMAG=DAMAG+DLAMD*CDAMA*DTIME
       ENDIF
       IF(DAMAG.GT.1.0) THEN
        DAMAG=1.0D0
        CALL RUNMEN('WARNING: DAMAG SET TO 1.0')
       ENDIF
       RETURN
C
C**** GURSON
C
   41  CONTINUE
       IF(IDAMG.EQ.1) THEN
c       DAMAG=DAMAG+DLAMD*CDAMA*DTIME
c       DAMAG=DAMAG+TLAMD*CDAMA
        DAMAG=DAMAG+CDAMA*EFFP1
       ENDIF
       IF(IDAMG.EQ.2.OR.IDAMG.EQ.3) THEN
        call runend('ERROR: idamg=2,3 not implemented in davari.f')
       ENDIF
       IF(IDAMG.EQ.4) THEN
        DAMAG=DAMAG+CCOB9/(CFRAC-CDAMA)*
     .              (DAMAP1+DAMAP2*(PMEAN/YIELD)*(PMEAN/YIELD))*
     .              (EFFPL**(2.0D00*CCOEB))*EFFP1
       ENDIF
       IF(DAMAG.GT.1.0) THEN
        DAMAG=1.0D0
        CALL RUNMEN('WARNING: DAMAG SET TO 1.0')
       ENDIF
       RETURN
C
      ELSE                ! damage models governed by a damage criterion
C
       GO TO (131,132) (IDAMG-20)
C
C**** CONCRETE DAMAGE MODEL
C
  131  CONTINUE
c      IF(IDAMG.EQ.1) THEN
c       DAMAG=DAMAG+DLAMD*CDAMA*DTIME
c      ENDIF
c      IF(DAMAG.GT.1.0) THEN
c       DAMAG=1.0D0
c       CALL RUNMEN('WARNING: DAMAG SET TO 1.0')
c      ENDIF
       RETURN
C
C**** not implemented !!
C
  132  call runend('idamg=22 not implemented')
       RETURN
C
      ENDIF                     ! idamg.le.20
C 
      END
