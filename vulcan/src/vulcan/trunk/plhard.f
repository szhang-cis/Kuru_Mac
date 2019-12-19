      SUBROUTINE PLHARD(HARDS,SIGMA,VECTF,VECTG,DMATP,DESIG,HARDE,BACKS,
     .                  NSTRS,NSTR1,NCRIT,
     .                  ISOTT,IKINE,IPORO,
     .                  CCOEF,A2COE,CCOEB,CCOEQ,CPOEF,CPOE1,CPOE2,
     .                  CKOEF,CKOEB,CKOEQ,
     .                  EFFPL,EFFP1,EFFP2,YIELD,PREYS,PREYA,PMEAN,POROS,
     .                  CCERO,DAMAG)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES "HARDS" (TERM OF THE "ABETA" MULTIPLICATOR)
C
C     Notes:
C
C     This computation is for plastic & viscoplastic models.
C
C     In general, this term is:
C     \sum_{i-1}^n_intF - \partial F / \partial \alpha_i * f_i
C     where n_intF is the number of internal variables \alpha_i
C     appearing in the yield function definition and f_i is the
C     corresponding flux (dot \alpha_i = \dot lambda f_i).
C     The particular cases are given below for each NCRIT.
C
C     The inclusion of damage must be carried out for every model!!!!!
C
C
C
C     A particular case of HARDS is:
C     - \partial F / \partial Cp * H_Cp
C     such that dot Cp = dot lambda * H_Cp
C
C     For the (temperature-dependent or not) Von Mises yield function, 
C     \partial F / \partial Cp = -1 (or A2COE=1-fpc^ls).
C
C     ISOTT defines the isotropic hardening model considered:
C
C     MODEL 1 (work hardening; i.e., thesis):
C             Knowing the evolution law of Cp (see plvari.f),
C             HARDS coincides with H_Cp = h_Cp \sigma:R,
C             where R=\partial G / \partial \sigma.
C             (In the improved model of thesis,
C             \partial F / \partial Cp = -(1-f_pc) =>
C             HARDS will be: (1-f_pc) H_Cp with the same H_cp.
C             This is taken into account by means of A2COE in the
C             version 2 of Gurson criterion. This is also valid for
C             the other MODELS).
C
C     MODEL 2 (linear strain hardening):
C             Knowing the evolution law of Cp (see plvari.f),
C             HARDS coincides with H_Cp = h_Cp \sqrt(2/3 R:R),
C             where R=\partial G / \partial \sigma.
C
C     MODEL 3 (version V2):
C             Once more, \partial F / \partial Cp = -1.
C             Knowing that the closed form of Cp is:
C             Cp = Q (1-exp(-b \barEp)) => 
C             \dot Cp = (HARDS) \dot lambda + A_Cp \dot T, where
C             HARDS is: Q b exp(-b \barEp) \sqrt(2/3 R:R).
C             Q & b can be temperature-dependent.
C             A_Cp involves temperature derivatives of Q & b.
C
C     MODEL 4 (linear version of MODEL 3):
C             HARDS coincides with H_Cp = h_Cp \sqrt(2/3 R:R)
C             of MODEL 2.
C
C     MODEL 5 (linear strain hardening):
C             Knowing the evolution law of Cp (see plvari.f),
C             HARDS coincides with H_Cp = h_Cp \sqrt(R:R),
C             where R=\partial G / \partial \sigma.
C             This model is not consistent with the Von Mises yield
C             function
C
C     MODEL 6 (non-linear strain hardening; non-linear version of
C             MODEL 2):
C             Once more, \partial F / \partial Cp = -1.
C             Knowing that the closed form of Cp is:
C             Cp = A \barEp^n =>
C             \dot Cp = (HARDS) \dot lambda + A_Cp \dot T, where
C             HARDS is: A n \barEp^(n-1) \sqrt(2/3 R:R).
C             A & n can be temperature-dependent.
C             A_Cp involves temperature derivatives of A & n.
C
C     MODEL 7 (non-linear strain hardening; similar to MODEL=6):
C             Once more, \partial F / \partial Cp = -1.
C             Knowing that the closed form of Cp is:
C             Cp = C = A (\barEp_o+\barEp)^n =>
C             \dot Cp = (HARDS) \dot lambda + A_Cp \dot T, where
C             HARDS is: A n (\barEp_o+\barEp)^(n-1) \sqrt(2/3 R:R).
C             A & n can be temperature-dependent.
C             A_Cp involves temperature derivatives of A & n.
C
C     MODEL 8 (model 3 + model 4):
C             Once more, \partial F / \partial Cp = -1.
C             Knowing that the closed form of Cp is:
C             Cp = h_Cp \barEp + Q (1-exp(-b \barEp)) =>
C             \dot Cp = (HARDS) \dot lambda + A_Cp \dot T, where
C             HARDS is: (h_Cp + Q b exp(-b \barEp)) \sqrt(2/3 R:R).
C             h_Cp, Q & b can be temperature-dependent.
C             A_Cp involves temperature derivatives of Q & b.
C
C     MODEL 9 Knowing that the closed form is C = f(\barEp,T), then
C             HARDS is dC/d\barEp \sqrt(2/3 R:R).
C
C     MODEL 10 Johnson & Cook hardening model (strain, strain rate &
C              temperature-dependent)
C              CCOEF, CCOEB & CCOEQ: H, n & A parameters
C              Once more, \partial F / \partial Cp = -1.
C              Knowing that the closed form of Cp is:
C              Cp = C = (Cth + A \barEp^n)*(1 + H ln \dot\barEp) =>
C              \dot Cp = (HARDS) \dot lambda + A_Cp \dot T + 
C                        A_Dp \dot\dot lambda, where
C              HARDS is: A n \barEp^(n-1) \sqrt(2/3 R:R)*
C                                                  (1 + H ln \dot\barEp)
C              The temperature-dependency is included in Cth & A.
C              H & n are normally temperature-independent.
C              A_Cp involves temperature derivatives of Cth & A.
C              A_Dp involves \dot lambda derivatives of Cp.
C
C     MODEL 11 Idem MODEL 10 with the strain hardening term of MODEL 7
C
C     MODEL 12 (extension of MODEL=7 with a critical effec. pl. strain):
C              Once more, \partial F / \partial Cp = -1.
C              Knowing that the closed form of Cp is:
C              Cp = C = A (\barEp_o+<\barEp-\barEp_c>)^n =>
C              \dot Cp = (HARDS) \dot lambda + A_Cp \dot T, where
C              HARDS is: A n (\barEp_o+<\barEp-\barEp_c>)^(n-1)
C              \sqrt(2/3 R:R).
C              A & n can be temperature-dependent.
C              A_Cp involves temperature derivatives of A & n.
C
C     MODEL 13 Idem MODEL 6 with the consideration of a critical
C              effective plastic strain in order to model the effect
C              of plastification without hardening (Idem MODEL 12)
C
C     When performing the product VECTG:VECTG, the fact that factor "2"
C     is included in VECTG must be taken into account.
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION SIGMA(*), VECTF(*), VECTG(*), DMATP(6,*), DESIG(*),
     .          HARDE(*), BACKS(*)
      DIMENSION UNOMA(6)
C
      DO ISTRS=1,NSTRS
       UNOMA(ISTRS)=0.0D0
      ENDDO
      UNOMA(1)=1.0D0
      UNOMA(2)=1.0D0
      UNOMA(4)=1.0D0
C
      DO ISTR1=1,NSTR1                                ! used in plrigi.f
       HARDE(ISTR1)=0.0D0
      ENDDO
C
C**** YIELD CRITERION CHOICE
C
      GO TO (31,32,33,34,35,36,37,38,39,40,41,42) (NCRIT-30)
C
C**** TRESCA
C
   31 CONTINUE
      CALL RUNEND('ERROR IN PLHARD       ')
      RETURN
C
C**** VON MISES
C
   32 CONTINUE
C
      HARDS=0.0D0
      IF(ISOTT.EQ.1) THEN
       DO ISTRS=1,NSTRS
        HARDS=HARDS+A2COE*CCOEF*SIGMA(ISTRS)*VECTG(ISTRS)
       ENDDO
      ENDIF
      IF(ISOTT.EQ.2.OR.ISOTT.EQ.4) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.3) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*EXP(-CCOEB*EFFPL)*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.5) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*HARDS
      ENDIF
      IF(ISOTT.EQ.6) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDX=0.0D0                               ! to avoid 0^n => NAN
       IF(EFFPL.GT.0.0D0) HARDX=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))
       HARDS=HARDX*DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                            DSQRT(2.0D0/3.0D0*HARDS)*(1.0D0-DAMAG)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*(CCOEF+CCOEQ*CCOEB*EXP(-CCOEB*EFFPL))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.9) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       CALL CCOEFET(CCOEFE,DCCOEFE,
     .              TEMPG,EFFPL,PROPS)
       HARDS=A2COE*DCCOEFE*DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.10) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDX=0.0D0                               ! to avoid 0^n => NAN
       IF(EFFPL.GT.0.0D0) HARDX=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))
       HARDS=HARDX*DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.11) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                            DSQRT(2.0D0/3.0D0*HARDS)*(1.0D0-DAMAG)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.12) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       IF(EFFPL.GE.CCOEF) THEN
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL-CCOEF)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ELSE                    ! approximation to improve convergence ??
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ENDIF
      ENDIF
      IF(ISOTT.EQ.13) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDX=0.0D0
c      IF(EFFPL.GT.CCOEF)                ! tangent: it does not converge
c    .  HARDX=A2COE*CCOEQ*CCOEB*((EFFPL-CCOEF)**(CCOEB-1.0D0))
       DEFO2=CCOEF                       ! secant: it converges!
       HARD2=0.0D0
       IF(EFFPL.GT.CCOEF) THEN
        DEFO2=EFFPL
        HARD2=A2COE*CCOEQ*(EFFPL-CCOEF)**CCOEB
       ENDIF
       DEFO1=CCOEF
       HARD1=0.0D0
       IF((EFFPL-EFFP1).GT.CCOEF) THEN
        DEFO1=EFFPL-EFFP1
        HARD1=A2COE*CCOEQ*(EFFPL-EFFP1-CCOEF)**CCOEB
       ENDIF
       IF(DEFO2.GT.DEFO1) THEN
        HARDX=(HARD2-HARD1)/(DEFO2-DEFO1)
       ELSE
        HARDX=CCOEB*CCERO/CCOEF   ! approximation to improve convergence
       ENDIF
       HARDS=HARDX*DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
C
      IF(IKINE.EQ.3) THEN
       DO ISTRS=1,NSTRS
        HARDS=HARDS+VECTF(ISTRS)*(CKOEB*VECTG(ISTRS)-CKOEQ*BACKS(ISTRS))
       ENDDO
      ENDIF
      RETURN
C
C**** MOHR-COULOMB & VERSION  WITH "TENSION CUT-OFF"
C
   33 CONTINUE
      CALL RUNEND('ERROR IN PLHARD       ')
      RETURN
C
C**** DRUCKER-PRAGER
C
   34 CONTINUE
      CALL RUNEND('ERROR IN PLHARD       ')
      RETURN
C
C**** J. LUBLINER'S THEORY
C
   35 CONTINUE
      CALL RUNEND('ERROR IN PLHARD       ')
      RETURN
C
C**** ABOUAF'S MODEL
C
   36 CONTINUE
      CALL RUNEND('ERROR IN PLHARD       ')
      RETURN
C
C**** WEBER BROWN'S MODEL
C
   37 CONTINUE
      CALL RUNEND('ERROR IN PLHARD       ')
      RETURN
C
C**** SG CAST IRON
C
   38 CONTINUE
C
      HARDS=0.0D0
      IF(ISOTT.EQ.1) THEN
       DO ISTRS=1,NSTRS
        HARDS=HARDS+A2COE*CCOEF*SIGMA(ISTRS)*VECTG(ISTRS)
       ENDDO
      ENDIF
      IF(ISOTT.EQ.2.OR.ISOTT.EQ.4) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.3) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*EXP(-CCOEB*EFFPL)*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.5) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*HARDS
      ENDIF
      IF(ISOTT.EQ.6) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*(CCOEF+CCOEQ*CCOEB*EXP(-CCOEB*EFFPL))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.10) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDX=0.0D0                               ! to avoid 0^n => NAN
       IF(EFFPL.GT.0.0D0) HARDX=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))
       HARDS=HARDX*DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.11) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.12) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       IF(EFFPL.GE.CCOEF) THEN
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL-CCOEF)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ELSE
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ENDIF
      ENDIF
      RETURN
C
C**** GREEN SAND
C
   39 CONTINUE
C
      HARDS=0.0D0
      IF(ISOTT.EQ.1) THEN
       DO ISTRS=1,NSTRS
        HARDS=HARDS+A2COE*CCOEF*SIGMA(ISTRS)*VECTG(ISTRS)
       ENDDO
      ENDIF
      IF(ISOTT.EQ.2.OR.ISOTT.EQ.4) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.3) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*EXP(-CCOEB*EFFPL)*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.5) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*HARDS
      ENDIF
      IF(ISOTT.EQ.6) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*(CCOEF+CCOEQ*CCOEB*EXP(-CCOEB*EFFPL))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.10) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDX=0.0D0                               ! to avoid 0^n => NAN
       IF(EFFPL.GT.0.0D0) HARDX=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))
       HARDS=HARDX*DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.11) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.12) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       IF(EFFPL.GE.CCOEF) THEN
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL-CCOEF)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ELSE
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ENDIF
      ENDIF
      RETURN
C
C**** GREEN SAND
C
   40 CONTINUE
C
      HARDS=0.0D0
      IF(ISOTT.EQ.1) THEN
       DO ISTRS=1,NSTRS
        HARDS=HARDS+A2COE*CCOEF*SIGMA(ISTRS)*VECTG(ISTRS)
       ENDDO
      ENDIF
      IF(ISOTT.EQ.2.OR.ISOTT.EQ.4) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.3) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*EXP(-CCOEB*EFFPL)*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.5) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*HARDS
      ENDIF
      IF(ISOTT.EQ.6) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*(CCOEF+CCOEQ*CCOEB*EXP(-CCOEB*EFFPL))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.10) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDX=0.0D0                               ! to avoid 0^n => NAN
       IF(EFFPL.GT.0.0D0) HARDX=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))
       HARDS=HARDX*DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.11) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.12) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       IF(EFFPL.GE.CCOEF) THEN
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL-CCOEF)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ELSE
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ENDIF
      ENDIF
      RETURN
C
C**** GURSON
C
   41 CONTINUE
C
C     For the (temperature-dependent or not) Gurson yield function, 
C     \partial F / \partial Cp = -1.
C
C     MODELS 1, 3 & 4 to be revised !!!!!!! (Feb/1997)
C
      HARDS=0.0D0
      IF(ISOTT.EQ.1) THEN
       DO ISTRS=1,NSTRS
        HARDS=HARDS+A2COE*CCOEF*SIGMA(ISTRS)*VECTG(ISTRS)
       ENDDO
      ENDIF
      IF(ISOTT.EQ.2.OR.ISOTT.EQ.4) THEN
       DO ISTRS=1,NSTRS
        HARDS=HARDS+A2COE*CCOEF/((1.0D0-POROS)*YIELD)*SIGMA(ISTRS)*
     .                                                      VECTG(ISTRS)
       ENDDO
      ENDIF
      IF(ISOTT.EQ.3) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*EXP(-CCOEB*EFFPL)*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.5) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*HARDS
      ENDIF
      IF(ISOTT.EQ.6) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*(CCOEF+CCOEQ*CCOEB*EXP(-CCOEB*EFFPL))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.10) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDX=0.0D0                               ! to avoid 0^n => NAN
       IF(EFFPL.GT.0.0D0) HARDX=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))
       HARDS=HARDX*DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.11) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.12) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       IF(EFFPL.GE.CCOEF) THEN
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL-CCOEF)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ELSE
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ENDIF
      ENDIF
C
C**** KINEMATIC HARDENING CONTRIBUTION IS LACKING !!! (Feb/11997)
C

C
C**** POROSITY CONTRIBUTION
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
C     Note: PREYS comes from pltoha.f
C           PMEAN comes from plinva.f
C           PREYA comes from plreco.f
C
      Q1PAM=1.0D+00
      Q2PAM=1.0D+00
      IF(IPORO.EQ.4) THEN
       Q1PAM=CPOE1
       Q2PAM=CPOE2
      ENDIF
      WWWAU=(3.0D+00*PMEAN*Q2PAM)/(2.0D+00*PREYS)
      WWWGG=1.0D+00+Q1PAM*Q1PAM*POROS*POROS-
     .                                  2.0D+00*Q1PAM*POROS*DCOSH(WWWAU)
C
C**** DERIVATIVE OF F WITH RESPECT TO POROS
C
      DFRFG=-0.5D0*YIELD/WWWGG*(2.0D0*POROS-2.0D0*DCOSH(WWWAU))
C
      IF(IPORO.EQ.1.OR.IPORO.EQ.3) THEN                     ! nucleation
       ALFAG=0.0D0
       BETAG=0.0D0
       HEAVI=0.0D0
       POAUX=1.0D0-POROS
       DO ISTRS=1,NSTRS
        HEAVI=HEAVI+UNOMA(ISTRS)*DESIG(ISTRS)
       ENDDO
       HEAVI=HEAVI/3.0D0+PREYS-PREYA
       IF(HEAVI.GT.0.0D0) THEN
        ALFAG=DFRFG*CPOEF/(3.0D0*PREYS)
        BETAG=DFRFG*A2COE*CCOEF*CPOEF/(PREYS*POAUX*YIELD)
       ENDIF
       DO ISTRS=1,NSTRS
        DO JSTRS=1,NSTRS
         HARDE(ISTRS)=HARDE(ISTRS)+ALFAG*UNOMA(JSTRS)*DMATP(JSTRS,ISTRS)
        ENDDO
        HARDS=HARDS+(HARDE(ISTRS)-BETAG*SIGMA(ISTRS))*VECTG(ISTRS)      
       ENDDO
      ENDIF
      IF(IPORO.EQ.4) THEN                                   ! nucleation


      ENDIF
      IF(IPORO.EQ.2.OR.IPORO.EQ.3.OR.IPORO.EQ.4) THEN       ! growth
       POAUX=1.0D0-POROS
       BETAG=0.0D0
       DO ISTR1=1,NSTR1
        BETAG=BETAG+POAUX*UNOMA(ISTR1)*VECTG(ISTR1)
       ENDDO
       HARDS=HARDS-DFRFG*BETAG
      ENDIF
      RETURN
C
C**** HILL 48
C
   42 CONTINUE
C
C     For the (temperature-dependent or not) Von Mises yield function, 
C     \partial F / \partial Cp = -1 (or A2COE=1-fpc^ls).
C
      HARDS=0.0D0
      IF(ISOTT.EQ.1) THEN
       DO ISTRS=1,NSTRS
        HARDS=HARDS+A2COE*CCOEF*SIGMA(ISTRS)*VECTG(ISTRS)
       ENDDO
      ENDIF
      IF(ISOTT.EQ.2.OR.ISOTT.EQ.4) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.3) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEQ*CCOEB*EXP(-CCOEB*EFFPL)*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.5) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*CCOEF*HARDS
      ENDIF
      IF(ISOTT.EQ.6) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDX=0.0D0                               ! to avoid 0^n => NAN
       IF(EFFPL.GT.0.0D0) HARDX=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))
       HARDS=HARDX*DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.7) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.8) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDS=A2COE*(CCOEF+CCOEQ*CCOEB*EXP(-CCOEB*EFFPL))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
      ENDIF
      IF(ISOTT.EQ.10) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       HARDX=0.0D0                               ! to avoid 0^n => NAN
       IF(EFFPL.GT.0.0D0) HARDX=A2COE*CCOEQ*CCOEB*(EFFPL**(CCOEB-1.0D0))
       HARDS=HARDX*DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.11) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       IF(EFFP2.GT.0.0D0) HARDX=1.0D0+CCOEF*DLOG(EFFP2)
       IF(HARDX.GT.1.0D0) HARDS=HARDS*HARDX      ! approx.; see plvari.f
      ENDIF
      IF(ISOTT.EQ.12) THEN
       DO ISTRS=1,NSTRS
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
         HARDS=HARDS+0.5D0*VECTG(ISTRS)*VECTG(ISTRS)
        ELSE
         HARDS=HARDS+VECTG(ISTRS)*VECTG(ISTRS)
        ENDIF
       ENDDO
       EFFPL0=(CCERO/CCOEQ)**(1.0D0/CCOEB)
       IF(EFFPL.GE.CCOEF) THEN
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0+EFFPL-CCOEF)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ELSE
        HARDS=A2COE*CCOEQ*CCOEB*((EFFPL0)**(CCOEB-1.0D0))*
     .                                          DSQRT(2.0D0/3.0D0*HARDS)
       ENDIF
      ENDIF
      RETURN
C
      END
