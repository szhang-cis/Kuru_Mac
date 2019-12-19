      SUBROUTINE CONPLA(SGTOT,DEVIA,VECTF,DESIG,DMATX,ABETA,DTEMG,
     .                  DCCER,NAUXI)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS SOME PLASTIC CONTROLS
C     
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION SGTOT(*),       DEVIA(6),
     .          VECTF(6),       DESIG(6), 
     .          DMATX(NAUXI,*)
C
C**** VERIFICATES THE PLASTIC DISSIPATION
C
C     In thesis: 4 sufficient conditions:
C                1) \sigma:\sigma' gt 0
C                2) dC^th / dT le 0 (thermal hardening)
C                3) H_Cp=h_Cp \sigma:R gt 0 => h_Cp gt 0 
C                   (H_Cp plastic harden. coefficient)
C                   (h_Cp plastic harden. modulus)
C                4) df_pc / dT gt 0
C
C     Note: only the first condition is verificated
C
      SIDEV=0.0         ! first condition: trivial in thesis formul.
      DO ISTRS=1,NSTRS
       IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
        SIDEV=SIDEV+2.0*SGTOT(ISTRS)*DEVIA(ISTRS)
       ELSE
        SIDEV=SIDEV+SGTOT(ISTRS)*DEVIA(ISTRS)
       ENDIF
      ENDDO
      IF(SIDEV.LT.0.0)
     . CALL RUNEND('ERROR: SIGMA:DEVIA LT 0            ')
C
C**** VERIFICATES THE "MAXIMUM PLASTIC DISSIPATION POSTULATE"
C     (FOR COUPLED & UNCOUPLED [with dtemg=dccer=0] PROBLEMS)
C
C     This postulate is derived from the normality rule for
C     the plastic variables and the Drucker's inequality. 
C     The Drucker's inequality is, at the same time, a consequence 
C     of the normality rule in the plasticity context (see Lubliner)
C
C     In the thesis, the Drucker's inequality reads: 
C         .       .        .   .
C     \sigma:\epsilon^p + T \eta^p ge 0. 
C
C
C     RELAI=measure of plastic dissipation
C     TOLIP=admissible out-of-balance plastic dissipation
C           (specially crucial for plane stress problems)
C
C     Note: In the MPDP, DESIG=\sigma-\sigma^* & DTEMG=T-T^*,
C           where \sigma^* & T^* must satisfy:
C           F(\sigma^*,\alpha_k,T^*) le 0
C                                          .                .
C     In thesis: MDPD=(\sigma-sigma^*):\epsilon^p+(T-T^*) \eta^p
C
C
C to be revised to take into account F(\sigma^*,\alpha_k,T^*) le 0
C
      RELAC=0.0
      RELAI=0.0
      DO ISTRS=1,NSTRS
       IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) THEN
        RELAC=RELAC+2.0*VECTF(ISTRS)*DESIG(ISTRS)
        RELAI=RELAI+2.0*VECTF(ISTRS)*SGTOT(ISTRS)
       ELSE
        RELAC=RELAC+VECTF(ISTRS)*DESIG(ISTRS)
        RELAI=RELAI+VECTF(ISTRS)*SGTOT(ISTRS)
       ENDIF
      ENDDO
      RELAC=RELAC-DTEMG*DCCER      ! dF/dT=-dC^th/dT
      TOLIP=-TOPLA*RELAI           ! (theoretically, tolip=0.0)
      IF(RELAC.LT.TOLIP) THEN
ctm       CALL RUNEND('ERROR: MPDP IS NOT VERIFICATED     ')
       CALL RUNMEN('WARNING: MPDP IS NOT VERIFICATED   ')
      ENDIF
C
C**** VERIFICATES ADDITIONAL PLASTIC CONDITIONS
C
      RRCRR=0.0          ! analyses the limit of hardening coeff. H_Cp
      DO ISTRS=1,NSTRS
       DO JSTRS=1,NSTRS
        RRCRR=RRCRR+VECTF(ISTRS)*DMATX(ISTRS,JSTRS)*VECTF(JSTRS)
       ENDDO
      ENDDO
      IF(RRCRR.LE.0.0)
     . CALL RUNEND('ERROR: R:C:R LE 0                  ')
C
      IF(ABETA.LE.0.0D0)                           ! 1/abeta=H_Cp+R:C:R
     . CALL RUNEND('ERROR: ABETA MULTIPLICATOR LE 0    ')
C
      RETURN
      END
