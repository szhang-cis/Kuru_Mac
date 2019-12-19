      SUBROUTINE PLHEAT(PWOEL,SHAPE,DSTRA,DEEPP,DEETH,DMAPL,SGTOT,TEMPG,
     .                  DBOLU,COUTD,COUTH)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE MECHANICAL COUPLING TERMS FOR THE
C     THERMAL PROBLEM (PLASTIC + VOLUME CHANGE WORKS)
C
C
C     COUFAC: COMPATIBILIZES ENERGY RATE UNITIES USED IN THERMAL AND
C             MECHANICAL PROBLEMS
C             (En. Ra. Unit. of Ther. Prob. / En. Ra. Unit. of Mech.
C             Prob.)
C
C     NITERC=0 >> STANDARD STAGGERED SCHEME
C           =1 >> IMPROVED STAGGERED SCHEME
C                 (isothermal split at Gauss points)
C           =2 >> IMPROVED STAGGERED SCHEME
C                 (isothermal split at nodes [smoothed])
C           =3 >> ADIABATIC STAGGERED SCHEME
C           =4 >> IMPROVED STAGGERED SCHEME (other version of 1)
C
C     PWOEL: MECHANICAL COUPLING TERM
C     SHAPE: SHAPE FUNCTION
C     DSTRA: INCREMENTAL TOTAL STRAIN
C     DEEPP: INCREMENTAL PLASTIC STRAIN
C     DEETH: INCREMENTAL THERMAL STRAIN
C     DMAPL: CONGUGATE OF THERMAL DILATATION TENSOR
C     SGTOT: STRESS TENSOR
C     TEMPG: TEMPERATURE
C     DBOLU: DIFFERENCIAL VOLUME
C     COUTD: COUPLING VARIABLE (T~\alpha_th:C:\alpha_th with 
C            \alpha_th: tangent thermal dilatation tensor)
C            Term used in the improved staggered scheme and 
C            transferred to the thermal problem in trasmeg.f (NITERC=1)
C     COUTH: ENERGY STORED AS HARDENING
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
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION PWOEL(NNODL), SHAPE(NNODL), DSTRA(*), SGTOT(*)
      DIMENSION DEEPP(6),     DMAPL(6),     DEETH(6)
C
      COUTE=0.0D+00           ! thermoelastic coupling term
      COUTA=0.0D+00           ! thermoplastic coupling term
C
C**** STANDARD & ADIABATIC STAGGERED SCHEMES
C
      IF(NITERC.EQ.0.OR.NITERC.EQ.3) THEN
       DO ISTRS=1,NSTRS
        COUTE=COUTE+DMAPL(ISTRS)*(DSTRA(ISTRS)-DEEPP(ISTRS))/DTIME
        COUTA=COUTA+SGTOT(ISTRS)*DEEPP(ISTRS)/DTIME
       END DO
       COUTA=COUTA*FPLCOU
       IF(NPLCOU.EQ.2) COUTA=COUTA+COUTH
       IF(NPLCOU.EQ.3) COUTE=0.0D0
C
       ISIMP=1                       ! better as input (see capcoft.f)
       IF(ISIMP.EQ.0) THEN           ! simplified thermoelasticity
        COUTE=COUTE*CENKEL
       ELSE                          ! general
        COUTE=COUTE*(TEMPG+CENKEL)
       ENDIF
      ENDIF           ! niterc.eq.0.or.niterc.eq.3
C
C**** IMPROVED STAGGERED SCHEME (ISOTHERMAL)
C
      IF(NITERC.EQ.1.OR.NITERC.EQ.2.OR.NITERC.EQ.4) THEN
       DO ISTRS=1,NSTRS
        COUTE=COUTE+DMAPL(ISTRS)*(DSTRA(ISTRS)-DEEPP(ISTRS)-
     .                            DEETH(ISTRS))/DTIME
        COUTA=COUTA+SGTOT(ISTRS)*DEEPP(ISTRS)/DTIME
       END DO
       COUTA=COUTA*FPLCOU
       IF(NPLCOU.EQ.2) COUTA=COUTA+COUTH
       IF(NPLCOU.EQ.3) COUTE=0.0D0
C
       ISIMP=1                       ! better as input (see capcoft.f)
       IF(ISIMP.EQ.0) THEN           ! simplified thermoelasticity
        COUTE=COUTE*CENKEL
        COUTD=COUTD*CENKEL/COUFAC
       ENDIF
       IF(ISIMP.EQ.1) THEN           ! general
        COUTE=COUTE*(TEMPG+CENKEL)
        IF((TEMPG+CENKEL).LT.0.0D0) COUTE=0.0D0            ! lower bound
        COUTD=COUTD*(TEMPG+CENKEL)/COUFAC
        IF((TEMPG+CENKEL).LT.0.0D0) COUTD=0.0D0            ! lower bound
       ENDIF
       IF(ISIMP.EQ.2) THEN    ! more general but inconsistent with coute
        COUTE=COUTE*(TEMPG+CENKEL)
        COUTD=COUTD/COUFAC
       ENDIF
      ENDIF           ! niterc.eq.1.or.niterc.eq.2
C
C**** PLASTIC DISSIPATION > 0
C
      IF(COUTA.LT.0.0D0) THEN
       CALL RUNMEN('PLASTIC DISSIPATION < 0 (SET TO ZERO)')
       COUTA=0.0D0    ! round-off errors can appear for small DEEPP
      ENDIF
C
      DO INODE=1,NNODL
       PWOEL(INODE)=PWOEL(INODE)+SHAPE(INODE)*(-COUTE+COUTA)*DBOLU/
     .              COUFAC
      END DO
C
#endif
      RETURN
      END
