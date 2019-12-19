      SUBROUTINE CALCST(SIGMA,DMATX,STRAN,STRAP,STRAT,DAMAX,
     .                  DETJM,BULKM,DISTM,CPENI,
     .                  FVOL1,STRRA,STRR4,
     .                  RCGTT,RCGTI,TRACC,TRACI,ENEYE,DETJC,SEINC,NAUXI,
     .                  FACTJ,STILD,ATILD,SGTOT,PRESG,DESIG,DSTRA,STRS0,
     .                  RCGPI,TRABE)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE STRESS TENSOR
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION SIGMA(*),       DMATX(NAUXI,*),
     .          STRAN(*),       STRAP(*),       STRAT(*)
      DIMENSION RCGTT(*),       RCGTI(*),
     .          ENEYE(*)
      DIMENSION STILD(*)
      DIMENSION SGTOT(*),       PRESG(*),       DESIG(*), DSTRA(*)
      DIMENSION STRS0(*)
      DIMENSION RCGPI(*)
C
      IF(LARGE.EQ.3) RETURN          ! to be revised !!
C
      ISTAN=2
      IF(LARGE.EQ.0) THEN
       ISTAN=1
      ELSE
#ifndef restricted
       IF(IFREN.EQ.1.OR.IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .    IFREN.EQ.4.OR.IFREN.EQ.5.OR.IFREN.EQ.6.OR.
     .    IFREN.EQ.9.OR.IFREN.EQ.10) ISTAN=1
#endif
      ENDIF
C
C**** INITIALISES STRESS
C
      DO ISTR1=1,NSTR1
       SIGMA(ISTR1)=0.0D0
      ENDDO
C
C**** TOTAL STRESS STATE PREDICTION
C
      IF(ISTAN.EQ.1) THEN          ! standard models
       DO ISTRS=1,NSTRS
        DO JSTRS=1,NSTRS
         SIGMA(ISTRS)=SIGMA(ISTRS)+DMATX(ISTRS,JSTRS)*
     .               (STRAN(JSTRS)-STRAP(JSTRS)-STRAT(JSTRS))
        ENDDO
        SIGMA(ISTRS)=SIGMA(ISTRS)+STRS0(ISTRS)
       ENDDO
       IF(NTYPE.EQ.1) STRR4=0.0D0
#ifndef restricted
      ELSE                         ! istan=2 (non-standard models)
       IF(IFREN.EQ.7) THEN         ! Simo's model
        DO ISTRS=1,NSTRS
         SIGMA(ISTRS)=DISTM*FACTJ*
     .                (RCGPI(ISTRS)-1.0D0/3.0D0*TRABE*RCGTI(ISTRS))+
     .                BULKM*DLOG(DETJM)*RCGTI(ISTRS)-
     .                3.0D0*BULKM*STRRA*RCGTI(ISTRS)+STRS0(ISTRS)
        ENDDO
        IF(NTYPE.EQ.1)             ! DETJC & TRACC redefined in fourth.f
     .   STRR4=0.5D0*(TRACC-RCGTT(1)-RCGTT(2)-1.0D0)
       ENDIF            ! ifren.eq.7
C
       IF(IFREN.EQ.8) THEN
        DO ISTRS=1,NSTRS
         SIGMA(ISTRS)=BULKM/2.0D0*(3.0D0-TRACI)*RCGTI(ISTRS)+
     .                DISTM/3.0D0*TRACI*RCGTI(ISTRS)-
     .                DISTM*ENEYE(ISTRS)-
     .                3.0D0*BULKM*STRRA*RCGTI(ISTRS)+STRS0(ISTRS)
        ENDDO
        IF(NTYPE.EQ.1) STRR4=0.0D0
       ENDIF            ! ifren.eq.8
C
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.
     .    IFREN.EQ.58) THEN        ! M-R, Y, O, D, H, M, G & K
        DO ISTRS=1,NSTRS           ! volumetric term
         SIGMA(ISTRS)=STILD(ISTRS)+
     .                2.0D0*CPENI*DETJC*RCGTI(ISTRS)*FVOL1+STRS0(ISTRS)
        ENDDO
        IF(NTYPE.EQ.1)      ! DETJC, SEINC & TRACC redefined in fourth.f
     .   STRR4=0.5D0*(TRACC-RCGTT(1)-RCGTT(2)-1.0D0)
       ENDIF            ! ifren.eq.51....
#endif
      ENDIF             ! istan.eq.1
C
C**** UPDATE TOTAL & INCREMENTAL STRESS
C
      DO ISTR1=1,NSTR1
       SGTOT(ISTR1)=SIGMA(ISTR1)*(1.0D0-DAMAX)
       DESIG(ISTR1)=SGTOT(ISTR1)-PRESG(ISTR1)   ! presg= ^t strsg
      ENDDO
C
C**** CHECKS CONVEXITY OF ENERGY FUNCTION
C
#ifndef restricted
      IF(ISTAN.EQ.2) THEN          ! istan=2 (non-standard models)
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.
     .    IFREN.EQ.58) THEN        ! M-R, Y, O, D, H, M, G & K
        ASEDE=0.0D0
        DO ISTR1=1,NSTR1
         ASEDE=ASEDE+DESIG(ISTR1)*DSTRA(ISTR1)
        ENDDO
        ASEDEI=-1.0D-18     ! to avoid initial round off negative values
        IF(ASEDE.LT.ASEDEI) THEN
         CALL RUNMEN('WARNING: dS:dE < 0 => LOSS OF CONVEXITY OF ENERGY 
     .FUNCTION                      ')
         WRITE(LURES,*) 'ELEMENT NUMBER=',IELEM,' dS:dE=',ASEDE
        ENDIF
       ENDIF            ! ifren.eq.51....
      ENDIF             ! istan.eq.2
C
#endif
      RETURN
      END
