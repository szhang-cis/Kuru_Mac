      SUBROUTINE CALCBT(DMAPL,COUTD,DMATX,UNOMA,BULKM,DAMAX,
     .                  STRAN,STRAP,STRAT,RCGTT,RCGTI,STRRD)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE TANGENT CONJUGATE OF THERMAL DILATATION
C     TENSOR (TANGENT BETA)
C
C     Notes:
C
C     This tensor is useful to calculate the coupling term for the
C     thermal: plastic + volumetric energy terms
C
C     For constant thermal dilatation coefficient and constitutive
C     tensor (i.e., non-temperature-dependent), the tangent BETA
C     coincides with the secant BETA. Even for temperature-dependent
C     thermal dilatation coefficient and constitutive tensor, a secant
C     BETA is computed in this routine
C
C     For LARGE=1, FACT2=1.0 for any LARGET as the coupling term is
C     integrated in the original configuration (see plheat.f) and
C     dV_o=rho/rho_o*dV=1/DETJM*dV
C
C     For LARGE=2, FACT2=1.0/DETJMU for any LARGET (DETJMU=rho_o/rho_u)
C
C     For LARGE=3, FACT2=1.0/DETJM for any LARGET (DETJM=rho_o/rho)
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
      DIMENSION DMAPL(*), DMATX(NSTRS,*), UBULK(6), UNOMA(*), UNOMX(6)
      DIMENSION STRAN(*), STRAP(*),       STRAT(*), RCGTT(*), RCGTI(*)
C
      IF(LARGE.EQ.3) RETURN
C
      ISTAN=2
      IF(LARGE.EQ.0) THEN
       ISTAN=1
      ELSE
       IF(IFREN.EQ.1.OR.IFREN.EQ.4.OR.IFREN.EQ.9.OR.IFREN.EQ.10) ISTAN=1
      ENDIF
C
      FACT2=1.0D0
C
      IF(ISTAN.EQ.1) THEN
       DO ISTRS=1,NSTRS
        UBULK(ISTRS)=0.0D0
        DO JSTRS=1,NSTRS
         UBULK(ISTRS)=UBULK(ISTRS)+(1.0D0-DAMAX)*DMATX(ISTRS,JSTRS)*
     .                                                      UNOMA(JSTRS)
        ENDDO
        UNOMX(ISTRS)=UNOMA(ISTRS)
       ENDDO
      ELSE               ! istan=2
       IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .    IFREN.EQ.5.OR.IFREN.EQ.6) THEN
        DO ISTRS=1,NSTRS
         UBULK(ISTRS)=0.0D0
         DO JSTRS=1,NSTRS
          UBULK(ISTRS)=UBULK(ISTRS)+(1.0D0-DAMAX)*DMATX(ISTRS,JSTRS)*
     .                                                      RCGTT(JSTRS)
         ENDDO
         UNOMX(ISTRS)=RCGTT(ISTRS)
        ENDDO
       ENDIF             ! ifren.eq.2.or.ifren.eq.3.or....
       IF(IFREN.EQ.7.OR.IFREN.EQ.8) THEN
        DO ISTRS=1,NSTRS
         UBULK(ISTRS)=3.0D0*(1.0D0-DAMAX)*BULKM*RCGTI(ISTRS)
         UNOMX(ISTRS)=RCGTT(ISTRS)
        ENDDO
       ENDIF             ! ifren.eq.7.or.ifren.eq.8
      ENDIF              ! istan.eq.1
C
C**** COMPUTES BETA TENSOR
C
      IF(LARGE.NE.0) THEN
       IF(LARGE.EQ.2) FACT2=1.0D0/DETJMU       ! DETJMU must be computed
       IF(LARGE.EQ.3) FACT2=1.0D0/DETJM        ! DETJM must be computed
      ENDIF                   ! large.ne.0
C
      DO ISTR1=1,NSTR1
       DMAPL(ISTR1)=0.0D0
      END DO
      DO ISTRS=1,NSTRS
       DMAPL(ISTRS)=STRRD*UBULK(ISTRS)                     ! BETA tensor
       DMAPL(ISTRS)=DMAPL(ISTRS)*FACT2     ! incorporates density change
      ENDDO
C
C**** COUPLING VARIABLE TO BE USED IN THE ISOTHERMAL IMPROVED SPLITS
C
C     Note: the density change for LARGE > 0 is included in capcoft.f
C
      IF(NITERC.EQ.1.OR.NITERC.EQ.2.OR.NITERC.EQ.4) THEN
       COUTD=0.0D0
       DO ISTRS=1,NSTRS
        COUTD=COUTD+DMAPL(ISTRS)*STRRD*UNOMX(ISTRS)/FACT2
        IF (NITERC.EQ.4) THEN
         COUTD=COUTD+DMAPL(ISTRS)*2.0D0*STRRD*
     .                    (STRAN(ISTRS)-STRAP(ISTRS)-STRAT(ISTRS))/FACT2
        ENDIF                ! niterc.eq.4
       END DO
      ENDIF                  ! niterc.eq.1.or.niterc.eq.2.or....
C
      RETURN
#endif
      END
