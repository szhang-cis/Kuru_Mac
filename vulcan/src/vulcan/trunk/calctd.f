      SUBROUTINE CALCTD(PROPS,UNOMA,STRAT,DEETH,STRRA,STRRD,
     .                  TEINI,TEMPO,TEMPG,DTEMG,
     .                  ALPHI,ALPHO,ALPHA,
     .                  STRAN,STRAP,RCGTT,YOUNG,POISM)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE THERMAL DEFORMATION
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
      DIMENSION PROPS(*), UNOMA(*), STRAT(*), DEETH(*), UNOMX(6)
      DIMENSION STRAN(*), STRAP(*), RCGTT(*)
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
      PATH1=1.0D0
      PATH2=1.0D0
C
C**** COMPUTES THE THERMAL DEFORMATION
C
      FTROX=1.0D0
      FTRAX=1.0D0
      STRRO=ALPHO*(TEMPO-TEREF)-PATH2/PATH1*ALPHI*(TEINI-TEREF)
      STRRA=ALPHA*(TEMPG-TEREF)-PATH2/PATH1*ALPHI*(TEINI-TEREF)
      STRRD=ALPHA
      DO ISTR1=1,NSTR1
       UNOMX(ISTR1)=UNOMA(ISTR1)
      ENDDO
C
      IF(ISTAN.EQ.1) THEN                       ! Green-Lagrange strains
       IF(LARGE.EQ.0) THEN
        IF(IFREN.EQ.1.OR.IFREN.EQ.9) THEN
         PATH1=YOUNG/(1.0D0-2.0D0*POISM)
         CALL YOUNGT(TEINI,PROPS,
     .               YOUNI,YOUNIF,YOUNIA)
         CALL POISMT(TEINI,PROPS,
     .               POISI,POISIF,POISIA)
         PATH2=YOUNI/(1.0D0-2.0D0*POISI)
         STRRO=ALPHO*(TEMPO-TEREF)-PATH2/PATH1*ALPHI*(TEINI-TEREF)
         STRRA=ALPHA*(TEMPG-TEREF)-PATH2/PATH1*ALPHI*(TEINI-TEREF)
        ENDIF
       ELSE                    ! large > 0
        IF(IFREN.EQ.4.OR.IFREN.EQ.10) THEN
         ALPHI=3.0D0*ALPHI                      ! volumetric coefficient
         ALPHO=3.0D0*ALPHO
         ALPHA=3.0D0*ALPHA
         STRRO=3.0D0*STRRO                      ! volumetric deformation
         STRRA=3.0D0*STRRA
         FTROX=1.0D0+STRRO
         FTRAX=1.0D0+STRRA
         STRRO=0.5D0*(FTROX**(2.0D0/3.0D0)-1.0D0)*FTRCX**(2.0D0/3.0D0)
         STRRA=0.5D0*(FTRAX**(2.0D0/3.0D0)-1.0D0)*FTRBX**(2.0D0/3.0D0)
         STRRD=1.0D0/3.0D0*FTRAX**(-1.0D0/3.0D0)*ALPHA
        ENDIF
       ENDIF
      ELSE                     ! istan=2         (Almansi strains)
       IF(IFREN.EQ.5.OR.IFREN.EQ.6) THEN
        ALPHI=3.0D0*ALPHI                       ! volumetric coefficient
        ALPHO=3.0D0*ALPHO
        ALPHA=3.0D0*ALPHA
        STRRO=3.0D0*STRRO                       ! volumetric deformation
        STRRA=3.0D0*STRRA
        FTROX=1.0D0-STRRO
        FTRAX=1.0D0-STRRA
        STRRO=0.5D0*(1.0D0-FTROX**(2.0D0/3.0D0))
        STRRA=0.5D0*(1.0D0-FTRAX**(2.0D0/3.0D0))
        STRRD=1.0D0/3.0D0*FTRAX**(-1.0D0/3.0D0)*ALPHA
       ENDIF
       DO ISTR1=1,NSTR1        ! transforms to Green-Lagrange strains
        UNOMX(ISTR1)=RCGTT(ISTR1)
       ENDDO
      ENDIF                    ! istan.eq.1
C
      DO ISTR1=1,NSTR1
       STRAT(ISTR1)=STRAT(ISTR1)+STRRA*UNOMX(ISTR1)
      ENDDO
C
C**** COMPUTES INCREMENTAL THERMAL DEFORMATION
C
C     Notes:
C
C     For ISTAN=2, DEETH only considers the temperature variation of
C     the thermal deformation, i.e., 
C     dot(E^th) = d STRRA/d T * C * dot T + STRRA * dot C,
C     where DEETH = d STRRA/d T * C * dot T = (STRRA-STRRO) / DTIME * C
C     (d STRRA/d T * dot T = (STRRA-STRRO) / DTIME = STRRD * dot T)
C     and the factor DTIME is included in plheat.f
C     This fact is considered in calcbt.f
C
C     The computation of STRRD (temperature derivative of STRRA) is
C     performed assuming a non-temperature dependent thermal dilatation
C     coefficient
C
      IF(ITERME.GT.0) THEN                          ! bidirect. problems
       IF(ITERMP.GT.0) THEN
        IF(NITERC.EQ.1.OR.NITERC.EQ.2.OR.NITERC.EQ.4) THEN
         DO ISTR1=1,NSTR1
          DEETH(ISTR1)=DEETH(ISTR1)+(STRRA-STRRO)*UNOMX(ISTR1)
          IF (NITERC.EQ.4) THEN
           DEETH(ISTR1)=DEETH(ISTR1)+2.0D0*STRRD*
     .                    (STRAN(ISTR1)-STRAP(ISTR1)-STRAT(ISTR1))*DTEMG
          ENDIF                ! niterc.eq.4
         ENDDO
        ENDIF                  ! niterc.eq.1.or.niterc.eq.2.or....
       ENDIF                   ! itermp.gt.0
      ENDIF                    ! iterme.gt.0
C
      RETURN
#endif
      END
