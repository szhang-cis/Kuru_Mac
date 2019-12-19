      SUBROUTINE CAEHIST(EHIST,EHISTT)
C***********************************************************************
C
C**** THIS ROUTINE TRANSFERS GAUSSIAN MECHANICAL COUPLING TERM FOR
C     THERMAL PROBLEM
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
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION EHIST(NHIST,*), EHISTT(NHISTT,*)
C
C**** LOOP OVER GAUSS POINTS
C
      IF(NGAUS.EQ.NGAUST) THEN
C
       IX=NDIMETO-1
       DO IGAUST=1,NGAUST
        IGAUS=IGAUST
        EHISTT(IX+7,IGAUST)=EHIST(IPLAS(6),IGAUS)
       END DO
C
      ELSE
C
       CALL RUNMENT
     .  ('WARNING: NGAUS NE NGAUST & NITERC=1,4=>APPROXIM. COMPUTATION')
C
       PICHU=0.0D+00
       DO IGAUS=1,NGAUS
        PICHU=PICHU+EHIST(IPLAS(6),IGAUS)
       END DO
       PICHU=PICHU/NGAUS
C
       DO IGAUST=1,NGAUST
        EHISTT(IX+7,IGAUST)=PICHU
       END DO
C
      ENDIF
C
      RETURN
      END
