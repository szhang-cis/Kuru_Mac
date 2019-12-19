      SUBROUTINE FRIN03(ELDIS,PROPS,BMSIG,TENOD,PWOEL)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING FORCES OF THE LINK 
C     (FREE SURFACE ELEMENT NO. 3)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION ELDIS(NDOFN,NNODL), PROPS(NPROP),
     .          BMSIG(NDOFN*NNODL)
      DIMENSION TENOD(NNODL),       DGAPS(3)
      DIMENSION PWOEL(NNODL)
C
C**** PROPS(1+IDOFN)=PENALTIES PARAMETERS IN X-Y-Z DIRECTIONS
C     PROPS(5)=LIQUIDUS TEMPERATURE
C
      TEMPL=PROPS(5)
C
C**** INITIALISE MECHANICAL COUPLING TERM
C
      IF(ITERME.GT.0) THEN
       DO INODL=1,NNODL
        PWOEL(INODL)=0.0
       END DO
      ENDIF
C
      IF(TENOD(1).GT.TEMPL.AND.TENOD(2).GT.TEMPL) THEN
C
C**** COMPUTE GAPS 
C
       DO IDOFN=1,NDOFN
c       DGAPS(IDOFN)=ELDIS(IDOFN,1)-ELDIS(IDOFN,2)
        DGAPS(IDOFN)=ELDIS(IDOFN,1)-ELDIS(IDOFN,2)-PROPS(6+IDOFN)
       END DO
C
C**** EVALUATE ELEMENT CONTRIBUTION 
C
       DO IDOFN=1,NDOFN
        BMSIG(IDOFN)      = PROPS(1+IDOFN)*DGAPS(IDOFN)
        BMSIG(IDOFN+NDOFN)=-PROPS(1+IDOFN)*DGAPS(IDOFN)
       END DO
C
      ENDIF
C
      RETURN
      END
