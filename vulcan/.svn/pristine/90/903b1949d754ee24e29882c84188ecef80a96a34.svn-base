      SUBROUTINE SMOMTXT(ELCODT,DVOLUT,SHAPET,EMASST,GPCODT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE ELEMENT SMOOTHING MATRIX AND INVERTS IT
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION ELCODT(NDIMET,*), EMASST(NNODST,*), GPCODT(NDIMET,*),
     .          SHAPET(NNODST,*), DVOLUT(*)
C
C**** FOR 'ONE-POINT' INTEGRATION RULES
C
      IF(NGAULT.EQ.1) RETURN
C
C**** "MASS" MATRIX
C
      RADIUS=1.0D0
      DO INODET=1,NNODST
       DO JNODET=INODET,NNODST
        EMASST(JNODET,INODET)=0.0D0
        DO IGAULT=1,NGAULT
         IF(NTYPET.EQ.3) THEN
          IF(ISMO1T.EQ.2) THEN
           GPCOA=0.0D0               ! average radial spatial coordinate
           DO IGAUXT=1,NGAULT
            GPCOA=GPCOA+GPCODT(1,IGAUXT)
           ENDDO
           GPCOA=GPCOA/NGAULT
           IF(GPCODT(1,IGAULT).GT.(1.0D-10*GPCOA)) THEN
            RADIUS=1.0D0/GPCODT(1,IGAULT)
           ELSE
            RADIUS=2.0D0/GPCOA
           ENDIF
          ENDIF
         ENDIF
         EMASST(JNODET,INODET)=EMASST(JNODET,INODET)+
     .                         SHAPET(INODET,IGAULT)*
     .                         SHAPET(JNODET,IGAULT)*DVOLUT(IGAULT)*
     .                         RADIUS
        ENDDO
       ENDDO
      ENDDO
C
      DO INODET=2,NNODST
       DO JNODET=1,INODET-1
        EMASST(JNODET,INODET)=EMASST(INODET,JNODET)
       ENDDO
      ENDDO
C
      IF(ISMO2T.EQ.1) CALL LUMPMMT(EMASST,NNODST,     0)
C
      CALL INVERT(EMASST,NNODST,NNODST)
C
      RETURN
      END
