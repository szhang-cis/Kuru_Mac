      SUBROUTINE SMOX05T(DVOLUT,GPCODT,EMASST,SHAPET,EHISTT,PROPST,
     .                   SFIPBT,SFIPAT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE ELEMENTAL CONTRIBUTION TO NODAL L*f_pc
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      COMMON/SMOGAUT/FACTAT
C
      DIMENSION DVOLUT(*),        SHAPET(NNODST,*), GPCODT(NDIMET,*),
     .          EHISTT(NHISTT,*), EMASST(NNODST,*),
     .          PROPST(*)
      DIMENSION SFIPBT(2,NNODST), SFIPAT(NNODST,2)
C
      IX=NDIMETO-1
C
C**** COMPUTE TOTAL VOLUME IF NECESSARY
C
      IF(IABS(KSGAUT).EQ.2) THEN       ! local smoothing
       FACTAT=0.D+00
       DO IGAUST=1,NGAULT
        FACTAT=FACTAT+DVOLUT(IGAUST)
       ENDDO
      ENDIF
C
C**** FOR 'ONE-POINT' INTEGRATION RULES
C
      IF(NGAULT.EQ.1) THEN
       DO INUINT=1,2
        IPUNST=3+INUINT+IX
        DO INODET=1,NNODST
         SFIPBT(INUINT,INODET)=EHISTT(IPUNST,1)
        ENDDO
       ENDDO
       RETURN
      ENDIF
C
C**** INTEGRAL PERFORMING
C
      RADIUS=1.0D0
      DO INUINT=1,2
       IPUNST=3+INUINT+IX
       DO INODET=1,NNODST
        DO IGAUST=1,NGAULT
         IF(NTYPET.EQ.3) THEN
          IF(ISMO1T.EQ.2) THEN
           GPCOA=0.0D0               ! average radial spatial coordinate
           DO IGAUXT=1,NGAULT
            GPCOA=GPCOA+GPCODT(1,IGAUXT)
           ENDDO
           GPCOA=GPCOA/NGAULT
           IF(GPCODT(1,IGAUST).GT.(1.0D-10*GPCOA)) THEN
            RADIUS=1.0D0/GPCODT(1,IGAUST)
           ELSE
            RADIUS=2.0D0/GPCOA
           ENDIF
          ENDIF
         ENDIF
         SFIPAT(INODET,INUINT)=SFIPAT(INODET,INUINT)+
     .                         SHAPET(INODET,IGAUST)*
     .                         EHISTT(IPUNST,IGAUST)*DVOLUT(IGAUST)*
     .                         RADIUS
        ENDDO
       ENDDO
      ENDDO
C
C**** PERFORMS MATRIX PRODUCT
C
      DO INUINT=1,2
       DO INODET=1,NNODST
        DO JNODET=1,NNODST
         SFIPBT(INUINT,INODET)=SFIPBT(INUINT,INODET)+
     .                       SFIPAT(JNODET,INUINT)*EMASST(INODET,JNODET)
        ENDDO
       ENDDO
      ENDDO
C
      RETURN
      END
