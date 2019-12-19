      SUBROUTINE SMOF05T(DVOLUT,EMASST,SHAPET,STRANT,PROPST,
     .                   SFISBT,SFISAT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE ELEMENTAL CONTRIBUTION TO NODAL
C     PRINCIPAL HEAT FLUXES OR TEMP. GRADIENT
C
C     Not implemented !!
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
      DIMENSION DVOLUT(*),             SHAPET(NNODST,*),
     .          STRANT(NSTR1T,*),      EMASST(NNODST,*),
     .          PROPST(*)
      DIMENSION SFISBT(NSTR1T,*),      SFISAT(NNODST,*)
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
       DO ISTR1T=1,NSTR1T
C$DIR SCALAR
        DO INODET=1,NNODST
         SFISBT(ISTR1T,INODET)=STRANT(ISTR1T,1)
        ENDDO
       ENDDO       ! istr1t=1,nstr1t
       RETURN
      ENDIF
C
C**** INTEGRAL PERFORMING
C
      DO ISTR1T=1,NSTR1T
C$DIR SCALAR
       DO INODET=1,NNODST
        DO IGAUST=1,NGAULT
         SFISAT(INODET,ISTR1T)=SFISAT(INODET,ISTR1T)+
     .                         SHAPET(INODET,IGAUST)*
     .                         STRANT(ISTR1T,IGAUST)*DVOLUT(IGAUST)
        ENDDO
       ENDDO
      ENDDO       !istr1t=1,nstr1t
C
C**** PERFORMS MATRIX PRODUCT
C
      DO ISTR1T=1,NSTR1T
C$DIR SCALAR
       DO INODET=1,NNODST
        DO JNODET=1,NNODST
         SFISBT(ISTR1T,INODET)=SFISBT(ISTR1T,INODET)+
     .                      SFISAT(JNODET,ISTR1T)*EMASST(INODET,JNODET)
        ENDDO
       ENDDO
      ENDDO       ! istr1t=1,nstr1t
C
      RETURN
      END
