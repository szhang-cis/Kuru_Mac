      SUBROUTINE SMOG05T(DVOLUT,EMASST,SHAPET,STRSGT,PROPST,
     .                   STREBT,STREAT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE ELEMENTAL CONTRIBUTION TO NODAL HEAT
C     FLUXES
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
     .          STRSGT(NSTR1T,*),      EMASST(NNODST,*),
     .          PROPST(*)
      DIMENSION STREBT(NSTR1T,*),      STREAT(NNODST,*)
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
         STREBT(ISTR1T,INODET)=STRSGT(ISTR1T,1)
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
         STREAT(INODET,ISTR1T)=STREAT(INODET,ISTR1T)+
     .                         SHAPET(INODET,IGAUST)*
     .                         STRSGT(ISTR1T,IGAUST)*DVOLUT(IGAUST)
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
         STREBT(ISTR1T,INODET)=STREBT(ISTR1T,INODET)+
     .                      STREAT(JNODET,ISTR1T)*EMASST(INODET,JNODET)
        ENDDO
       ENDDO
      ENDDO       ! istr1t=1,nstr1t
C
      RETURN
      END
