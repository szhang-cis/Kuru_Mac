      SUBROUTINE SMOG05S(DVOLUS,EMASSS,SHAPES,STRSGS,PROPSS,
     .                   STREBS,STREAS)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE ELEMENTAL CONTRIBUTION TO NODAL HEAT
C     FLUXES
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      COMMON/SMOGAUS/FACTAS
C
      DIMENSION DVOLUS(*),             SHAPES(NNODSS,*),
     .          STRSGS(NSTR1S,*),      EMASSS(NNODSS,*),
     .          PROPSS(*)
      DIMENSION STREBS(NSTR1S,*),      STREAS(NNODSS,*)
C
C**** COMPUTE TOTAL VOLUME IF NECESSARY
C
      IF(IABS(KSGAUS).EQ.2) THEN       ! local smoothing
       FACTAS=0.D+00
       DO IGAUSS=1,NGAULS
        FACTAS=FACTAS+DVOLUS(IGAUSS)
       ENDDO
      ENDIF
C
C**** FOR 'ONE-POINT' INTEGRATION RULES
C
      IF(NGAULS.EQ.1) THEN
       DO ISTR1S=1,NSTR1S
C$DIR SCALAR
        DO INODES=1,NNODSS
         STREBS(ISTR1S,INODES)=STRSGS(ISTR1S,1)
        ENDDO
       ENDDO       ! istr1s=1,nstr1s
       RETURN
      ENDIF
C
C**** INTEGRAL PERFORMING
C
      DO ISTR1S=1,NSTR1S
C$DIR SCALAR
       DO INODES=1,NNODSS
        DO IGAUSS=1,NGAULS
         STREAS(INODES,ISTR1S)=STREAS(INODES,ISTR1S)+
     .                         SHAPES(INODES,IGAUSS)*
     .                         STRSGS(ISTR1S,IGAUSS)*DVOLUS(IGAUSS)
        ENDDO
       ENDDO
      ENDDO       !istr1s=1,nstr1s
C
C**** PERFORMS MATRIX PRODUCT
C
      DO ISTR1S=1,NSTR1S
C$DIR SCALAR
       DO INODES=1,NNODSS
        DO JNODES=1,NNODSS
         STREBS(ISTR1S,INODES)=STREBS(ISTR1S,INODES)+
     .                      STREAS(JNODES,ISTR1S)*EMASSS(INODES,JNODES)
        ENDDO
       ENDDO
      ENDDO       ! istr1s=1,nstr1s
C
      RETURN
      END
