      SUBROUTINE SMOMTXS(ELCODS,DVOLUS,SHAPES,EMASSS,GPCODS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE ELEMENT SMOOTHING MATRIX AND INVERTS IT
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION ELCODS(NDIMES,*), EMASSS(NNODSS,*), GPCODS(NDIMES,*),
     .          SHAPES(NNODSS,*), DVOLUS(*)
C
C**** FOR 'ONE-POINT' INTEGRATION RULES
C
      IF(NGAULS.EQ.1) RETURN
C
C**** "MASS" MATRIX
C
      RADIUS=1.0D0
      DO INODES=1,NNODSS
       DO JNODES=INODES,NNODSS
        EMASSS(JNODES,INODES)=0.0D0
        DO IGAULS=1,NGAULS
         IF(NTYPES.EQ.3) THEN
          IF(ISMO1S.EQ.2) THEN
           GPCOA=0.0D0               ! average radial spatial coordinate
           DO IGAUXS=1,NGAULS
            GPCOA=GPCOA+GPCODS(1,IGAUXS)
           ENDDO
           GPCOA=GPCOA/NGAULS
           IF(GPCODS(1,IGAULS).GT.(1.0D-10*GPCOA)) THEN
            RADIUS=1.0D0/GPCODS(1,IGAULS)
           ELSE
            RADIUS=2.0D0/GPCOA
           ENDIF
          ENDIF
         ENDIF
         EMASSS(JNODES,INODES)=EMASSS(JNODES,INODES)+
     .                         SHAPES(INODES,IGAULS)*
     .                         SHAPES(JNODES,IGAULS)*DVOLUS(IGAULS)*
     .                         RADIUS
        ENDDO
       ENDDO
      ENDDO
C
      DO INODES=2,NNODSS
       DO JNODES=1,INODES-1
        EMASSS(JNODES,INODES)=EMASSS(INODES,JNODES)
       ENDDO
      ENDDO
C
      IF(ISMO2S.EQ.1) CALL LUMPMMT(EMASSS,NNODSS,     0)
C
      CALL INVERT(EMASSS,NNODSS,NNODSS)
C
      RETURN
      END
