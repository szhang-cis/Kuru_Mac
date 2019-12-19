      SUBROUTINE GATECA(TEMPG,DTEMG,TEINI,TENOD,DTENO,TENOI,SHAPE,KPART)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TEMPERATURES, TEMPERATURES RATES & 
C     AN INITIAL TEMPERATURES AT GAUSS POINTS
C
C     TENOD: nodal temperature
C     DTENO: nodal temperature rate
C     TENOI: nodal initial temperature
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
      DIMENSION SHAPE(NNODE), TENOD(*), DTENO(*), TENOI(*)
C
      TEMPG=0.0
      DTEMG=0.0
      TEINI=0.0
C
      IF(ITERME.LT.0) RETURN
C
      IF(INTERC.EQ.1) THEN
       CALL SMOMID(TENOD,NDIME,NNODL,NQUTR)
       CALL SMOMID(DTENO,NDIME,NNODL,NQUTR)
       CALL SMOMID(TENOI,NDIME,NNODL,NQUTR)
      ENDIF
C
      IF(KPART.EQ.0) THEN
       DO INODE=1,NNODL
        TEMPG=TEMPG+TENOD(INODE)/NNODL
        DTEMG=DTEMG+DTENO(INODE)/NNODL
        TEINI=TEINI+TENOI(INODE)/NNODL
       ENDDO
      ELSE
       DO INODE=1,NNODL
        TEMPG=TEMPG+SHAPE(INODE)*TENOD(INODE)
        DTEMG=DTEMG+SHAPE(INODE)*DTENO(INODE)
        TEINI=TEINI+SHAPE(INODE)*TENOI(INODE)
       ENDDO
      ENDIF
C
      RETURN
      END
