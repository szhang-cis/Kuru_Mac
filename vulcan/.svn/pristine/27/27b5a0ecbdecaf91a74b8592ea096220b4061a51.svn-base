      SUBROUTINE BDBCO1(BMATX,DMATX,EMATX,KSYMM,NCOLU,NRAWS,NRAW1,NCOL1,
     .                  NRAW2)
C***********************************************************************
C                                T 
C**** THIS ROUTINE PERFORMS E = B D B OF TWO GENERAL MATRICES
C                   ( upper triangle for the symmetric case)
C                   ( general product for the unsymmetric case)
C
C                  BMATX = NRAWS * NCOLU
C                  DMATX = NRAW1 * NRAW1
C                  EMATX = NCOLU * NCOLU
C
C                  NRAWS=NSTRE (from kmatri.f) or NSTR1 (from kmatrb.f)
C                  NRAW1=NSTRS
C                  NCOLU=NEVAB
C                  NCOL1=local NEVAB
C                  NRAW2=NSTRE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION BMATX(NRAWS,*), DMATX(NRAW1,*), EMATX(NCOLU,*)
C
      DO ICOLU=1,NCOL1
       INDEX=ICOLU                ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) INDEX=1     ! UNSYMMETRIC CASE
       DO JCOLU=INDEX,NCOL1
        EMATX(ICOLU,JCOLU)=0.
        DO KRAWS=1,NRAW2
         DO LRAWS=1,NRAW2
          EMATX(ICOLU,JCOLU)=EMATX(ICOLU,JCOLU)+
     .          BMATX(KRAWS,ICOLU)*DMATX(KRAWS,LRAWS)*BMATX(LRAWS,JCOLU)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C
      RETURN
      END
