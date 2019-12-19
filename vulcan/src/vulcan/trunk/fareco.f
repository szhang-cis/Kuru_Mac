      SUBROUTINE FARECO(FATIG,
     .                  FREDU,ANUMF,STENI,STEN1,STEN2,STMIN,STMAX,REVEF)
C***********************************************************************
C
C**** THIS SUBROUTINE RECOVERS FATIGUE (INTERNAL) VARIABLES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      COMMON/FATIPROE/PFATI,PCYCL,PREVF,REVCO
C
      DIMENSION FATIG(*)
C
      KFATI=INT(PFATI)
      IF(KFATI.EQ.0) RETURN
C
      IF(IFATI.EQ.0) RETURN
C
C**** FATIGUE VARIABLES
C
      IF(IFATM.EQ.1) THEN
       FREDU=FATIG(1)
       ANUMF=FATIG(2)
       STENI=FATIG(3)
       STEN1=FATIG(4)
       STEN2=FATIG(5)
       STMIN=FATIG(6)
       STMAX=FATIG(7)
       REVEF=FATIG(8)
      ENDIF
C
      RETURN
      END
