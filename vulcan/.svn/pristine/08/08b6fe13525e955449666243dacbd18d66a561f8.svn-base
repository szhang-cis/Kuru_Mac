      SUBROUTINE FASTOR(FATIG,
     .                  FREDU,ANUMF,STENI,STEN1,STEN2,STMIN,STMAX,REVEF)
C***********************************************************************
C
C**** THIS SUBROUTINE STORES FATIGUE (INTERNAL) VARIABLES
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
       FATIG(1)=FREDU
       FATIG(2)=ANUMF
       FATIG(3)=STENI
       FATIG(4)=STEN1
       FATIG(5)=STEN2
       FATIG(6)=STMIN
       FATIG(7)=STMAX
       FATIG(8)=REVEF
      ENDIF
C
      RETURN
      END
