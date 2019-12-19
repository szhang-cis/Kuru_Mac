      SUBROUTINE INIDISS(DISITS,IFFIXS,PRESCS,FICTOS,RLOADS)
C***********************************************************************
C
C**** THIS ROUTINE INITIALIZES ITERATIVE DISPLACEMENTS ACCORDING TO 
C     PRESCRIBED VALUES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION DISITS(*),        IFFIXS(*),
     .          RLOADS(NTOTVS)
      DIMENSION PRESCS(NTOTVS,2), FICTOS(NFUNCS)
C
      DO ITOTVS=1,NTOTVS
       DISITS(ITOTVS)=0.0
      END DO
C
      IF(IITERS.GT.1)                 RETURN
      IF(KARCLS.NE.0.AND.ISTEPS.GT.2) RETURN ! IITER = 1
      IF(ISTEPS.GT.1.AND.LACCES.NE.0) RETURN ! IITER = 1
C
      DO ITOTVS=1,NTOTVS
       IF(IFFIXS(ITOTVS).NE.0) THEN
        RLOADS(ITOTVS)=PRESCS(ITOTVS,1)
        IPRESC=INT(PRESCS(ITOTVS,2))
        DISITS(ITOTVS)=RLOADS(ITOTVS)*FICTOS(IPRESC)
       ENDIF
      END DO
C
      RETURN
      END
