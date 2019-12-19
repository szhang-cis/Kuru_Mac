      SUBROUTINE INIDIST(DISITT,IFFIXT,PRESCT,FICTOT,RLOADT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALIZES ITERATIVE DISPLACEMENTS ACCORDING TO 
C     PRESCRIBED VALUES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION DISITT(*),        IFFIXT(*),
     .          RLOADT(NTOTVT)
      DIMENSION PRESCT(NTOTVT,2), FICTOT(NFUNCT)
C
      DO ITOTVT=1,NTOTVT
       DISITT(ITOTVT)=0.0
      END DO
C
      IF(IITERT.GT.1)                 RETURN
      IF(KARCLT.NE.0.AND.ISTEPT.GT.2) RETURN ! IITER = 1
      IF(ISTEPT.GT.1.AND.LACCET.NE.0) RETURN ! IITER = 1
C
      DO ITOTVT=1,NTOTVT
       IF(IFFIXT(ITOTVT).NE.0) THEN
        RLOADT(ITOTVT)=PRESCT(ITOTVT,1)
        IPRESC=INT(PRESCT(ITOTVT,2))
        DISITT(ITOTVT)=RLOADT(ITOTVT)*FICTOT(IPRESC)
       ENDIF
      END DO
C
      RETURN
      END
