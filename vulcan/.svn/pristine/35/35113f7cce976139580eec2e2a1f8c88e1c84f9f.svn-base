      SUBROUTINE INIDIS(DISIT,IFFIX,PRESC,FICTO,RLOAD)
C***********************************************************************
C
C**** THIS ROUTINE INITIALIZES ITERATIVE DISPLACEMENTS ACCORDING TO 
C     PRESCRIBED VALUES
C    
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DISIT(*),       IFFIX(*),
     .          RLOAD(NTOTV)
      DIMENSION PRESC(NTOTV,2), FICTO(NFUNC)
C
      DO ITOTV=1,NTOTV
       DISIT(ITOTV)=0.0D0
      END DO
C
      IF(IITER.GT.1)                RETURN
      IF(KARCL.NE.0.AND.ISTEP.GT.2) RETURN ! IITER = 1
      IF(ISTEP.GT.1.AND.LACCE.NE.0) RETURN ! IITER = 1
C
      DO ITOTV=1,NTOTV
       IF(IFFIX(ITOTV).NE.0) THEN
        RLOAD(ITOTV)=PRESC(ITOTV,1)
        IPRESC=INT(PRESC(ITOTV,2))
        DISIT(ITOTV)=RLOAD(ITOTV)*FICTO(IPRESC)
       ENDIF
      END DO
C
      IF(KPORE.EQ.2) THEN
       DO ITOTV=NDOFC,NTOTV,NDOFC
        DISIT(ITOTV)=DISIT(ITOTV)/WLUMP
       END DO
      ENDIF
C
      RETURN
      END
