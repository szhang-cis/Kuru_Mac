      SUBROUTINE POSTSOS(DISITS,DSOLDS,REFORS,RFOLDS,VVECTS,WVECTS)
C***********************************************************************
C
C**** THIS ROUTINE MODIFIES THE ITERATIVE DISPLACEMENTS IF A 
C     SECANT-NEWTON OR QUASI-NEWTON METHOD IS USED
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION DISITS(*), DSOLDS(*), REFORS(*), RFOLDS(*)
      DIMENSION VVECTS(*), WVECTS(*)
C
C**** MODIFY ITERATIVE DISPLACEMENTS FOR SECANT-NEWTON METHOD
C
      CALL SECANTS(DISITT,DSOLDT,REFORT,RFOLDT)
C
C**** MODIFY ITERATIVE DISPLACEMENTS FOR QUASI-NEWTON METHOD ( BFGS )
C
c     IF((LACCET.EQ.3).AND.IITERT.GT.1.AND.KARCLT.EQ.0)
c    .  CALL BFGSUPT(DISITT,REFORT,VVECTT,WVECTT)
C
      IF(LACCES.LE.1) RETURN
C
C**** STORE DISIT & REFOR IN ALLOCATED SPACE
C
      DO ITOTVS=1,NTOTVS
        DSOLDS(ITOTVS)=DISITS(ITOTVS)
        RFOLDS(ITOTVS)=REFORS(ITOTVS)
      ENDDO
C
      RETURN
      END
