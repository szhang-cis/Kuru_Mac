      SUBROUTINE EQHEATC(BMSIG,SHAPE,DVOLU,CONVE,NNODE,NEVAB)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE EQUIVALENT NODAL "CONVECTIVE FORCES"
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION BMSIG(*), SHAPE(*)
C
      DO IEVAB=1,NNODE
       BMSIG(IEVAB)=BMSIG(IEVAB)+SHAPE(IEVAB)*CONVE*DVOLU
      ENDDO
C
      RETURN
      END
