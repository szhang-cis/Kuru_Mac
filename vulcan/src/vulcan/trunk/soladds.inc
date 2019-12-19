      INTERFACE 
        SUBROUTINE SOLADDS(IFFIXS,LNODSS,LPNTNS,PRESCS,
     .                     R1,R2,R3,R4,R5,R6,
     .                     R7,WORK1S)
        IMPLICIT REAL*8(A-H,O-Z)
        INCLUDE 'prob_oms.f'
        INCLUDE 'inte_oms.f'
        DIMENSION IFFIXS(NTOTVS,*), LPNTNS(*), LNODSS(NNODES,*),
     .            PRESCS(NTOTVS,2)
        REAL*8, ALLOCATABLE ,INTENT (INOUT) :: WORK1S(:)
       END SUBROUTINE
      END INTERFACE
