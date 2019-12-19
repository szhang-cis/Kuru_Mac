      INTERFACE

      SUBROUTINE SOLADDT(IFFIXT,LNODST,LPNTNT,PRESCT,
     .                   REAL1,REAL2,REAL3,REAL4,REAL5,REAL6,
     .                   REAL7,WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE ALLOCATES MEMORY FOR SOLVER
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE  :: WORK1T(:)

C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
C
      DIMENSION IFFIXT(NTOTVT,*), LPNTNT(*), LNODST(NNODET,*),
     .          PRESCT(NTOTVT,2)
      DIMENSION REAL1(*), REAL2(*), REAL3(*),
     .          REAL4(*), REAL5(*), REAL6(*),
     .          REAL7(*)
      END SUBROUTINE
      END INTERFACE
