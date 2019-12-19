      SUBROUTINE MODULSS(WORKS, IP, MPRIN,
     .                   WORKST,IPT,MPRINT,
     .                   WS,    IS, MPRINS,
     .                   W1)
C====================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: W1(:)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES (PERMANENT)
C
      DIMENSION WORKS(MPRIN)
      DIMENSION IP(50)                           ! see note in vul-*m*.f
C
C**** THERMAL VARIABLES (PERMANENT)
C
      DIMENSION WORKST(MPRINT), IPT(50)
C
C**** MICROSTRUCTURAL VARIABLES (PERMANENT)
C
      DIMENSION WS(MPRINS),     IS(50)
      END SUBROUTINE
