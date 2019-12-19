      INTERFACE
       SUBROUTINE MODULTS(WORKST,IPRINT,LPRINT,
     .                   WS,    IS,    LPRINS,
     .                   W1)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION WORKST(*),WS(*)
        DIMENSION IPRINT(50),IS(50)
        REAL*8, ALLOCATABLE, INTENT(INOUT) :: W1(:)
       END SUBROUTINE
      END INTERFACE

