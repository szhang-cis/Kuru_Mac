      INTERFACE
       SUBROUTINE MODUSS(WORKST,WORK1T,IPRINS,LPRINS)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION WORKST(*)
        DIMENSION IPRINS(50)
        REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1T(:)
       END SUBROUTINE
      END INTERFACE
