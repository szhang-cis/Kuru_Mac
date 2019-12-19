      SUBROUTINE ASSIFIC3
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS COUPLING FILES (VAX)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** COUPLING ASSIGN: INTERNAL FILE
C
C
C**** COUPLING ASSIGN: EXTERNAL FILE
C

      call runendt('error: input reading in vax not implemented')

C
      RETURN
      END
