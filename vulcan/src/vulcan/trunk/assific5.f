      SUBROUTINE ASSIFIC5
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS COUPLING FILES (PC)
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

      call runendt('error: input reading in pc not implemented')

C
      RETURN
      END
