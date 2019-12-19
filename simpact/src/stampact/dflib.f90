! THIS FILE PROVIDES AN INTERFACE TO GFORTRAN
! PORTABILITY LIBRARY ROUTINES.
!
! THESE ROUTINES PROVIDE THE FUNCTIONALITY OF MANY
! COMMON LIBRARY EXTENSIONS TO THE FORTRAN LANGUAGE.

MODULE DFLIB
  ! -----------------------------------------------------------------
  ! Keyboard and Speaker Routines
  ! -----------------------------------------------------------------
CONTAINS
  SUBROUTINE BEEP()
    IMPLICIT NONE
    INTEGER :: BEL=7
    WRITE(*,*) ACHAR(BEL)
  END SUBROUTINE BEEP
END MODULE DFLIB
