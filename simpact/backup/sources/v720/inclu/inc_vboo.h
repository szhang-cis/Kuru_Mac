SUBROUTINE inc_vboo(vboo,ninc,dflt)
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN),OPTIONAL:: ninc
  LOGICAL,INTENT(IN),OPTIONAL:: dflt
  LOGICAL,POINTER:: vboo(:)

END SUBROUTINE inc_vboo
