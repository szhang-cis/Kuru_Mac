      SUBROUTINE CHECKSM(values,ncoef,coefsum)
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION values(ncoef)
C
      coefsum=0.0
      do k=1,ncoef
         do l=1,4
            coefsum = coefsum + values(k)**l
         end do
      end do
C
      RETURN
      END
