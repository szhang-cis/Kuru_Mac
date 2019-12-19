SUBROUTINE readva ( )
!
!  Read Displacements, Velocities, Accelerations and Pressures
!

USE data_db
IMPLICIT NONE

REAL (kind=8), PARAMETER :: valor1=0.999D+99, valor2 =0.999D-99
INTEGER (kind=4) :: ipoin,idime,ndof   !node
REAL (kind=8) :: val(6)     !auxiliar values


!                    Read displacements
DO ipoin=1,npoin
  READ(17,end=100) displ(1:ndime,ipoin)
  DO idime=1,ndime
    IF(ABS(displ(idime,ipoin)) > valor1)    &
           displ(idime,ipoin) = SIGN(valor1,displ(idime,ipoin))
    IF(ABS(displ(idime,ipoin)) < valor2 .AND. displ(idime,ipoin) /= 0.0) &
           displ(idime,ipoin) = SIGN(valor2,displ(idime,ipoin))
  END DO
  coorf(:,ipoin) = coord(:,ipoin) + displ(:,ipoin)  !compute present coordinates
END DO

!                    Read nodal systems
IF (neulr > 0) THEN
  DO ipoin=1,npoin
    IF( weuler )THEN
      READ (17,end=100) euler(1:neulr,ipoin)
    ELSE
      READ (17,end=100)
    END IF
  END DO
  IF( addof ) THEN
    DO ipoin=1,npoin
      IF( waddof )THEN
        READ (17,end=100) psi(1:2,ipoin)
      ELSE
        READ (17,end=100)
      END IF
    END DO
  END IF
END IF

IF(idyna == 1) THEN
  !                    Read velocities
  ndof = ndime+nrotd
  val = 0d0
  DO ipoin=1,npoin
    READ(17,end=100) val(1:ndof)
    IF( wveloc ) veloc(1:ndime,ipoin) = val(1:ndime)
    IF( wangve ) anvel(1:nrotd,ipoin) = val(ndime+1:ndime+nrotd)
  END DO
  !                    Read accelerations
  DO ipoin=1,npoin
    READ(17,end=100) val(1:ndof)
    IF( waccel ) accel(1:ndime,ipoin) = val(1:ndime)
    IF( wangac ) anacc(1:nrotd,ipoin) = val(ndime+1:ndime+nrotd)
  END DO
END IF

IF( ndoft > 0 ) THEN
  !                    Read nodal temperatures
  DO ipoin=1,npoin
    IF( wtempe )THEN
      READ(17) tempe(:,ipoin)
    ELSE
      READ(17)
    END IF
  END DO
END IF


RETURN

100 fin = .TRUE.

END SUBROUTINE readva
