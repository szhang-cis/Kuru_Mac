      SUBROUTINE FUNLO9(FACTI,HTLOD,TTIME)
C***********************************************************************
C
C**** THIS ROUTINE MODELS AN ERROR FUNCTION
C
C            .
C        .       .
C       .         .
C     .             .
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION HTLOD(*)
C
      FACTI=0.0D0
      IF(TTIME.LE.HTLOD(2)) RETURN       ! NOT STARTED YET
      IF(TTIME.GT.HTLOD(3)) RETURN       ! ALREADY FINISHED
C
      TMEDI=(HTLOD(2)+HTLOD(3))*0.5D0    ! CENTERED TIME
      RTIME=HTLOD(4)*(TTIME-TMEDI)       ! RADIAL COORDINATE
      RADIL=(HTLOD(3)-HTLOD(2))*0.5D0*HTLOD(4)           ! RADIUS
      FACTI=DEXP(-2.0D0*(RTIME/RADIL)*(RTIME/RADIL))*HTLOD(5)
C
      RETURN
      END
