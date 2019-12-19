      SUBROUTINE CPUTIMT7(RTIMET)
C***********************************************************************
C
C**** THIS ROUTINE FINDS OUT THE CPU TIME IN SEC. (HP)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** FOR SILICON GRAPHICS (UNIX)
C
      REAL*4 STIMET,tarray,etime
      dimension tarray(2)
      SAVE 
C
      DATA IPASST/0/
C
      IF(IPASST.EQ.0) THEN
       RTIMET=0.0
       IPASST=1
      ELSE
       stimet=0
c       STIMET=etime(tarray)
c       RTIMET=tarray(1)
      ENDIF
C
      RETURN
      END
