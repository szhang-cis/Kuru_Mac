      SUBROUTINE CPUTIMT8(RTIMET)
C***********************************************************************
C
C**** THIS ROUTINE FINDS OUT THE CPU TIME IN SEC. (LINUX)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** FOR LINUX
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
       STIMET=etime(tarray)
       RTIMET=tarray(1)
c      rtimet=0.0
      ENDIF
C
      RETURN
      END
