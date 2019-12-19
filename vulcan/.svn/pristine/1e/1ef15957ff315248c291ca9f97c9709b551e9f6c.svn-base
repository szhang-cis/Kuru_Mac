      SUBROUTINE CPUTIM7(RTIME)
C***********************************************************************
C
C**** THIS ROUTINE FINDS OUT THE CPU TIME IN SEC. (HP)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** FOR SILICON GRAPHICS (UNIX)
C
      REAL*4 STIME,tarray,etime
      dimension tarray(2)
      SAVE
C
      DATA IPASS/0/
C
      IF(IPASS.EQ.0) THEN
       RTIME=0.0
       IPASS=1
      ELSE
       stime=0
c       STIME=etime(tarray)
c       RTIME=tarray(1)
      ENDIF
C
      RETURN
      END     
