      SUBROUTINE CPUTIM8(RTIME)
C***********************************************************************
C
C**** THIS ROUTINE FINDS OUT THE CPU TIME IN SEC. (LINUX)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** FOR LINUX
C
      REAL*4 STIME,tarray,etime
      dimension tarray(2)
      SAVE
C
      DATA IPASS/0/
C
      IF(IPASS.EQ.0) THEN
       STIME=omp_get_wtime()
       RTIME=0.0
       IPASS=1
      ELSE
c      STIME=etime(tarray)
c      RTIME=tarray(1)
       RTIME=omp_get_wtime()-STIME
c      rtime=0.0
      ENDIF
C
      RETURN
      END
