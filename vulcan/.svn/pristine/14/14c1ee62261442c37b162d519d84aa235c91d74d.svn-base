      SUBROUTINE EXPITET(DISIT,IFFIX,REFOR,
     .                   NTOTV,
     .                   CSTIT)
C***********************************************************************
C
C**** THIS ROUTINE SOLVES THE SYSTEM OF EQUATIONS USING  
C     AN EXPLICIT SOLVER
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DISIT(*), IFFIX(*), REFOR(*)
      DIMENSION CSTIT(NTOTV)
C
C**** OBTAIN ITERATIVE SOLUTION
C
      DO ITOTV=1,NTOTV
       IF(IFFIX(ITOTV).NE.0) REFOR(ITOTV)=0.0
       DISIT(ITOTV)=DISIT(ITOTV)+REFOR(ITOTV)/CSTIT(ITOTV)
      ENDDO
C
      RETURN
      END   
