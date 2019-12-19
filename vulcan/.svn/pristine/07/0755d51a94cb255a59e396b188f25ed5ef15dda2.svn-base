      SUBROUTINE STSNML(ALENG,OPENG,STRAN,SGTOT,YOUNG,conju,sssyy,templ,
     .                  teno3,teno4,conj1,SSSZZ,STRAZ,conj0,dstra)
C********************************************************************
C
C****THIS ROUTINE CALCULATES THE STRESS ACROSS THE JOINT
C
C   Input parameters:
C
C      ALENG      Length across the joint
C      OPENG      Initial opening gap
C      STRAN      Current total strain normal to joint
C      YOUNG      Young's modulus
C
C   Output parameters:
C
C      SGTOT      Current stress normal to joint
C
C********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C***CKECK IF THE JOINT IS OPEN
C
      sssyy=0.0
      SSSZZ=0.0
      taver=(teno3+teno4)/2.0
c
      if(taver.lt.templ) then
c
      THRES=-OPENG/ALENG
cccc      IF(STRAN.GE.THRES) THEN   
c
c si stran=thres, joint closed
c
ctm      IF(STRAN.GT.conj0) THEN   
      IF(STRAN.GT.0.0) THEN
       SGTOT=0.0
       conju=0.0
       RETURN
      ENDIF
C
C***IF THE JOINT IS CLOSED, EVALUATE NORMAL STRESS
C
      SOLID=ALENG-OPENG
c
ctm      CALL CONJUT(CONOD,CONOT,STRAN,CONJ1,CONJ0)
ctm      conju=conod+conot*stran
c
cc         if(taver.lt.conj0) then
cc          funcio=conj1
cc         else
cc          funcio=conj1-(taver-conj0)/(templ-conj0)
cc         endif
cc
cc         conju=1.0*funcio
         conju=1.0
c
ctm      SGTOT=conod*YOUNG*STRAN
cc         SGTOT=conju*YOUNG*STRAN
         SGTOT=YOUNG*STRAN
c
c
      else
c
C
C***JOINT IS CLOSED in liquid state
C
      conju=1.0
      SOLID=ALENG-OPENG
c
      SGTOT=YOUNG*STRAN
c
      endif
C
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cC
cC***CKECK IF THE JOINT IS OPEN
cC
c      sssyy=0.0
c      SSSZZ=0.0
c      taver=(teno3+teno4)/2.0
cc
c      if(taver.lt.templ) then
cctm
cctm      if(teno3.lt.templ.or.teno4.lt.templ) then
cctm
cc
c      THRES=-OPENG/ALENG
ccccc      IF(STRAN.GE.THRES) THEN   
cc
cc si stran=thres, joint closed
cc
cctm      IF(STRAN.GT.conj0) THEN   
c      IF(STRAN.GT.0.0) THEN
c       SGTOT=0.0
c       conju=0.0
c       RETURN
c      ENDIF
cC
cC***IF THE JOINT IS CLOSED, EVALUATE NORMAL STRESS
cC
c      SOLID=ALENG-OPENG
cc
cctm      CALL CONJUT(CONOD,CONOT,STRAN,CONJ1,CONJ0)
cctm      conju=conod+conot*stran
c         conju=1.0
cc
cctm      SGTOT=conod*YOUNG*STRAN
c      SGTOT=YOUNG*STRAN
cc
cctm
c      sgabb=dabs(sgtot)
c      if(sgabb.gt.conj0) sgtot=-conj0
cctm
cc
c      else
cc
cC
cC***JOINT IS CLOSED in liquid state
cC
ccccccccccccccccccccccccccccccc      conju=1.0
c      conju=1.0
c      SOLID=ALENG-OPENG
cc
cc
c      IF(STRAN.GT.0.0) THEN
c       SGTOT=0.0
c       conju=0.000001
c       RETURN
c      ENDIF
cc
cc
c      SGTOT=YOUNG*STRAN
cc
cc
cctm
c      sgabb=dabs(sgtot)
c      if(sgabb.gt.conj0) sgtot=-conj0
cctm
cc
cc
c      endif
cC
c      RETURN
c      END
