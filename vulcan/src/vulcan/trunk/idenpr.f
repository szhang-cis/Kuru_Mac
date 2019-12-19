      SUBROUTINE IDENPR(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE CHOOSE THE PROPER SUBROUTINE TO ORDER THE PROPERTIES
C     FOR THE DIFFERENT MODELS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*)
C
C**** INCOMPRESSIBILITY CONDITION
C
      IF(KPROB.EQ.5) RETURN
C
C**** LOOKS FOR MODEL DEFINITION
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IPEP1=100
      IF(ISOTR.EQ.1) IPEP1=200
      IPEP2=INT(PROPS(2))   ! 10=elastic; 20=plastic; 30=viscoplastic;
C                           ! 40=hyperelastic
      IPEP3=INT(PROPS(3))   ! 1=no temp.-dependent; 2=temp.-depend.
      IPEP4=INT(PROPS(4))   ! 0=standard; 1,...,n=non-standard
      IPEPE=IPEP1+IPEP2+IPEP3
C
      IPEP5=INT(PROPS(5))   ! 0=uni-compound; 1=multi-compound
C
      IPEP6=INT(PROPS(6))   ! 0=no damage; 1=damage     not used now!!!
C
      IF(IPEP5.EQ.0) THEN
       IF(IPEP4.EQ.0) THEN
        IF(IPEPE.EQ.111) THEN
         CALL IDEIEN(PROPS)
        ENDIF
#ifndef restricted
        IF(IPEPE.EQ.112) THEN
         CALL RUNEND('IDENPR: 112 not implemented yet')
        ENDIF
        IF(IPEPE.EQ.121) THEN
         CALL IDEIPN(PROPS)
        ENDIF
        IF(IPEPE.EQ.122) THEN
         CALL IDEIPT(PROPS)
        ENDIF
        IF(IPEPE.EQ.131) THEN
         CALL IDEIVN(PROPS)
        ENDIF
        IF(IPEPE.EQ.132) THEN
         CALL IDEIVT(PROPS)
        ENDIF
        IF(IPEPE.EQ.141) THEN
         CALL RUNEND('IDENPR: 141 not implemented yet')
        ENDIF
        IF(IPEPE.EQ.142) THEN
         CALL IDEPOL(PROPS)
        ENDIF
        IF(IPEPE.EQ.211) THEN
         CALL IDEOEN(PROPS)
        ENDIF
        IF(IPEPE.EQ.212) THEN
         CALL RUNEND('IDENPR: 212 not implemented yet')
        ENDIF
        IF(IPEPE.EQ.221) THEN
         CALL RUNEND('IDENPR: 221 not implemented yet')
        ENDIF
        IF(IPEPE.EQ.222) THEN
         CALL IDEOPT(PROPS)
        ENDIF
        IF(IPEPE.EQ.231) THEN
         CALL IDEOVN(PROPS)
        ENDIF
        IF(IPEPE.EQ.232) THEN
         CALL RUNEND('IDENPR: 232 not implemented yet')
        ENDIF
        IF(IPEPE.EQ.242) THEN
         CALL IDEPOL(PROPS)
        ENDIF
#endif
       ENDIF
#ifndef restricted
       IF(IPEP4.EQ.1) THEN
        CALL IDEIVT1(PROPS)
       ENDIF
#endif
       IF(IPEP4.EQ.2) THEN
        CALL RUNEND('IDENPR: 2 not implemented yet')
       ENDIF
#ifndef restricted
       IF(IPEP4.EQ.3) THEN
        CALL IDEIVT3(PROPS)
       ENDIF
#endif
      ELSE
       call runend('idenpr: composite not implemented')
      ENDIF
C
      IFATI=INT(PROPS(40))  ! 0=no fatigue; 1=fatigue
      IF(IFATI.EQ.1) THEN
       IFATM=INT(PROPS(41))
       IF(IFATM.EQ.1) THEN
        CALL IDEFA1(PROPS)
       ENDIF
c      IF(IFATM.EQ.2) THEN
c       CALL IDEFA2(PROPS)
c      ENDIF
      ENDIF
C
      RETURN
      END
