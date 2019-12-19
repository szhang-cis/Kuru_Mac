      SUBROUTINE INPLAS(PROPS,EBASE)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES SOME INTERNAL VARIABLES
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      DIMENSION PROPS(*), EBASE(*)
C
      NCRIT=INT(PROPS(36))
C
C**** YIELD CRITERION CHOICE
C
      GO TO (31,32,33,34,35,36,37,38,39,40,41,42) (NCRIT-30)
C
C**** TRESCA
C
   31 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** VON MISES
C
   32 CONTINUE
      IN=1
      IF(ISOTT.GT.0) THEN
       IFVER=INT(PROPS(35))
       IF(IFVER.EQ.4) THEN      ! rough estimation for temp.-dep. models
        TEMPG=0.0D+00
        DTEMG=0.0D+00
        CALL CCEROT(TEMPG,PROPS,DTEMG,
     .              CCERO,DCCER,CCEROF,CCEROA)
        EBASE(IN+1)=CCERO
       ENDIF
       IN=IN+1
      ENDIF
C
      IF(IKINE.GE.1.AND.IKINE.LE.4) THEN
       IN=IN+NSTR1
      ENDIF
      IF(IKINE.EQ.5) THEN                ! not implemented yet
       IN=IN+2*NSTR1
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       IN=IN+1
      ENDIF
C
      IF(IRECR.EQ.1) THEN
       EBASE(IN+1)=RECRY(7)
       EBASE(IN+2)=RECRY(8)
       IN=IN+2
      ENDIF
      RETURN
C
C**** MOHR-COULOMB & VERSION  WITH "TENSION CUT-OFF"
C
   33 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** DRUCKER-PRAGER
C
   34 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** J. LUBLINER'S THEORY
C
   35 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** ABOUAF'S MODEL
C
   36 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** WEBER BROWN'S MODEL
C
   37 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** SG CAST IRON
C
   38 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** GREEN SAND
C
   39 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** GREEN SAND
C
   40 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** GURSON
C
   41 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
C**** HILL 48
C
   42 CONTINUE
      CALL RUNEND('ERROR IN INPLAS       ')
      RETURN
C
      END
