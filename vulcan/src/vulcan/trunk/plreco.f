      SUBROUTINE PLRECO(EBASE,CAPAP,EFFPL,BACKA,BACKS,POROS,
     .                  PREYA,PROPS,CCERO,EFFP2,SRECR,ERECR,YIELO)
C***********************************************************************
C
C**** THIS SUBROUTINE RECOVERS INTERNAL VARIABLES & THEIR CONJUGATES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION EBASE(*), BACKA(*), BACKS(*), PROPS(*)
C
C**** FLOW POTENTIAL CHOICE
C
      GO TO (31,32,33,34,35,36,37,38,39,40,41,42) (NCRIP-30)
C
C**** TRESCA
C
   31 CONTINUE
      CALL RUNEND('ERROR IN PLRECO       ')
      RETURN
C
C**** VON MISES
C
   32 CONTINUE
      EFFPL=EBASE(1)                     ! effective plastic strain
      IN=1
C
      CAPAP=0.0D0
      IF(ISOTT.GT.0) THEN
       CAPAP=EBASE(IN+1)                 ! plastic hardening function
       IN=IN+1
      ENDIF
C
      DO ISTR1=1,NSTR1
       BACKS(ISTR1)=0.0D0
      ENDDO
      IF(IKINE.GE.1.AND.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        BACKS(ISTR1)=EBASE(IN+ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
      IF(IKINE.EQ.5) THEN                ! not implemented yet
       DO ISTR1=1,NSTR1                  ! to be revised; see pointe.f
        BACKA(ISTR1)=EBASE(IN+ISTR1)
        BACKS(ISTR1)=EBASE(IN+NSTR1+ISTR1)
       ENDDO
       IN=IN+2*NSTR1
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EFFP2=EBASE(IN+1)
       IN=IN+1
      ENDIF
C
      SRECR=0.0D0
      ERECR=0.0D0
      IF(IRECR.EQ.1) THEN
       SRECR=EBASE(IN+1)
       ERECR=EBASE(IN+2)
       IN=IN+2
      ENDIF
C
      IF(IVIFL.EQ.6) THEN                ! SMA
       YIELO=EBASE(IN+1)                 ! last converged equiv. stress
       IN=IN+1
      ENDIF
      RETURN
C
C**** MOHR-COULOMB & VERSION  WITH "TENSION CUT-OFF"
C
   33 CONTINUE
      CALL RUNEND('ERROR IN PLRECO       ')
      RETURN
C
C**** DRUCKER-PRAGER
C
   34 CONTINUE
      CALL RUNEND('ERROR IN PLRECO       ')
      RETURN
C
C**** J. LUBLINER'S THEORY
C
   35 CONTINUE
      CALL RUNEND('ERROR IN PLRECO       ')
      RETURN
C
C**** ABOUAF'S MODEL
C
   36 CONTINUE
      CALL RUNEND('ERROR IN PLRECO       ')
      RETURN
C
C**** WEBER BROWN'S MODEL
C
   37 CONTINUE
      CALL RUNEND('ERROR IN PLRECO       ')
      RETURN
C
C**** SG CAST IRON
C
   38 CONTINUE
      EFFPL=EBASE(1)   ! effective plastic strain
      IN=1
      CAPAP=0.0D0
      IF(ISOTT.GT.0) THEN
       IN=2
       CAPAP=EBASE(2)  ! plastic damage or plastic hard. function
      ENDIF
C
      DO ISTR1=1,NSTR1
       BACKS(ISTR1)=0.0D0
      ENDDO
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        BACKS(ISTR1)=EBASE(IN+ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EFFP2=EBASE(IN+1)
      ENDIF
      RETURN
C
C**** GREEN SAND
C
   39 CONTINUE
      EFFPL=EBASE(1)   ! effective plastic strain
      IN=1
      CAPAP=0.0D0
      IF(ISOTT.GT.0) THEN
       IN=2
       CAPAP=EBASE(2)  ! plastic damage or plastic hard. function
      ENDIF
C
      DO ISTR1=1,NSTR1
       BACKS(ISTR1)=0.0D0
      ENDDO
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        BACKS(ISTR1)=EBASE(IN+ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EFFP2=EBASE(IN+1)
      ENDIF
      RETURN
C
C**** GREEN SAND
C
   40 CONTINUE
      EFFPL=EBASE(1)   ! effective plastic strain
      IN=1
      CAPAP=0.0D0
      IF(ISOTT.GT.0) THEN
       IN=2
       CAPAP=EBASE(2)  ! plastic damage or plastic hard. function
      ENDIF
C
      DO ISTR1=1,NSTR1
       BACKS(ISTR1)=0.0D0
      ENDDO
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        BACKS(ISTR1)=EBASE(IN+ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EFFP2=EBASE(IN+1)
      ENDIF
      RETURN
C
C**** GURSON
C
   41 CONTINUE
      EFFPL=EBASE(1)   ! effective plastic strain
      IN=1
      CAPAP=0.0D0
      IF(ISOTT.GT.0) THEN
       IN=2
       CAPAP=EBASE(2)  ! plastic damage or plastic hard. function
      ENDIF
C
      DO ISTR1=1,NSTR1
       BACKS(ISTR1)=0.0D0
      ENDDO
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        BACKS(ISTR1)=EBASE(IN+ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
C
      POROS=0.0D0
      IF(IPORO.GT.0) THEN
       POROS=EBASE(IN+1)
       PREYA=EBASE(IN+2)
       IN=IN+2
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EFFP2=EBASE(IN+1)
      ENDIF
      RETURN
C
C**** HILL 48
C
   42 CONTINUE
      EFFPL=EBASE(1)   ! effective plastic strain
      IN=1
      CAPAP=0.0D0
      IF(ISOTT.GT.0) THEN
       IN=2
       CAPAP=EBASE(2)  ! plastic damage or plastic hard. function
      ENDIF
C
      DO ISTR1=1,NSTR1
       BACKS(ISTR1)=0.0D0
      ENDDO
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        BACKS(ISTR1)=EBASE(IN+ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EFFP2=EBASE(IN+1)
      ENDIF
      RETURN
C
      END
