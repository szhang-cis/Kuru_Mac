      SUBROUTINE PLSTOR(EBASE,CAPAP,EFFPL,BACKA,BACKS,POROS,
     .                  PREYS,STRAP,SPTOT,EFFP2,SRECR,ERECR,YIELD)
C***********************************************************************
C
C**** THIS SUBROUTINE STORES INTERNAL VARIABLES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION EBASE(*), BACKA(*), BACKS(*), STRAP(*), SPTOT(*)
C
C**** FLOW POTENTIAL CHOICE
C
      GO TO (31,32,33,34,35,36,37,38,39,40,41,42) (NCRIP-30)
C
C**** TRESCA
C
   31 CONTINUE
      CALL RUNEND('ERROR IN PLSTOR       ')
      RETURN
C
C**** VON MISES
C
   32 CONTINUE
      DO ISTR1=1,NSTR1
       STRAP(ISTR1)=SPTOT(ISTR1)
      ENDDO
C
      IN=1
      EBASE(1)=EFFPL                     ! effective plastic strain
C
      IF(ISOTT.GT.0) THEN
       EBASE(IN+1)=CAPAP                 ! plastic hardening function
       IN=IN+1
      ENDIF
C
      IF(IKINE.GE.1.AND.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        EBASE(IN+ISTR1)=BACKS(ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
      IF(IKINE.EQ.5) THEN                ! not implemented yet
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EBASE(IN+1)=EFFP2
       IN=IN+1
      ENDIF
C
      IF(IRECR.EQ.1) THEN
       EBASE(IN+1)=SRECR
       EBASE(IN+2)=ERECR
       IN=IN+2
      ENDIF
C
      IF(IVIFL.EQ.6) THEN                ! SMA
       EBASE(IN+1)=YIELD                 ! last converged equiv. stress
       IN=IN+1
      ENDIF
      RETURN
C
C**** MOHR-COULOMB & VERSION  WITH "TENSION CUT-OFF"
C
   33 CONTINUE
      CALL RUNEND('ERROR IN PLSTOR       ')
      RETURN
C
C**** DRUCKER-PRAGER
C
   34 CONTINUE
      CALL RUNEND('ERROR IN PLSTOR       ')
      RETURN
C
C**** J. LUBLINER'S THEORY
C
   35 CONTINUE
      CALL RUNEND('ERROR IN PLSTOR       ')
      RETURN
C
C**** ABOUAF'S MODEL
C
   36 CONTINUE
      CALL RUNEND('ERROR IN PLSTOR       ')
      RETURN
C
C**** WEBER BROWN'S MODEL
C
   37 CONTINUE
      CALL RUNEND('ERROR IN PLSTOR       ')
      RETURN
C
C**** SG CAST IRON
C
   38 CONTINUE
      DO ISTR1=1,NSTR1
       STRAP(ISTR1)=SPTOT(ISTR1)
      ENDDO
C
      IN=1
      EBASE(1)=EFFPL   ! effective plastic strain
      IF(ISOTT.GT.0) THEN
       IN=2
       EBASE(2)=CAPAP  ! plastic damage or plastic hard. function
      ENDIF
C
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        EBASE(IN+ISTR1)=BACKS(ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
      IF(IKINE.EQ.5) THEN                          ! not implemented yet
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EBASE(IN+1)=EFFP2
      ENDIF
      RETURN
C
C**** GREEN SAND
C
   39 CONTINUE
      DO ISTR1=1,NSTR1
       STRAP(ISTR1)=SPTOT(ISTR1)
      ENDDO
C
      IN=1
      EBASE(1)=EFFPL   ! effective plastic strain
      IF(ISOTT.GT.0) THEN
       IN=2
       EBASE(2)=CAPAP  ! plastic damage or plastic hard. function
      ENDIF
C
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        EBASE(IN+ISTR1)=BACKS(ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
      IF(IKINE.EQ.5) THEN                          ! not implemented yet
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EBASE(IN+1)=EFFP2
      ENDIF
      RETURN
C
C**** GREEN SAND
C
   40 CONTINUE
      DO ISTR1=1,NSTR1
       STRAP(ISTR1)=SPTOT(ISTR1)
      ENDDO
C
      IN=1
      EBASE(1)=EFFPL   ! effective plastic strain
      IF(ISOTT.GT.0) THEN
       IN=2
       EBASE(2)=CAPAP  ! plastic damage or plastic hard. function
      ENDIF
C
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        EBASE(IN+ISTR1)=BACKS(ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
      IF(IKINE.EQ.5) THEN                          ! not implemented yet
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EBASE(IN+1)=EFFP2
      ENDIF
      RETURN
C
C**** GURSON
C
   41 CONTINUE
      DO ISTR1=1,NSTR1
       STRAP(ISTR1)=SPTOT(ISTR1)
      ENDDO
C
      IN=1
      EBASE(1)=EFFPL   ! effective plastic strain
      IF(ISOTT.GT.0) THEN
       IN=2
       EBASE(2)=CAPAP  ! plastic damage or plastic hard. function
      ENDIF
C
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        EBASE(IN+ISTR1)=BACKS(ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
      IF(IKINE.EQ.5) THEN                          ! not implemented yet
      ENDIF
C
      IF(IPORO.GT.0) THEN
       EBASE(IN+1)=POROS
       EBASE(IN+2)=PREYS
       IN=IN+2
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EBASE(IN+1)=EFFP2
      ENDIF
      RETURN
C
C**** HILL 48
C
   42 CONTINUE
      DO ISTR1=1,NSTR1
       STRAP(ISTR1)=SPTOT(ISTR1)
      ENDDO
C
      IN=1
      EBASE(1)=EFFPL   ! effective plastic strain
      IF(ISOTT.GT.0) THEN
       IN=2
       EBASE(2)=CAPAP  ! plastic damage or plastic hard. function
      ENDIF
C
      IF(IKINE.GE.1.OR.IKINE.LE.4) THEN
       DO ISTR1=1,NSTR1
        EBASE(IN+ISTR1)=BACKS(ISTR1)
       ENDDO
       IN=IN+NSTR1
      ENDIF
      IF(IKINE.EQ.5) THEN                          ! not implemented yet
      ENDIF
C
      IF(ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       EBASE(IN+1)=EFFP2
      ENDIF
      RETURN
C
      END
