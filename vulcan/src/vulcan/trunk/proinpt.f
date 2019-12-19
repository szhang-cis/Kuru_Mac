      SUBROUTINE PROINPT
C***********************************************************************
C
C**** THIS ROUTINE REDEFINES THE DEFAULT PROBLEM PARAMETERS 
C
C     PARAMETERS         CONTROLLING WORDS
C
C     KDYNA              'DYNAM','TRANS'
C     KPORE              'COUPL'
C     KPROB              'SEEPA', 'STRUC', 'WAVES', 'THERM'
C     LARGE              'LARGE'
C     others             'DIMEN'
C     COMMUNICATION      'COMMU'
C     ICONVT             'ADVEC'
C     IGALET             'WEAK_'
C     NPOROT             'POROS'
C     NACTIT             'ACTIV'
C     NOCOIT             'CONTA'
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
      INCLUDE 'nuef_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      PARAMETER (MCOMMT=16)
      CHARACTER*5 COMMDT(MCOMMT)
C
      DATA COMMDT/'SEEPA','STRUC','WAVES','THERM','TRANS','DYNAM',
     .            'COUPL','LARGE','DIMEN','COMMU','ADVEC','WEAK_',
     .            'POROS','ACTIV','CONTA','END_P'/
C
C**** READ PROBLEM PARAMETERS
C
      NPRINT=0
      ITAPET=LUDATT
      CALL LISTENT('PROINPT',NPRINT,ITAPET)
      IF(WORDST(1).NE.'PROBL') GO TO 2000
 1000 CALL LISTENT('PROINPT',NPRINT,ITAPET)
C
C**** IDENTIFY COMMAND
C
      DO ICOMMT=1,MCOMMT
       IF(WORDST(1).EQ.COMMDT(ICOMMT)) GO TO 100
      END DO
      GO TO 2000 
C
C**** EXECUTE APROPRIATE COMMAND
C
  100 CONTINUE
      GO TO (1,1,1,1,2,2,3,4,5,6,7,8,9,10,11,12), ICOMMT
C
    1 IF(WORDST(1).EQ.'SEEPA')                          KPROBT=0
      IF(WORDST(1).EQ.'STRUC')                          KPROBT=1
      IF(WORDST(1).EQ.'STRUC'.AND.WORDST(2).EQ.'SHELL') KPROBT=2
      IF(WORDST(1).EQ.'WAVES')                          KPROBT=3
      IF(WORDST(1).EQ.'THERM')                          KPROBT=4
      GO TO 1000
C
    2 KDYNAT=1
      GO TO 1000
C
    3 KPORET=2
      IF(WORDST(2).EQ.'DISCR') KSMUST=-1
      IF(WORDST(2).EQ.'LOCAL') KSMUST=-2
      GO TO 1000
C
C**** LARGE DISPLACEMENTS/STRAINS ANALYSIS
C
C     Notes:
C
C     LARGET=1: Total Lagrangian Formulation
C     LARGET=2: Updated Lagrangian Formulation
C     LARGET=3: Eulerian Formulation
C
    4 LARGET=1                              ! default
C
      IF(ITERME.LT.0) THEN                  ! uncoupled thermal analysis
       LARGET=0
       CALL RUNMENT('WARNING: LARGET NE 0 IS NOT CONSIDERED')
      ELSE
       IF(WORDST(2).NE.'     ') THEN
        IF(WORDST(2).NE.'TOTAL'.AND.WORDST(2).NE.'UPDAT'.AND.
     .     WORDST(2).NE.'EULER')
     .   CALL RUNENDT('PROINPT: WRONG LARGE STRAINS INPUT DATA')
c       IF(WORDST(2).EQ.'TOTAL') LARGET=1
        IF(WORDST(2).EQ.'TOTAL')
     .   CALL RUNENDT('ERROR: LARGET=1 NOT IMPLEMENTED')
c       IF(WORDST(2).EQ.'UPDAT') LARGET=2
        IF(WORDST(3).EQ.'TOTAL')
     .   CALL RUNENDT('ERROR: LARGET=2 NOT IMPLEMENTED')
        IF(WORDST(2).EQ.'EULER') LARGET=3
       ENDIF                     ! wordst(2).ne....
      ENDIF                      ! iterme.lt.0
      GO TO 1000
C
    5 DWORDT(1)='NPOIN'
      DWORDT(2)='NELEM'
      DWORDT(3)='NDIME'
      DWORDT(4)='NNODE'
      DWORDT(5)='NGAUS'
      DWORDT(6)='NSETS'
      DWORDT(7)='NMATS'
      DWORDT(8)='NFUNC'
      DPARAT(1)=0.
      DPARAT(2)=0.
      DPARAT(3)=2.
      DPARAT(4)=4.
      DPARAT(5)=4.
      DPARAT(6)=1.
      DPARAT(7)=1.
      DPARAT(8)=1.
C
      CALL SORTERT(1)
C
      NPOINT=INT(DPARAT(1))
      NELEMT=INT(DPARAT(2))
      NDIMET=INT(DPARAT(3))
      NNODET=INT(DPARAT(4))
      NGAUST=INT(DPARAT(5))
      NGRUPT=INT(DPARAT(6))
      NMATST=INT(DPARAT(7))
      NFUNCT=INT(DPARAT(8))
C
      IF(NGRUPT.LT.NMATST)
     . CALL RUNENDT('NGRUPT CAN NOT BE LESS THAN NMATST')
      GO TO 1000
C
    6 CONTINUE
      GO TO 1000
C
C**** ADVECTIVE EFFECTS
C
C     Notes:
C
C     IEGFPC=0 => grad(f_pc)=df_pc/dT*grad(T)
C     IEGFPC=1 => grad(f_pc) directly computed from f_pc
C     (see forcint.f & elm005t.f)
C
C     ICUBIC=0 => only for linear f_pc, df_pc/dT is regularized with a 
C                 cubical polinomial to stabilize the solution
C     ICUBIC=1 => df_pc/dT is computed as it should be
C     (see capcoft.f)
C
    7 ICONVT=1
      IEGFPC=0
      ICUBIC=0
      IF(WORDST(2).EQ.'FPC_G') THEN
       IF(PARAMT(1).NE.0.0D0) IEGFPC=1
       IF(WORDST(3).EQ.'FPC_T') THEN
        IF(PARAMT(2).NE.0.0D0) ICUBIC=1
       ENDIF
      ENDIF
      IF(WORDST(2).EQ.'FPC_T') THEN
       IF(PARAMT(2).NE.0.0D0) ICUBIC=1
       IF(WORDST(2).EQ.'FPC_G') THEN
        IF(PARAMT(1).NE.0.0D0) IEGFPC=1
       ENDIF
      ENDIF
      IF(IEGFPC.EQ.1) THEN              ! control
       IF(NMEMO3.EQ.1.AND.KSGAUT.EQ.0)
     .  CALL RUNENDT('ERROR: KSGAUT=0 FOR IEGFPC=1 & NMEMO3=1')
      ENDIF
      GO TO 1000
C
    8 IF(WORDST(2).EQ.'GALER') IGALET=0 ! default value (see setdatt.f)
C
      IF(WORDST(2).EQ.'PETRO') THEN     ! upwinding techniques
       IGALET=1
       IX=1
       IF(WORDST(2+IX).EQ.'UPWIN') THEN
        IF(PARAMT(IX).NE.0.0D0) IUPWI=INT(PARAMT(IX)) ! upwind function
        IX=IX+1
       ENDIF
       IF(WORDST(2+IX).EQ.'PERTU') THEN
        IF(PARAMT(IX).NE.0.0D0) IPERT=INT(PARAMT(IX)) ! pertur. function
        IX=IX+1
       ENDIF
       IF(WORDST(2+IX).EQ.'PECLE') THEN
        IF(PARAMT(IX).NE.0.0D0) ISUPW=INT(PARAMT(IX)) ! Peclet def.
        IX=IX+1
       ENDIF
       IF(WORDST(2+IX).EQ.'FACTO') THEN
        IF(PARAMT(IX).NE.0.0D0) EXTUP=PARAMT(IX)      ! factor
        IX=IX+1
       ENDIF
C
       IF(WORDST(2+IX).EQ.'CROSS') THEN
        IF(PARAMT(IX).NE.0.0D0) IUPWIC=INT(PARAMT(IX)) ! crosswind func.
        IX=IX+1
        IF(IUPWIC.EQ.1) THEN
         IF(WORDST(2+IX).EQ.'PARAM') THEN
          IF(PARAMT(IX).NE.0.0D0) RAMON=PARAMT(IX)     ! Ramon's factor
          IX=IX+1
         ELSE
          RAMON=0.7D0                                  ! default
         ENDIF
        ENDIF
       ENDIF
       IF(WORDST(2+IX).EQ.'PERTU') THEN
        IF(PARAMT(IX).NE.0.0D0) IPERTC=INT(PARAMT(IX)) ! perturb. func.
        IX=IX+1
       ENDIF
       IF(WORDST(2+IX).EQ.'PECLE') THEN
        IF(PARAMT(IX).NE.0.0D0) ISUPWC=INT(PARAMT(IX)) ! Peclet def.
        IX=IX+1
       ENDIF
       IF(WORDST(2+IX).EQ.'FACTO') THEN
        IF(PARAMT(IX).NE.0.0D0) EXTUPC=PARAMT(IX)      ! factor
        IX=IX+1
       ENDIF
C
       IF(WORDST(2+IX).EQ.'TEMPO') THEN
        IF(PARAMT(IX).NE.0.0D0) IUPWIG=INT(PARAMT(IX)) ! temp. upwind f.
        IX=IX+1
       ENDIF
       IF(WORDST(2+IX).EQ.'PERTU') THEN
        IF(PARAMT(IX).NE.0.0D0) IPERTG=INT(PARAMT(IX)) ! pertur. f.
        IX=IX+1
       ENDIF
       IF(WORDST(2+IX).EQ.'FOURI') THEN
        IF(PARAMT(IX).NE.0.0D0) ISUPWG=INT(PARAMT(IX)) ! Fourier def.
        IX=IX+1
       ENDIF
       IF(WORDST(2+IX).EQ.'FACTO') THEN
        IF(PARAMT(IX).NE.0.0D0) EXTUPG=PARAMT(IX)      ! factor
        IX=IX+1
       ENDIF
C
       IF(WORDST(2+IX).EQ.'GALFA') THEN
        IF(PARAMT(IX).NE.0.0) THEN
         GALFA=PARAMT(IX)                              ! thermal tau
         IGALFA=1
        ENDIF
        IF(ITERMEF.LE.0) THEN
         CALL RUNMENT('WARNING: GALFA=0.0 IS CONSIDERED')
         GALFA=0.0D0                                   ! default
         IGALFA=0
        ENDIF
       ENDIF
C
C**** SOME CONTROLS
C
       IF(IUPWI.GT.0) THEN
        IF(IUPWI.GT.4)
     .   CALL RUNENDT('ERROR IN UPWIND FUNCTION CHOICE')
        IF(IUPWI.EQ.3)
     .   CALL RUNENDT('ERROR IN UPWIND FUNCTION CHOICE (not imp. yet)')
        IF(IPERT.LT.1.OR.IPERT.GT.3)
     .   CALL RUNENDT('ERROR IN PERTURBATION FUNCTION CHOICE')
        IF(ISUPW.LT.1.OR.ISUPW.GT.6)
     .   CALL RUNENDT('ERROR IN PECLET DEFINITION CHOICE')
        IF(EXTUP.LT.0.0D0)
     .   CALL RUNENDT('ERROR IN FACTOR CHOICE')
       ENDIF
C
       IF(IUPWIC.GT.0) THEN
        IF(IUPWIC.GT.1)
     .   CALL RUNENDT('ERROR IN CROSSWIND FUNCTION CHOICE')
        IF(IPERTC.LT.1.OR.IPERTC.GT.1)
     .   CALL RUNENDT('ERROR IN (CROSS) PERTURBATION FUNCTION CHOICE')
        IF(ISUPWC.LT.1.OR.ISUPWC.GT.4)
     .   CALL RUNENDT('ERROR IN (CROSS) PECLET DEFINITION CHOICE')
        IF(EXTUPC.LT.0.0D0)
     .   CALL RUNENDT('ERROR IN (CROSS) FACTOR CHOICE')
C
        IF(IUPWI.GT.0) THEN
         IF(IPERT.NE.3) THEN
          IPERT=3
          CALL RUNMENT('WARNING: IPERT SET TO 3 IN THIS CASE - 1')
         ENDIF
        ENDIF
       ENDIF
C
       IF(IUPWIG.GT.0) THEN
        IF(IUPWIG.GT.1)
     .   CALL RUNENDT('ERROR IN TEMPORAL UPWIND FUNCTION CHOICE')
        IF(IPERTG.LT.1.OR.IPERTG.GT.1)
     .   CALL RUNENDT('ERROR IN (TEMP) PERTURBATION FUNCTION CHOICE')
        IF(ISUPWG.LT.1.OR.ISUPWG.GT.1)
     .   CALL RUNENDT('ERROR IN (TEMP) FOURIER DEFINITION CHOICE')
        IF(EXTUPG.LT.0.0D0)
     .   CALL RUNENDT('ERROR IN (TEMP) FACTOR CHOICE')
C
        IF(IUPWI.GT.0) THEN
         IF(IPERT.NE.3) THEN
          IPERT=3
          CALL RUNMENT('WARNING: IPERT SET TO 3 IN THIS CASE - 2')
         ENDIF
        ENDIF
C
        IF(IUPWIC.GT.0) THEN
         IF(IPERTC.NE.1) THEN
          IPERTC=1
          CALL RUNMENT('WARNING: IPERTC SET TO 1 IN THIS CASE')
         ENDIF
        ENDIF
C
        IF(KDYNAT.EQ.0) THEN
         IUPWIG=0
         CALL RUNMENT('WARNING: IUPWIG SET TO ZERO FOR STEADY-STATE')
        ENDIF
       ENDIF
      ENDIF                              ! wordst(2).eq.'petro'
      GO TO 1000
C
    9 NPOROT=6                           ! porosity criteria
      IF(KSGAUT.EQ.0)
     . CALL RUNMENT('WARNING: SMOOTHING IS NECESSARY FOR POROSITY')
      GO TO 1000
C
   10 NACTIT=1                           ! active elements
      GO TO 1000
C
   11 IF(WORDST(2).EQ.'NON_C') NOCOIT=1  ! non-coincident mesh
      IF(WORDST(2).EQ.'BOTH_') NOCOIT=2  ! both coinc. & non-coinc. mesh
      IF(NOCOIT.GT.0) THEN
       IF(NMEMO4.NE.1)
     .  CALL RUNENDT('ERROR: NMEMO4 MUST BE EQUAL 1 FOR NON-COIN. MESH')
      ENDIF
      GO TO 1000
C
   12 CONTINUE
C
C**** ADDITIONAL CONTROLS
C
      IF(ICONVT.EQ.1.AND.LARGET.NE.0) THEN
       CALL RUNMENT('WARNING : ICONVT=1 & LARGET NE 0 => INCOMPATIBLE;
     . LARGET IS SET TO 0')
       LARGET=0
      ENDIF
      IF(ICONVT.EQ.0.AND.ITERMEF.EQ.1)
     . CALL RUNMENT('WARNING: ICONVT=0 & ITERMEF=1                 ')
      IF(IGALET.EQ.1) THEN
       IF(IUPWIC.EQ.1.AND.NDIMET.EQ.1) THEN
        IUPWIC=0
        CALL RUNMENT('WARNING: NO CROSSWIND FOR 1D PROBLEMS')
       ENDIF
      ENDIF
C
      RETURN 
 2000 CALL RUNENDT('PROINPT:ERROR IN PROBLEM DATA BLOCK')
      END
