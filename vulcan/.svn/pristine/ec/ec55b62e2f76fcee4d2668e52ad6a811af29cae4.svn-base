      SUBROUTINE PROPMIC6(ITAPET,PROPSS,NPOI2T,INUPM,NLINET)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MICROSTRUCTURAL MATERIAL PROPERTIES FOR
C     MODEL 6 (S.G. CAST IRON MODEL)
C                                                      
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'       ! thermal-microestructural
C
C**** MICROSCOPICAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      DIMENSION PROPSS(*)
C
C**** LOOK FOR 'VERSION' CARD
C
      NPRINT=0
      CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'VERSI') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=INT(PARAMS(1))
         WRITE(LURESS,829) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPMIC6: ERROR VERSION CARD NOT FOUND')
      ENDIF
C
      IVERSI=INT(PROPSS(NPOI2T))

c*****************************************************************
c     
c     VERSION 1 (SOLIDUS-SOLIDUS PHASE-CHANGE MODEL)
c
c*****************************************************************
      IF(IVERSI.EQ.1) THEN
    
C***  LOOK FOR 'NUCLEATION' CARD
         NPRINT=0
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'NUCLE') THEN
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! Ts
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! Kco
            WRITE(LURESS,830) PROPSS(NPOI2T-1),PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: NUCLEATION CARD NOT FOUND')
         ENDIF

C***  LOOKS FOR 'GROWTH' CARD
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'GROWT') THEN
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! kc1
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! kc2
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(3) ! nc
            WRITE(LURESS,831) PROPSS(NPOI2T-2),PROPSS(NPOI2T-1),
     .           PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: GROWTH CARD NOT FOUND')
         ENDIF

C***  LOOKS FOR 'LATENT HEAT' CARD
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'LATEN') THEN
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! L
            WRITE(LURESS,832) PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: LATENT HEAT CARD NOT FOUND')
         ENDIF
      ENDIF                     ! iversi.eq.1
C     


c*****************************************************************
c     
c     VERSION 2
c
c*****************************************************************
      IF(IVERSI.EQ.2) THEN
         
C***  LOOK FOR 'EUTECTOID_TEMPERATURES' CARD
         NPRINT=0
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'EUTEC') THEN
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! Tf
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! Tp
            WRITE(LURESS,833) PROPSS(NPOI2T-1),PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: EUTECTOID_TEMPERATURES
     .           CARD NOT FOUND')
         ENDIF
         
C***  LOOK FOR 'NUCLEATION' CARD
         NPRINT=0
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'NUCLE') THEN
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! Af
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! Ap
            WRITE(LURESS,834) PROPSS(NPOI2T-1),PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: NUCLEATION CARD NOT FOUND')
         ENDIF
         
C***  LOOKS FOR 'GROWTH' CARD
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'GROWT') THEN
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! Bf
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! Bp
            WRITE(LURESS,835) PROPSS(NPOI2T-1),PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: GROWTH CARD NOT FOUND')
         ENDIF
     
C**** LOOKS FOR 'LATENT HEAT' CARD
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'LATEN') THEN
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! L
            WRITE(LURESS,832) PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: LATENT HEAT CARD NOT FOUND')
         ENDIF
      ENDIF                     ! iversi.eq.2


c*****************************************************************
c     
c     VERSION 3 (Austempered Ductile Iron model with avrami's law)
c
c*****************************************************************
      IF(IVERSI.EQ.3) THEN
         NPRINT=0
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'INITI') THEN
c     number of phase change models that provide initial condition to
c     micros6/microm6 v3
            NCOPC=0
            DO I=1,NLINET
               IF(INT(PARAMS(I)).NE.0) NCOPC=NCOPC+1
            ENDDO
            IPLAC(INUPM,1,1)=NCOPC ! total number of ph-ch for ic
            
            IF(NCOPC.GT.1) THEN
               CALL RUNMENS('PROPMIC6(ERROR): NCOPC>1 initial')
               CALL RUNENDS('condition from only one model')
            ENDIF
            
            WRITE(LURESS,709) (INT(PARAMS(ICOPC)), ICOPC=1,NCOPC)
            WRITE(LURESS,899)
            
            DO ICOPC=1,NCOPC
               IPLAC(INUPM,ICOPC+1,1)= INT(PARAMS(ICOPC)) ! ph-ch numb. for ic
               IF(IPLAC(INUPM,ICOPC+1,1).EQ.INUPM)
     .              CALL RUNENDS('ERROR: SELF-COUPLING IS NOT POSSIBLE')
            ENDDO
         ELSE
            WRITE(LURESS,710)
            WRITE(LURESS,899)
c     chemical composition
            IF(WORDSS(1).NE.'CHEMI') THEN
               CALL RUNENDS('PROPMIC6: ERROR CHEMI COMPOSITION CARD
     .              NOT FOUND')
            ELSE
               WRITE(LURESS,711)
               WRITE(LURESS,899)
            ENDIF
         ENDIF

C         
C--   INITIAL CONDITION FROM INPUT FILE
C     
         IF (INT(IPLAC(INUPM,1,1)).EQ.0) THEN
c--   CHEMICAL COMPOSITION
c     CARBON
            NPRINT=0
            CALL LISTENS('PROPM14',NPRINT,ITAPET)
            IF(WORDSS(1).EQ.'CARBO') THEN
               NPOI2T=NPOI2T+1
               PROPSS(NPOI2T)=PARAMS(1) ! c
               WRITE(LURESS,712) PROPSS(NPOI2T)
               WRITE(LURESS,899)
            ELSE
               CALL RUNENDS('PROPM14: CARBON CARD NOT FOUND')
            ENDIF
C     silicon
            NPRINT=0
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            IF(WORDSS(1).EQ.'SILIC') THEN
               NPOI2T=NPOI2T+1
               PROPSS(NPOI2T)=PARAMS(1)
            ELSE
               CALL RUNENDS('PROPMIC6: ERROR NO SILICON CARD')
            ENDIF
            IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.4.0D0)
     .           CALL RUNENDS('PROPMIC6: ERROR WRONG SILICON
     .           COMPOSITION')
            WRITE(LURESS,713) PROPSS(NPOI2T)
            WRITE(LURESS,899)
C     copper
            NPRINT=0
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            IF(WORDSS(1).EQ.'COPPE') THEN
               NPOI2T=NPOI2T+1
               PROPSS(NPOI2T)=PARAMS(1)
            ELSE
               CALL RUNENDS('PROPMIC6: ERROR NO COPPER CARD')
            ENDIF
            IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.3.0D0)
     .           CALL RUNENDS('PROPMIC6: ERROR WRONG COPPER
     .           COMPOSITION')
            WRITE(LURESS,714) PROPSS(NPOI2T)
            WRITE(LURESS,899)
C     manganese
            NPRINT=0
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            IF(WORDSS(1).EQ.'MANGA') THEN
               NPOI2T=NPOI2T+1
               PROPSS(NPOI2T)=PARAMS(1)
            ELSE
               CALL RUNENDS('PROPMIC6: ERROR NO MANGANESE CARD')
            ENDIF
            IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.2.0D0)
     .           CALL RUNENDS('PROPMIC6: ERROR WRONG MANGANESE 
     .           COMPOSITION')
            WRITE(LURESS,715) PROPSS(NPOI2T)
            WRITE(LURESS,899)
C     nikel
            NPRINT=0
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            IF(WORDSS(1).EQ.'NICKE') THEN
               NPOI2T=NPOI2T+1
               PROPSS(NPOI2T)=PARAMS(1)
            ELSE
               CALL RUNENDS('PROPMIC6: ERROR NO NICKEL CARD')
            ENDIF
            IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.0.5D0)
     .           CALL RUNENDS('PROPMIC6: ERROR WRONG NICKEL
     .           COMPOSITION')
            WRITE(LURESS,716) PROPSS(NPOI2T)
            WRITE(LURESS,899)
C     molybdenum
            NPRINT=0
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            IF(WORDSS(1).EQ.'MOLYB') THEN
               NPOI2T=NPOI2T+1
               PROPSS(NPOI2T)=PARAMS(1)
            ELSE
               CALL RUNENDS('PROPMIC6: ERROR NO MOLYBDENUM CARD')
            ENDIF
            IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.0.5D0)
     .           CALL RUNENDS('PROPMIC6: ERROR WRONG MOLYBDENUM
     .           COMPOSITION')
            WRITE(LURESS,717) PROPSS(NPOI2T)
            WRITE(LURESS,899)
C     chromiun
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            IF(WORDSS(1).EQ.'CHROM') THEN
               NPOI2T=NPOI2T+1
               PROPSS(NPOI2T)=PARAMS(1)
            ELSE
               CALL RUNENDS('PROPMIC6: ERROR NO CHROMIUM CARD')
            ENDIF
            IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.0.5D0)
     .           CALL RUNENDS('PROPMIC6: ERROR WRONG CHROMIUM
     .           COMPOSITION')
            WRITE(LURESS,718) PROPSS(NPOI2T)
            WRITE(LURESS,899)
            
C--   LOOK FOR INITIAL AUSTENITE FRACTION CARD
            NPRINT=0
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            IF(WORDSS(1).EQ.'INITI') THEN
               NPOI2T=NPOI2T+1
               PROPSS(NPOI2T)=PARAMS(1) ! f_gamma_0
               WRITE(LURESS,719) PROPSS(NPOI2T)
               WRITE(LURESS,899)
            ELSE
               CALL RUNENDS('PROPMIC6: ERROR INITIAL AUSTENITE FR.
     .              CARD NOT FOUND')
            ENDIF

C--   LOOKS FOR AUSTENITIZATION TEMPERATURE CARD
            NPRINT=0
            CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
            IF(WORDSS(1).EQ.'AUSTE') THEN
               NPOI2T=NPOI2T+1
               PROPSS(NPOI2T)=PARAMS(1) ! f_gamma_0
               WRITE(LURESS,720) PROPSS(NPOI2T)
               WRITE(LURESS,899)
            ELSE
               CALL RUNENDS('PROPMIC6: ERROR AUSTENITIZATION 
     .              TEMPERATURE CARD NOT FOUND')
            ENDIF
         ENDIF

C
C--   LOOKS FOR BAINITE START TEMPERATURE CARD
C     
         NPRINT=0
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'BAINI') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1)
            WRITE(LURESS,721) PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: ERROR BAINITE START
     .           TEMPERATURE CARD NOT FOUND')
         ENDIF

C     
C--   LOOKS FOR BAINITE FINISH TEMPERATURE CARD
C     
         NPRINT=0
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'BAINI') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1)
            WRITE(LURESS,727) PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: ERROR BAINITE FINISH
     .           TEMPERATURE CARD NOT FOUND')
         ENDIF

C     
C--   LOOKS FOR 'CARBON CONCENTRATION IN BAINITIC FERRITE' CARD
C     
         NPRINT=0
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'CARBO') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! C_alpha
            WRITE(LURESS,722) PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: ERROR CARBON CONC. CARD NOT FOUND')
         ENDIF

C     
C--   LOOKS FOR 'AVRAMI PARAMETERS' CARD
C     
         NPRINT=0
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'AVRAM') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! k
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! m
            WRITE(LURESS,723) PROPSS(NPOI2T-1),PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: ERROR AVRAMI PARAMETERS CARD
     .           NOT FOUND')
         ENDIF

c     
c--   mechanicals parameters
c     
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET) ! lee una linea completa
         IF(WORDSS(1).NE.'MECHA') THEN
            CALL RUNENDS('PROPMIC6: MECHANICALS PARAM CARD NOT FOUND')
         ENDIF
C--   version mechanical model
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'MECHA') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) !mechanical model
            WRITE(LURESS,724) PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: MECHANICAL MODEL CARD NOT FOUND')
         ENDIF
C--   parameter of mechanical model 1
         IF (PARAMS(1).NE.1) CALL RUNENDS('PROPMIC6: ONLY MECHANICAL
     .        MODEL 1 IMPLEMENTED')
C     thermal dilatation ferrite
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'THERM') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) !alphaf
            WRITE(LURESS,725) PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: THERMAL DILATATION FERRITE
     .           CARD NOT FOUND')
         ENDIF
C     thermal dilatation austenite
         CALL LISTENS('PROPMIC6',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'THERM') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) !alphaa
            WRITE(LURESS,726) PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPMIC6: THERMAL DILATATION AUSTENITE 
     .           CARD NOT FOUND')
         ENDIF
         
      ENDIF
      
c      open(UNIT=12,FILE='salida.txt',STATUS='UNKNOWN')
c      write(12,*) NLINET
c      close(UNIT=12)
      
C     
C**** OUTPUT (see pointes.f)
C     
c     set to one
      KPLA6S=1
c     model number
      IPLLLS(INUPM)=6
C     
      RETURN
C     version: 3
 709  FORMAT(3X,'INITIAL CONDITION FROM PHASE-CHANGE NUMBER=',5I5)
 710  FORMAT(3X,'INITIAL CONDITION FROM INPUT FILE')
 711  FORMAT(3X,'INITIAL CHEMICAL COMPOSITION:')
 712  FORMAT(3X,'  CARBON CONTENT=',E15.6)
 713  FORMAT(3X,'  SILICON CONTENT=',E15.6)
 714  FORMAT(3X,'  COPPER CONTENT=',E15.6)
 715  FORMAT(3X,'  MANGANESE CONTENT=',E15.6)
 716  FORMAT(3X,'  NICKEL CONTENT=',E15.6)
 717  FORMAT(3X,'  MOLYBDENUM CONTENT=',E15.6)
 718  FORMAT(3X,'  CHROMIUM CONTENT=',E15.6)
 719  FORMAT(3X,'INITIAL AUSTENITE FRACTION=',E15.6)
 720  FORMAT(3X,'AUSTENITIZATION TEMPERATURE=',E15.6)
 721  FORMAT(3X,'BAINITE START TEMPERATURE=',E15.6)
 722  FORMAT(3X,'CARBON CONCENTRATION IN BAINITIE=',E15.6)
 723  FORMAT(3X,'AVRAMI PARAMETERS=',E15.6,E15.6)
 724  FORMAT(/,'MECHANICAL MODEL=',E15.6)
 725  FORMAT(3X,'THERMAL DILATATION FERRITE=',E15.6)
 726  FORMAT(3X,'THERMAL DILATATION AUSTENITE=',E15.6)
 727  FORMAT(3X,'BAINITE FINISH TEMPERATURE=',E15.6)
 728  FORMAT(3X,'NO INITIAL CONDITION FROM OTHER PHASE-CHANGE')
C     other versions:
 829  FORMAT(/,'MODEL - VERSION=',E15.6)
 830  FORMAT(/,'  EUTECTOID TEMPERATURE               =',E15.6,
     .     /,'  EUTECTOID NUCLEATION TIME           =',E15.6)
 831  FORMAT(/,'  EUTECTOID GROWTH LAW: CONSTANT 1    =',E15.6,
     .     /,'  EUTECTOID GROWTH LAW: CONSTANT 2    =',E15.6,
     .     /,'  EUTECTOID GROWTH LAW: EXPONENT      =',E15.6)
 832  FORMAT(/,'  LATENT HEAT VALUE                   =',E15.6)
 833  FORMAT(/,'  FERRITE EUTECTOID TEMPERATURE       =',E15.6,
     .     /,'  PEARLITE EUTECTOID TEMPERATURE      =',E15.6)
 834  FORMAT(/,'  FERRITE NUCLEATION CONSTANT         =',E15.6,
     .     /,'  PEARLITE NUCLEATION CONSTANT        =',E15.6)
 835  FORMAT(/,'  FERRITE GROWTH CONSTANT             =',E15.6,
     .     /,'  PEARLITE GROWTH CONSTANT            =',E15.6)
 899  FORMAT(/)
 900  FORMAT(//,'TOO MANY MATERIAL PROPERTIES TO READ:',/,
     .     20X,'NPROP =',I5,5X,'NUMBE =',I5,/)
      END
