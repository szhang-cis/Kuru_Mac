      SUBROUTINE PROPMIC12(ITAPET,PROPSS,NPOI2T,INUPM,NLINET)
C***********************************************************************
C     
C**** THIS ROUTINE READS THE MICROSTRUCTURAL MATERIAL PROPERTIES FOR
C     MODEL 12
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
C**** INITIAL CONDITION FROM OTHER PHASE-CHANGES
C     
C     Note:
C     In principle, the phase-changes that provide the initial condition
C     (ic) for the present model (eutectoid transformation) should be
C     saved in PROPSS and transferred to IPLAC in idepros.f.
C     However, this information is directly stored in IPLAC in the
C     present routine (propmic12.f) to avoid re-initialization of IPLAC
C     and, besides, to simplify the identification of the ALPHAM
C     pointers to be included in IPLAC (this last task in performed in
C     propmic.f).
C     
      NPRINT=0
      CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'INITI') THEN
         NCOPC=0
         DO I=1,NLINET
            IF(INT(PARAMS(I)).NE.0) NCOPC=NCOPC+1
         ENDDO
         IPLAC(INUPM,1,1)=NCOPC ! total number of ph-ch for ic
         IF(NCOPC.GT.1)
     .        CALL RUNMENS('WARNING: NCOPC > 1 - not implemented yet')
         WRITE(LURESS,710) (INT(PARAMS(ICOPC)), ICOPC=1,NCOPC)
         WRITE(LURESS,809)
         DO ICOPC=1,NCOPC
            IPLAC(INUPM,ICOPC+1,1)= INT(PARAMS(ICOPC)) ! ph-ch numb. for ic
            IF(IPLAC(INUPM,ICOPC+1,1).EQ.INUPM)
     .           CALL RUNENDS('ERROR: SELF-COUPLING IS NOT POSSIBLE')
         ENDDO
      ELSE
         WRITE(LURESS,711)
         WRITE(LURESS,809)
         IPLAC(INUPM,1,1)= 0    ! no initial condition from other ph-ch
      ENDIF
C     
C**** TEMPERATURES AND LATENT HEAT
C     
      CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'TEMPE') THEN
         CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)= PARAMS(1) ! stable temperature
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)= PARAMS(2) ! metastable temperature
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)= PARAMS(3) ! latent heat od ferrite
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)= PARAMS(4) ! latent heat od pearlite
         WRITE(LURESS,810) PROPSS(NPOI2T-3),PROPSS(NPOI2T-2),
     .        PROPSS(NPOI2T-1), PROPSS(NPOI2T)
         WRITE(LURESS,809)
      ELSE
         CALL RUNENDS('PROPMIC12: NO TEMPERATURE CARD')
      ENDIF
C     
C**** EUTECTOID GRAPHITE GROWTH PARAMETERS
C     
      CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'GROWT'.AND.WORDSS(2).EQ.'GRAPH') THEN
         CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)= PARAMS(1) ! kinetics params of graphite-ferrite
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)= PARAMS(2) ! ferrite difusion coeficient
         WRITE(LURESS,820) PROPSS(NPOI2T-1),PROPSS(NPOI2T)
         WRITE(LURESS,809)
      ELSE
         CALL RUNENDS('PROPMIC12: NO EUTECTOID_GRAPHITE GROWTH CARD')
      ENDIF
C     
C**** FERRITE GROWTH PARAMETERS
C     
      CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'GROWT'.AND.WORDSS(2).EQ.'FERRI'
     .     .AND.WORDSS(3).EQ.'MODEL') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! ferrite growing model
         IGRMFX=INT(PARAMS(1))
C     difusion from austenite to graphite and interface's equilibrium
C     difusion from austenite to austenite and graphite 
C     and interface's equilibrium
C     difusion from austenite to austenite and graphite, interface's
C     equilibrium and graphite/ferrite interfacial reaction
         IF(IGRMFX.EQ.1.OR.IGRMFX.EQ.2.OR.IGRMFX.EQ.3) THEN
            CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! ferrite density
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! austenite density
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(3) ! carbon diffusion coef. in austenite
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(4) ! latent heat of ferrite
            WRITE(LURESS,830) PROPSS(NPOI2T-3),PROPSS(NPOI2T-2), 
     .           PROPSS(NPOI2T-1), PROPSS(NPOI2T)
         ELSE
            CALL RUNENDS('PROPMIC12: INCORRECT FERRITE GROWTH 
     .           MODEL NUMBER ') 
         ENDIF
      ELSE
         CALL RUNENDS('PROPMIC12: NO FERRITE GROWTH CARD')
      ENDIF
      WRITE(LURESS,809)
C     
C**** PEARLITE NUCLEATION PARAMETERS
C     
      CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'NUCLE'.AND.WORDSS(2).EQ.'PEARL'
     .     .AND.WORDSS(3).EQ.'MODEL') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! nucleation model index
         INUCMPX=INT(PARAMS(1))
C     Lacaze or Own model (por el momento desordenado)
         IF(INUCMPX.EQ.1.OR.INUCMPX.EQ.2.OR.INUCMPX.EQ.3) THEN 
            CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! coefficient pearlite nucleation
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! exponent nucleation law
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(3) ! other value not used yet
            WRITE(LURESS,840) PROPSS(NPOI2T-2),PROPSS(NPOI2T-1), 
     .           PROPSS(NPOI2T)
         ELSE
            CALL RUNENDS('PROPMIC12: INCORRECT NUCLEATION MODEL NUMBER')
         ENDIF
      ELSE
         CALL RUNENDS('PROPMIC12: NO PEARLITE NUCLEATION CARD')
      ENDIF
      WRITE(LURESS,809)
C     
C**** ARREST CRITERION INDEX FOR PEARLITE NUCLEATION
C     
      CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'ARRES') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! nucleation arrest criterion index
         INUCAX=INT(PARAMS(1))
         IF(INUCAX.NE.1.AND.INUCAX.NE.2) ! Trecal & Tmin
     .        CALL RUNENDS('PROPMIC12: INCORRECT ARREST 
     .        CRITERION NUMBER')
      ELSE
         CALL RUNENDS('PROPMIC12: PEARLITE NO ARREST CRITERION CARD')
      ENDIF
C     
C**** PEARLITE GROWTH PARAMETERS
C     
      CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'GROWT'.AND.WORDSS(2).EQ.'PEARL'
     .     .AND.WORDSS(3).EQ.'MODEL') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! growth model index
         IGROMPX=INT(PARAMS(1))
C     lacaze or indues models
         IF(IGROMPX.EQ.1.OR.IGROMPX.EQ.2.OR.IGROMPX.EQ.3) THEN 
            CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! coefficient pearlite nucleation
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! exponent growth law
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(3) ! exponent nucleation law
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(4) ! pearlite model growtha 
            WRITE(LURESS,850) PROPSS(NPOI2T-3),PROPSS(NPOI2T-2),
     .           PROPSS(NPOI2T-1), PROPSS(NPOI2T)
         ELSE
            CALL RUNENDS('PROPMIC12: INCORRECT GROWTH MODEL NUMBER')
         ENDIF
      ELSE
         CALL RUNENDS('PROPMIC12: NO PEARLITE GROWTH CARD')
      ENDIF
      WRITE(LURESS,809)
C     
C**** MICROSTRUCTURE-DEPENDENT THERMAL PROPERTIES
C     
      CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'CONDU'.AND.WORDSS(2).EQ.'MODEL') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! micro-dependent conductivity index
         IKMICX=INT(PARAMS(1))
         IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
            CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) ! solid conductivity
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(2) ! mushy conductivity
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(3) ! liquid conductivity
            WRITE(LURESS,851) PROPSS(NPOI2T-3),PROPSS(NPOI2T-2),
     .           PROPSS(NPOI2T-1),PROPSS(NPOI2T)
         ELSE
            CALL RUNENDS('PROPMIC12: INCORRECT CONDUCTIVITY MODEL 
     .           NUMBER')
         ENDIF
         WRITE(LURESS,809)
      ELSE
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=0.0D0   ! conductivity does not depend on micro.
         GO TO 1
      ENDIF
C     
C**** NUMERICAL STRATEGY PARAMETERS
C     
      CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
    1 IF(WORDSS(1).EQ.'NUMER') THEN
         WRITE(LURESS,860)
         CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'TEMPE') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1)
            IFPCDT=INT(PARAMS(1))
            IF(IFPCDT.LT.1.OR.IFPCDT.GT.5)
     .           CALL RUNENDS('PROPMIC12: WRONG IFPCDT VALUE')
            WRITE(LURESS,861) IFPCDT
            WRITE(LURESS,809)
         ELSE
            CALL RUNENDS('PROPMIC12: NO TEMPERATURE DERIVATIVE FORM 
     .           OF FL')
         ENDIF
C     
         CALL LISTENS('PROPMIC12',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'FRACT') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1)
            IAFLOJ=INT(PARAMS(1))
            WRITE(LURESS,862) IAFLOJ
            WRITE(LURESS,809)
         ELSE
            CALL RUNENDS('PROPMIC12: NO FRACTION CORRECTION CARD')
         ENDIF
      ELSE
         CALL RUNENDS('PROPMIC12: NO NUMERICAL STRATEGY CARD')
      ENDIF
C     
C**** OUTPUT (see pointes.f)
C     
      KPLA12S=1
      IPLLLS(INUPM)=12
C     
      RETURN
C     
 710  FORMAT(3X,'INITIAL CONDITION FROM PHASE-CHANGE NUMBER=',5I5)
 711  FORMAT(3X,'NO INITIAL CONDITION FROM OTHER PHASE-CHANGE')
 809  FORMAT(/)
 810  FORMAT(3X,'TEMPERATURES AND LATENT HEAT:',/,
     .     3X,' STABLE TEMPERATURE       =',E15.6,/,
     .     3X,' METASTABLE TEMPERATURE   =',E15.6,/,
     .     3X,' LATENT HEAT OF FERRITE   =',E15.6,/,
     .     3X,' LATENT HEAT OF PEARLITE  =',E15.6)
 820  FORMAT(3X,'GROWTH GRAPHITE PARAMETERS:',/,
     .     3X,' KINETIC PARAM. OF GRA/FERR INTERF. REACTION= ',E15.6,/,
     .     3X,' DIFFUSION COEF. OF C IN FERRITE            = ',E15.6)
 830  FORMAT(3X,'GROWTH FERRITE PARAMETERS:',/,
     .     3X,' FERRITE   DENSITY                =',E15.6,/,
     .     3X,' AUSTENITE DENSITY                =',E15.6,/,
     .     3X,' DIFFUSION COEF. OF C IN AUSTENITE=',E15.6,/,
     .     3X,' FERRITE LATENT HEAT              =',E15.6)
 840  FORMAT(3X,'NUCLEATION PEARLITE PARAMETERS:',/,
     .     3X,' COEFICIENT PEARLITE NUCELATION=',E15.6,/,
     .     3X,' EXPONENT PEARLITE NUCLEATION  =',E15.6,/,
     .     3X,' OTHER VALUE NOT USED YET      =',E15.6)
 850  FORMAT(3X,'GROWTH PEARLITE PARAMETERS:',/,
     .     3X,' COEFICIENT PEARLITE GROWTH =',E15.6,/,
     .     3X,' EXPONENT OF PEARLITE GROWTH=',E15.6,/,
     .     3X,' PEARLITE LATENT HEAT       =',E15.6,/,
     .     3X,' OTHER VALUE NOT USED YET   =',E15.6)
 851  FORMAT(3X,'MICROSTRUCTURE-DEPENDENT CONDUCTIVITY:',/,
     .     3X,' MODEL              =',E15.6,/,
     .     3X,' SOLID CONDUCTIVITY =',E15.6,/,
     .     3X,' MUSHY CONDUCTIVITY =',E15.6,/,
     .     3X,' LIQUID CONDUCTIVITY=',E15.6)
 860  FORMAT(3X,'NUMERICAL STRATEGY PARAMETERS')
 861  FORMAT(3X,' TEMPERATURE DERIVATIVE FORM OF FL=',I5)
 862  FORMAT(3X,' FRACTION CORRECTION FORM         =',I5)
C     
      END
