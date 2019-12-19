      SUBROUTINE HACHEC(ELCODT,PROPST,ELDIST,VELCMT,HACHET,CENTRO,
     .                  ADVEMT,TEINIT,SHAPET,NSUPE, FPCHLT,
     .                  EHISTT,CARTDT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES SOME PARAMETERS OF THE WEIGHT FUNCTION
C     ACCORDING TO IUPWI, IPERT & ISUPW ( ELEMENT NO. 5 )
C
C
C     FOR PROBLEMS WITH ADVECTIVE EFFECTS:
C
C     IUPWI: DETERMINES THE UPWIND FUNCTION TO BE CONSIDERED
C
C        =0: NO UPWIND FUNCTION (DEFAULT)
C
C        =1: PROPERTIES & VELOCITY-DEPENDENT OPTIMAL FUNCTION
C            CONSIDERING ONE PECLET NUMBER
C
C        =2: PROPERTIES & VELOCITY-DEPENDENT CRITICAL FUNCTION
C            CONSIDERING ONE PECLET NUMBER
C
C        =3: PROPERTIES & VELOCITY-DEPENDENT OPTIMAL FUNCTION
C            CONSIDERING TWO PECLET NUMBERS
C
C        =4: PROPERTIES & VELOCITY-DEPENDENT CRITICAL FUNCTION
C            CONSIDERING TWO PECLET NUMBERS
C
C
C     IPERT: PERTURBATION FUNCTION TYPE TO BE CONSIDERED
C            (see whape1.f, whape2.f & whape3.f)
C
C        =1: "BUBBLE" FUNCTION WITHOUT COMPUTING K_p
C
C        =2: "BUBBLE" FUNCTION COMPUTING K_p
C
C        =3: CARTESIAN DERIVATIVE OF SHAPE FUNCTIONS (COMPUTING K_p)
C
C
C     ISUPW: DETERMINES THE PECLET DEFINITION TO BE CONSIDERED
C    
C        =1: STANDARD
C   
C        =2: STANDARD (averaged velocity & temperature)
C
C        =3: STANDARD (averaged velocity & properties)
C
C        =4: STANDARD + PHASE-CHANGE
C
C        =5: STANDARD + PHASE-CHANGE (averaged velocity, temperature &
C                                     df_pc/dT)
C
C        =6: STANDARD + PHASE-CHANGE (df_pc/dT at time t)
C
C
C     EXTUP: FACTOR MULTIPLYING THE PERTURBATION FUNCTION
C
C
C
C     FOR PROBLEMS WITH ADVECTIVE EFFECTS (CROSSWIND DISSIPATION):
C
C     IUPWIC: DETERMINES THE CROSSWIND FUNCTION TO BE CONSIDERED
C
C         =0: NO CROSSWIND FUNCTION (DEFAULT)
C
C         =1: Ramon Codina's proposal
C
C
C     IPERTC: CROSSWIND PERTURBATION FUNCTION TYPE TO BE CONSIDERED
C             (see whape1.f, whape2.f & whape3.f)
C
C         =1: CARTESIAN DERIVATIVE OF SHAPE FUNCTIONS (COMPUTING K_p)
C
C
C     ISUPWC: DETERMINES THE CROSSWIND PECLET DEFINITION TO BE
C             CONSIDERED
C
C         =1: STANDARD CONSIDERING u_||
C
C         =2: STANDARD + PHASE-CHANGE CONSIDERING u_||
C
C         =3: LINEARIZED VERSION OF 1
C
C         =4: LINEARIZED VERSION OF 2
C
C
C     EXTUPC: FACTOR MULTIPLYING THE CROSSWIND PERTURBATION FUNCTION
C
C
C
C     FOR PROBLEMS WITH HIGH TEMPERATURES GRADIENTS (TG):
C
C     IUPWIG: DETERMINES THE TEMPORAL UPWIND FUNCTION TO BE CONSIDERED
C
C        =0: NO TEMPORAL UPWIND FUNCTION (DEFAULT)
C
C        =1: DJC's proposal
C
C
C     IPERTG: PERTURBATION FUNCTION TYPE TO BE CONSIDERED
C             (see whape1.f, whape2.f & whape3.f)
C
C        =1: CARTESIAN DERIVATIVE OF SHAPE FUNCTIONS (COMPUTING K_p)
C
C
C     ISUPWG: DETERMINES THE FOURIER DEFINITION TO BE CONSIDERED
C    
C        =1: STANDARD
C   
C        =2: STANDARD + PHASE-CHANGE
C
C
C     EXTUPG: FACTOR MULTIPLYING THE PERTURBATION FUNCTION
C
C
C
C     NSUPE: number of element subdivisions
C
C     Warning: Isotropic Conductivity or Conductivity x (orthotropic
C              material) is only considered.
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
      INCLUDE 'nuec_om.f'        ! thermal-mechanical
      INCLUDE 'nued_om.f'        ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION ELCODT(NDIMET,*), PROPST(*),
     .          ELDIST(NDOFCT,*)
      DIMENSION VELCMT(*),        HACHET(NNODLT)
      DIMENSION CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          SHAPET(*)
      DIMENSION VELOCT(3),        VELOCP(3),        VELOCX(3)
      DIMENSION BASKK(9)
      DIMENSION EHISTT(NHISTT),   CARTDT(NDIMET,*)
C
      IX=NDIMETO-1
      VELTO=1.0D-10
C
C**** INITIALIZATION
C
      DO INODLT=1,NNODLT
       HACHET(INODLT)=0.0D0
      ENDDO
C
      IF(ICONVT.EQ.0) THEN                ! non-advection effects
       IF(IUPWIG.EQ.0) RETURN             ! no temp. grad. stabilization
      ENDIF
C
C**** COMPUTES CROSSWIND DISSIPATION
C
      IF(IUPWI.GT.0) THEN
C
C**** COMPUTES THE PECLET NUMBER ACCORDING TO ISUPW
C
       GO TO (11,12,13,14,15,16) ISUPW
C
   11  CONTINUE
C
C**** COMPUTES VELOCITY
C
       DO IDIMET=1,NDIMET
        VELOCT(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCT(IDIMET)=VELOCT(IDIMET)+
     .                  SHAPET(INODLT)*ADVEMT(IDIMET,INODLT)
        ENDDO
       ENDDO
C
C**** COMPUTES TEMPERATURE & TEMPERATURE RATE
C
       TGAUST=0.0D0
       DTEMPT=0.0D0
       TGAUIT=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        TGAUST=TGAUST+SHAPET(INODLT)*ELDIST(1,INODLT)
        IF(KDYNAT.EQ.1)
     .   DTEMPT=DTEMPT+SHAPET(INODLT)*VELCMT(INODLT)
        TGAUIT=TGAUIT+SHAPET(INODLT)*TEINIT(1,INODLT)
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO=PSEUDO+SHAPET(INODLT)*FPCHLT(IPSEU,INODLT)
        ENDIF
       ENDDO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
       TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** 1) COMPUTES DENSITY & SPECIFIC HEAT
C
       IF(IMICR.EQ.0) THEN
        ILAHET=0
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        BASMM=EHISTT(1)
        BASCC=EHISTT(2)
       ENDIF
C
C**** 2) COMPUTES CONDUCTIVITY
C
       IF(IMICR.EQ.0) THEN
        CALL CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
       ELSE
        BASKK(1)=EHISTT(3)
       ENDIF
C
C**** 3) COMPUTES THERMAL DIFFUSIVITY
C
       IPECN1=0
       THEDI1=0.0D0
       BASMC=BASMM*BASCC
       IF(BASMC.GT.0.0D0) THEN
        IPECN1=1
        THEDI1=BASKK(1)/BASMC
       ENDIF
       IPECN2=0
       THEDI2=0.0D0
       GO TO 100
C
   12  CONTINUE
C
C**** COMPUTES AVERAGE VELOCITY
C
       DO IDIMET=1,NDIMET
        VELOCT(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCT(IDIMET)=VELOCT(IDIMET)+ADVEMT(IDIMET,INODLT)
        ENDDO
        VELOCT(IDIMET)=VELOCT(IDIMET)/NNODLT
       END DO
C
C**** COMPUTES AVERAGE TEMPERATURE & TEMPERATURE RATE
C
       TAVER=0.0D0
       DTVER=0.0D0
       TAVEI=0.0D0
       PSAVER=0.0D0
       DO INODLT=1,NNODLT
        TAVER=TAVER+ELDIST(1,INODLT)
        IF(KDYNAT.EQ.1)
     .   DTVER=DTVER+VELCMT(INODLT)
        TAVEI=TAVEI+TEINIT(1,INODLT)
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSAVER=PSAVER+FPCHLT(IPSEU,INODLT)
        ENDIF
       ENDDO
       TAVER=TAVER/NNODLT
       DTVER=DTVER/NNODLT
       TAVEI=TAVEI/NNODLT
       PSAVER=PSAVER/NNODLT
C
       TGAUST=TAVER
       DTEMPT=DTVER
       TGAUIT=TAVEI
       PSEUDO=PSAVER
C
       TGAUSX=TGAUST
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** 1) COMPUTES DENSITY & SPECIFIC HEAT
C
       ILAHET=0
       CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .              DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .              DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
C
C**** 2) COMPUTES CONDUCTIVITY
C
       CALL CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
C
C**** 3) COMPUTES THERMAL DIFFUSIVITY
C
       THEDI=0.0D0
       BASMC=BASMM*BASCC
       IF(BASMC.GT.0.0D0) THEDI=BASKK(1)/BASMC
       GO TO 100
C
   13  CONTINUE
C
C**** COMPUTES AVERAGE VELOCITY
C
       DO IDIMET=1,NDIMET
        VELOCT(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCT(IDIMET)=VELOCT(IDIMET)+ADVEMT(IDIMET,INODLT)
        ENDDO
        VELOCT(IDIMET)=VELOCT(IDIMET)/NNODLT
       END DO
C
C**** COMPUTES TEMPERATURE & TEMPERATURE RATE
C
       TAVER=0.0D0
       DTVER=0.0D0
       TAVEI=0.0D0
       PSAVER=0.0D0

       basm1=0.0d0
       basc1=0.0d0
       bask1=0.0d0
       sour1=0.0d0
       sour2=0.0d0
C
       DO INODLT=1,NNODLT
        TAVER=ELDIST(1,INODLT)
        IF(KDYNAT.EQ.1)
     .   DTVER=VELCMT(INODLT)
        TAVEI=TEINIT(1,INODLT)
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSAVER=FPCHLT(IPSEU,INODLT)
        ENDIF
C
C**** COMPUTES PARAMETERS:
C
        TGAUST=TAVER
        DTEMPT=DTVER
        TGAUIT=TAVEI
        PSEUDO=PSAVER
C
        TGAUSX=TGAUST
        IF(KDYNAT.EQ.1) THEN
         IF(KINTET.EQ.1)              ! Euler method
     .    TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
        ENDIF
C
C**** 1) COMPUTES DENSITY & SPECIFIC HEAT
C
        ILAHET=0
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
        basm1=basm1+basmm
        basc1=basc1+bascc
C
C**** 2) COMPUTES L*(df_pc/dT)
C
        ILAHET=2
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
        sour1=sour1+sour1t
        sour2=sour2+sour2t
C
C**** 3) COMPUTES CONDUCTIVITY
C
        CALL CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
        bask1=bask1+baskk(1)
C
       ENDDO    ! inodlt=1,nnodlt

       basmm    =basm1/nnodlt
       bascc    =basc1/nnodlt
       baskk(1) =bask1/nnodlt
       sour1t   =sour1/nnodlt
       sour2t   =sour2/nnodlt
C
C**** 4) COMPUTES THERMAL DIFFUSIVITY
C
C     3 OPTIONS: a) STANDARD
C                b) STANDARD + PHASE-CHANGE
C                c) STANDARD + PHASE-CHANGE (at time t)
C
       THEDI=0.0D0
       IF(BASMM*BASCC.GT.0.0D0) THEN
        IF(ISUPW.EQ.1) THEDI=BASKK(1)/(BASMM*BASCC)
        IF(ISUPW.EQ.2) THEDI=BASKK(1)/(BASMM*(BASCC+SOUR2T))
        IF(ISUPW.EQ.3) THEDI=BASKK(1)/(BASMM*(BASCC+SOUR1T))
       ENDIF
       GO TO 100
C
   14 CONTINUE
C
C**** COMPUTES VELOCITY
C
       DO IDIMET=1,NDIMET
        VELOCT(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCT(IDIMET)=VELOCT(IDIMET)+
     .                  SHAPET(INODLT)*ADVEMT(IDIMET,INODLT)
        ENDDO
       ENDDO
C
C**** COMPUTES TEMPERATURE & TEMPERATURE RATE
C
       TGAUST=0.0D0
       DTEMPT=0.0D0
       TGAUIT=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        TGAUST=TGAUST+SHAPET(INODLT)*ELDIST(1,INODLT)
        IF(KDYNAT.EQ.1)
     .   DTEMPT=DTEMPT+SHAPET(INODLT)*VELCMT(INODLT)
        TGAUIT=TGAUIT+SHAPET(INODLT)*TEINIT(1,INODLT)
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO=PSEUDO+SHAPET(INODLT)*FPCHLT(IPSEU,INODLT)
        ENDIF
       ENDDO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
       TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** 1) COMPUTES DENSITY & SPECIFIC HEAT
C
       IF(IMICR.EQ.0) THEN
        ILAHET=0
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        BASMM=EHISTT(1)
        BASCC=EHISTT(2)
       ENDIF
C
C**** 2) COMPUTES L*(df_pc/dT)
C
       IF(IMICR.EQ.0) THEN
        ILAHET=2
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        SOUR2T=EHISTT(8+IX)
       ENDIF
C
C**** 3) COMPUTES CONDUCTIVITY
C
       IF(IMICR.EQ.0) THEN
        CALL CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
       ELSE
        BASKK(1)=EHISTT(3)
       ENDIF
C
C**** 4) COMPUTES THERMAL DIFFUSIVITY
C
       IPECN1=0
       THEDI1=0.0D0
       BASMC=BASMM*BASCC
       IF(BASMC.GT.0.0D0) THEN
        IPECN1=1
        THEDI1=BASKK(1)/BASMC
       ENDIF
       IPECN2=0
       THEDI2=0.0D0
       BASMC=BASMM*SOUR2T
       IF(BASMC.GT.0.0D0) THEN
        IPECN2=1
        THEDI2=BASKK(1)/BASMC
       ENDIF
       GO TO 100
C
   15  CONTINUE
C
C**** OTHER MODEL
C
       GO TO 100
C
   16  CONTINUE
C
C**** COMPUTES VELOCITY
C
       DO IDIMET=1,NDIMET
        VELOCT(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCT(IDIMET)=VELOCT(IDIMET)+
     .                  SHAPET(INODLT)*ADVEMT(IDIMET,INODLT)
        ENDDO
       ENDDO
C
C**** COMPUTES TEMPERATURE & TEMPERATURE RATE
C
       TGAUST=0.0D0
       DTEMPT=0.0D0
       TGAUIT=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        TGAUST=TGAUST+SHAPET(INODLT)*ELDIST(1,INODLT)
        IF(KDYNAT.EQ.1)
     .   DTEMPT=DTEMPT+SHAPET(INODLT)*VELCMT(INODLT)
         TGAUIT=TGAUIT+SHAPET(INODLT)*TEINIT(1,INODLT)
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO=PSEUDO+SHAPET(INODLT)*FPCHLT(IPSEU,INODLT)
        ENDIF
       ENDDO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
       TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** 1) COMPUTES DENSITY & SPECIFIC HEAT
C
       ILAHET=0
       CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .              DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .              DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
C
C**** 2) COMPUTES L*(df_pc/dT)
C
       ILAHET=2
       CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .              DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .              DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
C
C**** 3) COMPUTES CONDUCTIVITY
C
       CALL CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
C
C**** 4) COMPUTES THERMAL DIFFUSIVITY
C
       THEDI1=0.0D0
       BASMC=BASMM*BASCC
       IF(BASMC.GT.0.0D0) THEDI1=BASKK(1)/BASMC
       THEDI2=0.0D0
       BASMC=BASMM*SOUR1T
       IF(BASMC.GT.0.0D0) THEDI2=BASKK(1)/BASMC
       GO TO 100
C
  100  CONTINUE
C
C**** COMPUTES ELEMENTAL CHARACTERISTIC LENGTH
C
       CALL CHALEI(ELCODT,NDIMET,NNODLT,CHALE)
c      CHALE=CHALE/NSUPE
C
C**** COMPUTES VELOCITY NORM
C
       VELNO=0.0D0
       DO IDIMET=1,NDIMET
        VELNO=VELNO+VELOCT(IDIMET)*VELOCT(IDIMET)
       ENDDO
       IF(VELNO.GT.VELTO) THEN
        VELNO=DSQRT(VELNO)
       ELSE
        VELNO=0.0D0
       ENDIF
C
C**** COMPUTES ELEMENTAL PECLET NUMBER
C
       PECNU1=0.0D0
       IF(IPECN1.EQ.1) THEN
        PECNU1=VELNO*CHALE/(2.0D0*THEDI1)
        IPECN1=0
        IF(PECNU1.GT.0.0D0) IPECN1=1
       ENDIF
       PECNU2=0.0D0
       IF(IPECN2.EQ.1) THEN
        PECNU2=VELNO*CHALE/(2.0D0*THEDI2)
        IPECN2=0
        IF(PECNU2.GT.0.0D0) IPECN2=1
       ENDIF
C
C**** COMPUTES THE UPWIND FUNCTION ACCORDING TO IUPWI
C
       GO TO (21,22,23,24) IUPWI
C
   21  CONTINUE
C
C**** COMPUTES THE OPTIMUM UPWING FUNCTION
C
       UPCOE=0.0D0
       IF(IPECN1.EQ.1.OR.IPECN2.EQ.1) THEN
        PECNU=PECNU1+PECNU2
        UPCOE=1.0D0/(DTANH(PECNU))-1.0D0/PECNU
       ENDIF
       GO TO 200
C
   22  CONTINUE
C
C**** COMPUTES THE CRITICAL UPWING FUNCTION
C
       UPCOE=0.0D0
       IF(IPECN1.EQ.1.OR.IPECN2.EQ.1) THEN
        PECNU=PECNU1+PECNU2
        UPCOE=1.0D0
        UPCOX=PECNU/3.0D0
        IF(UPCOX.LT.UPCOE) UPCOE=UPCOX
       ENDIF
       GO TO 200
C
   23  CONTINUE
C
C**** COMPUTES THE OPTIMUM UPWING FUNCTION
C
       UPCOE1=0.0D0
       IF(IPECN1.EQ.1) THEN
        UPCOE1=1.0D0/(DTANH(PECNU1))-1.0D0/PECNU1
       ENDIF
       UPCOE2=0.0D0
       IF(IPECN2.EQ.1) THEN
        UPCOE2=1.0D0/(DTANH(PECNU2))-1.0D0/PECNU2
       ENDIF
c      UPCOE=UPCOE1+UPCOE2
       UPCOE=UPCOE1+UPCOE2*EXTUP
       GO TO 200
C
   24  CONTINUE
C
C**** COMPUTES THE CRITICAL UPWING FUNCTION
C
       UPCOE1=0.0D0
       IF(IPECN1.EQ.1) THEN
        UPCOE1=1.0D0
        UPCOX=PECNU1/3.0D0
        IF(UPCOX.LT.UPCOE1) UPCOE1=UPCOX
       ENDIF
       UPCOE2=0.0D0
       IF(IPECN2.EQ.1) THEN
        UPCOE2=1.0
        UPCOX=PECNU2/3.0D0
        IF(UPCOX.LT.UPCOE2) UPCOE2=UPCOX
       ENDIF
       UPCOE=UPCOE1+UPCOE2
       GO TO 200
C
  200  CONTINUE
C
C**** COMPUTES THE MULTIPLICATOR OF THE PERTURBATION FUNCTION ACCORDING
C     TO IPERT
C
       GO TO (31,31,33) IPERT
C
   31  CONTINUE
C
C**** COMPUTES THE VELOCITY-DEPENDENT PARAMETERS
C
       CALL CENTROI(ELCODT,NDIMET,NNODLT,CENTRO)
C
       DO INODLT=1,NNODLT
        HACHET(INODLT)=0.0D0
        DO IDIMET=1,NDIMET
         HACHET(INODLT)=HACHET(INODLT)+VELOCT(IDIMET)*
     .                  CENTRO(IDIMET,INODLT)
        ENDDO
C
        CENNO=0.0D0
        DO IDIMET=1,NDIMET
         CENNO=CENNO+CENTRO(IDIMET,INODLT)*CENTRO(IDIMET,INODLT)
        ENDDO
        CENTO=1.0D-10
        IF(CENNO.GT.CENTO) THEN
         CENNO=DSQRT(CENNO)
        ELSE
         CENNO=0.0D0
        ENDIF
C
        VELCE=VELNO*CENNO
        IF(VELCE.GT.VELTO) THEN
c        HACHET(INODLT)=UPCOE*EXTUP*HACHET(INODLT)/(VELNO*CENNO)
         HACHET(INODLT)=UPCOE*HACHET(INODLT)/(VELNO*CENNO)
        ELSE
         HACHET(INODLT)=0.0D0
        ENDIF
       ENDDO           ! inodlt=1,nnodlt
C
       DO INODLT=1,NNODLT
        HACHET(INODLT)=-HACHET(INODLT)
       ENDDO
       GO TO 300
C
   33  CONTINUE
C
C**** NOTES:
C
C     IN THIS CASE, HACHET(NDIMET) (HACHET has been dimensioned with 
C     NNODET, but NNODET > NDIMET always)
C
C     tau=EXTUP*UPCOE*CHALE/(2.0*VELNO)
C
       VELCE=VELNO
       IF(VELCE.GT.VELTO) THEN
        DO IDIMET=1,NDIMET
         HACHET(IDIMET)=EXTUP*VELOCT(IDIMET)*UPCOE*CHALE/(2.0D0*VELNO)
c        HACHET(IDIMET)=VELOCT(IDIMET)*UPCOE*CHALE/(2.0D0*VELNO)
        ENDDO
       ELSE
        DO IDIMET=1,NDIMET
         HACHET(IDIMET)=0.0D0
        ENDDO
       ENDIF
       GO TO 300
C
  300  CONTINUE
      ENDIF               ! iupwi.gt.0
C
C**** COMPUTES CROSSWIND DISSIPATION
C
      IF(IUPWIC.GT.0) THEN
C
C**** COMPUTES THE PECLET NUMBER ACCORDING TO ISUPWC
C
       GO TO (41,42,41,42,43) ISUPWC
C
   41  CONTINUE
C
C**** COMPUTES VELOCITY & VELOCITY ||
C
       VELOXX=0.0D0
       DO IDIMET=1,NDIMET
        VELOCT(IDIMET)=0.0D0
        VELOCP(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCT(IDIMET)=VELOCT(IDIMET)+
     .                  SHAPET(INODLT)*ADVEMT(IDIMET,INODLT)
         VELOCP(IDIMET)=VELOCP(IDIMET)+                    ! grad T
     .                  CARTDT(IDIMET,INODLT)*ELDIST(1,INODLT)
         IF(ISUPWC.EQ.3)
     .    VELOCP(IDIMET)=VELOCP(IDIMET)-                   ! grad T at t
     .                   CARTDT(IDIMET,INODLT)*VELCMT(INODLT)*DTIMET
        ENDDO
        VELOXX=VELOXX+VELOCT(IDIMET)*VELOCP(IDIMET)
       ENDDO
C
       TEMNO=0.0D0
       DO IDIMET=1,NDIMET
        TEMNO=TEMNO+VELOCP(IDIMET)*VELOCP(IDIMET)
       ENDDO
       IF(TEMNO.GT.VELTO) THEN
        DO IDIMET=1,NDIMET
         VELOCP(IDIMET)=VELOCP(IDIMET)*VELOXX/TEMNO
        ENDDO
       ELSE
        DO IDIMET=1,NDIMET
         VELOCP(IDIMET)=0.0D0
        ENDDO
       ENDIF
C
C**** COMPUTES TEMPERATURE & TEMPERATURE RATE
C
       TGAUST=0.0D0
       DTEMPT=0.0D0
       TGAUIT=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        TGAUST=TGAUST+SHAPET(INODLT)*ELDIST(1,INODLT)
        IF(KDYNAT.EQ.1)
     .   DTEMPT=DTEMPT+SHAPET(INODLT)*VELCMT(INODLT)
         TGAUIT=TGAUIT+SHAPET(INODLT)*TEINIT(1,INODLT)
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO=PSEUDO+SHAPET(INODLT)*FPCHLT(IPSEU,INODLT)
        ENDIF
       ENDDO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
       TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** 1) COMPUTES DENSITY & SPECIFIC HEAT
C
       IF(IMICR.EQ.0) THEN
        ILAHET=0
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        BASMM=EHISTT(1)
        BASCC=EHISTT(2)
       ENDIF
C
C**** 2) COMPUTES CONDUCTIVITY
C
       IF(IMICR.EQ.0) THEN
        CALL CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
       ELSE
        BASKK(1)=EHISTT(3)
       ENDIF
C
C**** 3) COMPUTES THERMAL DIFFUSIVITY
C
       IPECN1=0
       THEDI1=0.0D0
       BASMC=BASMM*BASCC
       IF(BASMC.GT.0.0D0) THEN
        IPECN1=1
        THEDI1=BASKK(1)/BASMC
       ENDIF
       IPECN2=0
       THEDI2=0.0D0
       GO TO 400
C
   42  CONTINUE
C
C**** COMPUTES VELOCITY & VELOCITY ||
C
       VELOXX=0.0D0
       DO IDIMET=1,NDIMET
        VELOCT(IDIMET)=0.0D0
        VELOCP(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCT(IDIMET)=VELOCT(IDIMET)+
     .                  SHAPET(INODLT)*ADVEMT(IDIMET,INODLT)
         VELOCP(IDIMET)=VELOCP(IDIMET)+                    ! grad T
     .                  CARTDT(IDIMET,INODLT)*ELDIST(1,INODLT)
         IF(ISUPWC.EQ.4)
     .    VELOCP(IDIMET)=VELOCP(IDIMET)-                   ! grad T at t
     .                   CARTDT(IDIMET,INODLT)*VELCMT(INODLT)*DTIMET
        ENDDO
        VELOXX=VELOXX+VELOCT(IDIMET)*VELOCP(IDIMET)
       ENDDO
C
       TEMNO=0.0D0
       DO IDIMET=1,NDIMET
        TEMNO=TEMNO+VELOCP(IDIMET)*VELOCP(IDIMET)
       ENDDO
       IF(TEMNO.GT.VELTO) THEN
        DO IDIMET=1,NDIMET
         VELOCP(IDIMET)=VELOCP(IDIMET)*VELOXX/TEMNO
        ENDDO
       ELSE
        DO IDIMET=1,NDIMET
         VELOCP(IDIMET)=0.0D0
        ENDDO
       ENDIF
C
C**** COMPUTES TEMPERATURE & TEMPERATURE RATE
C
       TGAUST=0.0D0
       DTEMPT=0.0D0
       TGAUIT=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        TGAUST=TGAUST+SHAPET(INODLT)*ELDIST(1,INODLT)
        IF(KDYNAT.EQ.1)
     .   DTEMPT=DTEMPT+SHAPET(INODLT)*VELCMT(INODLT)
        TGAUIT=TGAUIT+SHAPET(INODLT)*TEINIT(1,INODLT)
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO=PSEUDO+SHAPET(INODLT)*FPCHLT(IPSEU,INODLT)
        ENDIF
       ENDDO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
       TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** 1) COMPUTES DENSITY & SPECIFIC HEAT
C
       IF(IMICR.EQ.0) THEN
        ILAHET=0
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        BASMM=EHISTT(1)
        BASCC=EHISTT(2)
       ENDIF
C
C**** 2) COMPUTES L*(df_pc/dT)
C
       IF(IMICR.EQ.0) THEN
        ILAHET=2
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        SOUR2T=EHISTT(8+IX)
       ENDIF
C
C**** 3) COMPUTES CONDUCTIVITY
C
       IF(IMICR.EQ.0) THEN
        CALL CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
       ELSE
        BASKK(1)=EHISTT(3)
       ENDIF
C
C**** 4) COMPUTES THERMAL DIFFUSIVITY
C
       IPECN1=0
       THEDI1=0.0D0
       BASMC=BASMM*BASCC
       IF(BASMC.GT.0.0D0) THEN
        IPECN1=1
        THEDI1=BASKK(1)/BASMC
       ENDIF
       IPECN2=0
       THEDI2=0.0D0
       BASMC=BASMM*SOUR2T
       IF(BASMC.GT.0.0D0) THEN
        IPECN2=1
        THEDI2=BASKK(1)/BASMC
       ENDIF
       GO TO 400
C
   43  CONTINUE
C
C**** OTHER MODEL
C
       GO TO 400
C
  400  CONTINUE
C
C**** COMPUTES ELEMENTAL CHARACTERISTIC LENGTH
C
       CALL CHALEI(ELCODT,NDIMET,NNODLT,CHALE)
c      CHALE=CHALE/NSUPE
C
C**** COMPUTES VELOCITY & VELOCITY || NORMS
C
       VELNO=0.0D0
       VELNP=0.0D0
       DO IDIMET=1,NDIMET
        VELNO=VELNO+VELOCT(IDIMET)*VELOCT(IDIMET)
        VELNP=VELNP+VELOCP(IDIMET)*VELOCP(IDIMET)
       ENDDO
       IF(VELNO.GT.VELTO) THEN
        VELNO=DSQRT(VELNO)
       ELSE
        VELNO=0.0D0
       ENDIF
       IF(VELNP.GT.VELTO) THEN
        VELNP=DSQRT(VELNP)
       ELSE
        VELNP=0.0D0
       ENDIF
C
C**** COMPUTES ELEMENTAL PECLET NUMBER
C
       PECNU1C=0.0D0
       IF(IPECN1.EQ.1) THEN
        PECNU1C=VELNP*CHALE/(2.0D0*THEDI1)
        IPECN1=0
        IF(PECNU1C.GT.0.0D0) IPECN1=1
       ENDIF
       PECNU2C=0.0D0
       IF(IPECN2.EQ.1) THEN
        PECNU2C=VELNP*CHALE/(2.0D0*THEDI2)
        IPECN2=0
        IF(PECNU2C.GT.0.0D0) IPECN2=1
       ENDIF
C
C**** COMPUTES THE UPWIND FUNCTION ACCORDING TO IUPWIC
C
       GO TO (51,52) IUPWIC
C
   51  CONTINUE
C
C**** COMPUTES THE CROSSWIND FUNCTION
C
       UPCOEC=0.0D0
       IF(IPECN1.EQ.1.OR.IPECN2.EQ.1) THEN
        PECNUC=PECNU1C+PECNU2C
        UPCOEC=0.0D0
        UPCOXC=RAMON-1.0D0/PECNUC
        IF(UPCOXC.GT.UPCOEC) UPCOEC=UPCOXC
       ENDIF
       GO TO 500
C
   52  CONTINUE
C
C**** OTHER MODEL
C
       GO TO 500
C
  500  CONTINUE
C
C**** COMPUTES THE MULTIPLICATOR OF THE PERTURBATION FUNCTION ACCORDING
C     TO IPERTC
C
       GO TO (61,62) IPERTC
C
   61  CONTINUE
C
C**** NOTES:
C
C     IN THIS CASE, HACHET(NDIMET) (HACHET has been dimensioned with 
C     NNODET, but NNODET > NDIMET always)
C
C     tau=EXTUP*UPCOE*CHALE/(2.0*VELNP)
C
       VELOCX(1)=0.0D0
       VELOCX(2)=0.0D0
       VELOCX(3)=0.0D0
C
       VELCE=VELNO
       IF(VELCE.GT.VELTO) THEN
        IF(NDIMET.EQ.1) THEN
         VELOCX(1)=0.0D0                           ! no crosswind for 1D
        ELSE
         VELN2=VELNO*VELNO
         VELOCX(1)=VELOCP(1)*(1.0D0-VELOCT(1)*VELOCT(1)/VELN2)-
     .             VELOCP(2)*VELOCT(1)*VELOCT(2)/VELN2
         VELOCX(2)=VELOCP(2)*(1.0D0-VELOCT(2)*VELOCT(2)/VELN2)-
     .             VELOCP(1)*VELOCT(1)*VELOCT(2)/VELN2
         IF(NDIMET.EQ.3) THEN
          VELOCX(1)=VELOCX(1)-
     .              VELOCP(3)*VELOCT(1)*VELOCT(3)/VELN2
          VELOCX(2)=VELOCX(2)-
     .              VELOCP(3)*VELOCT(2)*VELOCT(3)/VELN2
          VELOCX(3)=VELOCP(3)*(1.0D0-VELOCT(3)*VELOCT(3)/VELN2)-
     .              VELOCP(1)*VELOCT(1)*VELOCT(3)/VELN2-
     .              VELOCP(2)*VELOCT(2)*VELOCT(3)/VELN2
         ENDIF
        ENDIF
       ENDIF
C
       VELCE=VELNP
       IF(VELCE.GT.VELTO) THEN
        DO IDIMET=1,NDIMET
         HACHET(IDIMET)=HACHET(IDIMET)+
     .                  EXTUPC*VELOCX(IDIMET)*UPCOEC*CHALE/(2.0D0*VELNP)
        ENDDO
       ENDIF
       GO TO 600
C
   62  CONTINUE
C
C**** OTHER MODEL
C
       GO TO 600
C
  600  CONTINUE
      ENDIF               ! iupwic.gt.0
C
C**** COMPUTES TEMPERATURE GRADIENT STABILIZATION (TEMPORAL UPWIND)
C
      IF(IUPWIG.GT.0) THEN
C
C**** COMPUTES THE FOURIER NUMBER ACCORDING TO ISUPWG
C
       GO TO (71,72) ISUPWG
C
   71  CONTINUE
C
C**** COMPUTES TEMPERATURE & TEMPERATURE RATE
C
       TGAUST=0.0D0
       DTEMPT=0.0D0
       TGAUIT=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        TGAUST=TGAUST+SHAPET(INODLT)*ELDIST(1,INODLT)
        IF(KDYNAT.EQ.1)
     .   DTEMPT=DTEMPT+SHAPET(INODLT)*VELCMT(INODLT)
        TGAUIT=TGAUIT+SHAPET(INODLT)*TEINIT(1,INODLT)
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO=PSEUDO+SHAPET(INODLT)*FPCHLT(IPSEU,INODLT)
        ENDIF
       ENDDO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
       TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** 1) COMPUTES DENSITY & SPECIFIC HEAT
C
       IF(IMICR.EQ.0) THEN
        ILAHET=0
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        BASMM=EHISTT(1)
        BASCC=EHISTT(2)
       ENDIF
C
C**** 2) COMPUTES CONDUCTIVITY
C
       IF(IMICR.EQ.0) THEN
        CALL CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
       ELSE
        BASKK(1)=EHISTT(3)
       ENDIF
C
C**** 3) COMPUTES THERMAL DIFFUSIVITY
C
       IFOUN1=0
       THEDI1=0.0D0
       BASMC=BASMM*BASCC
       IF(BASMC.GT.0.0D0) THEN
        IFOUN1=1
        THEDI1=BASKK(1)/BASMC
       ENDIF
       IFOUN2=0
       THEDI2=0.0D0
       GO TO 700
C
   72  CONTINUE
C
C**** COMPUTES TEMPERATURE & TEMPERATURE RATE
C
       TGAUST=0.0D0
       DTEMPT=0.0D0
       TGAUIT=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        TGAUST=TGAUST+SHAPET(INODLT)*ELDIST(1,INODLT)
        IF(KDYNAT.EQ.1)
     .   DTEMPT=DTEMPT+SHAPET(INODLT)*VELCMT(INODLT)
        TGAUIT=TGAUIT+SHAPET(INODLT)*TEINIT(1,INODLT)
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO=PSEUDO+SHAPET(INODLT)*FPCHLT(IPSEU,INODLT)
        ENDIF
       ENDDO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
       TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** 1) COMPUTES DENSITY & SPECIFIC HEAT
C
       IF(IMICR.EQ.0) THEN
        ILAHET=0
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        BASMM=EHISTT(1)
        BASCC=EHISTT(2)
       ENDIF
C
C**** 2) COMPUTES L*(df_pc/dT)
C
       IF(IMICR.EQ.0) THEN
        ILAHET=2
        CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        SOUR2T=EHISTT(8+IX)
       ENDIF
C
C**** 3) COMPUTES CONDUCTIVITY
C
       IF(IMICR.EQ.0) THEN
        CALL CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
       ELSE
        BASKK(1)=EHISTT(3)
       ENDIF
C
C**** 4) COMPUTES THERMAL DIFFUSIVITY
C
       IFOUN1=0
       THEDI1=0.0D0
       BASMC=BASMM*BASCC
       IF(BASMC.GT.0.0D0) THEN
        IFOUN1=1
        THEDI1=BASKK(1)/BASMC
       ENDIF
       IFOUN2=0
       THEDI2=0.0D0
       BASMC=BASMM*SOUR2T
       IF(BASMC.GT.0.0D0) THEN
        IFOUN2=1
        THEDI2=BASKK(1)/BASMC
       ENDIF
       GO TO 700
C
  700  CONTINUE
C
C**** COMPUTES ELEMENTAL CHARACTERISTIC LENGTH
C
       CALL CHALEI(ELCODT,NDIMET,NNODLT,CHALE)
c      CHALE=CHALE/NSUPE
C
C**** COMPUTES TEMPERATURE GRADIENT & ITS NORM
C
       TEMNO=0.0D0
       DO IDIMET=1,NDIMET
        VELOCP(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCP(IDIMET)=VELOCP(IDIMET)+                    ! grad T
     .                  CARTDT(IDIMET,INODLT)*ELDIST(1,INODLT)
        ENDDO
        TEMNO=TEMNO+VELOCP(IDIMET)*VELOCP(IDIMET)
       ENDDO
       IF(TEMNO.GT.VELTO) THEN
        TEMNO=DSQRT(TEMNO)
       ELSE
        TEMNO=0.0D0
       ENDIF
C
C**** COMPUTES ELEMENTAL FOURIER NUMBER
C
       FOUNU1=0.0D0
       IF(IFOUN1.EQ.1) THEN
        FOUNU1=THEDI1*DTIMET/(CHALE*CHALE)
        IFOUN1=0
        IF(FOUNU1.GT.0.0D0) IFOUN1=1
       ENDIF
       FOUNU2=0.0D0
       IF(IFOUN2.EQ.1) THEN
        FOUNU2=THEDI2*DTIMET/(CHALE*CHALE)
        IFOUN2=0
        IF(FOUNU2.GT.0.0D0) IFOUN2=1
       ENDIF
C
C**** COMPUTES THE UPWIND FUNCTION ACCORDING TO IUPWIG
C
       GO TO (81,82) IUPWIG
C
   81  CONTINUE
C
C**** COMPUTES THE DJC TEMPORAL UPWING FUNCTION
C
       UPCOEG=0.0D0
       IF(IFOUN1.EQ.1.OR.IFOUN2.EQ.1) THEN
        FOUNU=FOUNU1+FOUNU2
        UPCOEG=1.0D0/FOUNU
        IF(UPCOEG.GT.1.0D0) UPCOEG=1.0D0             ! to be revised!
       ENDIF
       GO TO 800
C
   82  CONTINUE
C
C**** OTHER MODEL
C
       GO TO 800
C
  800  CONTINUE
C
C**** COMPUTES THE MULTIPLICATOR OF THE PERTURBATION FUNCTION ACCORDING
C     TO IPERTG
C
       GO TO (91,92) IPERTG
C
   91  CONTINUE
C
C**** NOTES:
C
C     IN THIS CASE, HACHET(NDIMET) (HACHET has been dimensioned with 
C     NNODET, but NNODET > NDIMET always)
C
C     tau=EXTUPG*UPCOEG*CHALE*SIGN(DTEMPT)/(2.0*VELNO)
C
       VELCE=TEMNO
       IF(VELCE.GT.VELTO) THEN
        DTEMPS=1.0D0
        IF(DTEMPT.LT.0.0D0) DTEMPS=-1.0D0
        DO IDIMET=1,NDIMET
         HACHET(IDIMET)=HACHET(IDIMET)+
     .           EXTUPG*VELOCP(IDIMET)*UPCOEG*CHALE*DTEMPS/(2.0D0*TEMNO)
        ENDDO
       ENDIF
       IPERT=3            ! considered in whafun.f (to be revised!)
       GO TO 900
C
   92  CONTINUE
C
C**** OTHER MODEL
C
       GO TO 900
C
  900  CONTINUE
      ENDIF               ! iupwig.gt.0
C
      RETURN
      END
