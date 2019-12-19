C***********************************************************************
C
C**** GENERAL DIMENSIONS OF PROGRAM VULCAN (thermal variables)
C
C     IF ANY CHANGE IS PRODUCED, COMPILE:
C
C       adddatt.f addelmt.f addprit.f addsolt.f addwort.f
C       datbast.f memodtt.f (comp-t)
C
C     IF "segmentation fault" is produced just starting VULCAN, check
C     maximum dimensions
C
C     Command "max" works for every MMACHI except for MMACHI=5 & 8
C     (see also para_omt.f)
C
C***********************************************************************
C
C     PROBLEM: crankshaft
C
C***********************************************************************
C
C**** WRITE MAXIMUM DIMENSIONS (0=WRITE, 1=WRITE & STOP, 2=WRITE MORE)
C
      PARAMETER(
     .   MWRITE=2 )
C
C**** GENERAL DIMENSIONS (change if necessary; see input data)
C
      PARAMETER(
     .   MPOINT=5450,
     .   MELEMT=5450,
     .   MDIMET=3,
     .   MNODET=20,
     .   MGAUST=8,
     .   MGRUPT=20,
     .   MMATST=20,
     .   MFUNCT=500 )
C
      PARAMETER(
     .   MDYNAT=1,       ! 0:steady-state; 1:transient
     .   MSMOMT=1,       ! smoothing (0:no; 1:yes)
     .   MRENUT=0,       ! renumbering (0:no; 1:yes)
     .   MHOURT=0,       ! hourglass control (0:no; 1:yes)
     .   MPORET=0 )      ! pore water pressure prob. (0:no; 2:yes)
C
      PARAMETER(
     .   MCONVT=1,       ! 0:no convection effects; 1:convec. effec.
     .   MGALET=1,       ! 0:no upwinding technique; 1:upwin. tech.
     .   MPOROT=6,       ! porosity criteria
     .   MACTIT=1,       ! 0: no active elements; 1: active elements
     .   MOCOIT=1 )      ! 0: coincident contact mesh; 1: non-coin. c.m.
C
C**** COUPLING (thermal-mechanical) DIMENSIONS
C     (change if necessary; see input data)
C
      PARAMETER(
     .   MTERME=1,  ! -1: pure ther. 0:no mech. coupling; 1:mech. coupl.
     .   MFPCH=8,   ! 2*NNUPT+NNUPO+NFILL+IGALFA+2*ICONVT*IEGFPC
C                                                        (see setdatt.f)
     .   MITERC=0 ) ! 0:standard stag. scheme; 1:improved stag. scheme
C
C**** COUPLING (thermal-microstructural) DIMENSIONS
C     (change if necessary; see input data)
C
      PARAMETER(
     .   MMICR=1,        ! 0:no microstructural coupling; 1:micr. coupl.
     .   MMICO=0 )       ! 0:weak micr. coupling; 1:full micr. coupl.
C
C**** OTHER DIMENSIONS (in general, do not change; however, if any 
C                       change is produced, each dimesion is set to 
C                       its corresponding maximum in addelmt.f)
C
C     NUMBER OF MATERIAL PROPERTIES
C
C     MPROPT:        thermal & microstructural properties
C     MPROPM:        thermal properties
C     MPROPT-MPROPM: microstructural properties
C
C     In a thermal analyisis: MPROPT=MPROPM
C
C     Recommended values (see below):
C     MPROPT=1000 & MPROPM=500 (thermal & microstructural analysis)
C     MPROPT=500 (thermal analysis only)
C
C
C     NUMBER OF INTERNAL VARIABLES
C
C     MHISTT:        thermal & microstructural variables
C     MHISTM:        microstructural variables
C     MHISTT-MHISTM: thermal variables
C
C     In a thermal analysis: MHISTM=0
C
C     Recommended values (see below):
C     MHISTT=310 & MHISTM=300 (thermal & microstructural analysis)
C     MHISTT=10 (thermal analysis only)
C
C     In all cases, MHISTT depends on MMEMO3 given below
C
      PARAMETER(
     .  MSTR1T=MDIMET,         ! nstr1t=ndimet; see consett.f
     .  MPROPT=500+MMICR*500,  ! see setdatt.f & setdats.f
     .  MDOFCT=1,              ! always; see consett.f
     .  MHLODT=17,             ! see setdatt.f
     .  MSUBFT=50,             ! see setdatt.f
     .  MPRELT=11,             ! see setdatt.f
     .  MNUM4T=300*MMICR,      ! see setdatd.f
c    .  MNUM4TM=(11*MNUM4T+46)*MMICR,  ! estimated; see setdatd.f
     .  MNUM4TM=(3*MNUM4T+6)*MMICR,  ! estimated; see setdatd.f
     .  MHIST1=10,             ! see consett.f
     .  MHISTT=(MHIST1+MNUM4TM)*MMICR, ! simplification; see consett.f
c    .  MHISTT=MHIST1+MNUM4TM*MMICR,   ! simplification; see consett.f
     .  MNUINT=MNUM4TM )       ! comp. of int.var.to print;see setdatt.f
C
C**** DEGREES OF FREEDOMS (do not change)
C
      PARAMETER(
     .   MEVABT=MNODET*MDOFCT,
     .   MTOTVT=MPOINT*MDOFCT )
C
      PARAMETER(                       ! mechanical dimensions in solids
     .   MDOFCM=MDOFCT*MDIMET,
     .   MTOTVM=MTOTVT*MDIMET )
C
C***********************************************************************
C
C**** ADDITIONAL PARAMETERS OF PROGRAM VULCAN (thermal problem)
C
C     it should be read in the input data !!!  (see check0t.f) ctm
C
C***********************************************************************
C
C**** ADDITIONAL PARAMETERS
C
C               MDISKD=0: database is out of core (datbast.f)
C               MDISKD=1: database is in core (datbast.f)
C
C               MFURES=0: a future restart will not be made (datrstt.f)
C               MFURES>0: a future restart will be made (datrstt.f)
C
C               MMACHI=1: CONVEX
C               MMACHI=2: SILICON GRAPHICS
C               MMACHI=3: VAX COMPUTER (not implemented)
C               MMACHI=4: SUN
C               MMACHI=5: PERSONAL COMPUTER
C               MMACHI=6: SILICON GRAPHICS POWER CHALLENGE
C               MMACHI=7: HEWLETT PACKARD
C               MMACHI=8: LINUX
C
C     THEORETICAL DEFAULTS: MDISKD=0, MMURES=1, MMACHI=1
C     REAL DEFAULTS:        MDISKD=1, MFURES=0, MMACHI=2
C
C***********************************************************************
C
      PARAMETER(
     .   MDISKD=1,
     .   MFURES=1,
     .   MMACHI=8 )
C
      PARAMETER(
c    .   MCHA1T=max((1+(MMACHI/6)*6-MMACHI),0),   ! sg (1 for MMACHI=6)
c    .   MCHA2T=max((1+(MMACHI/10)*10-MMACHI),0), ! sg (1 for MMACHI=10)
c    .   MCHAAT=max(MCHA1T,MCHA2T),               ! sg
c    .   MCHALT=max(4,8*MCHAAT) )                 ! sg
     .   MCHALT=4 )                               ! PC & linux
C
      PARAMETER(
c    .   MWINDT=1 )                               ! PC
     .   MWINDT=0 )                               ! linux
C
C***********************************************************************
C
C**** FULL OR PARTIAL MEMORY (MMEMO=1: full; MMEMO=0: partial)
C
C     MMEMO concerns to: EPMTX (elastic const. tensor)
C                        RMAT1 (anisotropic const. tensor)
C                        STRA0 (initial temperature gradient)
C                        STRS0 (initial heat flux)
C                        TEMPC (temperature)
C
C     THEORETICAL DEFAULT: MMEMO=1
C
C     (change if necessary; see input data [now in setdatt.f])
C
C***********************************************************************
C
      PARAMETER(
     .   MMEMO=0 )
C
C***********************************************************************
C
C**** ADDITIONAL MEMORY PARAMETERS
C
C     MMEMO1 concerns to: COORDT(NDIMET,NNODET) as a global array
C
C               MMEMO1=0: COORDT is an elemental array (ELCODT)
C               MMEMO1=1: COORDT is a global array
C
C     MMEMA1=1-MMEMO1 (do not change)
C
C     THEORETICAL DEFAULT: MMEMO1=0
C
C
C     MMEMO2 concerns to: shape functions, cartesian derivatives, 
C                         etc. computed initially or every time
C                         when needed
C
C               MMEMO2=0: SHAPE, CARTD, etc. computed initially
C               MMEMO2=1: SHAPE, CARTD, etc. computed every time
C
C     MMEMA1=1-MMEMO1 (do not change)
C
C     THEORETICAL DEFAULT: MMEMO2=0
C
C
C     MMEMO3 concerns to: history-dependent material thermal properties
C
C               MMEMO3=0: no history-dependent thermal properties
C               MMEMO3=1: history-dependent material properties 
C
C     THEORETICAL DEFAULT: MMEMO3=0
C
C
C     MMEMO4 concerns to: history-dependent heat flux
C
C               MMEMO4=0: no history-dependent heat flux
C               MMEMO4=1: history-dependent heat flux
C
C     THEORETICAL DEFAULT: MMEMO4=0
C
C
C     MMEMO5 concerns to: ELDIST in/out of ELVART
C
C               MMEMO5=0: ELDIST in ELVART
C               MMEMO5=1: ELDIST out of ELVART
C
C     THEORETICAL DEFAULT: MMEMO5=0
C
C
C     MMEMO6 concerns to: elemental assembly process in other/same
C                         NELEMT loop as the evaluation of each
C                         contribution of the jacobian matrix
C
C               MMEMO6=0: elemental assembly process in other loop
C               MMEMO6=1: elemental assembly process in the same loop
C
C     THEORETICAL DEFAULT: MMEMO6=0
C
C
C     MMEMO7 concerns to: the jacobian matrix is (not) evaluated in
C                         the NELEMT loop in the solver routines. 
C
C               MMEMO7=0: jacobian not evaluated in NELEMT solver loop
C               MMEMO7=1: jacobian evaluated in NELEMT solver loop. For
C                         this case, MMEMO6 must be equals 1
C
C     THEORETICAL DEFAULT: MMEMO7=0
C
C
C     MMEMO8 concerns to: the second derivative respect to time of the
C                         temperature
C
C               MMEMO8=0: d2T/dt2 not used
C               MMEMO8=1: d2T/dt2 used
C
C     THEORETICAL DEFAULT: MMEMO8=0
C
C
C     MMEMO9 concerns to: the second component of DISITT
C
C               MMEMO9=0: second component of DISITT not used
C               MMEMO9=1: second component of DISITT used
C
C     THEORETICAL DEFAULT: MMEMO9=0
C
C
C     MMEMO10 concerns to: the temperature-dependency of the density
C                          changes the geometry (volume and boundary)
C                          MMEMO10 is not compatible with ITERMD=1
C                          (see couinc.f)
C
C               MMEMO10=0: rho(T) does not change the geometry
C               MMEMO10=1: rho(T) changes the geometry
C
C     THEORETICAL DEFAULT: MMEMO10=0
C
C
C     MMEMO11 concerns to: normal gap & normal pressure to be considered
C                          in the heat transfer value computed in
C                          thermal or mechanical problems
C
C               MMEMO11=0: g_n & p_n computed in mechanical problem
C               MMEMO11=1: g_n & p_n computed in thermal problem
C
C     THEORETICAL DEFAULT: MMEMO11=0
C
C
C     (change if necessary; see input data [now in setdatt.f])
C
C***********************************************************************
C
      PARAMETER(
     .   MMEMO1=1,
     .   MMEMO2=1,
c    .   MMEMO3=max(0,MMICR),                       ! sg
     .   MMEMO3=MMICR,                              ! PC & linux
     .   MMEMO4=1,
     .   MMEMO5=1,
     .   MMEMO6A=1,
     .   MMEMO7A=0,
     .   MMEMO8=0,
     .   MMEMO9=0,
     .   MMEMO10=1,
     .   MMEMO11=1 )
C
C***********************************************************************
C
C     SOLVER VARIABLES
C
C     1) DEFAULTS OR CHOSEN VARIABLES
C     2) SOLVER TO BE USED
C     3) SYMMETRIC OR UNSYMMETRIC CASE
C
C     (change if necessary; see below)
C
C***********************************************************************
C
C**** 1) DEFAULTS (MDEFA1=1 & MDEFA2=0) OR
C        CHOSEN VARIABLES (MDEFA1=0 & MDEFA2=1)
C
C        (change if necessary; see input data)
C
C***********************************************************************
C
      PARAMETER(
     .   MDEFA1=0,
     .   MDEFA2=1 )
C
C***********************************************************************
C
C**** 2) SOLVER TO BE USED:
C
C      MSOLV1=1, MSOLV2=0, MSOLV3=0, MSOLV4=0 & MSOLV5=0: SKYLINE SOLVER
C      MSOLV1=0, MSOLV2=1, MSOLV3=0, MSOLV4=0 & MSOLV5=0: FRONTAL SOLVER
C      MSOLV1=0, MSOLV2=0, MSOLV3=1, MSOLV4=0 & MSOLV5=0: PCG SOLVER
C      MSOLV1=0, MSOLV2=0, MSOLV3=0, MSOLV4=1 & MSOLV5=0: GMRES SOLVER
C      MSOLV1=0, MSOLV2=0, MSOLV3=0, MSOLV4=0 & MSOLV5=1: EXPLICIT
C
C        (change if necessary; see input data)
C
C***********************************************************************
C
      PARAMETER(
     .   MSOLV1=1,
     .   MSOLV2=0,
     .   MSOLV3=0,
     .   MSOLV4=0,
     .   MSOLV5=0 )
C
C***********************************************************************
C
C**** 3) SYMMETRIC CASE (MUNSY1=1 & MUNSY2=0) (skyline & frontal)
C        OR
C        UNSYMMETRIC CASE (MUNSY1=0 & MUNSY2=1) (skyline & frontal)
C
C        Notes:
C        pcg: always symmetric case
C        explicit & gmres: always unsymmetric case
C        Therefore, for these solvers MUNSY1 and MUNSY2 are not used
C
C        (change if necessary; see input data)
C
C***********************************************************************
C
      PARAMETER(
     .   MUNSY1=0,
     .   MUNSY2=1 )
C
C***********************************************************************
C
C     (from here onwards change only for chosen variables of the 
C      chosen solver)
C
C***********************************************************************
C
C**** SKYLINE SOLVER
C
C     VARIABLES:
C     MEQNST: number of equations
C     MLASTT: bandwidth
C
      PARAMETER(
     .   MEQNST1=MDEFA1*MTOTVT,                 ! default
     .   MEQNST2=MDEFA2*90000 )                 ! chosen variable
C
      PARAMETER(
c    .   MEQNST=max(MEQNST1,MEQNST2) )          ! sg
     .   MEQNST=(MEQNST1+MEQNST2) )             ! PC & linux
C
      PARAMETER(
     .   MLASTT1=MDEFA1*MEQNST*(MEQNST+1)/2,    ! default
     .   MLASTT2=MDEFA2*30000000 )               ! chosen variable
C
      PARAMETER(
c    .   MLASTT=max(MLASTT1,MLASTT2) )          ! sg
     .   MLASTT=(MLASTT1+MLASTT2) )             ! PC & linux
C
      PARAMETER(         ! mfitet=0: direct solver
     .   MFITET=0 )      ! mfitet=1: semidirect solver
C
C**** FRONTAL SOLVER
C
C     VARIABLES:
C     MFRONT: frontwidth
C     MBUFAT: buffer size
C
      PARAMETER(
     .   MFRONT1=MDEFA1*MTOTVT,                 ! default
     .   MFRONT2=MDEFA2*2000 )                  ! chosen variable
C
      PARAMETER(
c    .   MFRONT=max(MFRONT1,MFRONT2) )          ! sg
     .   MFRONT=(MFRONT1+MFRONT2) )             ! PC & linux
C
      PARAMETER(
     .   MBUFAT1=MDEFA1*(MFRONT+1),             ! default
     .   MBUFAT2=MDEFA2*70 )                    ! chosen variable
C
      PARAMETER(
c    .   MBUFAT=max(MBUFAT1,MBUFAT2) )          ! sg
     .   MBUFAT=(MBUFAT1+MBUFAT2) )             ! PC & linux
C
      PARAMETER(
     .   MSTIFT1=MUNSY1*MFRONT*(MFRONT+1)/2,    ! symmetric case
     .   MSTIFT2=MUNSY2*MFRONT*MFRONT  )        ! unsymmetric case 
C
      PARAMETER(
c    .   MSTIFT=max(MSTIFT1,MSTIFT2) )          ! sg
     .   MSTIFT=(MSTIFT1+MSTIFT2) )             ! PC & linux
C
C**** PCG SOLVER
C
      PARAMETER(
     .   MWIDTT1=0,                             ! default
     .   MWIDTT2=MDEFA2 )                       ! chosen variable
C
      PARAMETER(
c    .   MWIDTT=max(MWIDTT1,MWIDTT2) )          ! sg
     .   MWIDTT=(MWIDTT1+MWIDTT2) )             ! PC & linux
C
      PARAMETER(
     .   MSIZET1=MDOFCT,
     .   MSIZET2=MDOFCT*MDOFCT*MWIDTT )
C
      PARAMETER(
c    .   MSIZET=max(MSIZET1,MSIZET2) )          ! sg
     .   MSIZET=(MSIZET1+MSIZET2) )             ! PC & linux
C
C**** GMRES SOLVER
C
      PARAMETER(
     .   MMKRYLT1=MDEFA1*(MEQNST/100+1),        ! default
     .   MMKRYLT2=MDEFA2*1000 )                 ! chosen variable
C
      PARAMETER(
c    .   MMKRYLT=max(MMKRYLT1,MMKRYLT2) )       ! sg
     .   MMKRYLT=(MMKRYLT1+MMKRYLT2) )          ! PC & linux
C
      PARAMETER(
     .   MWORKGT=MEQNST*(MMKRYLT+1)+(MMKRYLT*(MMKRYLT+1))/2+
     .           4*MMKRYLT+2 )
C
      PARAMETER(
     .   MCEROT=0 )
C
C***********************************************************************
C
C**** ADDITIONAL MEMORY PARAMETERS (do not change)
C
C***********************************************************************
C
      PARAMETER(
c    .   MMEMO7B=max(MMEMO7A,MSOLV5),           ! sg
c    .   MMEMO7=max(MMEMO7A,MMEMO7B) )          ! sg
     .   MMEMO7=0 )                             ! PC & linux
C
      PARAMETER(
c    .   MMEMO6B=max(MMEMO7,0),                 ! sg
c    .   MMEMO6=max(MMEMO6A,MMEMO6B) )          ! sg
     .   MMEMO6=1 )                             ! PC & linux
C
      PARAMETER(
     .   MMEMA1=1-MMEMO1,
     .   MMEMA2=1-MMEMO2,
     .   MMEMA3=1-MMEMO3,
     .   MMEMA4=1-MMEMO4,
     .   MMEMA5=1-MMEMO5,
     .   MMEMA6=1-MMEMO6,
     .   MMEMA7=1-MMEMO7,
     .   MMEMA8=1-MMEMO8,
     .   MMEMA9=1-MMEMO9,
     .   MMEMA10=1-MMEMO10,
     .   MMEMA11=1-MMEMO11 )
C
C***********************************************************************
C
C**** MATRICES DIMENSIONS (do not change)
C
C***********************************************************************
C
      PARAMETER(
     .   MEVACT=MEVABT,
     .   MKOVAT1=MUNSY1*MEVABT*(MEVABT+1)/2,    ! symmetric case
     .   MKOVAT2=MUNSY2*MEVABT*MEVABT )         ! unsymmetric case
C
      PARAMETER(
c    .   MKOVAT=max(MKOVAT1,MKOVAT2) )          ! sg
     .   MKOVAT=(MKOVAT1+MKOVAT2) )             ! PC & linux
C
      PARAMETER(
     .   MKONDT=MKOVAT,
     .   MKOSTT=(MSTR1T+1)*MSTR1T/2 )
C
