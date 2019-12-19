C***********************************************************************
C
C**** GENERAL DIMENSIONS OF PROGRAM VULCAN (microscopical variables)
C
C     IF ANY CHANGE IS PRODUCED, COMPILE:
C
C       adddats.f addelms.f addpris.f addsols.f addwors.f
C       datbass.f memodts.f (comp-s)
C
C     IF "segmentation fault" is produced just starting VULCAN, check
C     maximum dimensions
C
C     Command "max" works for every MMACHI except for MMACHI=5 & 8
C     (see also para_oms.f & vul-ts.f)
C
C***********************************************************************
C
C     PROBLEM: crankshaft
C
C***********************************************************************
C
C**** WRITE MAXIMUM DIMENSIONS (0=write, 1=write & stop, 2=write more)
C
      PARAMETER(
     .   MWRITES=2 )
C
C**** EVOLUTION (MEVFI=0) OR FIELD (MEVFI=1) EQUATIONS
C
      PARAMETER(
     .   MEVFI=1 ) 
C
C**** GENERAL DIMENSIONS (change if necessary; see input data)
C
      PARAMETER(
     .   MPOINS=MEVFI*100,
     .   MELEMS=MEVFI*100,
     .   MDIMES=MEVFI*2,
     .   MNODES=MEVFI*4,
     .   MGAUSS=MEVFI*4,
     .   MGRUPS=20,
     .   MMATSS=20,
     .   MFUNCS=MEVFI*2 )
C
      PARAMETER(
     .   MDYNAS=1,             ! 0:steady-state; 1:transient
     .   MSMOMS=0,             ! smoothing (0:no; 1:yes)
     .   MRENUS=0,             ! renumbering (0:no; 1:yes)
     .   MHOURS=0,             ! hourglass control (0:no; 1:yes)
     .   MPORES=0 )            ! pore water pressure prob. (0:no; 1:yes)
C
      PARAMETER(
     .   MCONVS=1,       ! 0:no convection effects; 1:convec. effec.
     .   MGALES=0,       ! 0:no upwinding technique; 1:upwin. tech.
     .   MPOROS=0,       ! porosity criteria
     .   MACTIS=0 )      ! 0: no active elements; 1: active elements
C
C**** COUPLING (thermal-mechanical) DIMENSIONS
C     (change if necessary; see input data)
C
      PARAMETER(
     .   MTERMES=0,  !-1: pure ther. 0:no mech. coupling; 1:mech. coupl.
     .   MFPCHS=0,   ! number of phase-changes (see setdats.f)
     .   MITERCS=0 ) ! 0:standard stag. scheme; 1:improved stag. scheme
C
C**** COUPLING (thermal-microstructural) DIMENSIONS
C     (change if necessary; see input data)
C
      PARAMETER(
     .   MMICRS=0*MEVFI, ! 0:no microstructural coupling; 1:micr. coupl.
     .   MMICOS=0 )      ! 0:weak micr. coupling; 1:full micr. coupl.
C
C**** OTHER DIMENSIONS (in general, do not change; however, if any 
C                       change is produced, each dimesion is set to 
C                       its corresponding maximum in addelmt.f)
C
C     NUMBER OF MATERIAL PROPERTIES
C
C     MPROPS:        microstructural properties
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
C     NUMBER OF INTERNAL VARIABLES     to be revised !!!!!!
C
C     MHISTS:         thermal & microstructural variables
C     MHISTMS:        thermal variables
C     MHISTS-MHISTMS: microstructural variables
C
C     In a thermal analysis: MHISTS=MHISTMS
C
C     Recommended values (see below):
C     MHISTS=60 & MHISTMS=10 (thermal & microstructural analysis)
C     MHISTS=6 (thermal analysis only)
C     MHISTS=10 (thermal analysis with advective effects)
C
C     In all cases, MHISTT depends on MMEMO3 given below
C
      PARAMETER(
     . MSTR1S=MDIMES,          ! nstr1s=ndimes; see consets.f
c    . MPROPS=1000+MMICRS*500, ! see setdatt.f & setdats.f
     . MPROPS=1000,            ! see setdatt.f & setdats.f
     . MDOFCS=1,               ! to be improved; see consett.f
     . MHLODS=5,               ! see setdats.f
     . MSUBFS=5,               ! see setdats.f
     . MPRELS=7,               ! see setdats.f
     . MNUM4S=200*(1-MMICRS),  ! see setdatd.f
     . MNUM4SM=(11*MNUM4S+46)*(1-MMICRS),   ! see setdatd.f
     . MHIST1S=MMICRS*50,    ! see consets.f & setdats.f (to be revised)
     . MNUINS=MMICRS*20 )    ! int. var. to print; see setdats.f
C
C**** DEGREES OF FREEDOMS (do not change)
C
      PARAMETER(
     .   MEVABS=MNODES*MDOFCS,
     .   MTOTVS=MPOINS*MDOFCS )
C
      PARAMETER(                       ! mechanical dimensions in solids
     .   MDOFCMS=MDOFCS*MDIMES,
     .   MTOTVMS=MTOTVS*MDIMES )
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
C               MFURES=0: a future restart will be made (datrstt.f)
C               MFURES=1: a future restart will not be made (datrstt.f)
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
C     THEORETICAL DEFAULTS: MDISKD=0, MMURES=0, MMACHI=1
C     REAL DEFAULTS:        MDISKD=1, MFURES=1, MMACHI=2
C
C***********************************************************************
C
      PARAMETER(
     .   MDISKDS=1,
     .   MFURESS=1,
     .   MMACHIS=8 )
C
      PARAMETER(
c    .   MCHA1S=max((1+(MMACHIS/6)*6-MMACHIS),0),      ! 1 for MMACHI=6
c    .   MCHA2S=max((1+(MMACHIS/10)*10-MMACHIS),0),    ! 1 for MMACHI=10
c    .   MCHAAS=max(MCHA1S,MCHA2S),
c    .   MCHALS=max(4,8*MCHAAS) )
     .   MCHALS=4 )                                    ! PC & linux
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
     .   MMEMOS=0 )
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
C     (change if necessary; see input data (now in setdatt.f])
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
     .   MMEMO1S=1,
     .   MMEMO2S=1,
c    .   MMEMO3S=max(0,MMICRS),                   ! sg
     .   MMEMO3S=MMICRS,                          ! PC & linux
     .   MMEMO4S=0,
     .   MMEMO5S=1,
     .   MMEMO6AS=1,
     .   MMEMO7AS=0,
     .   MMEMO8S=0,
     .   MMEMO9S=0,
     .   MMEMO10S=1,
     .   MMEMO11S=1 )
C
C**** DEFINTION OF MHISTS
C
      PARAMETER(
     .   MHISTS=MMEMO3S*10+MHIST1S )
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
     .   MDEFA1S=0,
     .   MDEFA2S=1 )
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
     .   MSOLV1S=0,
     .   MSOLV2S=1,
     .   MSOLV3S=0,
     .   MSOLV4S=0,
     .   MSOLV5S=0 )
C
C***********************************************************************
C
C**** 3) SYMMETRIC CASE (MUNSY1=1 & MUNSY2=0) (skyline, frontal & pcg)
C        OR
C        UNSYMMETRIC CASE (MUNSY1=0 & MUNSY2=1) (skyline & frontal)
C
C        (change if necessary; see input data)
C
C***********************************************************************
C
      PARAMETER(
     .   MUNSY1S=0,
     .   MUNSY2S=1 )
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
     .   MEQNSS1=MDEFA1S*MTOTVS,                 ! default
     .   MEQNSS2=MDEFA2S*1 )                     ! chosen variable
C
      PARAMETER(
c    .   MEQNSS=max(MEQNSS1,MEQNSS2) )           ! sg
     .   MEQNSS=(MEQNSS1+MEQNSS2) )              ! PC & linux
C
      PARAMETER(
     .   MLASTS1=MDEFA1S*MEQNSS*(MEQNSS+1)/2,    ! default
     .   MLASTS2=MDEFA2S*1 )                     ! chosen variable
C
      PARAMETER(
c    .   MLASTS=max(MLASTS1,MLASTS2) )           ! sg
     .   MLASTS=(MLASTS1+MLASTS2) )              ! PC & linux
C
      PARAMETER(         ! mfitet=0: direct solver
     .   MFITES=0 )      ! mfitet=1: semidirect solver
C
C**** FRONTAL SOLVER
C
C     VARIABLES:
C     MFRONS: frontwidth
C     MBUFAS: buffer size
C
      PARAMETER(
     .   MFRONS1=MDEFA1S*MTOTVS,                 ! default
     .   MFRONS2=MDEFA2S*200 )                   ! chosen variable
C
      PARAMETER(
c    .   MFRONS=max(MFRONS1,MFRONS2) )           ! sg
     .   MFRONS=(MFRONS1+MFRONS2) )              ! PC & linux
C
      PARAMETER(
     .   MBUFAS1=MDEFA1S*(MFRONS+1),             ! default
     .   MBUFAS2=MDEFA2S*7000 )                  ! chosen variable
C
      PARAMETER(
c    .   MBUFAS=max(MBUFAS1,MBUFAS2) )           ! sg
     .   MBUFAS=(MBUFAS1+MBUFAS2) )              ! PC & linux
C
      PARAMETER(
     .   MSTIFS1=MUNSY1S*MFRONS*(MFRONS+1)/2,    ! symmetric case
     .   MSTIFS2=MUNSY2S*MFRONS*MFRONS  )        ! unsymmetric case 
C
      PARAMETER(
c    .   MSTIFS=max(MSTIFS1,MSTIFS2) )           ! sg
     .   MSTIFS=(MSTIFS1+MSTIFS2) )              ! PC & linux
C
C**** PCG SOLVER
C
      PARAMETER(
     .   MWIDTS1=0,                             ! default
     .   MWIDTS2=MDEFA2S )                      ! chosen variable
C
      PARAMETER(
c    .   MWIDTS=max(MWIDTS1,MWIDTS2) )          ! sg
     .   MWIDTS=(MWIDTS1+MWIDTS2) )             ! PC & linux
C
      PARAMETER(
     .   MSIZES1=MDOFCS,
     .   MSIZES2=MDOFCS*MDOFCS*MWIDTS )
C
      PARAMETER(
c    .   MSIZES=max(MSIZES1,MSIZES2) )          ! sg
     .   MSIZES=(MSIZES1+MSIZES2) )             ! PC & linux
C
C***********************************************************************
C
C**** ADDITIONAL MEMORY PARAMETERS (do not change)
C
C***********************************************************************
C
      PARAMETER(
c    .   MMEMO7BS=max(MMEMO7AS,MSOLV5S),         ! sg
c    .   MMEMO7S=max(MMEMO7AS,MMEMO7BS) )        ! sg
     .   MMEMO7S=0 )                             ! PC & linux
C
      PARAMETER(
c    .   MMEMO6BS=max(MMEMO7S,0),               ! sg
c    .   MMEMO6S=max(MMEMO6AS,MMEMO6BS) )       ! sg
     .   MMEMO6S=1 )                            ! PC & linux
C
      PARAMETER(
     .   MMEMA1S=1-MMEMO1S,
     .   MMEMA2S=1-MMEMO2S,
     .   MMEMA3S=1-MMEMO3S,
     .   MMEMA4S=1-MMEMO4S,
     .   MMEMA5S=1-MMEMO5S,
     .   MMEMA6S=1-MMEMO6S,
     .   MMEMA7S=1-MMEMO7S,
     .   MMEMA8S=1-MMEMO8S,
     .   MMEMA9S=1-MMEMO9S,
     .   MMEMA10S=1-MMEMO10S,
     .   MMEMA11S=1-MMEMO11S )
C
C***********************************************************************
C
C**** MATRICES DIMENSIONS (do not change)
C
C***********************************************************************
C
      PARAMETER(
     .   MEVACS=MEVABS,
     .   MKOVAS1=MUNSY1S*MEVABS*(MEVABS+1)/2,    ! symmetric case
     .   MKOVAS2=MUNSY2S*MEVABS*MEVABS )         ! unsymmetric case
C
      PARAMETER(
c    .   MKOVAS=max(MKOVAS1,MKOVAS2) )           ! sg
     .   MKOVAS=(MKOVAS1+MKOVAS2) )              ! PC & linux
C
      PARAMETER(
     .   MKONDS=MKOVAS,
     .   MKOSTS=(MSTR1S+1)*MSTR1S/2 )
C
