      SUBROUTINE SETDAT
C***********************************************************************
C
C**** THIS ROUTINE SETS DATA 
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'       ! thermal-mechanical
      INCLUDE 'nued_om.f'       ! thermal-microestructural
      INCLUDE 'nuee_om.f'       ! mechanical-microestructural
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
C**** ADDDAT
C
      IREL1=0
C
C**** ADDSOL
C
      IRELE=0
C
C**** ASELMT
C
      NSKEW=0
C
C**** CONINP
C
      KPOST=0
      KSGAU=0
      NBLIM=11
      NHOUR=0
      TLIMT=1.0E+20
C
      KRENU=0
C
      KSOLV=0
      KSYMM=1
      NWIDT=50000
      MITCG=0
      NBUFA=0
      TOLCG=0.01
      TOLC1=1.0D-15
      NPRIR=0
      MITGM=0
      MKRYL=0
      TOLGM=1.0D-15
      IPGMR=0
C
C**** CPUTIM
C
      CPUIN=0.0D+00
      CPUST=0.0D+00
      CPUSF=0.0D+00
      CPUAS=0.0D+00
      CPUSO=0.0D+00
      CPURE=0.0D+00
      CPURS=0.0D+00
      CPUOU=0.0D+00
      CPUDA=0.0D+00
C
C**** PROINP
C
      KDYNA=0
      NDISR=1
      NDISO=1
      KPORE=0
      KPROB=1
      KSMUS=0
      KTEMP=0
      LARGE=0
      NHLOD=5       ! number of parameters that define the load function
      NSUBF=50      ! number of subfunctions of a load function
      NPREL=22
      NPROP=1000                                         ! see setdatt.f
      IF(IMICR.EQ.0.AND.IMICRM.EQ.0) NPROP=500
      NPOIC=0
      IAUGM=0
      IAUG3=0
      NNODC=0
      TOLERC=0.01
      ICONC=0
C
      NKOST=0
C
C**** NNUIN: NUMBER OF COMPONENTS OF THE NUMBER OF INTERNAL VARIABLES TO
C            PRINT; see pointe.f
C
      NNUIN=23
      NNUNO=17
C
      KPLA1=0
      KPLA2=0
      KPLA3=0
      KPLA4=0
      KPLA5=0
      KPLA6=0
      KPLA7=0
      KPLA8=0
      KPLA9=0
      KPLA10=0
      KPLA11=0
      KPLA12=0
C
      KVLA1=0       ! not used now
      KVLA2=0
      KVLA3=0
      KVLA4=0
      KVLA5=0
      KVLA6=0
      KVLA7=0
      KVLA8=0
C
C**** INITIALISES POINTERS (see auxl_om.f)
C
      DO I=1,50       ! check pointe.f
       IPLAS(I)=0
       IPLAO(I)=0
       IPLAN(I)=0
       IPLAM(I)=0
C
       IVLAS(I)=0     ! not used now
       IVLAO(I)=0
       IVLAN(I)=0
       IVLAM(I)=0
      ENDDO
C
C**** PRESSURE + NON-STANDARD INITIAL CONDITIONS
C
      NPRE1=1
      NPRE2=0
      NPRE3=0
      NPRE4=0
      NPRE5=0
      NPREA=NPRE1+NPRE2+NPRE3+NPRE4+NPRE5
C
C**** ACTIVE ELEMENTS
C
      NACTI=0
C
C**** CONTACT - NON-COINCIDENT MESH
C
      NOCOI=0
      LARGC=0
      LICOI=0                             ! linearized computation of n
      GAPPC=10.0D+10                      ! infinite value
      NSKIC=1                             ! compute n at every time step
C
C**** DEFORMATION-DEPENDENT FACE LOAD
C
      NLDSF=0                             ! only to dimension
      ILDSF=0                             ! global index
      LLDSF=0                             ! linearized computation
C
C**** ANISOTROPIC CONSTITUTIVE MODELS
C
      NANIS=0                             ! only to dimension
      IANIS=0                             ! global index
C
C**** OUTPUT OF INITIAL CONDITIONS
C
      IPRCO=0
C
C**** OPEN EXTERNAL FILES
C
      IOFIX=0
      IOLOA=0
      IOACT=0
C
      IOINI=0
C
C**** STRINP
C
      KALGO=1
      KARCL=0
      KINTE=0
      KOPTI=0
      KCONV=1
      KSAVE=0
      LACCE=0
      LAUTO=0
      LINES=0
      MITER=50
      NALGO=1
      NBACK=0
C
      DO I=1,50        ! see inte_om.f
       NOUTP(I)=0
      ENDDO
C
      STIFI=0
      TOLER=0.01
      XTIME=1.0
      WLUMP=1.0
C
      NCOLD=0
      NCURV=0
      DO I=1,MMCUR
       NPONT(I)=0
      ENDDO
C
      KPPCG=0       ! no Gaussian variables to postprocess
      KPPCN=1       ! nodal variables to postprocess
C
      NFORZ=1
      MKONT=200     ! check dimensions of VCONV in plastc*.f & viscoc*.f
      TOPLA=1.0D-06
      ALFAP=1.0D0
      MSUBP=1
      EXCTP=0.0D0
      ITEPL=0
      ICOCO=0
      NALGP=1       ! updates elasto-plastic constitutive tensor Cep
C
      TALFA=0.0D0
      TBETA=0.0D0
      TGAMA=1.0D0
C
C**** INPDAT
C
      GRAVY=1.0D+00
C
      INITP=0           ! initial condition index for internal variables
C
C**** PLOINP
C
      DO IPLOT=1,MMCUR
       DO IPLO2= 1,2
        DO IPLO3= 1,3
         MPLOT(IPLOT,IPLO2,IPLO3)=0
        ENDDO
       ENDDO
      ENDDO
C
      ISUVOAC=0           ! useful when GRAVI are used with ACTIVE_ELEM.
C
C***********************************************************************
C
C**** ADDITIONAL PARAMETERS OF PROGRAM VULCAN (mechanical problem)
C
C     it should be read in the input data !!! (see check0.f) ctm
C
C***********************************************************************
C
C**** SOLVER ADDITIONAL PARAMETERS
C
C               NDISKDM=0: database is out of core (datbas.f)
C               NDISKDM=1: database is in core (datbas.f)
C
C               NFURESM=0: a future restart will not be made (datrst.f)
C               NFURESM>0: a future restart will be made (datrst.f)
C
C               NMACHIM=1: CONVEX
C               NMACHIM=2: SILICON GRAPHICS
C               NMACHIM=3: VAX COMPUTER (not implemented)
C               NMACHIM=4: SUN
C               NMACHIM=5: PERSONAL COMPUTER
C               NMACHIM=6: SILICON GRAPHICS POWER CHALLENGE
C               NMACHIM=7: HEWLETT PACKARD
C               NMACHIM=8: LINUX
C
C     THEORETICAL DEFAULTS: NDISKDM=0, NFURESM=1, NMACHIM=1
C     REAL DEFAULTS:        NDISKDM=1, NFURESM=0, NMACHIM=2
C
C***********************************************************************
C
      NDISKDM=1
      NFURESM=0
      NMACHIM=8
C
      ICHAL=4
      IF(NMACHIM.EQ.6) ICHAL=8
C
C***********************************************************************
C
C**** MEMORY ADDTIONAL PARAMETER
C
C     NMEMO concerns to: EPTMX (elastic const. tensor)
C                        RMAT1 (anisotropic const. tensor)
C                        STRA0 (initial strains)
C                        STRS0 (initial stress)
C                        TEMPC (temperature)
C
C               NMEMO=0: partial memory
C               NMEMO=1: full memory
C
C     THEORETICAL DEFAULT: NMEMO=1
C
C     Note: STRS0 (initial stress) is used for INITV=1 (non-standard
C           initial conditions)
C
C***********************************************************************
C
      NMEMOM=0
C
C***********************************************************************
C
C**** MEMORY ADDTIONAL PARAMETERS
C
C     NMEMO1M concerns to: COORD(NDIME,NNODE) as a global array
C
C               NMEMO1M=0: COORD is an elemental array (ELCOD)
C               NMEMO1M=1: COORD is a global array
C
C     THEORETICAL DEFAULT: NMEMO1M=0
C
C
C     NMEMO2M concerns to: shape functions, cartesian derivatives,
C                          etc. computed initially or every time
C                          when needed
C
C               NMEMO2M=0: SHAPE, CARTD, etc. computed initially
C               NMEMO2M=1: SHAPE, CARTD, etc. computed every time
C
C     THEORETICAL DEFAULT: NMEMO2M=0
C
C
C     NMEMO3M concerns to: ---
C
C
C     NMEMO4M concerns to: ---
C
C
C     NMEMO5M concerns to: ELDIS in/out of ELVAR
C
C               NMEMO5M=0: ELDIS in ELVAR
C               NMEMO5M=1: ELDIS out of ELVAR
C
C     THEORETICAL DEFAULT: NMEMO5M=0
C
C
C     NMEMO6M concerns to: elemental assembly process in other/same
C                          NELEM loop as the evaluation of each
C                          contribution of the jacobian matrix
C
C               NMEMO6M=0: elemental assembly process in other loop
C               NMEMO6M=1: elemental assembly process in the same loop
C
C     THEORETICAL DEFAULT: NMEMO6=0
C
C
C     NMEMO7M concerns to: the jacobian matrix is (not) evaluated in
C                          the NELEM loop in the solver routines. 
C
C               NMEMO7M=0: jacobian not evaluated in NELEM solver loop
C               NMEMO7M=1: jacobian evaluated in NELEM solver loop. For
C                         this case, NMEMO6M must be equals 1
C
C     With NMEMO7M=1, the jacobian matrix is computed:
C
C     Skyline solver: 1) skyass.f (kresl=1)
C                     2) skyrhs.f (inhrs=1)
C                     3) skyite.f => skypro.f
C                     4) skyrhs.f (ipass=2)
C
C     Frontal solver: 1) froass.f (kresl=1)
C
C     PCG solver:     1) pcgass.f (kresl=1)
C                     2) skyrhs.f (inhrs=1)
C                     3) pcgite.f => skyrhs.f
C                     4) skyrhs.f (ipass=2)
C
C     GMRES solver:   not implemented yet !! (see gmresa.f)
C
C     Comment: NMEMO7M=1 implies KRESL=KSTIF=1 (see algors.f).
C
C     Conclusion: NMEMO7M=1 is very convenient for the frontal solver
C                 but can be most time consuming for the skyline and
C                 pcg solvers.
C
C     THEORETICAL DEFAULT: NMEMO7M=0
C
C
C     NMEMO8M concerns to: the mass matrix is computed in other/same
C                          NELEM loop as the stiffness matrix
C
C               NMEMO8M=0: mass matrix computed in other loop (NMEMO6M
C                          must be equals 0)
C               NMEMO8M=1: mass matrix computed in the same (stiffness)
C                          loop
C
C     THEORETICAL DEFAULT: NMEMO8M=0
C
C
C     NMEMO9M concerns to: the second component of DISIT
C
C               NMEMO9M=0: second component of DISIT not used
C               NMEMO9M=1: second component of DISIT used
C
C     THEORETICAL DEFAULT: NMEMO9M=0
C
C
C     NMEMO10M concerns to: ---
C
C
C     NMEMO11M concerns to: normal gap & normal pressure vectors
C                           considered in mechanical problem are
C                           necessary for printing (average values for
C                           contact elements ITYPE=4 or ITYPE=32)
C
C     => NMEMO11M is not used
C
C***********************************************************************
C
      NMEMO1M=1
      NMEMO2M=1
      NMEMO5M=1
      NMEMO6M=1
      NMEMO7M=0
      IF(NMEMO7M.EQ.1) NMEMO6M=1
      NMEMO8M=1
      IF(NMEMO8M.EQ.0) NMEMO6M=0
      NMEMO9M=0
C
C***********************************************************************
C
C**** LOGICAL UNITS (defaults)
C
C***********************************************************************
C
      LUDTS=  101
      LUSOL=  102
      LUFRO=  103
      LUFRH=  104
      LUDAT=  105
      LUPRI=  106
      LURES=  107
      LUSO2=  108
      LUFR2=  109
      LUPOS=  110
      LURST=  111
      LUBFG=  112
      LUPIP=  113
      LUPAN=  114
      LUINF=  502
      LUFAN=  139
C
      LUGEO=  140
      LUSET=  141
      LUMAT=  142
      LUINI=  143
      LULOA=  144
      LUFIX=  145
      LUIN1=  146
      LUTUN=  147
      LUCON=  148
      LUACT=  149
C
      LUCU1=  191
      LUCU2=  192
      LUCU3=  193
      LUCU4=  194
      LUCU5=  195
      LUCU6=  196
      LUCU7=  197
      LUCU8=  198
      LUCU9=  199
      LUC10=  200
C
      RETURN
      END
