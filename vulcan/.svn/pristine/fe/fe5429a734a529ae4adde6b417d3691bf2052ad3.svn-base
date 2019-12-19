      SUBROUTINE SETDATS
C***********************************************************************
C
C**** THIS ROUTINE SETS DATA FOR THE MICROSTRUCTURAL PROBLEM
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'       ! thermal-mechanical
      INCLUDE 'nued_om.f'       ! thermal-microestructural
      INCLUDE 'nuee_om.f'       ! mechanical-microestructural
      INCLUDE 'nuef_om.f'       ! thermal-flow
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'prob_oms.f'
C
C**** ADDDATS
C
      IREL1S=0
C
C**** ADDSOLS
C
      IRELES=0
C
C**** ADDPRIS
C
C     NFPCH should be greater (or equal) than:
C     - 2*number of phase-changes (f_pc & dot f_pc) +
C     - other microst. variables to transfer to mechanical problem by
C     means of FPHCAT (see outmic.f) or to use in thermally-coupled
C     problems to evaluate the convective term associated with the
C     material time evolution of a microstructural variable (see
C     pointes.f for both cases) +
C     - 1 (only for pseudo-concentration for filling prob.; NFILL=1) +
C     - 1 (only for thermal-flow problems)
C
C     => NFPCH ge 2*NNUPT + NNUPO + NFILL + IGALFA
C
c     NFPCH=8                 ! see setdatt.f
C
C**** CONINPS
C
      KPOSTS=0
      KSGAUS=0
      NBLIMS=11
      NHOURS=0
      TLIMTS=1.0E+20
C
      KRENUS=0
C
      KSOLVS=0
      KSYMMS=1
      NWIDTS=50000
      MITCGS=0
      NBUFAS=0
      TOLCGS=0.01
      TOLC1S=1.0D-15
      NPRIRS=0
C
C**** CPUTIMS
C
      CPUINS=0.0D+00
      CPUSTS=0.0D+00
      CPUSFS=0.0D+00
      CPUASS=0.0D+00
      CPUSOS=0.0D+00
      CPURES=0.0D+00
      CPURSS=0.0D+00
      CPUOUS=0.0D+00
      CPUDAS=0.0D+00
C
C**** PROINPS
C
      KDYNAS=0
      KPORES=0
      KPROBS=5      ! species
      KSMUSS=0
      KTEMPS=0
      LARGES=0
      NHLODS=5      ! number of parameters that define the load function
      NSUBFS=5      ! number of subfunctions of a load function
      NPRELS=7
      NPROPS=1000   ! number of microstructural/species properties
C                   !   1: 500=microstructural properties; see setdatd.f
C                   ! 501:1000=species properties
      IF(IMICR.EQ.0) NPROPM=500   ! microst. prop.; see setdatd.f
C
C**** INDEXES FOR MICROSTRUCTURAL VARIABLES; see pointes.f
C
      IF(IMICR.EQ.1) THEN
       NNUM4S=NNUM4TS                         ! see setdatd.f
       NNUM4SM=NNUM4TMS
       NNUINS=NNUM4TMS
       NNUNOS=NNUM4TMS
C
       KPLA1S=0
       KPLA2S=0
       KPLA3S=0
       KPLA4S=0
       KPLA5S=0
       KPLA6S=0
       KPLA7S=0
       KPLA8S=0
       KPLA9S=0
       KPLA10S=0
       KPLA10S=0
       KPLA11S=0
       KPLA12S=0
       KPLA13S=0
       KPLA14S=0
C
       KVLA1S=0           ! not used now
       KVLA2S=0
       KVLA3S=0
       KVLA4S=0
       KVLA5S=0
C
C**** INITIALISES POINTERS (see auxl_oms.f & pointes.f)
C
       DO I=1,NNUM4SM
        IPLASS(I)=0
        IPLAOS(I)=0
        IPLANS(I)=0
        IPLAMS(I)=0
        IPLUAS(I)=0    ! simplification: dimension should be NNUPM
        IPLUOS(I)=0    ! simplification: dimension should be max(INDEX3)
        IPLLLS(I)=0    ! simplification: dimension should be NNUPM
        IPLXXS(I)=0    ! simplification: dimension should be max(INDEX3)
       ENDDO
      ENDIF
C
      NNUPCS=0         ! number macroscopic phase-changes
      NNUPTS=0         ! total number of phase-changes
C
      ICONVS=0
C
      IGALES=0
      IUPWIS=0
      IPERTS=0
      ISUPWS=0
      EXTUPS=0.0
C
C**** POROSITY
C
      NPREAS=0
      NPOROS=0
      NPRE1S=0
C
C**** ACTIVE ELEMENTS
C
      NACTIS=0
C
C**** FILLING AND/OR FLUID PROBLEM
C
      NFILLS=0
      IGALFAS=0
C
C**** ISOTROPIC OR ORTHOTROPIC MATERIAL
C
      NSOTRS=0
C
C**** STRINPS
C
      KALGOS=1
      KARCLS=0
      KINTES=0
      KOPTIS=0
      KCONVS=1
      KSAVES=0
      LACCES=0
      LAUTOS=0
      LINESS=0
      MITERS=50
      NALGOS=1
      NBACKS=0
C
      DO I=1,50          ! see inte_oms.f
       NOUTPS(I)=0
      ENDDO
C
      STIFIS=0
      TOLERS=0.01
      XTIMES=1.0
      WLUMPS=1.0
C
      NCOLDS=0
      NCURVS=0
      DO I=1,MMCURS
       NPONTS(I)=0
      ENDDO
C
      KPPCGS=0       ! no Gaussian variables to postprocess
      KPPCNS=1       ! nodal variables to postprocess
C
C**** OUTPUT OF INITIAL CONDITIONS
C
      IPRCOS=0
C
C**** STEADY STATE: TALFAT=1.0 ----- TRANSIENT: TALFAT=TALFAT
C
      TALFAS=1.0
C
C**** INPDATS
C
      GRAVYS=1.0D+00
C
C**** PLOINPS
C
      DO IPLOTS= 1,MMCURS
       DO IPLO2S= 1,2
        DO IPLO3S= 1,3
         MPLOTS(IPLOTS,IPLO2S,IPLO3S)=0
        ENDDO
       ENDDO
      ENDDO
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
C               NDISKD=0: database is out of core (datbast.f)
C               NDISKD=1: database is in core (datbast.f)
C
C               NFURES=0: a future restart will not be made (datrstt.f)
C               NFURES>0: a future restart will be made (datrstt.f)
C
C               NMACHI=1: CONVEX
C               NMACHI=2: SILICON GRAPHICS
C               NMACHI=3: VAX COMPUTER (not implemented)
C               NMACHI=4: SUN
C               NMACHI=5: PERSONAL COMPUTER
C               NMACHI=6: SILICON GRAPHICS POWER CHALLENGE
C               NMACHI=7: HEWLETT PACKARD
C               NMACHI=8: LINUX
C
C     THEORETICAL DEFAULTS: NDISKD=0, NFURES=1, NMACHI=1
C     REAL DEFAULTS:        NDISKD=1, NFURES=0, NMACHI=2
C
C***********************************************************************
C
      NDISKDS=1
      NFURESS=0
      NMACHIS=8
C
      ICHALS=4
      IF(NMACHIS.EQ.6) ICHALS=8
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
C***********************************************************************
C
      NMEMOS=0
C
C***********************************************************************
C
C**** MEMORY ADDTIONAL PARAMETERS
C
C     NMEMO1 concerns to: COORDT(NDIMET,NNODET) as a global array
C
C               NMEMO1=0: COORDT is an elemental array (ELCODT)
C               NMEMO1=1: COORDT is a global array
C
C     THEORETICAL DEFAULT: NMEMO1=0
C
C
C     NMEMO2 concerns to: shape functions, cartesian derivatives,
C                         etc. computed initially or every time
C                         when needed
C
C               NMEMO2=0: SHAPE, CARTD, etc. computed initially
C               NMEMO2=1: SHAPE, CARTD, etc. computed every time
C
C     THEORETICAL DEFAULT: NMEMO2=0
C
C
C     NMEMO3 concerns to: history-dependent material thermal properties
C
C               NMEMO3=0: no history-dependent thermal properties
C               NMEMO3=1: history-dependent material properties 
C
C     THEORETICAL DEFAULT: NMEMO3=0
C
C
C     NMEMO4 concerns to: history-dependent heat flux
C
C               NMEMO4=0: no history-dependent heat flux
C               NMEMO4=1: history-dependent heat flux
C
C     THEORETICAL DEFAULT: NMEMO4=0
C
C
C     NMEMO5 concerns to: ELDIST in/out of ELVART
C
C               NMEMO5=0: ELDIST in ELVART
C               NMEMO5=1: ELDIST out of ELVART
C
C     THEORETICAL DEFAULT: NMEMO5=0
C
C
C     NMEMO6 concerns to: elemental assembly process in other/same
C                         NELEMT loop as the evaluation of each
C                         contribution of the jacobian matrix
C
C               NMEMO6=0: elemental assembly process in other loop
C               NMEMO6=1: elemental assembly process in the same loop
C
C     THEORETICAL DEFAULT: NMEMO6=0
C
C
C     NMEMO7 concerns to: the jacobian matrix is (not) evaluated in
C                         the NELEMT loop in the solver routines. 
C
C               NMEMO7=0: jacobian not evaluated in NELEMT solver loop
C               NMEMO7=1: jacobian evaluated in NELEMT solver loop. For
C                         this case, NMEMO6 must be equals 1
C
C     With NMEMO7=1, the jacobian matrix is computed:
C
C     Skyline solver: 1) skyasst.f (kreslt=1)
C                     2) skyrhst.f (inhrst=1)
C                     3) skyitet.f => skyprot.f
C                     4) skyrhst.f (ipasst=2)
C
C     Frontal solver: 1) froasst.f (kreslt=1)
C
C     PCG solver:     1) pcgasst.f (kreslt=1)
C                     2) skyrhst.f (inhrst=1)
C                     3) pcgitet.f => skyrhst.f
C                     4) skyrhst.f (ipasst=2)
C
C     Comment: NMEMO7=1 implies KRESLT=KSTIFT=1 (see algorst.f).
C              In phase-change problems, if IITERT > NSOL1, KSTIFT=0 but
C              KRESLT=1. Therefore, NMEMO7=1 is not recommended for this
C              class of problems.
C
C     THEORETICAL DEFAULT: NMEMO7=0
C
C
C     NMEMO8 concerns to: the second derivative respect to time of the
C                         temperature
C
C               NMEMO8=0: d2T/dt2 not used
C               NMEMO8=1: d2T/dt2 used
C
C     THEORETICAL DEFAULT: NMEMO8=0
C
C
C     NMEMO9 concerns to: the second component of DISITT
C
C               NMEMO9=0: second component of DISITT not used
C               NMEMO9=1: second component of DISITT used
C
C     THEORETICAL DEFAULT: NMEMO9=0
C
C
C     NMEMO10 concerns to: the temperature-dependency of the density
C                          changes the geometry (volume and boundary)
C                          NMEMO10 is not compatible with ITERMD=1
C                          (see couinc.f)
C
C               NMEMO10=0: rho(T) does not change the geometry
C               NMEMO10=1: rho(T) changes the geometry
C
C     THEORETICAL DEFAULT: NMEMO10=0
C
C
C     NMEMO11 concerns to: normal gap & normal pressure to be considered
C                          in the heat transfer value computed in
C                          thermal or mechanical problems
C
C               NMEMO11=0: g_n & p_n computed in mechanical problem
C               NMEMO11=1: g_n & p_n computed in thermal problem
C
C     THEORETICAL DEFAULT: NMEMO11=0
C
C***********************************************************************
C
      NMEMO1S=1
      NMEMO2S=1
      NMEMO3S=0
c     IF(IMICR.GT.0) NMEMO3S=1           ! to be improved
      NMEMO4S=0
      NMEMO5S=1
      NMEMO6S=1
      NMEMO7S=0
      IF(NMEMO7S.EQ.1) NMEMO6S=1
      NMEMO8S=0
      NMEMO9S=0
      NMEMO10S=0
      NMEMO11S=1
C
C***********************************************************************
C
C**** LOGICAL UNITS (defaults)
C
C***********************************************************************
C
c     LUDTSS=201
c     LUSOLS=202
c     LUFROS=203
c     LUFRHS=204
c     LUDATS=205
c     LUPRIS=206
c     LURESS=207
c     LUSO2S=208
c     LUFR2S=209
c     LUPOSS=210          ! see psopens.f
c     LURSTS=211
c     LUBFGS=212
c     LUPIPS=213
c     LUPANS=214
c     LUFANS=215
C
c     LUGEOS=240
c     LUSETS=241
c     LUMATS=242
c     LUINIS=243
c     LULOAS=244
c     LUFIXS=245
c     LUADVS=246
C
      LUINFS=502
c     LUCU1S=151
c     LUCU2S=152
c     LUCU3S=153
c     LUCU4S=154
c     LUCU5S=155
c     LUCU6S=156
c     LUCU7S=157
c     LUCU8S=158
c     LUCU9S=159
c     LUC10S=160
C
      RETURN
      END
