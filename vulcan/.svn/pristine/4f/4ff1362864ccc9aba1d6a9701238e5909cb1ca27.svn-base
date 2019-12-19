      SUBROUTINE SETDATT
C***********************************************************************
C
C**** THIS ROUTINE SETS DATA 
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'       ! thermal-mechanical
      INCLUDE 'nued_om.f'       ! thermal-microestructural
      INCLUDE 'nuee_om.f'       ! mechanical-microestructural
      INCLUDE 'nuef_om.f'       ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'prob_omt.f'
C
C**** ADDDATT
C
      IREL1T=0
C
C**** ADDSOLT
C
      IRELET=0
C
C**** ADDPRI
C
C     NFPCH should be greater (or equal) than:
C     - 2*number of (macro & micro) phase-changes (f_pc & dot f_pc) + 
C     - other microst. variables to transfer to mechanical problem by
C     means of FPHCAT (see outmic.f) or to use in thermally-coupled
C     problems to evaluate the convective term associated with the
C     material time evolution of a microstructural variable (see
C     pointes.f for both cases) +
C     - 1 (only for pseudo-concentration for filling prob.; NFILL=1) +
C     - 1 (only for thermal-flow problems)
C     - 2 (only for alternative option for (macro or micro) advective
C     problems; ICONVT=1 & IEGFPC=1
C
C     => NFPCH ge 2*NNUPT + NNUPO + NFILL + IGALFA + 2*ICONVT*IEGFPC
C
      NFPCH=8
C
C**** CONINPT
C
      KPOSTT=0
      KSGAUT=0
      NBLIMT=11
      NHOURT=0
      TLIMTT=1.0E+20
C
      KRENUT=0
C
      KSOLVT=0
      KSYMMT=1
      NWIDTT=50000
      MITCGT=0
      NBUFAT=0
      TOLCGT=0.01
      TOLC1T=1.0D-15
      NPRIRT=0
      MITGMT=0
      MKRYLT=0
      TOLGMT=1.0D-15
      IPGMRT=0
C
C**** CPUTIMT
C
      CPUINT=0.0D+00
      CPUSTT=0.0D+00
      CPUSFT=0.0D+00
      CPUAST=0.0D+00
      CPUSOT=0.0D+00
      CPURET=0.0D+00
      CPURST=0.0D+00
      CPUOUT=0.0D+00
      CPUDAT=0.0D+00
C
C**** PROINPT
C
      KDYNAT=0
      KPORET=0
      KPROBT=4      ! thermal
      KSMUST=0
      KTEMPT=0
      LARGET=0
      NHLODT=17     ! number of parameters that define the load function
C                   ! =5 for time-dependent only heats
C                   ! =12 to 17 for space-time-dependent heats
      NSUBFT=50     ! number of subfunctions of a load function
      ILDV1T=0      ! index for time-dependent only heats
      ILDV2T=0      ! index for space-time-dependent only heats
      NPRELT=11
      NPROPT=500
      IF(IMICR.EQ.1) NPROPT=NPROPT+NPROPM     ! see setdatd.f
C
C**** INDEXES FOR MICROSTRUCTURAL VARIABLES; see setdatd.f & pointes.f
C
      IF(IMICR.EQ.0) THEN
       NNUM4T=0
       NNUM4TM=0
       NNUINT=0
       NNUNOT=0
      ELSE
       NNUINT=NNUM4TMS                        ! see setdatd.f
       NNUNOT=NNUM4TMS
      ENDIF
C
      KPLA1T=0
      KPLA2T=0
      KPLA3T=0
      KPLA4T=0
      KPLA5T=0
      KPLA6T=0
      KPLA7T=0
      KPLA8T=0
      KPLA9T=0
      KPLA10T=0
C
      KVLA1T=0          ! not used now
      KVLA2T=0
      KVLA3T=0
      KVLA4T=0
      KVLA5T=0
C
C**** INITIALISES POINTERS (see auxl_omt.f)
C
      IF(NNUM4TM.GT.0) THEN
       DO I=1,NNUM4TM
        IPLAST(I)=0
        IPLAOT(I)=0
        IPLANT(I)=0
        IPLAMT(I)=0
        IPLUAT(I)=0
        IPLUOT(I)=0
        IPLLLT(I)=0
        IPLXXT(I)=0       ! auxiliar array
        APLUOT(I)=0       ! auxiliar array for advective problems
       ENDDO
      ENDIF
C
      NNUPC=0             ! number macroscopic phase-changes
      NNUPT=0             ! total number of phase-changes
C
      ICONVT=0
C
C**** PETROV-GALERKIN STABILIZATION
C
C     Recommended i(standard) values:
C     Streamline dissipation: IUPWI=1, IPERT=1, ISUPW=1, EXTUP=1.0
C     Crosswind dissipation: IUPWIC=1, RAMON=0.7, IPERTC=1, ISUPWC=1,
C                            EXTUPC=1.0
C     Temporal dissipation: IUPWIG=1, IPERTG=1, ISUPWG=1, EXTUPG=1.0
C
      IGALET=0
      IUPWI=0            ! default values
      IPERT=0
      ISUPW=0
      EXTUP=0.0D0
C
      IUPWIC=0
      RAMON=0.0D0
      IPERTC=0
      ISUPWC=0
      EXTUPC=0.0D0
C
      IUPWIG=0
      IPERTG=0
      ISUPWG=0
      EXTUPG=0.0D0
C
      IGALFA=0
      GALFA=0.0D0
C
C**** POROSITY
C
      NPREAT=0
      NPOROT=0
      NPRE1T=0
C
C**** ACTIVE ELEMENTS
C
      NACTIT=0
C
C**** CONTACT - NON-COINCIDENT MESH
C
      NOCOIT=0
C
C**** FILLING PROBLEM
C
      NFILL=0
C
C**** ISOTROPIC OR ORTHOTROPIC MATERIAL
C
      NSOTRT=0
C
C**** OUTPUT OF INITIAL CONDITIONS
C
      IPRCOT=0
C
C**** OPEN EXTERNAL FILES
C
      IOFIXT=0
      IOLOAT=0
      IOADVT=0
      IOACTT=0
      IOSTRT=0
C
      IOINIT=0
C
C**** STRINPT
C
      KALGOT=1
      KARCLT=0
      KINTET=0
      KOPTIT=0
      KCONVT=1
      KSAVET=0
      LACCET=0
      LAUTOT=0
      LINEST=0
      MITERT=50
      NALGOT=1
      NBACKT=0
C
      DO I=1,50          ! see inte_omt.f
       NOUTPT(I)=0
      ENDDO
C
      STIFIT=0
      TOLERT=0.01
      XTIMET=1.0
      WLUMPT=1.0
C
      NCOLDT=0
      NCURVT=0
      DO I=1,MMCURT
       NPONTT(I)=0
      ENDDO
C
      KPPCGT=0       ! no Gaussian variables to postprocess
      KPPCNT=1       ! nodal variables to postprocess
C
C**** STEADY STATE: TALFAT=1.0 ----- TRANSIENT: TALFAT=TALFAT
C
      TALFAT=1.0D0
      TBETAT=0.0D0               ! not used now !!
C
C**** INPDATT
C
      GRAVYT=1.0D+00
C
C**** PLOINPT
C
      DO IPLOTT= 1,MMCURT
       DO IPLO2T= 1,2
        DO IPLO3T= 1,3
         MPLOTT(IPLOTT,IPLO2T,IPLO3T)=0
        ENDDO
       ENDDO
      ENDDO
C
      ISUVOACT=0 ! useful when SURFA or VOLUM are used with ACTIVE_ELEM.
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
      NDISKD=1
      NFURES=0
      NMACHI=8
C
      ICHALT=4
      IF(NMACHI.EQ.6) ICHALT=8
C
      IWINDT=0
      IF(NMACHI.EQ.5) IWINDT=1  ! needed to dimension not used variables
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
      NMEMO=0
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
C     GMRES solver:   not implemented yet !! (see gmresat.f)
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
      NMEMO1=1
      NMEMO2=1
      NMEMO3=1
      IF(IMICR.EQ.0) NMEMO3=0
      NMEMO4=1
      NMEMO5=1
      NMEMO6=1
      NMEMO7=0
      IF(NMEMO7.EQ.1) NMEMO6=1
      NMEMO8=0
      if(nmemo8.eq.1)                      ! see addprit.f
     . call runendt('error: DISTOT(NTOTVT,NCOMPT) must be implemented')
      NMEMO9=0
      NMEMO10=1
      NMEMO11=1
C
C***********************************************************************
C
C**** LOGICAL UNITS (defaults)
C
C***********************************************************************
C
      LUDTST=201
      LUSOLT=202
      LUFROT=203
      LUFRHT=204
      LUDATT=205
      LUPRIT=206
      LUREST=207
      LUSO2T=208
      LUFR2T=209
      LUPOST=210
      LURSTT=211
      LUBFGT=212
      LUPIPT=213
      LUPANT=214
      LUFANT=215
C
      LUGEOT=240
      LUSETT=241
      LUMATT=242
      LUINIT=243
      LULOAT=244
      LUFIXT=245
      LUADVT=246
      LUACTT=247
      LUSTRT=248
C
      LUINFT=502

      LUCU1T=151
      LUCU2T=152
      LUCU3T=153
      LUCU4T=154
      LUCU5T=155
      LUCU6T=156
      LUCU7T=157
      LUCU8T=158
      LUCU9T=159
      LUC10T=160
C
      RETURN
      END
