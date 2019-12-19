      SUBROUTINE ELM030(LNODS,PROEL,PROPS,INFRI,COFRI,NOPRF,PRESF,VANIS,
     .                  WORK1,
     .                  ELCOD,CARTD,DVOLU,GPCOD,SHAPE,EPMTX,RMAT1,
     .                  EMASS,STIFH,BSBAR,                        !ELDAT
     .                  STRA0,STRS0,TEMPC,                        !ELPRE
     .                  ELDIS,EHIST,STRAN,STRSG,                  !ELVAR
     .                  CSTIF,ESTIF,WSTIF,HSTIF,PSTIF,QSTIF,      !ELMAT
     .                                                      ITASK)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS COMPUTATIONS AT ELEMENT LEVEL FOR
C     ELEMENT NO. 1 :
C
C     1D/2D/3D THERMOELASTO-PLASTIC ELEMENT
C
C     FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  2   Evaluate Mass Matrix
C       ITASK =  3   Evaluate Stiffness ( & Coupling ) Matrices
C       ITASK =  4   Evaluate Internal Resisting Forces
C       ITASK =  5   ------
C       ITASK =  6   Evaluate Internal Resisting Dynamic Forces
C       ITASK =  7   Evaluate equivalent Volume Forces
C       ITASK =  8   Evaluate equivalent Surface Forces
C       ITASK =  9   Check correctness of integration rule
C       ITASK = 10   ------ (Increment Non-Tensional Strains)
C       ITASK = 11   ------ (Read Initial State Variables at Gauss P.)
C       ITASK = 12   Output Gaussian Variables
C       ITASK = 13   Evaluate contribution for Nodal stresses
C       ITASK = 14   Evaluate contribution for Nodal strains
C       ITASK = 15   Evaluate contribution for Nodal internal variables
C       ITASK = 16   Evaluate variables for outpos.f
C       ITASK = 17   ------
C       ITASK = 18   Nothing (see outsmo.f)
C       ITASK = 19   ------
C
C
C
C     TREATMENT OF SMOOTHING OPERATION (IF NEEDED)
C
C     ISMO1=0 => use standard integration rule for smoothing
C          =1 => use Lobatto type integration rule for smoothing
C          =2 => idem ISMO1=0 but, for axisymmetric problems, computes
C                the 2D DVOLU (warning: this option is not compatible
C                with nodal integration rules)
C
C     ISMO2=0 => do not diagonalise mass matrix for smoothing
C          =1 => diagonalise mass matrix for smoothing
C
C     Notes:
C
C     The option ISMO1=0 and ISMO2=0 does not work well
C     The option ISMO1=1 (for any ISMO2) works well but it is not
C     absolutely correct (see rulius.f)
C     The option ISMO1=0 and ISMO2=1 works well and it is more correct
C     than the former
C     Both parameters are input in the data file (see coninp.f)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
      COMMON/FATIPROE/PFATI,PCYCL,PREVF,REVCO
      COMMON/NEGATITE/NEGATT                           ! like in thermal
C
      DIMENSION LNODS(*), PROEL(*), PROPS(*), INFRI(*), COFRI(*),
     .          NOPRF(*), PRESF(*), VANIS(*)
      DIMENSION ELCOD(*), CARTD(*), DVOLU(*), GPCOD(*), SHAPE(*), 
     .          EPMTX(*), RMAT1(*), EMASS(*), STIFH(*), BSBAR(*)
      DIMENSION STRA0(*), STRS0(*), TEMPC(*)
      DIMENSION ELDIS(*), EHIST(*), STRAN(*), STRSG(*)
      DIMENSION CSTIF(*), ESTIF(*), WSTIF(*), HSTIF(*), PSTIF(*),
     .          QSTIF(*)
      DIMENSION WORK1(*)
      DIMENSION NDATO(2,5)
C
      DATA NDATO/3,3,   ! 2D PLANE STRESS
     .           3,4,   ! 2D PLANE STRAIN
     .           4,4,   ! 2D AXISYMMETRIC
     .           6,6,   ! 3D
     .           1,1/   ! 1D
C
C**** SELECT ELEMENT PROPERTIES
C
      NNODL=INT(PROEL( 2))
      NRULE=INT(PROEL( 3))
      NGAUL=INT(PROEL( 4))
      NTYPE=INT(PROEL( 6))
      THICK=    PROEL( 7)
C
      PPART=    PROEL( 8)
      PPARB=    PROEL( 9)
      PPARI=    PROEL(10)
      ESTAB=    PROEL(11)
C
      PFATI=    PROEL(15)              ! fatigue parameters
      PCYCL=    PROEL(16)
      PREVF=    PROEL(17)
      REVCO=    PROEL(18)
C
      NCRIT=INT(PROPS(36))
C
      NSTRE=NDATO(1,NTYPE)
      NSTRS=NDATO(2,NTYPE)
      IF(KPROB.EQ.5) THEN      ! incompressibility condition
       NSTRE=1
       NSTRS=1
      ENDIF
C
      NNODS=NNODL              ! NNODS used for smoothing operations
C
      NEGATT=1                 ! always 1 except for mass matrix
C
C**** CONSIDERS DIFFERENT ELEMENTS WITH THE SAME NUMBER OF NODES
C
      NQUTR=1                  ! to be completed !!!!!!!!!!!
      IF(NDIME.EQ.2) THEN
c      IF(NNODL.EQ.3) THEN     ! to be revised
c       IF(NRULE.EQ.1.OR.NRULE.EQ.5) NQUTR=2
c      ENDIF
       IF(NNODL.EQ.4) THEN
        IF(NRULE.EQ.3.OR.NRULE.EQ.6.OR.NRULE.EQ.7) NQUTR=2
       ENDIF
      ENDIF
      IF(NDIME.EQ.3) THEN
c      IF(NNODL.EQ.4) THEN     ! to be revised
c       IF(NRULE.EQ.1.OR.NRULE.EQ.5) NQUTR=2
c      ENDIF
      ENDIF
C
C**** CONTROLS THAT NO B-BAR SHOULD BE USED FOR PLANE STRESS CASES
C
      IF(NTYPE.EQ.1) THEN
       KPARB=INT(PPARB)
c      IF(KPARB.NE.0) THEN                               ! to be revised
       IF(KPARB.EQ.1.OR.KPARB.EQ.4.OR.KPARB.EQ.7) THEN
        CALL RUNMEN('WARNING: KPARB IS SET TO 0 FOR PLANE STRESS PROB.')
        PROEL(9)=0.0D0
        PPARB=PROEL( 9)
       ENDIF
      ENDIF
C
C**** CONTROLS FOR HYPERELASTIC MATERIALS FOR WHICH B-BAR IS NOT
C     AVAILABLE   this is not true!!! to be revised!!!
C
      IPEP2=INT(PROPS(2))   ! 10=elastic; 20=plastic; 30=viscoplastic;
C                           ! 40=hyperelastic
#ifndef restricted
      IF(IPEP2.EQ.40) THEN
       KPARB=INT(PPARB)
       IF(KPARB.NE.0) THEN
c       CALL RUNMEN('WARNING: POLYMERS INCOMP. WITH KPARB NE 0; elm030')
c       PROEL(9)=0.0D0
c       PPARB=PROEL( 9)
       ENDIF
      ENDIF
#endif
C
C**** DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO ( 10, 20, 30, 40,  1, 60, 70, 80, 90,100,
     .       110,120,130,140,150,160,  1,  1,  1), ITASK
    1 RETURN ! Nothing
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
   10 CONTINUE
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(ISETM(11)),WORK1(ISETM(12)),       ! CARTD,DVOLU
     .            WORK1(ISETM(13)),WORK1(ISETM(14)),       ! GPCOD,SHAPE
     .            WORK1(ISETM(15)),                        ! BSBAR
     .            WORK1(ISETM( 1)),WORK1(ISETM( 2)),WORK1(ISETM( 3)),
     .            WORK1(ISETM( 4)),WORK1(ISETM( 5)),WORK1(ISETM( 6)),
     .            WORK1(ISETM( 7)),WORK1(ISETM( 8)),WORK1(ISETM( 9)),
     .            WORK1(ISETM(10)),WORK1(ISETM(16)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(ISETM(17)),                        ! ELDI1
     .            WORK1(ISETM(18)),WORK1(ISETM(19)),       ! CARTS,SHAPS
     .            WORK1(ISETM(20)),                        ! GPCOS
     .            DUMMY,
     .            DUMMY,DUMMY,
     .            WORK1(ISETM(32)),WORK1(ISETM(33)),       ! STREI,VARII
     .            STRS0,EHIST,
     .            ITASK)
C
      IF(INITP.EQ.1)
     . CALL INPL30(PROPS,EHIST)
      RETURN
C
C**** EVALUATES MASS MATRIX (ONLY FOR NMEMO8M=0) &
C     ASSIGNS               (ONLY FOR NMEMO6M=0)
C
   20 CONTINUE
      IF(NMEMO8M.EQ.1) RETURN
C
      NEGATT=2
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(IFORD( 6)),WORK1(IFORD( 7)),       ! CARTD,DVOLU
     .            WORK1(IFORD( 8)),WORK1(IFORD( 9)),       ! GPCOD,SHAPE
     .            WORK1(IFORD(10)),                        ! BSBAR
     .            WORK1(IFORD(11)),WORK1(IFORD(12)),WORK1(IFORD(13)),
     .            WORK1(IFORD(14)),WORK1(IFORD(15)),WORK1(IFORD(16)),
     .            WORK1(IFORD(17)),WORK1(IFORD(18)),WORK1(IFORD(19)),
     .            WORK1(IFORD(20)),WORK1(IFORD(21)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(IFORD(22)),                        ! ELDI1
     .            WORK1(IFORD(23)),WORK1(IFORD(24)),       ! CARTS,SHAPS
     .            WORK1(IFORD(25)),                        ! GPCOS
     .            DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
C
      CALL MASS01(WORK1(IFORD(7)),                         ! DVOLU
     .            PROPS,
     .            WORK1(IFORD(9)),                         ! SHAPE
     .            WORK1(IFORD(4)))                         ! WSTIF
C
      CALL ASAL00(WSTIF,WSTIF,
     .            WORK1(IFORD( 4)),                        ! ESTIF
     .            WORK1(IFORD( 4)),                        ! WSTIF
     .            NKOVA,KDYNA,NMEMO6M)
      RETURN
C
C**** EVALUATES STIFFNESS ( & COUPLING MATRICES ),
C     EVALUATES MASS MATRIX FOR DYNAMIC PROBLEMS (ONLY FOR NMEMO8M=1),
C     PERFORMES A LUMPED MATRIX (ONLY FOR NCETA=1) &
C     PERFORMS THE ASSEMBLY PROCESS (ONLY FOR NMEMO6M=1)
C
   30 CONTINUE
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(ISTIF(10)),WORK1(ISTIF(11)),       ! CARTD,DVOLU
     .            WORK1(ISTIF(12)),WORK1(ISTIF(13)),       ! GPCOD,SHAPE
     .            WORK1(ISTIF(14)),                        ! BSBAR
     .            WORK1(ISTIF(15)),WORK1(ISTIF(16)),WORK1(ISTIF(17)),
     .            WORK1(ISTIF(18)),WORK1(ISTIF(19)),WORK1(ISTIF(20)),
     .            WORK1(ISTIF(21)),WORK1(ISTIF(22)),WORK1(ISTIF(23)),
     .            WORK1(ISTIF( 9)),WORK1(ISTIF(24)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(ISTIF(27)),                        ! ELDI1
     .            WORK1(ISTIF(29)),WORK1(ISTIF(30)),       ! CARTS,SHAPS
     .            WORK1(ISTIF(31)),                        ! GPCOS
     .            WORK1(ISTIF(37)),                        ! ELCO2
     .            WORK1(ISTIF(40)),WORK1(ISTIF( 7)),       ! DISI2,DISPL
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
C
      CALL ASIL00(WORK1(ISTIF(25)),WORK1(ISTIF(26)),       ! ESTIF,WSTIF
     .            NKOVA,KDYNA,
     .            WORK1(ISTIF(34)),                        ! EMATX
     .            NEVAB,    1)
C
      CALL STIF30(WORK1(ISTIF(10)),WORK1(ISTIF(11)),       ! CARTD,DVOLU
     .            EHIST,
     .            WORK1(ISTIF( 9)),                        ! ELCOD
     .            WORK1(ISTIF(27)),                        ! ELDIS
     .            EPMTX,
     .            WORK1(ISTIF(12)),                        ! GPCOD
     .            PROPS,
     .            WORK1(ISTIF(13)),                        ! SHAPE
     .            STRSG,
     .            WORK1(ISTIF(25)),                        ! ESTIF
     .            HSTIF,
     .            WORK1(ISTIF(24)),                        ! EMASS
     .            WORK1(ISTIF( 1)),WORK1(ISTIF( 2)),
     .            WORK1(ISTIF( 3)),WORK1(ISTIF( 4)),
     .            WORK1(ISTIF( 5)),
     .            WORK1(ISTIF(14)),                        ! BSBAR
     .            WORK1(ISTIF( 6)),                        ! DTENO
     .            WORK1(ISTIF(18)),WORK1(ISTIF(34)))       ! XJACN,EMATX
C
      IF(KDYNA.EQ.1) THEN
       IF(NMEMO8M.EQ.1) THEN
        NEGATT=2
        CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .              RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .              WORK1(ISTIF(10)),WORK1(ISTIF(11)),     ! CARTD,DVOLU
     .              WORK1(ISTIF(12)),WORK1(ISTIF(13)),     ! GPCOD,SHAPE
     .              WORK1(ISTIF(14)),                      ! BSBAR
     .              WORK1(ISTIF(15)),WORK1(ISTIF(16)),WORK1(ISTIF(17)),
     .              WORK1(ISTIF(18)),WORK1(ISTIF(19)),WORK1(ISTIF(20)),
     .              WORK1(ISTIF(21)),WORK1(ISTIF(22)),WORK1(ISTIF(23)),
     .              WORK1(ISTIF( 9)),WORK1(ISTIF(24)),     ! ELCO1,EMASS
     .              BSBAR,
     .              WORK1(ISTIF(27)),                      ! ELDI1
     .              WORK1(ISTIF(29)),WORK1(ISTIF(30)),     ! CARTS,SHAPS
     .              WORK1(ISTIF(31)),                      ! GPCOS
     .              DUMMY,
     .              DUMMY,DUMMY,
     .              DUMMY,DUMMY,
     .              DUMMY,DUMMY,
     .              ITASK)
C
        CALL MASS01(WORK1(ISTIF(11)),                      ! DVOLU
     .              PROPS,
     .              WORK1(ISTIF(13)),                      ! SHAPE
     .              WORK1(ISTIF(26)))                      ! WSTIF
       ENDIF
C
       IF(NCETA.EQ.1) CALL LUMPMT(WORK1(ISTIF(26)),        ! WSTIF
     .                            NKOVA,NNODL,NDOFN,KSYMM,NEVAB)
C
      ENDIF                         ! kdyna.eq.1
C
      CALL ASAL00(ESTIF,WSTIF,
     .            WORK1(ISTIF(25)),                        ! ESTIF
     .            WORK1(ISTIF(26)),                        ! WSTIF
     .            NKOVA,KDYNA,NMEMO6M)
C
      CALL ASEL30(CSTIF,
     .            WORK1(ISTIF(25)),                        ! ESTIF
     .            WORK1(ISTIF(26)),                        ! WSTIF
     .            PROPS,
     .            WORK1(ISTIF(28)),                        ! CSTI1
     .            LNODS,
     .            WORK1(ISTIF(32)),WORK1(ISTIF(33)),       ! AUXMA,TRAMA
     .            INFRI,COFRI)
C
      RETURN
C
C**** EVALUATES INTERNAL RESISTING FORCES &
C     EVALUATES DEFORMATION-DEPENDENT SURFACE LOADS (ONLY FOR ILDSF=1)
C
   40 CONTINUE
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(IFORC(21)),WORK1(IFORC(22)),       ! CARTD,DVOLU
     .            WORK1(IFORC(23)),WORK1(IFORC(24)),       ! GPCOD,SHAPE
     .            WORK1(IFORC(25)),                        ! BSBAR
     .            WORK1(IFORC(26)),WORK1(IFORC(27)),WORK1(IFORC(28)),
     .            WORK1(IFORC(29)),WORK1(IFORC(30)),WORK1(IFORC(31)),
     .            WORK1(IFORC(32)),WORK1(IFORC(33)),WORK1(IFORC(34)),
     .            WORK1(IFORC(19)),WORK1(IFORC(35)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(IFORC(37)),                        ! ELDI1
     .            WORK1(IFORC(39)),WORK1(IFORC(40)),       ! CARTS,SHAPS
     .            WORK1(IFORC(41)),                        ! GPCOS
     .            WORK1(IFORC(43)),                        ! ELCO2
     .            WORK1(IFORC(38)),WORK1(IFORC(15)),       ! DISI2,DISPL
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
C
      CALL FRIN30(WORK1(IFORC(21)),WORK1(IFORC(22)),       ! CARTD,DVOLU
     .            EHIST,
     .            WORK1(IFORC(19)),                        ! ELCOD
     .            WORK1(IFORC(37)),                        ! ELDIS
     .            WORK1(IFORC(35)),                        ! EMASS
     .            EPMTX,
     .            WORK1(IFORC(23)),                        ! GPCOD
     .            LNODS,PROPS,RMAT1,
     .            WORK1(IFORC(24)),                        ! SHAPE
     .            STRAN,STRSG,STRA0,STRS0,THICK,VANIS,
     .            WORK1(IFORC( 1)),WORK1(IFORC( 2)),WORK1(IFORC( 3)),
     .            WORK1(IFORC( 4)),WORK1(IFORC( 5)),WORK1(IFORC( 6)),
     .            WORK1(IFORC( 7)),WORK1(IFORC( 8)),WORK1(IFORC( 9)),
     .            WORK1(IFORC(10)),WORK1(IFORC(11)),WORK1(IFORC(12)),
     .            WORK1(IFORC(14)),                        ! PWOEL
     .            WORK1(IFORC(16)),WORK1(IFORC(17)),       ! PREAL,TGAPL
     .                             WORK1(IFORC(20)),       ! TENOI
     .            WORK1(IFORC(25)),                        ! BSBAR
     .            WORK1(IFORC(29)),                        ! XJACN
     .                             WORK1(IFORC(36)))       ! FPCHL
C
      IF(NLDSF.EQ.1) THEN
       IF(ILDSF.EQ.1) THEN
        CALL LDSF01(WORK1(IFORC(43)),                      ! ELCO2
     .              LNODS,PROPS,THICK,
     .              WORK1(IFORC( 1)),WORK1(IFORC(26)),WORK1(IFORC(46)),
     .                               WORK1(IFORC(18)),NOPRF,
     .              WORK1(IFORC(27)),PRESF,           WORK1(IFORC(34)),
     .              WORK1(IFORC(24)),                 WORK1(IFORC(28)),
     .              WORK1(IFORC(29)),    1,-1.0D0)
       ENDIF                         ! ildsf.eq.1
      ENDIF                          ! nldsf.eq.1
C
      IF(IGALE.EQ.1)
     . CALL SMDG30(WORK1(IFORC(22)),                             ! DVOLU
     .             WORK1(IFORC(35)),                             ! EMASS
     .             WORK1(IFORC(24)),                             ! SHAPE
     .             STRSG,EHIST,
     .             WORK1(IFORC(48)),WORK1(IFORC(50)),      ! STDGL,STDGA
     .             WORK1(IFORC(49)),WORK1(IFORC(51)),      ! CTDGL,CTDGA
     .             WORK1(IFORC(21)),WORK1(IFORC(37)),      ! CARTD,ELDIS
     .             STRAN,
     .             WORK1(IFORC(10)),WORK1(IFORC(23)),      ! XJACM,GPCOD
     .             WORK1(IFORC(25)))                       ! BSBAR
      RETURN
C
C**** EVALUATE INTERNAL DYNAMIC FORCES
C
   60 CONTINUE
      NEGATT=2
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(IFORD( 6)),WORK1(IFORD( 7)),       ! CARTD,DVOLU
     .            WORK1(IFORD( 8)),WORK1(IFORD( 9)),       ! GPCOD,SHAPE
     .            WORK1(IFORD(10)),                        ! BSBAR
     .            WORK1(IFORD(11)),WORK1(IFORD(12)),WORK1(IFORD(13)),
     .            WORK1(IFORD(14)),WORK1(IFORD(15)),WORK1(IFORD(16)),
     .            WORK1(IFORD(17)),WORK1(IFORD(18)),WORK1(IFORD(19)),
     .            WORK1(IFORD(20)),WORK1(IFORD(21)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(IFORD(22)),                        ! ELDI1
     .            WORK1(IFORD(23)),WORK1(IFORD(24)),       ! CARTS,SHAPS
     .            WORK1(IFORD(25)),                        ! GPCOS
     .            DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
C
      CALL FRDY30(PROPS,
     .            WORK1(IFORD(1)),                         ! ELELM
     .            WORK1(IFORD(2)),                         ! ACELM
     .            WORK1(IFORD(7)),                         ! DVOLU
     .            WORK1(IFORD(9)))                         ! SHAPE
      RETURN
C
C**** EVALUATE EQUIVALENT GRAVITY LOADS
C
   70 CONTINUE
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(ILOSU(17)),WORK1(ILOSU(18)),       ! CARTD,DVOLU
     .            WORK1(ILOSU(19)),WORK1(ILOSU(20)),       ! GPCOD,SHAPE
     .            WORK1(ILOSU(21)),                        ! BSBAR
     .            WORK1(ILOSU(22)),WORK1(ILOSU(23)),WORK1(ILOSU(24)),
     .            WORK1(ILOSU(25)),WORK1(ILOSU(26)),WORK1(ILOSU(27)),
     .            WORK1(ILOSU(28)),WORK1(ILOSU(29)),WORK1(ILOSU(30)),
     .            WORK1(ILOSU(16)),WORK1(ILOSU(31)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(ILOSU(32)),                        ! ELDI1
     .            WORK1(ILOSU(33)),WORK1(ILOSU(34)),       ! CARTS,SHAPS
     .            WORK1(ILOSU(35)),                        ! GPCOS
     .            DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
C
      CALL LDGR01(WORK1(ILOSU(18)),                        ! DVOLU
     .            PROPS,
     .            WORK1(ILOSU(20)),                        ! SHAPE
     .            WORK1(ILOSU( 1)))
      RETURN
C
C**** EVALUATE EQUIVALENT SURFACE LOADS
C
   80 CONTINUE
      IF(KLDSF.EQ.0) THEN            ! deformation-independent face load
       CALL LDSF01(WORK1(ILOSU(16)),                        ! ELCOD
     .             LNODS,PROPS,THICK,
     .             WORK1(ILOSU( 1)),WORK1(ILOSU( 2)),WORK1(ILOSU( 3)),
     .                              WORK1(ILOSU( 5)),WORK1(ILOSU( 6)),
     .             WORK1(ILOSU( 7)),WORK1(ILOSU( 8)),WORK1(ILOSU( 9)),
     .             WORK1(ILOSU(10)),                 WORK1(ILOSU(12)),
     .             WORK1(ILOSU(13)),    0, 1.0D0)
      ELSE                           ! deformation-dependent face load
       CALL LDSF01(WORK1(ILOSU(16)),                        ! ELCOD
     .             LNODS,PROPS,THICK,
     .             WORK1(ILOSU( 1)),WORK1(ILOSU( 2)),WORK1(ILOSU( 3)),
     .                              WORK1(ILOSU( 5)),NOPRF,
     .             WORK1(ILOSU( 7)),PRESF,           WORK1(ILOSU( 9)),
     .             WORK1(ILOSU(10)),                 WORK1(ILOSU(12)),
     .             WORK1(ILOSU(13)),    0, 0.0D0)
      ENDIF                          ! kldsf.eq.0
      RETURN
C
C**** CHECK CORRECTNESS OF INTEGRATION RULES
C
   90 CONTINUE
      CALL CHEK01(NDIME,IELEM,NNODL,NRULE,NGAUL,NQUTR)
      RETURN
C
C**** INCREMENT NON-TENSIONAL STRAINS
C
  100 CONTINUE
      CALL INCS01(EHIST,PROPS,STRA0,TEMPC)
      RETURN
C
C**** READ INITIAL STATE VARIABLES AT GAUSS POINTS
C
  110 CONTINUE
      CALL PRGA01(PROPS,EHIST,STRA0,STRS0,STRSG)
      RETURN
C
C**** OUTPUT GAUSSIAN STATE VARIABLES
C
  120 CONTINUE
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(IGSMO(11)),WORK1(IGSMO(12)),       ! CARTD,DVOLU
     .            WORK1(IGSMO(13)),WORK1(IGSMO(14)),       ! GPCOD,SHAPE
     .            WORK1(IGSMO(15)),                        ! BSBAR
     .            WORK1(IGSMO(16)),WORK1(IGSMO(17)),WORK1(IGSMO(18)),
     .            WORK1(IGSMO(19)),WORK1(IGSMO(20)),WORK1(IGSMO(21)),
     .            WORK1(IGSMO(22)),WORK1(IGSMO(23)),WORK1(IGSMO(24)),
     .            WORK1(IGSMO(25)),WORK1(IGSMO(26)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(IGSMO(27)),                        ! ELDI1
     .            WORK1(IGSMO(28)),WORK1(IGSMO(29)),       ! CARTS,SHAPS
     .            WORK1(IGSMO(30)),                        ! GPCOS
     .            DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
C
      IF(LARGE.NE.0)
     . CALL OUTCAU(WORK1(IGSMO(11)),                       ! CARTD
     .             WORK1(IGSMO(27)),                       ! ELDIS
     .             STRSG,STRAN,
     .             WORK1(IGSMO(19)),                       ! XJACM
     .             WORK1(IGSMO(13)),WORK1(IGSMO(14)),      ! GPCOD,SHAPE
     .             WORK1(IGSMO(31)),WORK1(IGSMO(15)))      ! XJACI,BSBAR
C
      KPARI=INT(PPARI)
C
      IF(KPARI.EQ.0) THEN
       CALL OUTG01(EHIST,STRAN,STRSG)
      ELSE
       CALL OUTG30(EHIST,STRAN,STRSG,
     .             WORK1(IGSMO(12)),                             ! DVOLU
     .             KPARI)
      ENDIF
      RETURN
C
C**** COMPUTE ELEMENTAL CONTRIBUTION FOR NODAL STRESSES
C
  130 CONTINUE
      IF(ISMO1.EQ.1) THEN
       NRUAU=NRULE
       NRULE=5
      ENDIF
      CALL SMONOD(NDIME,NNODL,NNODS,NQUTR)
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(IGSMO(11)),WORK1(IGSMO(12)),       ! CARTD,DVOLU
     .            WORK1(IGSMO(13)),WORK1(IGSMO(14)),       ! GPCOD,SHAPE
     .            WORK1(IGSMO(15)),                        ! BSBAR
     .            WORK1(IGSMO(16)),WORK1(IGSMO(17)),WORK1(IGSMO(18)),
     .            WORK1(IGSMO(19)),WORK1(IGSMO(20)),WORK1(IGSMO(21)),
     .            WORK1(IGSMO(22)),WORK1(IGSMO(23)),WORK1(IGSMO(24)),
     .            WORK1(IGSMO(25)),WORK1(IGSMO(26)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(IGSMO(27)),                        ! ELDI1
     .            WORK1(IGSMO(28)),WORK1(IGSMO(29)),       ! CARTS,SHAPS
     .            WORK1(IGSMO(30)),                        ! GPCOS
     .            DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
C
      IF(LARGE.NE.0)
     . CALL OUTCAU(WORK1(IGSMO(11)),                       ! CARTD
     .             WORK1(IGSMO(27)),                       ! ELDIS
     .             STRSG,STRAN,
     .             WORK1(IGSMO(19)),                       ! XJACM
     .             WORK1(IGSMO(13)),WORK1(IGSMO(14)),      ! GPCOD,SHAPE
     .             WORK1(IGSMO(31)),WORK1(IGSMO(15)))      ! XJACI,BSBAR
C
      KPARI=INT(PPARI)
C
      IF(KPARI.EQ.0) THEN
       CALL SMOG01(WORK1(IGSMO(12)),                             ! DVOLU
     .             WORK1(IGSMO(26)),                             ! EMASS
     .             WORK1(IGSMO(14)),                             ! SHAPE
     .             STRSG,
     .             WORK1(IGSMO(5)),WORK1(IGSMO(4)))        ! STREB,STREA
      ELSE
       CALL SMOG30(WORK1(IGSMO(12)),                             ! DVOLU
     .             WORK1(IGSMO(26)),                             ! EMASS
     .             WORK1(IGSMO(14)),                             ! SHAPE
     .             STRSG,
     .             WORK1(IGSMO(5)),WORK1(IGSMO(4)))        ! STREB,STREA
      ENDIF
      IF(ISMO1.EQ.1) NRULE=NRUAU
      RETURN
C
C**** COMPUTE ELEMENTAL CONTRIBUTION FOR NODAL STRAINS
C
  140 CONTINUE
      IF(ISMO1.EQ.1) THEN
       NRUAU=NRULE
       NRULE=5
      ENDIF
      CALL SMONOD(NDIME,NNODL,NNODS,NQUTR)
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(IGSMO(11)),WORK1(IGSMO(12)),       ! CARTD,DVOLU
     .            WORK1(IGSMO(13)),WORK1(IGSMO(14)),       ! GPCOD,SHAPE
     .            WORK1(IGSMO(15)),                        ! BSBAR
     .            WORK1(IGSMO(16)),WORK1(IGSMO(17)),WORK1(IGSMO(18)),
     .            WORK1(IGSMO(19)),WORK1(IGSMO(20)),WORK1(IGSMO(21)),
     .            WORK1(IGSMO(22)),WORK1(IGSMO(23)),WORK1(IGSMO(24)),
     .            WORK1(IGSMO(25)),WORK1(IGSMO(26)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(IGSMO(27)),                        ! ELDI1
     .            WORK1(IGSMO(28)),WORK1(IGSMO(29)),       ! CARTS,SHAPS
     .            WORK1(IGSMO(30)),                        ! GPCOS
     .            DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
C
      CALL SMOF30(WORK1(IGSMO(12)),                              ! DVOLU
     .            WORK1(IGSMO(26)),                              ! EMASS
     .            WORK1(IGSMO(14)),                              ! SHAPE
     .            STRAN,
     .            WORK1(IGSMO(7)),WORK1(IGSMO(6)))         ! SFISB,SFISA
      IF(ISMO1.EQ.1) NRULE=NRUAU
      RETURN
C
C**** COMPUTE ELEMENTAL CONTRIBUTION FOR NODAL INTERNAL VARIABLES
C
  150 CONTINUE
      IF(ISMO1.EQ.1) THEN
       NRUAU=NRULE
       NRULE=5
      ENDIF
      CALL SMONOD(NDIME,NNODL,NNODS,NQUTR)
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(IGSMO(11)),WORK1(IGSMO(12)),       ! CARTD,DVOLU
     .            WORK1(IGSMO(13)),WORK1(IGSMO(14)),       ! GPCOD,SHAPE
     .            WORK1(IGSMO(15)),                        ! BSBAR
     .            WORK1(IGSMO(16)),WORK1(IGSMO(17)),WORK1(IGSMO(18)),
     .            WORK1(IGSMO(19)),WORK1(IGSMO(20)),WORK1(IGSMO(21)),
     .            WORK1(IGSMO(22)),WORK1(IGSMO(23)),WORK1(IGSMO(24)),
     .            WORK1(IGSMO(25)),WORK1(IGSMO(26)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(IGSMO(27)),                        ! ELDI1
     .            WORK1(IGSMO(28)),WORK1(IGSMO(29)),       ! CARTS,SHAPS
     .            WORK1(IGSMO(30)),                        ! GPCOS
     .            DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
C
      CALL SMOI30(WORK1(IGSMO(12)),                              ! DVOLU
     .            WORK1(IGSMO(26)),                              ! EMASS
     .            WORK1(IGSMO(14)),                              ! SHAPE
     .            EHIST,PROPS,
     .            WORK1(IGSMO(10)),WORK1(IGSMO(9)))        ! SFIPB,SFIPA
      IF(ISMO1.EQ.1) NRULE=NRUAU
      RETURN
C
C**** OUTPOST
C
  160 CONTINUE
      CALL SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .            WORK1(ISTAR( 3)),WORK1(ISTAR( 4)),       ! CARTD,DVOLU
     .            WORK1(ISTAR( 5)),WORK1(ISTAR( 6)),       ! GPCOD,SHAPE
     .            WORK1(ISTAR( 7)),                        ! BSBAR
     .            WORK1(ISTAR( 8)),WORK1(ISTAR( 9)),WORK1(ISTAR(10)),
     .            WORK1(ISTAR(11)),WORK1(ISTAR(12)),WORK1(ISTAR(13)),
     .            WORK1(ISTAR(14)),WORK1(ISTAR(15)),WORK1(ISTAR(16)),
     .            WORK1(ISTAR(17)),WORK1(ISTAR(18)),       ! ELCO1,EMASS
     .            BSBAR,
     .            WORK1(ISTAR(19)),                        ! ELDI1
     .            WORK1(ISTAR(20)),WORK1(ISTAR(21)),       ! CARTS,SHAPS
     .            WORK1(ISTAR(22)),                        ! GPCOS
     .            DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            DUMMY,DUMMY,
     .            ITASK)
      RETURN
C
      END
