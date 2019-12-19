      SUBROUTINE ELM005S(LNODSS,PROELS,PROPSS,WORK1S,
     .                   ELCODS,CARTDS,DVOLUS,GPCODS,SHAPES,EPMTXS,
     .                   RMAT1S,EMASSS,                         !ELDAT
     .                   STRA0S,STRS0S,TEMPCS,                  !ELPRE
     .                   ELDISS,EHISTS,STRANS,STRSGS,           !ELVAR
     .                   CSTIFS,ESTIFS,WSTIFS,HSTIFS,PSTIFS,
     .                   QSTIFS,                                !ELMAT
     .                                               ITASKS)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS COMPUTATIONS AT ELEMENT LEVEL FOR
C     ELEMENT NO. 5 :
C
C     1D/2D/3D THERMAL ELEMENT
C
C     FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  2   -------
C       ITASK =  3   Evaluate : Conductivity Matrix
C                               Heat Capacity Matrix (only for
C                               Transient Problems)
C       ITASK =  4   Evaluate Internal Resisting Heat
C       ITASK =  5   -------
C       ITASK =  6   Evaluate : Internal Dynamic Resisting Heat
C                               (only for Transient Problems)
C                               Internal Dynamic Phase Change Heat
C       ITASK =  7   Evaluate equivalent Volume Heat (Internal Heat)
C       ITASK =  8   ------- (Evaluate equivalent Surface Forces)
C       ITASK =  9   Check correctness of integration rule
C       ITASK = 10   ------- (Increment Non-Tensional Strains)
C       ITASK = 11   ------- (Read Initial State Variables at Gauss P.)
C       ITASK = 12   Output Gaussian variables
C       ITASK = 13   Evaluate contribution for nodal heat fluxes
C       ITASK = 14   ------- (Evaluate contribution for nodal strains)
C       ITASK = 15   Evaluate contribution for nodal internal variables
C       ITASK = 16   Evaluate variables for outpost.f
C       ITASK = 17   Evaluate contribution for nodal porosity criteria
C       ITASK = 18   Nothing (see outsmot.f)
C       ITASK = 19   -------
C
C
C
C     TREATMENT OF SMOOTHING OPERATION (IF NEEDED)
C
C     ISMO1S=0 => use standard integration rule for smoothing
C           =1 => use Lobatto type integration rule for smoothing
C
C     ISMO2S=0 => do not diagonalise mass matrix for smoothing
C           =1 => diagonalise mass matrix for smoothing
C           =2 => idem ISMO1S=0 but, for axisymmetric problems, computes
C                 the 2D DVOLU (warning: this option is not compatible
C                 with nodal integration rules)
C
C
C     Notes:
C
C     The option ISMO1S=0 and ISMO2S=0 does not work well
C     The option ISMO1S=1 (for any ISMO2S) works well but it is not
C     absolutely correct (see rulius.f)
C     The option ISMO1S=0 and ISMO2S=1 works well and it is more correct
C     than the former
C     Both parameters are input in the data file (see coninps.f)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      COMMON/THICKNESSS/THICKS
      COMMON/MULPHASES/MULPHS
      COMMON/NEGATITES/NEGATS
C
      DIMENSION LNODSS(*), PROELS(*), PROPSS(*)
      DIMENSION ELCODS(*), CARTDS(*), DVOLUS(*), GPCODS(*), SHAPES(*), 
     .          EPMTXS(*), RMAT1S(*), EMASSS(*)
      DIMENSION STRA0S(*), STRS0S(*), TEMPCS(*)
      DIMENSION ELDISS(*), EHISTS(*), STRANS(*), STRSGS(*)
      DIMENSION CSTIFS(*), ESTIFS(*), WSTIFS(*), HSTIFS(*), PSTIFS(*),
     .          QSTIFS(*)
      DIMENSION WORK1S(*)
C
C**** SELECTS ELEMENT PROPERTIES
C
      NNODLS=INT(PROELS( 2))
      NRULES=INT(PROELS( 3))
      NGAULS=INT(PROELS( 4))
      NTYPES=INT(PROELS( 6))
      THICKS=    PROELS( 7)
C
      NNODSS=NNODLS                ! NNODS used for smoothing operations
C
C**** SELECTS TREATMENT OF SOURCE TERM (IF IT EXISTS)
C
      ILDGRS=1                     ! better as input
C
C**** CONSIDERS DIFFERENT ELEMENTS WITH THE SAME NUMBER OF NODES
C
      NQUTRS=1                     ! to be completed !!!!!! not used yet
      IF(NDIMES.EQ.2) THEN
       IF(NNODLS.EQ.4) THEN
        IF(NRULES.EQ.3.OR.NRULES.EQ.6.OR.NRULES.EQ.7) NQUTRS=2
       ENDIF
      ENDIF
c     IF(NDIMES.EQ.3) THEN
c      CALL RUNENDS('ELM005S: 3D NOT IMPLEMENTED YET')
c     ENDIF
C
C**** DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO ( 10, 20, 30, 40,  1, 60, 70, 80, 90,100,
     .       110,120,130,140,150,160,170,  1,  1), ITASKS
    1 RETURN ! Nothing
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
   10 CONTINUE
      NEGATS=1
c     IF(IMICR.EQ.1) NEGATS=2
      CALL SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,PROPSS,
     .             SHAPES,THICKS,ELDISS,
     .             WORK1S(ISETMS( 7)),WORK1S(ISETMS( 8)),  ! cartd,dvolu
     .             WORK1S(ISETMS( 9)),WORK1S(ISETMS(10)),  ! gpcod,shape
     .             WORK1S(ISETMS( 1)),WORK1S(ISETMS( 2)),  ! deriv,posgp
     .             WORK1S(ISETMS( 3)),WORK1S(ISETMS( 4)),  ! weigp,xjacm
     .             WORK1S(ISETMS( 6)),WORK1S(ISETMS(11)),  ! elco1,emas1
     .                                WORK1S(ISETMS(12)),  ! eldi1
     .             WORK1S(ISETMS(16)),WORK1S(ISETMS(17)),  ! dvoli,elcoi
     .             WORK1S(ISETMS(18)),                     ! dispt
     .             0,ITASKS)
C
c     IF(IMICR.EQ.1)
c    . CALL INMI05T(PROPST,EHISTT)
C
c     CALL FPCN05T(WORK1T(ISETMT(13)),                     ! eldis=teini
c    .             WORK1T(ISETMT(12)),                     ! not used
c    .             WORK1T(ISETMT(14)),     0)              ! FPCHL
      RETURN
C
C**** NOTHING
C
   20 CONTINUE
      RETURN
C
C**** EVALUATES CONDUCTIVITY MATRIX. IN TRANSIENT ANALYSIS, EVALUATES:
C     1) HEAT CAPACITY MATRIX & 2) "PHASE-CHANGE" MATRIX
C
   30 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATS=1
c     IF(IMICR.EQ.1) NEGATS=2
      CALL SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,PROPSS,
     .             SHAPES,THICKS,ELDISS,
     .             WORK1S(ISTIFS(26)),WORK1S(ISTIFS(27)),  ! cartd,dvolu
     .             WORK1S(ISTIFS(28)),WORK1S(ISTIFS(29)),  ! gpcod,shape
     .             WORK1S(ISTIFS(30)),WORK1S(ISTIFS(31)),  ! deriv,posgp
     .             WORK1S(ISTIFS(32)),WORK1S(ISTIFS(33)),  ! weigp,xjacm
     .             WORK1S(ISTIFS(18)),WORK1S(ISTIFS(34)),  ! elco1,emas1
     .                                WORK1S(ISTIFS(39)),  ! eldi1
     .             WORK1S(ISTIFS(46)),WORK1S(ISTIFS(47)),  ! dvoli,elcoi
     .             WORK1S(ISTIFS(44)),                     ! dispt
     .             0,ITASKS)
C
C**** EVALUATES WEIGHTING FUNCTIONS & ITS CARTESIAN DERIVATIVES
C
      CALL CART05S(WORK1S(ISTIFS(18)),                     ! ELCOD
     .             PROPSS,
     .             WORK1S(ISTIFS(39)),                     ! ELDIS
     .             WORK1S(ISTIFS(29)),WORK1S(ISTIFS(26)),  ! SHAPE,CARTD
     .             WORK1S(ISTIFS( 7)),WORK1S(ISTIFS( 8)),  ! posgp,weigp
     .             WORK1S(ISTIFS(12)),WORK1S(ISTIFS( 5)),  ! gpcdd,shapd
     .             WORK1S(ISTIFS( 9)),WORK1S(ISTIFS(10)),  ! derid,xjacm
     .             WORK1S(ISTIFS(11)),                     ! cardd
     .             WORK1S(ISTIFS(13)),                     ! velcm
     .             WORK1S(ISTIFS(19)),WORK1S(ISTIFS(20)),  ! whape,hache
     .             WORK1S(ISTIFS(21)),WORK1S(ISTIFS(22)),  ! wartd,weriv
     .                                WORK1S(ISTIFS(24)),  ! centr
     .             WORK1S(ISTIFS(25)),                     ! advem
     .             WORK1S(ISTIFS(40)),WORK1S(ISTIFS(45)),  ! teini,fpchl
     .             EHISTS)
C
      CALL ASIL00(WORK1S(ISTIFS(37)),WORK1S(ISTIFS(38)),   ! ESTIF,WSTIF
     .            NKOVAS,KDYNAS,
     .            WORK1S(ISTIFS(49)),                      ! EMATX
     .            NEVABS,     1)
C
C**** EVALUATES CONDUCTIVITY MATRIX
C
      CALL STIF05S(WORK1S(ISTIFS(26)),WORK1S(ISTIFS(27)),  ! CARTD,DVOLU
     .             EHISTS,
     .             WORK1S(ISTIFS(39)),                     ! ELDIS
     .             EPMTXS,
     .             WORK1S(ISTIFS(28)),                     ! GPCOD
     .             PROPSS,
     .             WORK1S(ISTIFS(29)),                     ! SHAPE
     .             STRSGS,
     .             WORK1S(ISTIFS(37)),                     ! ESTIF
     .             HSTIFS,
     .             WORK1S(ISTIFS(18)),                     ! ELCOD
     .             WORK1S(ISTIFS( 1)),WORK1S(ISTIFS( 2)),
     .             WORK1S(ISTIFS( 3)),WORK1S(ISTIFS( 4)),
     .             WORK1S(ISTIFS(13)),                     ! velcm
     .             WORK1S(ISTIFS(21)),                     ! wartd
     .             WORK1S(ISTIFS(40)),WORK1S(ISTIFS(45)),  ! teini,fpchl
     .             WORK1S(ISTIFS(46)),WORK1S(ISTIFS(49)))  ! dvoli,EMATX
C
      IF(KDYNAS.EQ.1) THEN
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
       NEGATS=2
       CALL SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,PROPSS,
     .              SHAPES,THICKS,ELDISS,
     .              WORK1S(ISTIFS(26)),WORK1S(ISTIFS(27)), ! cartd,dvolu
     .              WORK1S(ISTIFS(28)),WORK1S(ISTIFS(29)), ! gpcod,shape
     .              WORK1S(ISTIFS(30)),WORK1S(ISTIFS(31)), ! deriv,posgp
     .              WORK1S(ISTIFS(32)),WORK1S(ISTIFS(33)), ! weigp,xjacm
     .              WORK1S(ISTIFS(18)),WORK1S(ISTIFS(34)), ! elco1,emas1
     .                                 WORK1S(ISTIFS(39)), ! eldi1
     .              WORK1S(ISTIFS(46)),WORK1S(ISTIFS(47)), ! dvoli,elcoi
     .              WORK1S(ISTIFS(44)),                    ! dispt
     .              1,ITASKS)
C
C**** EVALUATES WEIGHTING FUNCTIONS & ITS CARTESIAN DERIVATIVES
C
       CALL CART05S(WORK1S(ISTIFS(18)),                    ! ELCOD
     .              PROPSS,
     .              WORK1S(ISTIFS(39)),                    ! ELDIS
     .              WORK1S(ISTIFS(29)),WORK1S(ISTIFS(26)), ! SHAPE,CARTD
     .              WORK1S(ISTIFS( 7)),WORK1S(ISTIFS( 8)), ! posgp,weigp
     .              WORK1S(ISTIFS(12)),WORK1S(ISTIFS( 5)), ! gpcdd,shapd
     .              WORK1S(ISTIFS( 9)),WORK1S(ISTIFS(10)), ! derid,xjacm
     .              WORK1S(ISTIFS(11)),                    ! cardd
     .              WORK1S(ISTIFS(13)),                    ! velcm
     .              WORK1S(ISTIFS(19)),WORK1S(ISTIFS(20)), ! whape,hache
     .              WORK1S(ISTIFS(21)),WORK1S(ISTIFS(22)), ! wartd,weriv
     .                                 WORK1S(ISTIFS(24)), ! centr
     .              WORK1S(ISTIFS(25)),                    ! advem
     .              WORK1S(ISTIFS(40)),WORK1S(ISTIFS(45)), ! teini,fpchl
     .              EHISTS)
C
       CALL MASS05S(WORK1S(ISTIFS(27)),                    ! DVOLU
     .              PROPSS,
     .              WORK1S(ISTIFS(29)),                    ! SHAPE
     .              WORK1S(ISTIFS(38)),                    ! WSTIF
     .              EHISTS,
     .              WORK1S(ISTIFS(19)),                    ! whape
     .              WORK1S(ISTIFS(39)),                    ! ELDIS
     .              WORK1S(ISTIFS(13)),                    ! velcm
     .              WORK1S(ISTIFS(40)),WORK1S(ISTIFS(45)), ! teini,fpchl
     .              WORK1S(ISTIFS(46)))                    ! dvoli
C
cc     IF(IITERT.NE.1) THEN
c       CALL CERO00T(WORK1T(ISTIFT(15)),WORK1T(ISTIFT(14)),NKOVAT)
c       CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
c    .               WORK1T(ISTIFT(39)),                   ! ELDIS
c    .               PROPST,WORK1T(ISTIFT(13)),1,TEMTLT,
c    .               TEMDLT)
C
c       IF(MULPHT.EQ.1) THEN
c        CALL MAUF05T(WORK1T(ISTIFT(27)),                  ! DVOLU
c    .                PROPST,
c    .                WORK1T(ISTIFT(29)),                  ! SHAPE
c    .                WORK1T(ISTIFT(38)),                  ! WSTIF
c    .                EHISTT,
c    .                WORK1T(ISTIFT(13)),WORK1T(ISTIFT(14)),
c    .                WORK1T(ISTIFT(15)),
c    .                WORK1T(ISTIFT(39)),                  ! ELDIS
c    .                WORK1T(ISTIFT(17)),                  ! disit
c    .                WORK1T(ISTIFT(19)),                  ! whape
c    .                WORK1T(ISTIFT(40)),                  ! teini
c    .                WORK1T(ISTIFT(45)),                  ! fpchl
c    .                WORK1T(ISTIFT(46)))                  ! dvoli
c       ELSE
c        CALL MAMF05T(WORK1T(ISTIFT(27)),                  ! DVOLU
c    .                PROPST,
c    .                WORK1T(ISTIFT(29)),                  ! SHAPE
c    .                WORK1T(ISTIFT(38)),                  ! WSTIF
c    .                EHISTT,
c    .                WORK1T(ISTIFT(5)),WORK1T(ISTIFT(6)),
c    .                WORK1T(ISTIFT(7)),
c    .                WORK1T(ISTIFT(39)),                  ! ELDIS
c    .                WORK1T(ISTIFT(18)),                  ! ELCOD
c    .                WORK1T(ISTIFT(8)),WORK1T(ISTIFT(9)),
c    .                WORK1T(ISTIFT(10)),WORK1T(ISTIFT(11)),
c    .                WORK1T(ISTIFT(12)),WORK1T(ISTIFT(13)),
c    .                WORK1T(ISTIFT(14)),WORK1T(ISTIFT(15)),
c    .                LNODST,
c    .                WORK1T(ISTIFT(17)),
c    .               WORK1T(ISTIFT(20)),WORK1T(ISTIFT(23)), !hache,whade
c    .               WORK1T(ISTIFT(24)),WORK1T(ISTIFT(25)), !centr,advem
c    .               WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)), !teini,fpchl
c    .               WORK1T(ISTIFT(46)))                    !dvoli
c       ENDIF            ! mulpht.eq.1
C
C**** CHECK THE CORRECTNESS OF THE PHASE CHANGE MATRIX AND ADDS ITS 
C     CONTRIBUTION TO THE "JACOBIAN MATRIX"
C
c       CALL VERIPCT(WORK1T(ISTIFT(15)),WORK1T(ISTIFT(14)),NKOVAT,
c    .               IVEPCT)
c       IF(IVEPCT.NE.0) WRITE(LUREST,2000) IVEPCT
C
c2000   FORMAT(' *** JACOBIAN MATRIX WARNING *** IVEPCT ELEMENTS OF',
c    1         ' PHASE-CHANGE MATRIX ARE NEGATIVE',I5)
C
c       CALL SUMMATT(WORK1T(ISTIFT(38)),                  ! WSTIF
c    .               WORK1T(ISTIFT(15)),WORK1T(ISTIFT(14)),
c    .               NKOVAT)
C
cc     ENDIF                        ! iitert.ne.1
C
C**** PERFORM A LUMPED HEAT CAPACITY + PHASE CHANGE MATRICES
C
c      IF(NCETAT.EQ.1) CALL LUMPMTT(WORK1T(ISTIFT(38)),          ! WSTIF
c    .                              NKOVAT,NNODLT,KSYMMT,NEVABT)
C
      ENDIF                         ! kdynas.eq.1
C
C**** EVALUATES JACOBIAN MATRIX DUE TO CONVECTION (ADVECTIVE) EFFECTS
C
      IF(ICONVS.EQ.1) THEN
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
       NEGATS=1
c      IF(IMICR.EQ.1) NEGATS=2
       CALL SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,PROPSS,
     .              SHAPES,THICKS,ELDISS,
     .              WORK1S(ISTIFS(26)),WORK1S(ISTIFS(27)), ! cartd,dvolu
     .              WORK1S(ISTIFS(28)),WORK1S(ISTIFS(29)), ! gpcod,shape
     .              WORK1S(ISTIFS(30)),WORK1S(ISTIFS(31)), ! deriv,posgp
     .              WORK1S(ISTIFS(32)),WORK1S(ISTIFS(33)), ! weigp,xjacm
     .              WORK1S(ISTIFS(18)),WORK1S(ISTIFS(34)), ! elco1,emas1
     .                                 WORK1S(ISTIFS(39)), ! eldi1
     .              WORK1S(ISTIFS(46)),WORK1S(ISTIFS(47)), ! dvoli,elcoi
     .              WORK1S(ISTIFS(44)),                    ! dispt
     .              1,ITASKS)
C
C**** EVALUATES WEIGHTING FUNCTIONS & ITS CARTESIAN DERIVATIVES
C
       CALL CART05S(WORK1S(ISTIFS(18)),                    ! ELCOD
     .              PROPSS,
     .              WORK1S(ISTIFS(39)),                    ! ELDIS
     .              WORK1S(ISTIFS(29)),WORK1S(ISTIFS(26)), ! SHAPE,CARTD
     .              WORK1S(ISTIFS( 7)),WORK1S(ISTIFS( 8)), ! posgp,weigp
     .              WORK1S(ISTIFS(12)),WORK1S(ISTIFS( 5)), ! gpcdd,shapd
     .              WORK1S(ISTIFS( 9)),WORK1S(ISTIFS(10)), ! derid,xjacm
     .              WORK1S(ISTIFS(11)),                    ! cardd
     .              WORK1S(ISTIFS(13)),                    ! velcm
     .              WORK1S(ISTIFS(19)),WORK1S(ISTIFS(20)), ! whape,hache
     .              WORK1S(ISTIFS(21)),WORK1S(ISTIFS(22)), ! wartd,weriv
     .                                 WORK1S(ISTIFS(24)), ! centr
     .              WORK1S(ISTIFS(25)),                    ! advem
     .              WORK1S(ISTIFS(40)),WORK1S(ISTIFS(45)), ! teini,fpchl
     .              EHISTS)
C
       CALL STIC05S(WORK1S(ISTIFS(26)),WORK1S(ISTIFS(27)), ! CARTD,DVOLU
     .              EHISTS,
     .              WORK1S(ISTIFS(39)),                    ! ELDIS
     .              EPMTXS,
     .              WORK1S(ISTIFS(28)),                    ! GPCOD
     .              PROPSS,
     .              WORK1S(ISTIFS(29)),                    ! SHAPE
     .              STRSGS,
     .              WORK1S(ISTIFS(37)),                    ! ESTIF
     .              HSTIFS,
     .              WORK1S(ISTIFS( 1)),WORK1S(ISTIFS( 2)),
     .              WORK1S(ISTIFS( 3)),WORK1S(ISTIFS( 4)),
     .              WORK1S(ISTIFS(19)),                    ! whape
     .              WORK1S(ISTIFS(25)),                    ! advem
     .              WORK1S(ISTIFS(13)),                    ! velcm
     .              WORK1S(ISTIFS(40)),WORK1S(ISTIFS(45))) ! teini,fpchl
C
C**** EVALUATES JACOBIAN MATRIX DUE TO CONVECTION PHASE-CHANGE EFFECTS
C
c      CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
c    .              WORK1T(ISTIFT(39)),                    ! ELDIS
c    .              PROPST,WORK1T(ISTIFT(13)),       1,TEMTLT,
c    .              TEMDLT)
C
c      MULPHT=MULPDT
c      IF(MULPHT.EQ.1) THEN
c       CALL STUC05T(WORK1T(ISTIFT(26)),WORK1T(ISTIFT(27)), !CARTD,DVOLU
c    .               EHISTT,
c    .               WORK1T(ISTIFT(39)),                    !ELDIS
c    .               EPMTXT,
c    .               WORK1T(ISTIFT(28)),                    !GPCOD
c    .               PROPST,
c    .               WORK1T(ISTIFT(29)),                    !SHAPE
c    .               STRSGT,
c    .               WORK1T(ISTIFT(37)),                    !ESTIF
c    .               HSTIFT,
c    .               WORK1T(ISTIFT( 1)),WORK1T(ISTIFT( 2)),
c    .               WORK1T(ISTIFT( 3)),WORK1T(ISTIFT( 4)),
c    .               WORK1T(ISTIFT(13)),WORK1T(ISTIFT(14)), !velcm,wsti1
c    .               WORK1T(ISTIFT(15)),                    !wsti2
c    .               WORK1T(ISTIFT(19)),                    !whape
c    .               WORK1T(ISTIFT(25)),                    !advem
c    .               WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45))) !teini,fpchl
c      ELSE
c       CALL STMC05T(WORK1T(ISTIFT(26)),WORK1T(ISTIFT(27)), !CARTD,DVOLU
c    .               EHISTT,
c    .               WORK1T(ISTIFT(39)),                    !ELDIS
c    .               EPMTXT,WORK1T(ISTIFT(28)),             !GPCOD
c    .               PROPST,WORK1T(ISTIFT(29)),             !SHAPE
c    .               STRSGT,WORK1T(ISTIFT(37)),             !ESTIF
c    .               HSTIFT,
c    .               WORK1T(ISTIFT(1)), WORK1T(ISTIFT(2)),
c    .               WORK1T(ISTIFT(3)), WORK1T(ISTIFT(4)),
c    .               WORK1T(ISTIFT(5)), WORK1T(ISTIFT(6)),  !shapd,dvold
c    .               WORK1T(ISTIFT(7)), WORK1T(ISTIFT(8)),  !posgp,weigp
c    .               WORK1T(ISTIFT(9)),                     !derid
c    .               WORK1T(ISTIFT(11)),WORK1T(ISTIFT(12)), !cardd,gpcdd
c    .               WORK1T(ISTIFT(13)),WORK1T(ISTIFT(14)), !velcm,wsti1
c    .               WORK1T(ISTIFT(15)),WORK1T(ISTIFT(18)), !wsti2,ELCOD
c    .               WORK1T(ISTIFT(20)),WORK1T(ISTIFT(23)), !hache,whade
c    .               WORK1T(ISTIFT(24)),WORK1T(ISTIFT(25)), !centr,advem
c    .               WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)), !teini,fpchl
c    .               LNODST)
c      ENDIF                        ! mulpdt.eq.1
C
      ENDIF                         ! iconvs.eq.1
C
C**** ASSEMBLY PROCESS (NMEMO6=1)
C
      CALL ASAL00(ESTIFS,WSTIFS,
     .            WORK1S(ISTIFS(37)),                       ! ESTIF
     .            WORK1S(ISTIFS(38)),                       ! WSTIF
     .            NKOVAS,KDYNAS,NMEMO6S)
C
      CALL ASEL05S(CSTIFS,
     .             WORK1S(ISTIFS(37)),                      ! ESTIF
     .             WORK1S(ISTIFS(38)),                      ! WSTIF
     .             WORK1S(ISTIFS(41)))                      ! CSTI1
C
      RETURN
C
C**** EVALUATE INTERNAL RESISTING HEATS
C
   40 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATS=1
c     IF(IMICR.EQ.1) NEGATS=2
      CALL SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,PROPSS,
     .             SHAPES,THICKS,ELDISS,
     .             WORK1S(IFORCS(35)),WORK1S(IFORCS(36)),  ! cartd,dvolu
     .             WORK1S(IFORCS(37)),WORK1S(IFORCS(38)),  ! gpcod,shape
     .             WORK1S(IFORCS(39)),WORK1S(IFORCS(40)),  ! deriv,posgp
     .             WORK1S(IFORCS(41)),WORK1S(IFORCS(42)),  ! weigp,xjacm
     .             WORK1S(IFORCS(16)),WORK1S(IFORCS(43)),  ! elco1,emas1
     .                                WORK1S(IFORCS(44)),  ! eldi1
     .             WORK1S(IFORCS(49)),WORK1S(IFORCS(50)),  ! dvoli,elcoi
     .             WORK1S(IFORCS(13)),                     ! dispt
     .             0,ITASKS)
C
C**** EVALUATES WEIGHTING FUNCTIONS & ITS CARTESIAN DERIVATIVES
C
      CALL CART05S(WORK1S(IFORCS(16)),                    ! ELCOD
     .             PROPSS,
     .             WORK1S(IFORCS(44)),                    ! ELDIS
     .             WORK1S(IFORCS(38)),WORK1S(IFORCS(35)), ! SHAPE,CARTD
     .             WORK1S(IFORCS(19)),WORK1S(IFORCS(20)), ! posgp,weigp
     .             WORK1S(IFORCS(24)),WORK1S(IFORCS(17)), ! gpcdd,shapd
     .             WORK1S(IFORCS(21)),WORK1S(IFORCS(22)), ! derid,xjacm
     .             WORK1S(IFORCS(23)),                    ! cardd
     .             WORK1S(IFORCS(12)),                    ! velcm
     .             WORK1S(IFORCS(27)),WORK1S(IFORCS(28)), ! whape,hache
     .             WORK1S(IFORCS(29)),WORK1S(IFORCS(30)), ! wartd,weriv
     .                                WORK1S(IFORCS(32)), ! centr
     .                                WORK1S(IFORCS(34)), ! advem
     .             WORK1S(IFORCS(45)),WORK1S(IFORCS(48)), ! teini,fpchl
     .             EHISTS)
C
C**** TAKES INTO ACCOUNT THE MICROSTRUCTURAL MODEL
C
c     IF(IMICR.EQ.1) THEN
c      CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
c    .              WORK1T(IFORCT(44)),                    ! ELDIS
c    .              PROPST,
c    .              WORK1T(IFORCT(12)),
c    .                   2,TEMTLT,TEMDLT)
c      CALL MICR05T(WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), ! CARTD,DVOLU
c    .              EHISTT,
c    .              WORK1T(IFORCT(16)),                    ! ELCOD
c    .              WORK1T(IFORCT(44)),                    ! ELDIS
c    .              EPMTXT,
c    .              WORK1T(IFORCT(37)),                    ! GPCOD
c    .              LNODST,PROPST,RMAT1T,
c    .              WORK1T(IFORCT(38)),                    ! SHAPE
c    .              STRANT,STRSGT,STRA0T,STRS0T,
c    .              WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
c    .              WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
c    .              WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
c    .              WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
c    .              WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
c    .              WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
c    .              WORK1T(IFORCT(33)),                    ! vel1m
c    .              WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) ! teini,fpchl
c     ENDIF
C
C**** EVALUATES INTERNAL RESISTING HEAT DUE TO CONDUCTIVITY
C
      CALL FRIN05S(WORK1S(IFORCS(35)),WORK1S(IFORCS(36)),  ! CARTD,DVOLU
     .             EHISTS,
     .             WORK1S(IFORCS(16)),                     ! ELCOD
     .             WORK1S(IFORCS(44)),                     ! ELDIS
     .             EPMTXS,
     .             WORK1S(IFORCS(37)),                     ! GPCOD
     .             LNODSS,PROPSS,RMAT1S,
     .             WORK1S(IFORCS(38)),                     ! SHAPE
     .             STRANS,STRSGS,STRA0S,STRS0S,
     .             WORK1S(IFORCS( 1)),WORK1S(IFORCS( 2)),
     .             WORK1S(IFORCS( 3)),WORK1S(IFORCS( 4)),
     .             WORK1S(IFORCS( 5)),WORK1S(IFORCS( 6)),
     .             WORK1S(IFORCS( 7)),WORK1S(IFORCS( 8)),
     .             WORK1S(IFORCS( 9)),WORK1S(IFORCS(10)),
     .             WORK1S(IFORCS(11)),WORK1S(IFORCS(12)),
     .             WORK1S(IFORCS(29)),                     ! wartd
     .             WORK1S(IFORCS(45)),WORK1S(IFORCS(48)),  ! teini,fpchl
     .             WORK1S(IFORCS(49)))                     ! dvoli
C
C**** COMPUTES NODAL BOUNDARY CHANGES DUE TO DENSITY VARIATION
C
c     CALL BOUC05T(WORK1T(IFORCT(45)),                     ! TEINI
c    .             WORK1T(IFORCT(44)),                     ! ELDIS
c    .             WORK1T(IFORCT(46)))                     ! BOUCH
C
C**** EVALUATES INTERNAL RESISTING HEATS DUE TO CONVECTION EFFECTS
C
      IF(ICONVS.EQ.1) THEN
       CALL FRIC05S(WORK1S(IFORCS(35)),WORK1S(IFORCS(36)), ! CARTD,DVOLU
     .              EHISTS,
     .              WORK1S(IFORCS(16)),                    ! ELCOD
     .              WORK1S(IFORCS(44)),                    ! ELDIS
     .              EPMTXS,
     .              WORK1S(IFORCS(37)),                    ! GPCOD
     .              LNODSS,PROPSS,RMAT1S,
     .              WORK1S(IFORCS(38)),                    ! SHAPE
     .              STRANS,STRSGS,STRA0S,STRS0S,
     .              WORK1S(IFORCS( 1)),WORK1S(IFORCS( 2)),
     .              WORK1S(IFORCS( 3)),WORK1S(IFORCS( 4)),
     .              WORK1S(IFORCS( 5)),WORK1S(IFORCS( 6)),
     .              WORK1S(IFORCS( 7)),WORK1S(IFORCS( 8)),
     .              WORK1S(IFORCS( 9)),WORK1S(IFORCS(10)),
     .              WORK1S(IFORCS(11)),WORK1S(IFORCS(12)),
     .              WORK1S(IFORCS(27)),                          ! whape
     .                                 WORK1S(IFORCS(34)),       ! advem
     .              WORK1S(IFORCS(45)),WORK1S(IFORCS(48))) ! teini,fpchl
C
c      CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
c    .              WORK1T(IFORCT(44)),                    ! ELDIS
c    .              PROPST,
c    .              WORK1T(IFORCT(12)),
c    .                   2,TEMTLT,TEMDLT)
C
c      MULPHT=MULPDT
c      NPHCHT=NPHCDT
c      NISOTT=NISDTT
c      TEMPLT=TEMDLT
C
c      IF(MULPHT.EQ.1) THEN
c       CALL FRUC05T(WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), !CARTD,DVOLU
c    .               EHISTT,
c    .               WORK1T(IFORCT(16)),                    !ELCOD
c    .               WORK1T(IFORCT(44)),                    !ELDIS
c    .               EPMTXT,
c    .               WORK1T(IFORCT(37)),                    !GPCOD
c    .               LNODST,PROPST,RMAT1T,
c    .               WORK1T(IFORCT(38)),                    !SHAPE
c    .               STRANT,STRSGT,STRA0T,STRS0T,
c    .               WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
c    .               WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
c    .               WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
c    .               WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
c    .               WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
c    .               WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
c    .               WORK1T(IFORCT(27)),                         ! whape
c    .                                  WORK1T(IFORCT(34)),      ! advem
c    .               WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) !teini,fpchl
c      ELSE
C
c       IF(NISOTT.EQ.1) 
c    .   CALL RUNENDT('ERROR IN ELM005T-NISOTT=1') ! not implemented yet
C
c       IF(NISOTT.EQ.2) 
c    .   CALL FRMC05T(
c    .               WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), !CARTD,DVOLU
c    .               EHISTT,
c    .               WORK1T(IFORCT(16)),WORK1T(IFORCT(44)), !ELCOD,ELDIS
c    .               EPMTXT,WORK1T(IFORCT(37)),             !GPCOD
c    .               LNODST,PROPST,RMAT1T,WORK1T(IFORCT(38)), !SHAPE
c    .               STRANT,STRSGT,STRA0T,STRS0T,
c    .  WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),WORK1T(IFORCT( 3)),
c    .  WORK1T(IFORCT( 4)),WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
c    .  WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),WORK1T(IFORCT( 9)),
c    .  WORK1T(IFORCT(10)),WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
c    .               WORK1T(IFORCT(17)),WORK1T(IFORCT(18)), !shapd,dvold
c    .               WORK1T(IFORCT(19)),WORK1T(IFORCT(20)), !posgp,weigp
c    .               WORK1T(IFORCT(21)),                    !derid
c    .               WORK1T(IFORCT(23)),WORK1T(IFORCT(24)), !cardd,gpcdd
c    .               WORK1T(IFORCT(25)),WORK1T(IFORCT(26)), !elel1,elel2
c    .               WORK1T(IFORCT(28)),WORK1T(IFORCT(31)), !hache,whade
c    .               WORK1T(IFORCT(32)),WORK1T(IFORCT(34)), !centr,advem
c    .               WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) !teini,fpchl
c      ENDIF                        ! mulpht.eq.1
C
      ENDIF                         ! iconvs.eq.1
C
C**** EVALUATE EQUIVALENT VOLUME HEAT (NOT CONSTANT, i.e. AS INT. HEAT)
C
      IF(ILDGRS.EQ.1)
     . CALL LDGR05S(WORK1S(IFORCS(36)),       ! DVOLU
     .              PROPSS,
     .              WORK1S(IFORCS(38)),       ! SHAPE
     .              WORK1S(IFORCS( 1)),
     .              WORK1S(IFORCS(44)),       ! ELDIS
     .              WORK1S(IFORCS(12)),       ! VELCM
     .              EHISTS,ILDGRS,
     .              WORK1S(IFORCS(27)),WORK1S(IFORCS(45)), ! whape,teini
     .              WORK1S(IFORCS(48)),WORK1S(IFORCS(49)), ! fpchl,dvoli
     .              WORK1S(IFORCS(34)),WORK1S(IFORCS(35)), ! advem,cartd
     .              WORK1S(IFORCS(10)))                    ! XJACM
      RETURN
C
C**** EVALUATE INTERNAL TRANSIENT RESISTING HEATS
C
   60 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATS=2
      CALL SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,PROPSS,
     .             SHAPES,THICKS,ELDISS,
     .             WORK1S(IFORDS(23)),WORK1S(IFORDS(24)),  ! cartd,dvolu
     .             WORK1S(IFORDS(25)),WORK1S(IFORDS(26)),  ! gpcod,shape
     .             WORK1S(IFORDS(27)),WORK1S(IFORDS(28)),  ! deriv,posgp
     .             WORK1S(IFORDS(29)),WORK1S(IFORDS(30)),  ! weigp,xjacm
     .             WORK1S(IFORDS(17)),WORK1S(IFORDS(31)),  ! elco1,emas1
     .                                WORK1S(IFORDS(32)),  ! eldi1
     .             WORK1S(IFORDS(35)),WORK1S(IFORDS(36)),  ! dvoli,elcoi
     .             WORK1S(IFORDS(37)),                     ! dispt
     .             0,ITASKS)
C
C**** EVALUATES WEIGHTING FUNCTIONS
C
      CALL PETR05S(WORK1S(IFORDS(17)),                     ! ELCOD
     .             PROPSS,
     .             WORK1S(IFORDS(32)),                     ! ELDIS
     .             WORK1S(IFORDS(26)),WORK1S(IFORDS(23)),  ! SHAPE,CARTD
     .             WORK1S(IFORDS( 9)),WORK1S(IFORDS(10)),  ! posgp,weigp
     .             WORK1S(IFORDS(14)),WORK1S(IFORDS( 7)),  ! gpcod,shapd
     .             WORK1S(IFORDS(11)),WORK1S(IFORDS(12)),  ! derid,xjacm
     .             WORK1S(IFORDS(13)),                     ! cardd
     .             WORK1S(IFORDS( 2)),                     ! velcm
     .             WORK1S(IFORDS(18)),WORK1S(IFORDS(19)),  ! whape,hache
     .             WORK1S(IFORDS(11)),                     ! werid=derid
     .             WORK1S(IFORDS(21)),WORK1S(IFORDS(22)),  ! centr,advem
     .             WORK1S(IFORDS(33)),WORK1S(IFORDS(34)),  ! teini,fpchl
     .             EHISTS)
C
C**** EVALUATES INTERNAL TRANSIENT HEAT DUE TO SPECIFIC HEAT
C
      CALL FRDY05S(PROPSS,WORK1S(IFORDS(2)),WORK1S(IFORDS(1)),
     .             WORK1S(IFORDS(24)),                     ! DVOLU
     .             WORK1S(IFORDS(26)),                     ! SHAPE
     .             EHISTS,
     .             WORK1S(IFORDS(32)),                     ! ELDIS
     .             WORK1S(IFORDS( 3)),WORK1S(IFORDS( 6)),
     .                                WORK1S(IFORDS(18)),  ! whape
     .             WORK1S(IFORDS(33)),WORK1S(IFORDS(34)),  ! teini,fpchl
     .             WORK1S(IFORDS(35)))                     ! dvoli
C
C**** EVALUATE PHASE CHANGE INTERNAL DYNAMIC RESISTING HEATS
C
c     CALL CERO00T(WORK1T(IFORDT(16)),WORK1T(IFORDT(15)),NEVABT)
c     CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
c    .             WORK1T(IFORDT(32)),                     ! ELDIS
c    .             PROPST,WORK1T(IFORDT(2)), 2,TEMTLT,
c    .             TEMDLT)
C
c     DO INPCCT=1,2
c      IF(INPCCT.EQ.1) THEN
c       MULPHT=MULPTT
c       NPHCHT=NPHCTT
c       NISOTT=NISTTT
c       TEMPLT=TEMTLT
c      ELSE
c       MULPHT=MULPDT
c       NPHCHT=NPHCDT
c       NISOTT=NISDTT
c       TEMPLT=TEMDLT
c      ENDIF
C
c      IF(MULPHT.EQ.1) THEN
c       CALL FRUF05T(PROPST,WORK1T(IFORDT(2)),WORK1T(IFORDT(1)),
c    .               WORK1T(IFORDT(24)),                    ! DVOLU
c    .               WORK1T(IFORDT(26)),                    ! SHAPE
c    .               EHISTT,
c    .               WORK1T(IFORDT(32)),                    ! ELDIS
c    .               WORK1T(IFORDT( 3)),
c    .                                  WORK1T(IFORDT( 6)),
c    .               WORK1T(IFORDT(15)),WORK1T(IFORDT(16)),
c    .                                  WORK1T(IFORDT(18)), ! whape
c    .               INPCCT,
c    .              WORK1T(IFORDT(33)),WORK1T(IFORDT(34)),  !teini,fpchl
c    .              WORK1T(IFORDT(35)))                     !dvoli
c      ELSE
c       IF(NISOTT.EQ.1) 
c    .   CALL FRTF05T(PROPST,WORK1T(IFORDT(2)),WORK1T(IFORDT(1)),
c    .                WORK1T(IFORDT(24)),                        ! DVOLU
c    .                WORK1T(IFORDT(26)),                        ! SHAPE
c    .                EHISTT,
c    .                WORK1T(IFORDT(32)),                        ! ELDIS
c    .                WORK1T(IFORDT( 3)),
c    .                WORK1T(IFORDT( 6)),WORK1T(IFORDT( 7)),
c    .                WORK1T(IFORDT( 8)),
c    .                WORK1T(IFORDT(25)),                        ! GPCOD
c    .                NPHCHT,
c    .                WORK1T(IFORDT( 9)),
c    .                WORK1T(IFORDT(17)),                        ! ELCOD
c    .                WORK1T(IFORDT(10)),
c    .                WORK1T(IFORDT(11)),WORK1T(IFORDT(12)),
c    .                WORK1T(IFORDT(13)),WORK1T(IFORDT(14)),
c    .                WORK1T(IFORDT(15)),WORK1T(IFORDT(16)),
c    .               WORK1T(IFORDT(19)),WORK1T(IFORDT(20)), !hache,whade
c    .               WORK1T(IFORDT(21)),WORK1T(IFORDT(22)), !centr,advem
c    .               WORK1T(IFORDT(33)),WORK1T(IFORDT(34)), !teini,fpchl
c    .               WORK1T(IFORDT(35)),                    !dvoli
c    .                LNODST,INPCCT,TEMPLT)
c       IF(NISOTT.EQ.2) 
c    .   CALL FRMF05T(PROPST,WORK1T(IFORDT(2)),WORK1T(IFORDT(1)),
c    .                WORK1T(IFORDT(24)),                        ! DVOLU
c    .                WORK1T(IFORDT(26)),                        ! SHAPE
c    .                EHISTT,
c    .                WORK1T(IFORDT(32)),                        ! ELDIS
c    .                WORK1T(IFORDT(3)),
c    .                WORK1T(IFORDT(6)),
c    .                WORK1T(IFORDT(7)),WORK1T(IFORDT(8)),
c    .                WORK1T(IFORDT(25)),                        ! GPCOD
c    .                WORK1T(IFORDT(9)),
c    .                WORK1T(IFORDT(17)),                        ! ELCOD
c    .                WORK1T(IFORDT(10)),
c    .                WORK1T(IFORDT(11)),WORK1T(IFORDT(12)),
c    .                WORK1T(IFORDT(13)),WORK1T(IFORDT(14)),
c    .                WORK1T(IFORDT(15)),WORK1T(IFORDT(16)),
c    .               WORK1T(IFORDT(19)),WORK1T(IFORDT(20)), !hache,whade
c    .               WORK1T(IFORDT(21)),WORK1T(IFORDT(22)), !centr,advem
c    .               WORK1T(IFORDT(33)),WORK1T(IFORDT(34)), !teini,fpchl
c    .               WORK1T(IFORDT(35)),                    !dvoli
c    .                LNODST,INPCCT)
c      ENDIF                   ! mulpht.eq.1
c     ENDDO                    ! inpcct=1,2
C
c     CALL SUMMATT(WORK1T(IFORDT(1)),WORK1T(IFORDT(16)),
c    .             WORK1T(IFORDT(15)),NEVABT)
C
C**** COMPUTES NODAL PHASE-CHANGE FUNCTION AT TIME t+dt
C     (useful in the mechanical problem)
C
c     CALL FPCN05T(WORK1T(IFORDT(32)),                     ! ELDIS
c    .             WORK1T(IFORDT(2)),                      ! VELCM
c    .             WORK1T(IFORDT(34)),     1)              ! FPCHL
      RETURN
C
C**** EVALUATE EQUIVALENT VOLUME HEAT (CONSTANT !!!, i.e. AS EXT. HEAT)
C
   70 CONTINUE
      IF(ILDGRS.EQ.1) RETURN
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATS=1
c     IF(IMICR.EQ.1) NEGATS=2
      CALL SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,PROPSS,
     .             SHAPES,THICKS,ELDISS,
     .             WORK1S(ILOSUS(21)),WORK1S(ILOSUS(22)),  ! cartd,dvolu
     .             WORK1S(ILOSUS(23)),WORK1S(ILOSUS(24)),  ! gpcod,shape
     .             WORK1S(ILOSUS(25)),WORK1S(ILOSUS(26)),  ! deriv,posgp
     .             WORK1S(ILOSUS(27)),WORK1S(ILOSUS(28)),  ! weigp,xjacm
     .             WORK1S(ILOSUS(17)),WORK1S(ILOSUS(29)),  ! elco1,emas1
     .                                WORK1S(ILOSUS(31)),  ! eldi1
     .             WORK1S(ILOSUS(35)),WORK1S(ILOSUS(36)),  ! dvoli,elcoi
     .             WORK1S(ILOSUS(37)),                     ! dispt
     .             0,ITASKS)
C
C**** EVALUATES WEIGHTING FUNCTIONS
C
      CALL PETR05S(WORK1S(ILOSUS(17)),                     ! ELCOD
     .             PROPSS,
     .             WORK1S(ILOSUS(31)),       ! TEINI=ELDIS for this case
     .             WORK1S(ILOSUS(24)),WORK1S(ILOSUS(21)),  ! SHAPE,CARTD
     .             WORK1S(ILOSUS(26)),WORK1S(ILOSUS(27)),  ! posgp,weigp
     .             WORK1S(ILOSUS( 4)),WORK1S(ILOSUS(10)),  ! gpcod,shapd
     .             WORK1S(ILOSUS( 2)),WORK1S(ILOSUS(13)),  ! derid,xjacm
     .             WORK1S(ILOSUS(14)),                     ! cardd
     .             WORK1S(ILOSUS(16)),                     ! velcm
     .             WORK1S(ILOSUS(18)),WORK1S(ILOSUS(19)),  ! whape,hache
     .             WORK1S(ILOSUS(25)),                     ! weriv=derid
     .             WORK1S(ILOSUS(20)),WORK1S(ILOSUS(32)),  ! centr,advem
     .             WORK1S(ILOSUS(31)),WORK1S(ILOSUS(34)),  ! teini,fpchl
     .             EHISTS)
C
C**** EVALUATE HEAT
C
      CALL LDGR05S(WORK1S(ILOSUS(22)),                    ! DVOLU
     .             PROPSS,
     .             WORK1S(ILOSUS(24)),                    ! SHAPE
     .             WORK1S(ILOSUs( 1)),
     .             WORK1S(ILOSUS(31)),      ! TEINI=ELDIS for this case
     .             WORK1S(ILOSUS(16)),                    ! VELCM
     .             EHISTS,ILDGRS,
     .             WORK1S(ILOSUS(18)),WORK1S(ILOSUS(31)), ! whape,teini
     .             WORK1S(ILOSUS(34)),WORK1S(ILOSUS(35)), ! fpchl,dvoli
     .             WORK1S(ILOSUS(32)),WORK1S(ILOSUS(21)), ! advem,cartd
     .             WORK1S(IFORCS(13)))                    ! XJACM
      RETURN
C
C**** NOTHING
C
   80 CONTINUE
      RETURN
C
C**** CHECK CORRECTNESS OF INTEGRATION RULES
C
   90 CONTINUE
      CALL CHEK01S(NDIMES,IELEMS,NNODLS,NRULES,NGAULS,     5)
      RETURN
C
C**** NOTHING
C
  100 CONTINUE
      RETURN
C
C**** NOTHING
C
  110 CONTINUE
      RETURN
C
C**** OUTPUT GAUSSIAN VARIABLES
C
  120 CONTINUE
c     CALL OUTG05S(EHISTS,STRANS,STRSGS)
      RETURN
C
C**** OUTPUT NODAL (SMOOTHED) HEAT FLUXES
C
  130 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATS=1
      IF(IMICR.EQ.1) NEGATS=2
      IF(ISMO1S.EQ.1) THEN
       NRUAUX=NRULES
       NRULES=5
      ENDIF
      CALL SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,PROPSS,
     .             SHAPES,THICKS,ELDISS,
     .             WORK1S(IGSMOS(11)),WORK1S(IGSMOS(12)),  ! cartd,dvolu
     .             WORK1S(IGSMOS(13)),WORK1S(IGSMOS(14)),  ! gpcod,shape
     .             WORK1S(IGSMOS(15)),WORK1S(IGSMOS(16)),  ! deriv,posgp
     .             WORK1S(IGSMOS(17)),WORK1S(IGSMOS(18)),  ! weigp,xjacm
     .             WORK1S(IGSMOS(19)),WORK1S(IGSMOS(20)),  ! elco1,emas1
     .                                WORK1S(IGSMOS(21)),  ! eldi1
     .             WORK1S(IGSMOS(28)),WORK1S(IGSMOS(29)),  ! dvoli,elcoi
     .             WORK1S(IGSMOS(30)),                     ! dispt
     .             0,ITASKS)
C
      CALL SMOG05S(WORK1S(IGSMOS(12)),                     ! DVOLU
     .             WORK1S(IGSMOS(20)),                     ! EMASS
     .             WORK1S(IGSMOS(14)),                     ! SHAPE
     .             STRSGS,PROPSS,
     .             WORK1S(IGSMOS( 5)),WORK1S(IGSMOS(4)))   ! STREB,STREA
      IF(ISMO1S.EQ.1) NRULES=NRUAUX
      RETURN
C
C**** OUTPUT NODAL (SMOOTHED) TEMP.GRADIENTS
C
  140 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
c     NEGATT=1
c     IF(IMICR.EQ.1) NEGATT=2
c     IF(ISMO1T.EQ.1) THEN
c      NRUAUX=NRULET
c      NRULET=5
c     ENDIF
c     CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
c    .             SHAPET,THICKT,ELDIST,
c    .             WORK1T(IGSMOT(11)),WORK1T(IGSMOT(12)),  ! cartd,dvolu
c    .             WORK1T(IGSMOT(13)),WORK1T(IGSMOT(14)),  ! gpcod,shape
c    .             WORK1T(IGSMOT(15)),WORK1T(IGSMOT(16)),  ! deriv,posgp
c    .             WORK1T(IGSMOT(17)),WORK1T(IGSMOT(18)),  ! weigp,xjacm
c    .             WORK1T(IGSMOT(19)),WORK1T(IGSMOT(20)),  ! elco1,emas1
c    .                                WORK1T(IGSMOT(21)),  ! eldi1
c    .             WORK1T(IGSMOT(28)),WORK1T(IGSMOT(29)),  ! dvoli,elcoi
c    .             WORK1T(IGSMOT(30)),                     ! dispt
c    .             0,ITASKT)
C
c     CALL SMOF05T(WORK1T(IGSMOT(12)),                     ! DVOLU
c    .             WORK1T(IGSMOT(20)),                     ! EMASS
c    .             WORK1T(IGSMOT(14)),                     ! SHAPE
c    .             STRANT,PROPST,
c    .             WORK1T(IGSMOT( 7)),WORK1T(IGSMOT(6)))   ! SFISB,SFISA
c     IF(ISMO1T.EQ.1) NRULET=NRUAUX
      RETURN
C
C**** OUTPUT NODAL (SMOOTHED) INTERNAL VARIABLES
C
  150 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
c     NEGATT=1
c     IF(IMICR.EQ.1) NEGATT=2
c     IF(ISMO1T.EQ.1) THEN
c      NRUAUX=NRULET
c      NRULET=5
c     ENDIF
c     CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
c    .             SHAPET,THICKT,ELDIST,
c    .             WORK1T(IGSMOT(11)),WORK1T(IGSMOT(12)),  ! cartd,dvolu
c    .             WORK1T(IGSMOT(13)),WORK1T(IGSMOT(14)),  ! gpcod,shape
c    .             WORK1T(IGSMOT(15)),WORK1T(IGSMOT(16)),  ! deriv,posgp
c    .             WORK1T(IGSMOT(17)),WORK1T(IGSMOT(18)),  ! weigp,xjacm
c    .             WORK1T(IGSMOT(19)),WORK1T(IGSMOT(20)),  ! elco1,emas1
c    .                                WORK1T(IGSMOT(21)),  ! eldi1
c    .             WORK1T(IGSMOT(28)),WORK1T(IGSMOT(29)),  ! dvoli,elcoi
c    .             WORK1T(IGSMOT(30)),                     ! dispt
c    .             0,ITASKT)
C
c     CALL SMOI05T(WORK1T(IGSMOT(12)),                     ! DVOLU
c    .             WORK1T(IGSMOT(20)),                     ! EMASS
c    .             WORK1T(IGSMOT(14)),                     ! SHAPE
c    .             EHISTT,PROPST,
c    .             WORK1T(IGSMOT(10)),WORK1T(IGSMOT(9)))   ! SFIPB,SFIPA
c     IF(ISMO1T.EQ.1) NRULET=NRUAUX
      RETURN
C
C**** OUTPOST
C
  160 CONTINUE
      NEGATS=1
c     IF(IMICR.EQ.1) NEGATS=2
      CALL SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,PROPSS,
     .             SHAPES,THICKS,ELDISS,
     .             WORK1S(ISTARS( 3)),WORK1S(ISTARS( 4)),  ! cartd,dvolu
     .             WORK1S(ISTARS( 5)),WORK1S(ISTARS( 6)),  ! gpcod,shape
     .             WORK1S(ISTARS( 7)),WORK1S(ISTARS( 8)),  ! deriv,posgp
     .             WORK1S(ISTARS( 9)),WORK1S(ISTARS(10)),  ! weigp,xjacm
     .             WORK1S(ISTARS(11)),WORK1S(ISTARS(12)),  ! elco1,emas1
     .                                WORK1S(ISTARS(13)),  ! eldi1
     .             WORK1S(ISTARS(15)),WORK1S(ISTARS(16)),  ! dvoli,elcoi
     .             WORK1S(ISTARS(17)),                     ! dispt
     .             0,ITASKS)
      RETURN
C
C**** OUTPUT NODAL (SMOOTHED) POROSITY CRITERIA
C
  170 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
c     IF(ISMO1T.EQ.1) THEN
c      NRUAUX=NRULET
c      NRULET=5
c     ENDIF
c     CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
c    .             SHAPET,THICKT,ELDIST,
c    .             WORK1T(IGSMOT(11)),WORK1T(IGSMOT(12)),  ! cartd,dvolu
c    .             WORK1T(IGSMOT(13)),WORK1T(IGSMOT(14)),  ! gpcod,shape
c    .             WORK1T(IGSMOT(15)),WORK1T(IGSMOT(16)),  ! deriv,posgp
c    .             WORK1T(IGSMOT(17)),WORK1T(IGSMOT(18)),  ! weigp,xjacm
c    .             WORK1T(IGSMOT(19)),WORK1T(IGSMOT(20)),  ! elco1,emas1
c    .                                WORK1T(IGSMOT(21)),  ! eldi1
c    .             WORK1T(IGSMOT(28)),WORK1T(IGSMOT(29)),  ! dvoli,elcoi
c    .             WORK1T(IGSMOT(30)),                     ! dispt
c    .             0,ITASKT)
C
c     CALL PORO05T(WORK1T(IGSMOT(11)),WORK1T(IGSMOT(12)),  ! CARTD,DVOLU
c    .             WORK1T(IGSMOT(20)),                     ! EMASS
c    .             WORK1T(IGSMOT(14)),                     ! SHAPE
c    .             PROPST,
c    .             WORK1T(IGSMOT(24)),WORK1T(IGSMOT(23)),  ! SPIPB,SPIPA
c    .             WORK1T(IGSMOT(25)),                     ! SPIPC
c    .             WORK1T(IGSMOT(21)),WORK1T(IGSMOT(26)),  ! ELDI1,FPCHL
c    .             WORK1T(IGSMOT(27)))                     ! VELCM
c     IF(ISMO1T.EQ.1) NRULET=NRUAUX
      RETURN
C
      END
