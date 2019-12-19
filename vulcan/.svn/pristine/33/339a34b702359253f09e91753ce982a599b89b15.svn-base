      SUBROUTINE ELM005T(LNODST,PROELT,PROPST,WORK1T,
     .                   ELCODT,CARTDT,DVOLUT,GPCODT,SHAPET,EPMTXT,
     .                   RMAT1T,EMASST,                         !ELDAT
     .                   STRA0T,STRS0T,TEMPCT,                  !ELPRE
     .                   ELDIST,EHISTT,STRANT,STRSGT,           !ELVAR
     .                   CSTIFT,ESTIFT,WSTIFT,HSTIFT,PSTIFT,
     .                   QSTIFT,                                !ELMAT
     .                                               ITASKT)
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
C                               Phase Change Matrix
C       ITASK =  4   Evaluate Internal Resisting Heat
C       ITASK =  5   Evaluate Internal Resisting Heat
C                    (alternative option for micro advective effects)
C       ITASK =  6   Evaluate : Internal Dynamic Resisting Heat
C                               (only for Transient Problems)
C                               Internal Dynamic Phase Change Heat
C       ITASK =  7   Evaluate equivalent Volume Heat (Internal Heat)
C       ITASK =  8   -------
C       ITASK =  9   Check correctness of integration rule
C       ITASK = 10   -------
C       ITASK = 11   -------
C       ITASK = 12   Output Gaussian variables
C       ITASK = 13   Evaluate contribution for nodal heat fluxes
C       ITASK = 14   -------
C       ITASK = 15   Evaluate contribution for nodal internal variables
C       ITASK = 16   Evaluate variables for outpost.f
C       ITASK = 17   Evaluate contribution for nodal porosity criteria
C       ITASK = 18   Nothing (see outsmot.f)
C       ITASK = 19   Evaluate contribution for nodal L*f_pc
C
C
C
C     TREATMENT OF SMOOTHING OPERATION (IF NEEDED)
C
C     ISMO1T=0 => use standard integration rule for smoothing
C           =1 => use Lobatto type integration rule for smoothing
C           =2 => idem ISMO1T=0 but, for axisymmetric problems, computes
C                 the 2D DVOLU (warning: this option is not compatible
C                 with nodal integration rules)
C
C     ISMO2T=0 => do not diagonalise mass matrix for smoothing
C           =1 => diagonalise mass matrix for smoothing
C
C
C     Notes:
C
C     The option ISMO1T=0 and ISMO2T=0 does not work well
C     The option ISMO1T=1 (for any ISMO2T) works well but it is not
C     absolutely correct (see rulius.f)
C     The option ISMO1T=0 and ISMO2T=1 works well and it is more correct
C     than the former
C     Both parameters are input in the data file (see coninpt.f)
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
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
      INCLUDE 'nuef_om.f'   ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/THICKNESS/THICKT
      COMMON/MULPHASE/MULPHT
      COMMON/NEGATITE/NEGATT
C
      DIMENSION LNODST(*), PROELT(*), PROPST(*)
      DIMENSION ELCODT(*), CARTDT(*), DVOLUT(*), GPCODT(*), SHAPET(*), 
     .          EPMTXT(*), RMAT1T(*), EMASST(*)
      DIMENSION STRA0T(*), STRS0T(*), TEMPCT(*)
      DIMENSION ELDIST(*), EHISTT(*), STRANT(*), STRSGT(*)
      DIMENSION CSTIFT(*), ESTIFT(*), WSTIFT(*), HSTIFT(*), PSTIFT(*),
     .          QSTIFT(*)
      DIMENSION WORK1T(*)
C
C**** SELECTS ELEMENT PROPERTIES
C
      NNODLT=INT(PROELT( 2))
      NRULET=INT(PROELT( 3))
      NGAULT=INT(PROELT( 4))
      NTYPET=INT(PROELT( 6))
      THICKT=    PROELT( 7)
C
      NNODST=NNODLT                ! NNODS used for smoothing operations
C
C**** SELECTS TREATMENT OF SOURCE TERM (IF IT EXISTS)
C
      ILDGRT=1                     ! better as input
C
C**** CONSIDERS DIFFERENT ELEMENTS WITH THE SAME NUMBER OF NODES
C
      NQUTRT=1                     ! to be completed !!!!!! not used yet
      IF(NDIMET.EQ.2) THEN
       IF(NNODLT.EQ.4) THEN
        IF(NRULET.EQ.3.OR.NRULET.EQ.6.OR.NRULET.EQ.7) NQUTRT=2
       ENDIF
      ENDIF
c     IF(NDIMET.EQ.3) THEN
c      CALL RUNENDT('ELM005T: 3D NOT IMPLEMENTED YET')
c     ENDIF
C
C**** DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90,100,
     .       110,120,130,140,150,160,170,  1,190), ITASKT
    1 RETURN ! Nothing
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
   10 CONTINUE
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(ISETMT( 7)),WORK1T(ISETMT( 8)),  ! cartd,dvolu
     .             WORK1T(ISETMT( 9)),WORK1T(ISETMT(10)),  ! gpcod,shape
     .             WORK1T(ISETMT( 1)),WORK1T(ISETMT( 2)),  ! deriv,posgp
     .             WORK1T(ISETMT( 3)),WORK1T(ISETMT( 4)),  ! weigp,xjacm
     .             WORK1T(ISETMT( 6)),WORK1T(ISETMT(11)),  ! elco1,emas1
     .                                WORK1T(ISETMT(12)),  ! eldi1
     .             WORK1T(ISETMT(16)),WORK1T(ISETMT(17)),  ! dvoli,elcoi
     .             WORK1T(ISETMT(18)),                     ! dispt
     .             0,ITASKT)
C
      IF(IMICR.EQ.1)
     . CALL INMI05T(PROPST,EHISTT)
C
      CALL FPCN05T(WORK1T(ISETMT(13)),                     ! eldis=teini
     .             WORK1T(ISETMT(12)),                     ! not used
     .             PROPST,                                 ! not used
     .             WORK1T(ISETMT(14)),     0,     0)       ! FPCHL
      RETURN
C
C**** NOTHING
C
   20 CONTINUE
      RETURN
C
C**** EVALUATES CONDUCTIVITY MATRIX. IN TRANSIENT ANALYSIS, EVALUATES:
C     1) HEAT CAPACITY MATRIX, 2) "PHASE-CHANGE" MATRIX & 3) ADVECTIVE
C     MATRICES
C
   30 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(ISTIFT(26)),WORK1T(ISTIFT(27)),  ! cartd,dvolu
     .             WORK1T(ISTIFT(28)),WORK1T(ISTIFT(29)),  ! gpcod,shape
     .             WORK1T(ISTIFT(30)),WORK1T(ISTIFT(31)),  ! deriv,posgp
     .             WORK1T(ISTIFT(32)),WORK1T(ISTIFT(33)),  ! weigp,xjacm
     .             WORK1T(ISTIFT(18)),WORK1T(ISTIFT(34)),  ! elco1,emas1
     .                                WORK1T(ISTIFT(39)),  ! eldi1
     .             WORK1T(ISTIFT(46)),WORK1T(ISTIFT(47)),  ! dvoli,elcoi
     .             WORK1T(ISTIFT(44)),                     ! dispt
     .             0,ITASKT)
C
C**** EVALUATES WEIGHTING FUNCTIONS & THEIR CARTESIAN DERIVATIVES
C
      CALL CART05T(WORK1T(ISTIFT(18)),                     ! ELCOD
     .             PROPST,
     .             WORK1T(ISTIFT(39)),                     ! ELDIS
     .             WORK1T(ISTIFT(29)),WORK1T(ISTIFT(26)),  ! SHAPE,CARTD
     .             WORK1T(ISTIFT( 7)),WORK1T(ISTIFT( 8)),  ! posgp,weigp
     .             WORK1T(ISTIFT(12)),WORK1T(ISTIFT( 5)),  ! gpcdd,shapd
     .             WORK1T(ISTIFT( 9)),WORK1T(ISTIFT(10)),  ! derid,xjacm
     .             WORK1T(ISTIFT(11)),                     ! cardd
     .             WORK1T(ISTIFT(13)),                     ! velcm
     .             WORK1T(ISTIFT(19)),WORK1T(ISTIFT(20)),  ! whape,hache
     .             WORK1T(ISTIFT(21)),WORK1T(ISTIFT(22)),  ! wartd,weriv
     .                                WORK1T(ISTIFT(24)),  ! centr
     .             WORK1T(ISTIFT(25)),                     ! advem
     .             WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)),  ! teini,fpchl
     .             EHISTT)
C
      CALL ASIL00(WORK1T(ISTIFT(37)),WORK1T(ISTIFT(38)),   ! ESTIF,WSTIF
     .            NKOVAT,KDYNAT,
     .            WORK1T(ISTIFT(49)),                      ! EMATX
     .            NEVABT,     1)
C
C**** EVALUATES CONDUCTIVITY MATRIX
C
      CALL STIF05T(WORK1T(ISTIFT(26)),WORK1T(ISTIFT(27)),  ! CARTD,DVOLU
     .             EHISTT,
     .             WORK1T(ISTIFT(39)),                     ! ELDIS
     .             EPMTXT,
     .             WORK1T(ISTIFT(28)),                     ! GPCOD
     .             PROPST,
     .             WORK1T(ISTIFT(29)),                     ! SHAPE
     .             STRSGT,
     .             WORK1T(ISTIFT(37)),                     ! ESTIF
     .             HSTIFT,
     .             WORK1T(ISTIFT(18)),                     ! ELCOD
     .             WORK1T(ISTIFT( 1)),WORK1T(ISTIFT( 2)),
     .             WORK1T(ISTIFT( 3)),WORK1T(ISTIFT( 4)),
     .             WORK1T(ISTIFT(13)),                     ! velcm
     .             WORK1T(ISTIFT(21)),                     ! wartd
     .             WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)),  ! teini,fpchl
     .             WORK1T(ISTIFT(46)),WORK1T(ISTIFT(49)))  ! dvoli,EMATX
C
      IF(KDYNAT.EQ.1) THEN
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
       NEGATT=2
       CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .              SHAPET,THICKT,ELDIST,
     .              WORK1T(ISTIFT(26)),WORK1T(ISTIFT(27)), ! cartd,dvolu
     .              WORK1T(ISTIFT(28)),WORK1T(ISTIFT(29)), ! gpcod,shape
     .              WORK1T(ISTIFT(30)),WORK1T(ISTIFT(31)), ! deriv,posgp
     .              WORK1T(ISTIFT(32)),WORK1T(ISTIFT(33)), ! weigp,xjacm
     .              WORK1T(ISTIFT(18)),WORK1T(ISTIFT(34)), ! elco1,emas1
     .                                 WORK1T(ISTIFT(39)), ! eldi1
     .              WORK1T(ISTIFT(46)),WORK1T(ISTIFT(47)), ! dvoli,elcoi
     .              WORK1T(ISTIFT(44)),                    ! dispt
     .              1,ITASKT)
C
C**** EVALUATES WEIGHTING FUNCTIONS & THEIR CARTESIAN DERIVATIVES
C
       CALL CART05T(WORK1T(ISTIFT(18)),                    ! ELCOD
     .              PROPST,
     .              WORK1T(ISTIFT(39)),                    ! ELDIS
     .              WORK1T(ISTIFT(29)),WORK1T(ISTIFT(26)), ! SHAPE,CARTD
     .              WORK1T(ISTIFT( 7)),WORK1T(ISTIFT( 8)), ! posgp,weigp
     .              WORK1T(ISTIFT(12)),WORK1T(ISTIFT( 5)), ! gpcdd,shapd
     .              WORK1T(ISTIFT( 9)),WORK1T(ISTIFT(10)), ! derid,xjacm
     .              WORK1T(ISTIFT(11)),                    ! cardd
     .              WORK1T(ISTIFT(13)),                    ! velcm
     .              WORK1T(ISTIFT(19)),WORK1T(ISTIFT(20)), ! whape,hache
     .              WORK1T(ISTIFT(21)),WORK1T(ISTIFT(22)), ! wartd,weriv
     .                                 WORK1T(ISTIFT(24)), ! centr
     .              WORK1T(ISTIFT(25)),                    ! advem
     .              WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)), ! teini,fpchl
     .              EHISTT)
C
C**** EVALUATES CAPACITY MATRIX
C
       CALL MASS05T(WORK1T(ISTIFT(27)),                    ! DVOLU
     .              PROPST,
     .              WORK1T(ISTIFT(29)),                    ! SHAPE
     .              WORK1T(ISTIFT(38)),                    ! WSTIF
     .              EHISTT,
     .              WORK1T(ISTIFT(19)),                    ! whape
     .              WORK1T(ISTIFT(39)),                    ! ELDIS
     .              WORK1T(ISTIFT(13)),                    ! velcm
     .              WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)), ! teini,fpchl
     .              WORK1T(ISTIFT(46)),WORK1T(ISTIFT(48))) ! dvoli,EMATX
C
       CALL CERO00T(WORK1T(ISTIFT(15)),WORK1T(ISTIFT(14)),NKOVAT)
       CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
     .              WORK1T(ISTIFT(39)),                    ! ELDIS
     .              PROPST,WORK1T(ISTIFT(13)),1,TEMTLT,
     .              TEMDLT)
C
C**** EVALUATES PHASE-CHANGE MATRIX
C
       IF(MULPHT.EQ.1) THEN
        CALL MAUF05T(WORK1T(ISTIFT(27)),                   ! DVOLU
     .               PROPST,
     .               WORK1T(ISTIFT(29)),                   ! SHAPE
     .               WORK1T(ISTIFT(38)),                   ! WSTIF
     .               EHISTT,
     .               WORK1T(ISTIFT(13)),WORK1T(ISTIFT(14)),
     .               WORK1T(ISTIFT(15)),
     .               WORK1T(ISTIFT(39)),                   ! ELDIS
     .               WORK1T(ISTIFT(17)),                   ! disit
     .               WORK1T(ISTIFT(19)),                   ! whape
     .               WORK1T(ISTIFT(40)),                   ! teini
     .               WORK1T(ISTIFT(45)),                   ! fpchl
     .               WORK1T(ISTIFT(46)))                   ! dvoli
       ELSE
        CALL MAMF05T(WORK1T(ISTIFT(27)),                   ! DVOLU
     .               PROPST,
     .               WORK1T(ISTIFT(29)),                   ! SHAPE
     .               WORK1T(ISTIFT(38)),                   ! WSTIF
     .               EHISTT,
     .               WORK1T(ISTIFT(5)),WORK1T(ISTIFT(6)),
     .               WORK1T(ISTIFT(7)),
     .               WORK1T(ISTIFT(39)),                   ! ELDIS
     .               WORK1T(ISTIFT(18)),                   ! ELCOD
     .               WORK1T(ISTIFT(8)),WORK1T(ISTIFT(9)),
     .               WORK1T(ISTIFT(10)),WORK1T(ISTIFT(11)),
     .               WORK1T(ISTIFT(12)),WORK1T(ISTIFT(13)),
     .               WORK1T(ISTIFT(14)),WORK1T(ISTIFT(15)),
     .               LNODST,
     .               WORK1T(ISTIFT(17)),
     .               WORK1T(ISTIFT(20)),WORK1T(ISTIFT(23)), !hache,whade
     .               WORK1T(ISTIFT(24)),WORK1T(ISTIFT(25)), !centr,advem
     .               WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)), !teini,fpchl
     .               WORK1T(ISTIFT(46)))                    !dvoli
       ENDIF            ! mulpht.eq.1
C
C**** CHECK THE CORRECTNESS OF THE PHASE CHANGE MATRIX AND ADDS ITS 
C     CONTRIBUTION TO THE "JACOBIAN MATRIX"
C
       CALL VERIPCT(WORK1T(ISTIFT(15)),WORK1T(ISTIFT(14)),NKOVAT,
     .              IVEPCT)
       IF(IVEPCT.NE.0) WRITE(LUREST,2000) IVEPCT
C
 2000  FORMAT(' *** JACOBIAN MATRIX WARNING *** IVEPCT ELEMENTS OF',
     1        ' PHASE-CHANGE MATRIX ARE NEGATIVE',I5)
C
       CALL SUMMATT(WORK1T(ISTIFT(38)),                  ! WSTIF
     .              WORK1T(ISTIFT(15)),WORK1T(ISTIFT(14)),
     .              NKOVAT)
C
C**** PERFORM A LUMPED HEAT CAPACITY + PHASE CHANGE MATRICES
C
       IF(NCETAT.EQ.1) CALL LUMPMTT(WORK1T(ISTIFT(38)),          ! WSTIF
     .                              NKOVAT,NNODLT,KSYMMT,NEVABT)
C
      ENDIF                         ! kdynat.eq.1
C
      IF(ICONVT.EQ.1) THEN
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
       NEGATT=1
       IF(IMICR.EQ.1) NEGATT=2
       CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .              SHAPET,THICKT,ELDIST,
     .              WORK1T(ISTIFT(26)),WORK1T(ISTIFT(27)), ! cartd,dvolu
     .              WORK1T(ISTIFT(28)),WORK1T(ISTIFT(29)), ! gpcod,shape
     .              WORK1T(ISTIFT(30)),WORK1T(ISTIFT(31)), ! deriv,posgp
     .              WORK1T(ISTIFT(32)),WORK1T(ISTIFT(33)), ! weigp,xjacm
     .              WORK1T(ISTIFT(18)),WORK1T(ISTIFT(34)), ! elco1,emas1
     .                                 WORK1T(ISTIFT(39)), ! eldi1
     .              WORK1T(ISTIFT(46)),WORK1T(ISTIFT(47)), ! dvoli,elcoi
     .              WORK1T(ISTIFT(44)),                    ! dispt
     .              1,ITASKT)
C
C**** EVALUATES WEIGHTING FUNCTIONS & THEIR CARTESIAN DERIVATIVES
C
       CALL CART05T(WORK1T(ISTIFT(18)),                    ! ELCOD
     .              PROPST,
     .              WORK1T(ISTIFT(39)),                    ! ELDIS
     .              WORK1T(ISTIFT(29)),WORK1T(ISTIFT(26)), ! SHAPE,CARTD
     .              WORK1T(ISTIFT( 7)),WORK1T(ISTIFT( 8)), ! posgp,weigp
     .              WORK1T(ISTIFT(12)),WORK1T(ISTIFT( 5)), ! gpcdd,shapd
     .              WORK1T(ISTIFT( 9)),WORK1T(ISTIFT(10)), ! derid,xjacm
     .              WORK1T(ISTIFT(11)),                    ! cardd
     .              WORK1T(ISTIFT(13)),                    ! velcm
     .              WORK1T(ISTIFT(19)),WORK1T(ISTIFT(20)), ! whape,hache
     .              WORK1T(ISTIFT(21)),WORK1T(ISTIFT(22)), ! wartd,weriv
     .                                 WORK1T(ISTIFT(24)), ! centr
     .              WORK1T(ISTIFT(25)),                    ! advem
     .              WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)), ! teini,fpchl
     .              EHISTT)
C
C**** EVALUATES JACOBIAN MATRIX DUE TO CONVECTION (ADVECTIVE) EFFECTS
C
       CALL STIC05T(WORK1T(ISTIFT(26)),WORK1T(ISTIFT(27)), ! CARTD,DVOLU
     .              EHISTT,
     .              WORK1T(ISTIFT(39)),                    ! ELDIS
     .              EPMTXT,
     .              WORK1T(ISTIFT(28)),                    ! GPCOD
     .              PROPST,
     .              WORK1T(ISTIFT(29)),                    ! SHAPE
     .              STRSGT,
     .              WORK1T(ISTIFT(37)),                    ! ESTIF
     .              HSTIFT,
     .              WORK1T(ISTIFT( 1)),WORK1T(ISTIFT( 2)),
     .              WORK1T(ISTIFT( 3)),WORK1T(ISTIFT( 4)),
     .              WORK1T(ISTIFT(19)),                    ! whape
     .              WORK1T(ISTIFT(25)),                    ! advem
     .              WORK1T(ISTIFT(13)),                    ! velcm
     .              WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45))) ! teini,fpchl
C
       CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
     .              WORK1T(ISTIFT(39)),                    ! ELDIS
     .              PROPST,WORK1T(ISTIFT(13)),       1,TEMTLT,
     .              TEMDLT)
C
C**** EVALUATES JACOBIAN MATRIX DUE TO CONVECTION PHASE-CHANGE EFFECTS
C
       MULPHT=MULPDT
       IF(MULPHT.EQ.1) THEN
        CALL STUC05T(WORK1T(ISTIFT(26)),WORK1T(ISTIFT(27)), !CARTD,DVOLU
     .               EHISTT,
     .               WORK1T(ISTIFT(39)),                    !ELDIS
     .               EPMTXT,
     .               WORK1T(ISTIFT(28)),                    !GPCOD
     .               PROPST,
     .               WORK1T(ISTIFT(29)),                    !SHAPE
     .               STRSGT,
     .               WORK1T(ISTIFT(37)),                    !ESTIF
     .               HSTIFT,
     .               WORK1T(ISTIFT( 1)),WORK1T(ISTIFT( 2)),
     .               WORK1T(ISTIFT( 3)),WORK1T(ISTIFT( 4)),
     .               WORK1T(ISTIFT(13)),WORK1T(ISTIFT(14)), !velcm,wsti1
     .               WORK1T(ISTIFT(15)),                    !wsti2
     .               WORK1T(ISTIFT(19)),                    !whape
     .               WORK1T(ISTIFT(25)),                    !advem
     .               WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45))) !teini,fpchl
       ELSE
        CALL STMC05T(WORK1T(ISTIFT(26)),WORK1T(ISTIFT(27)), !CARTD,DVOLU
     .               EHISTT,
     .               WORK1T(ISTIFT(39)),                    !ELDIS
     .               EPMTXT,WORK1T(ISTIFT(28)),             !GPCOD
     .               PROPST,WORK1T(ISTIFT(29)),             !SHAPE
     .               STRSGT,WORK1T(ISTIFT(37)),             !ESTIF
     .               HSTIFT,
     .               WORK1T(ISTIFT(1)), WORK1T(ISTIFT(2)),
     .               WORK1T(ISTIFT(3)), WORK1T(ISTIFT(4)),
     .               WORK1T(ISTIFT(5)), WORK1T(ISTIFT(6)),  !shapd,dvold
     .               WORK1T(ISTIFT(7)), WORK1T(ISTIFT(8)),  !posgp,weigp
     .               WORK1T(ISTIFT(9)),                     !derid
     .               WORK1T(ISTIFT(11)),WORK1T(ISTIFT(12)), !cardd,gpcdd
     .               WORK1T(ISTIFT(13)),WORK1T(ISTIFT(14)), !velcm,wsti1
     .               WORK1T(ISTIFT(15)),WORK1T(ISTIFT(18)), !wsti2,ELCOD
     .               WORK1T(ISTIFT(20)),WORK1T(ISTIFT(23)), !hache,whade
     .               WORK1T(ISTIFT(24)),WORK1T(ISTIFT(25)), !centr,advem
     .               WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)), !teini,fpchl
     .               LNODST)
       ENDIF                        ! mulpdt.eq.1
C
C**** EVALUATES JACOBIAN MATRIX DUE TO IMPROVED THERMO-FLUID ALGORITHM
C
       CALL MARA05T(WORK1T(ISTIFT(27)),                    ! DVOLU
     .              PROPST,
     .              WORK1T(ISTIFT(29)),                    ! SHAPE
     .              WORK1T(ISTIFT(37)),                    ! ESTIF
     .              EHISTT,
     .              WORK1T(ISTIFT(19)),                    ! whape
     .              WORK1T(ISTIFT(39)),                    ! ELDIS
     .              WORK1T(ISTIFT(13)),                    ! velcm
     .              WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)), ! teini,fpchl
     .              WORK1T(ISTIFT(46)))                    ! dvoli
C
C**** EVALUATES JACOBIAN MATRIX DUE TO IMPROVED THERMO-FLUID ALGORITHM
C
       CALL MARX05T(WORK1T(ISTIFT(27)),                    ! DVOLU
     .              PROPST,
     .              WORK1T(ISTIFT(29)),                    ! SHAPE
     .              WORK1T(ISTIFT(37)),                    ! ESTIF
     .              EHISTT,
     .              WORK1T(ISTIFT(19)),                    ! whape
     .              WORK1T(ISTIFT(39)),                    ! ELDIS
     .              WORK1T(ISTIFT(13)),                    ! velcm
     .              WORK1T(ISTIFT(40)),WORK1T(ISTIFT(45)), ! teini,fpchl
     .              WORK1T(ISTIFT(46)),WORK1T(ISTIFT(48))) ! dvoli,EMATX
      ENDIF                         ! iconvt.eq.1
C
C**** ASSEMBLY PROCESS (NMEMO6=1)
C
      CALL ASAL00(ESTIFT,WSTIFT,
     .            WORK1T(ISTIFT(37)),                       ! ESTIF
     .            WORK1T(ISTIFT(38)),                       ! WSTIF
     .            NKOVAT,KDYNAT,NMEMO6)
C
      CALL ASEL05T(CSTIFT,
     .             WORK1T(ISTIFT(37)),                      ! ESTIF
     .             WORK1T(ISTIFT(38)),                      ! WSTIF
     .             WORK1T(ISTIFT(41)))                      ! CSTI1
C
      RETURN
C
C**** EVALUATE INTERNAL RESISTING HEATS
C
   40 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(IFORCT(35)),WORK1T(IFORCT(36)),  ! cartd,dvolu
     .             WORK1T(IFORCT(37)),WORK1T(IFORCT(38)),  ! gpcod,shape
     .             WORK1T(IFORCT(39)),WORK1T(IFORCT(40)),  ! deriv,posgp
     .             WORK1T(IFORCT(41)),WORK1T(IFORCT(42)),  ! weigp,xjacm
     .             WORK1T(IFORCT(16)),WORK1T(IFORCT(43)),  ! elco1,emas1
     .                                WORK1T(IFORCT(44)),  ! eldi1
     .             WORK1T(IFORCT(49)),WORK1T(IFORCT(50)),  ! dvoli,elcoi
     .             WORK1T(IFORCT(13)),                     ! dispt
     .             0,ITASKT)
C
C**** TAKES INTO ACCOUNT THE MICROSTRUCTURAL MODEL
C
      IF(IMICR.EQ.1) THEN
       CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
     .              WORK1T(IFORCT(44)),                    ! ELDIS
     .              PROPST,
     .              WORK1T(IFORCT(12)),
     .                   2,TEMTLT,TEMDLT)
C
       CALL MICR05T(WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), ! CARTD,DVOLU
     .              EHISTT,
     .              WORK1T(IFORCT(16)),                    ! ELCOD
     .              WORK1T(IFORCT(44)),                    ! ELDIS
     .              EPMTXT,
     .              WORK1T(IFORCT(37)),                    ! GPCOD
     .              LNODST,PROPST,RMAT1T,
     .              WORK1T(IFORCT(38)),                    ! SHAPE
     .              STRANT,STRSGT,STRA0T,STRS0T,
     .              WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .              WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .              WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .              WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .              WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .              WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .              WORK1T(IFORCT(33)),WORK1T(IFORCT(34)), ! vel1m,advem
     .              WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) ! teini,fpchl
      ENDIF                          ! imicr.eq.1
C
C**** EVALUATES WEIGHTING FUNCTIONS & THEIR CARTESIAN DERIVATIVES
C
C     Note: for microstructural problems, these weighting functions
C           depend on microstructural variables.
C
      CALL CART05T(WORK1T(IFORCT(16)),                    ! ELCOD
     .             PROPST,
     .             WORK1T(IFORCT(44)),                    ! ELDIS
     .             WORK1T(IFORCT(38)),WORK1T(IFORCT(35)), ! SHAPE,CARTD
     .             WORK1T(IFORCT(19)),WORK1T(IFORCT(20)), ! posgp,weigp
     .             WORK1T(IFORCT(24)),WORK1T(IFORCT(17)), ! gpcdd,shapd
     .             WORK1T(IFORCT(21)),WORK1T(IFORCT(22)), ! derid,xjacm
     .             WORK1T(IFORCT(23)),                    ! cardd
     .             WORK1T(IFORCT(12)),                    ! velcm
     .             WORK1T(IFORCT(27)),WORK1T(IFORCT(28)), ! whape,hache
     .             WORK1T(IFORCT(29)),WORK1T(IFORCT(30)), ! wartd,weriv
     .                                WORK1T(IFORCT(32)), ! centr
     .                                WORK1T(IFORCT(34)), ! advem
     .             WORK1T(IFORCT(45)),WORK1T(IFORCT(48)), ! teini,fpchl
     .             EHISTT)
C
C**** EVALUATES INTERNAL RESISTING HEAT DUE TO CONDUCTIVITY
C
      CALL FRIN05T(WORK1T(IFORCT(35)),WORK1T(IFORCT(36)),  ! CARTD,DVOLU
     .             EHISTT,
     .             WORK1T(IFORCT(16)),                     ! ELCOD
     .             WORK1T(IFORCT(44)),                     ! ELDIS
     .             EPMTXT,
     .             WORK1T(IFORCT(37)),                     ! GPCOD
     .             LNODST,PROPST,RMAT1T,
     .             WORK1T(IFORCT(38)),                     ! SHAPE
     .             STRANT,STRSGT,STRA0T,STRS0T,
     .             WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .             WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .             WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .             WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .             WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .             WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .             WORK1T(IFORCT(29)),                     ! wartd
     .             WORK1T(IFORCT(45)),WORK1T(IFORCT(48)),  ! teini,fpchl
     .             WORK1T(IFORCT(49)))                     ! dvoli
C
C**** COMPUTES NODAL BOUNDARY CHANGES DUE TO DENSITY VARIATION
C
      CALL BOUC05T(WORK1T(IFORCT(45)),                     ! TEINI
     .             WORK1T(IFORCT(44)),                     ! ELDIS
     .             WORK1T(IFORCT(46)))                     ! BOUCH
C
C**** EVALUATES INTERNAL RESISTING HEATS DUE TO CONVECTION EFFECTS
C
      IF(ICONVT.EQ.1) THEN
       CALL FRIC05T(WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), ! CARTD,DVOLU
     .              EHISTT,
     .              WORK1T(IFORCT(16)),                    ! ELCOD
     .              WORK1T(IFORCT(44)),                    ! ELDIS
     .              EPMTXT,
     .              WORK1T(IFORCT(37)),                    ! GPCOD
     .              LNODST,PROPST,RMAT1T,
     .              WORK1T(IFORCT(38)),                    ! SHAPE
     .              STRANT,STRSGT,STRA0T,STRS0T,
     .              WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .              WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .              WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .              WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .              WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .              WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .              WORK1T(IFORCT(27)),                          ! whape
     .                                 WORK1T(IFORCT(34)),       ! advem
     .              WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) ! teini,fpchl
C
       CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
     .              WORK1T(IFORCT(44)),                   ! ELDIS
     .              PROPST,
     .              WORK1T(IFORCT(12)),
     .                   2,TEMTLT,TEMDLT)
C
       MULPHT=MULPDT
       NPHCHT=NPHCDT
       NISOTT=NISDTT
       TEMPLT=TEMDLT
C
       IF(IEGFPC.EQ.0) THEN
        IF(MULPHT.EQ.1) THEN
         CALL FRUC05T(
     .               WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), !CARTD,DVOLU
     .               EHISTT,
     .               WORK1T(IFORCT(16)),                    !ELCOD
     .               WORK1T(IFORCT(44)),                    !ELDIS
     .               EPMTXT,
     .               WORK1T(IFORCT(37)),                    !GPCOD
     .               LNODST,PROPST,RMAT1T,
     .               WORK1T(IFORCT(38)),                    !SHAPE
     .               STRANT,STRSGT,STRA0T,STRS0T,
     .               WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .               WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .               WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .               WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .               WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .               WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .               WORK1T(IFORCT(27)),                         ! whape
     .                                  WORK1T(IFORCT(34)),      ! advem
     .               WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) !teini,fpchl
        ELSE
         IF(NISOTT.EQ.1)
     .    CALL FRMC05T(                         ! it should be frtc05t.f
     .               WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), !CARTD,DVOLU
     .               EHISTT,
     .               WORK1T(IFORCT(16)),WORK1T(IFORCT(44)), !ELCOD,ELDIS
     .               EPMTXT,WORK1T(IFORCT(37)),             !GPCOD
     .               LNODST,PROPST,RMAT1T,WORK1T(IFORCT(38)), !SHAPE
     .               STRANT,STRSGT,STRA0T,STRS0T,
     .               WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .               WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .               WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .               WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .               WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .               WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .               WORK1T(IFORCT(17)),WORK1T(IFORCT(18)), !shapd,dvold
     .               WORK1T(IFORCT(19)),WORK1T(IFORCT(20)), !posgp,weigp
     .               WORK1T(IFORCT(21)),                    !derid
     .               WORK1T(IFORCT(23)),WORK1T(IFORCT(24)), !cardd,gpcdd
     .               WORK1T(IFORCT(25)),WORK1T(IFORCT(26)), !elel1,elel2
     .               WORK1T(IFORCT(28)),WORK1T(IFORCT(31)), !hache,whade
     .               WORK1T(IFORCT(32)),WORK1T(IFORCT(34)), !centr,advem
     .               WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) !teini,fpchl
C
         IF(NISOTT.EQ.2) 
     .    CALL FRMC05T(
     .               WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), !CARTD,DVOLU
     .               EHISTT,
     .               WORK1T(IFORCT(16)),WORK1T(IFORCT(44)), !ELCOD,ELDIS
     .               EPMTXT,WORK1T(IFORCT(37)),             !GPCOD
     .               LNODST,PROPST,RMAT1T,WORK1T(IFORCT(38)), !SHAPE
     .               STRANT,STRSGT,STRA0T,STRS0T,
     .               WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .               WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .               WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .               WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .               WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .               WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .               WORK1T(IFORCT(17)),WORK1T(IFORCT(18)), !shapd,dvold
     .               WORK1T(IFORCT(19)),WORK1T(IFORCT(20)), !posgp,weigp
     .               WORK1T(IFORCT(21)),                    !derid
     .               WORK1T(IFORCT(23)),WORK1T(IFORCT(24)), !cardd,gpcdd
     .               WORK1T(IFORCT(25)),WORK1T(IFORCT(26)), !elel1,elel2
     .               WORK1T(IFORCT(28)),WORK1T(IFORCT(31)), !hache,whade
     .               WORK1T(IFORCT(32)),WORK1T(IFORCT(34)), !centr,advem
     .               WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) !teini,fpchl
        ENDIF                        ! mulpht.eq.1
       ENDIF                         ! iegfpc.eq.0
       IF(IEGFPC.EQ.1.AND.NMEMO3.EQ.0) THEN
        CALL FPCN05T(WORK1T(IFORCT(44)),                    ! ELDIS
     .               WORK1T(IFORCT(12)),                    ! not used
     .               PROPST,                                ! not used
     .               WORK1T(IFORCT(48)),     0,     1)      ! FPCHL
C
        IF(MULPHT.EQ.1) THEN
         CALL FRUX05T(
     .               WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), !CARTD,DVOLU
     .               EHISTT,
     .               WORK1T(IFORCT(16)),                    !ELCOD
     .               WORK1T(IFORCT(44)),                    !ELDIS
     .               EPMTXT,
     .               WORK1T(IFORCT(37)),                    !GPCOD
     .               LNODST,PROPST,RMAT1T,
     .               WORK1T(IFORCT(38)),                    !SHAPE
     .               STRANT,STRSGT,STRA0T,STRS0T,
     .               WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .               WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .               WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .               WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .               WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .               WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .               WORK1T(IFORCT(27)),                         ! whape
     .                                  WORK1T(IFORCT(34)),      ! advem
     .               WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) !teini,fpchl
        ELSE
         IF(NISOTT.EQ.1)
     .    CALL FRMX05T(                         ! it should be frtx05t.f
     .               WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), !CARTD,DVOLU
     .               EHISTT,
     .               WORK1T(IFORCT(16)),WORK1T(IFORCT(44)), !ELCOD,ELDIS
     .               EPMTXT,WORK1T(IFORCT(37)),             !GPCOD
     .               LNODST,PROPST,RMAT1T,WORK1T(IFORCT(38)), !SHAPE
     .               STRANT,STRSGT,STRA0T,STRS0T,
     .               WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .               WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .               WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .               WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .               WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .               WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .               WORK1T(IFORCT(17)),WORK1T(IFORCT(18)), !shapd,dvold
     .               WORK1T(IFORCT(19)),WORK1T(IFORCT(20)), !posgp,weigp
     .               WORK1T(IFORCT(21)),                    !derid
     .               WORK1T(IFORCT(23)),WORK1T(IFORCT(24)), !cardd,gpcdd
     .               WORK1T(IFORCT(25)),WORK1T(IFORCT(26)), !elel1,elel2
     .               WORK1T(IFORCT(28)),WORK1T(IFORCT(31)), !hache,whade
     .               WORK1T(IFORCT(32)),WORK1T(IFORCT(34)), !centr,advem
     .               WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) !teini,fpchl
C
         IF(NISOTT.EQ.2) 
     .    CALL FRMX05T(
     .               WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), !CARTD,DVOLU
     .               EHISTT,
     .               WORK1T(IFORCT(16)),WORK1T(IFORCT(44)), !ELCOD,ELDIS
     .               EPMTXT,WORK1T(IFORCT(37)),             !GPCOD
     .               LNODST,PROPST,RMAT1T,WORK1T(IFORCT(38)), !SHAPE
     .               STRANT,STRSGT,STRA0T,STRS0T,
     .               WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .               WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .               WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .               WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .               WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .               WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .               WORK1T(IFORCT(17)),WORK1T(IFORCT(18)), !shapd,dvold
     .               WORK1T(IFORCT(19)),WORK1T(IFORCT(20)), !posgp,weigp
     .               WORK1T(IFORCT(21)),                    !derid
     .               WORK1T(IFORCT(23)),WORK1T(IFORCT(24)), !cardd,gpcdd
     .               WORK1T(IFORCT(25)),WORK1T(IFORCT(26)), !elel1,elel2
     .               WORK1T(IFORCT(28)),WORK1T(IFORCT(31)), !hache,whade
     .               WORK1T(IFORCT(32)),WORK1T(IFORCT(34)), !centr,advem
     .               WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) !teini,fpchl
        ENDIF                        ! mulpht.eq.1
       ENDIF                         ! iegfpc.eq.1.and.nmemo3.eq.0
C
       IF(KDYNAT.EQ.0)
     .  CALL FPCN05T(WORK1T(IFORCT(44)),                    ! ELDIS
     .               WORK1T(IFORCT(12)),                    ! not used
     .               PROPST,                                ! not used
     .               WORK1T(IFORCT(48)),     0,     0)      ! FPCHL
C
C**** EVALUATES INTERNAL HEAT FOR IMPROVED THERMO-FLUID ALGORITHM
C
       CALL FRRA05T(PROPST,
     .              WORK1T(IFORCT(12)),WORK1T(IFORCT(1)),  ! VELCM,ELELM
     .              WORK1T(IFORCT(36)),                    ! DVOLU
     .              WORK1T(IFORCT(38)),                    ! SHAPE
     .              EHISTT,
     .              WORK1T(IFORCT(44)),                    ! ELDIS
     .                                 WORK1T(IFORCT(27)), ! whape
     .              WORK1T(IFORCT(45)),WORK1T(IFORCT(48)), ! teini,fpchl
     .              WORK1T(IFORCT(49)),WORK1T(IFORCT(34)), ! dvoli,advem
     .              WORK1T(IFORCT(35)),WORK1T(IFORCT(10))) ! cartd,xjacm
C
C**** EVALUATES INTERNAL HEAT FOR IMPROVED THERMO-FLUID ALGORITHM
C
       CALL FRRX05T(PROPST,
     .              WORK1T(IFORCT(12)),WORK1T(IFORCT(1)),  ! VELCM,ELELM
     .              WORK1T(IFORCT(36)),                    ! DVOLU
     .              WORK1T(IFORCT(38)),                    ! SHAPE
     .              EHISTT,
     .              WORK1T(IFORCT(44)),                    ! ELDIS
     .                                 WORK1T(IFORCT(27)), ! whape
     .              WORK1T(IFORCT(45)),WORK1T(IFORCT(48)), ! teini,fpchl
     .              WORK1T(IFORCT(49)),WORK1T(IFORCT(55))) ! dvoli,eldiq
      ENDIF                         ! iconvt.eq.1
C
C**** EVALUATE EQUIVALENT VOLUME HEAT (NOT CONSTANT, i.e. AS INT. HEAT)
C
      IF(ILDGRT.EQ.1)
     . CALL LDGR05T(WORK1T(IFORCT(36)),       ! DVOLU
     .              PROPST,
     .              WORK1T(IFORCT(38)),       ! SHAPE
     .              WORK1T(IFORCT( 1)),
     .              WORK1T(IFORCT(44)),       ! ELDIS
     .              WORK1T(IFORCT(12)),       ! VELCM
     .              EHISTT,ILDGRT,
     .              WORK1T(IFORCT(27)),WORK1T(IFORCT(45)), ! whape,teini
     .              WORK1T(IFORCT(48)),WORK1T(IFORCT(49)), ! fpchl,dvoli
     .              WORK1T(IFORCT(34)),WORK1T(IFORCT(35)), ! advem,cartd
     .              WORK1T(IFORCT(10)))                    ! XJACM
      RETURN
C
C**** EVALUATE INTERNAL RESISTING HEATS (ALTERNATIVE OPTION FOR MICRO
C     ADVECTIVE EFFECTS, THAT IS, IEGFPC=1 & NMEMO3=1; see forcint.f)
C
   50 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(IFORCT(35)),WORK1T(IFORCT(36)),  ! cartd,dvolu
     .             WORK1T(IFORCT(37)),WORK1T(IFORCT(38)),  ! gpcod,shape
     .             WORK1T(IFORCT(39)),WORK1T(IFORCT(40)),  ! deriv,posgp
     .             WORK1T(IFORCT(41)),WORK1T(IFORCT(42)),  ! weigp,xjacm
     .             WORK1T(IFORCT(16)),WORK1T(IFORCT(43)),  ! elco1,emas1
     .                                WORK1T(IFORCT(44)),  ! eldi1
     .             WORK1T(IFORCT(49)),WORK1T(IFORCT(50)),  ! dvoli,elcoi
     .             WORK1T(IFORCT(13)),                     ! dispt
     .             0,ITASKT)
C
C**** EVALUATES WEIGHTING FUNCTIONS & THEIR CARTESIAN DERIVATIVES
C
      CALL CART05T(WORK1T(IFORCT(16)),                    ! ELCOD
     .             PROPST,
     .             WORK1T(IFORCT(44)),                    ! ELDIS
     .             WORK1T(IFORCT(38)),WORK1T(IFORCT(35)), ! SHAPE,CARTD
     .             WORK1T(IFORCT(19)),WORK1T(IFORCT(20)), ! posgp,weigp
     .             WORK1T(IFORCT(24)),WORK1T(IFORCT(17)), ! gpcdd,shapd
     .             WORK1T(IFORCT(21)),WORK1T(IFORCT(22)), ! derid,xjacm
     .             WORK1T(IFORCT(23)),                    ! cardd
     .             WORK1T(IFORCT(12)),                    ! velcm
     .             WORK1T(IFORCT(27)),WORK1T(IFORCT(28)), ! whape,hache
     .             WORK1T(IFORCT(29)),WORK1T(IFORCT(30)), ! wartd,weriv
     .                                WORK1T(IFORCT(32)), ! centr
     .                                WORK1T(IFORCT(34)), ! advem
     .             WORK1T(IFORCT(45)),WORK1T(IFORCT(48)), ! teini,fpchl
     .             EHISTT)
C
C**** EVALUATES INTERNAL RESISTING HEATS DUE TO PC CONVECTION EFFECTS
C
      CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
     .             WORK1T(IFORCT(44)),                    ! ELDIS
     .             PROPST,
     .             WORK1T(IFORCT(12)),
     .                  2,TEMTLT,TEMDLT)
C
      CALL FRUX05T(
     .               WORK1T(IFORCT(35)),WORK1T(IFORCT(36)), !CARTD,DVOLU
     .               EHISTT,
     .               WORK1T(IFORCT(16)),                    !ELCOD
     .               WORK1T(IFORCT(44)),                    !ELDIS
     .               EPMTXT,
     .               WORK1T(IFORCT(37)),                    !GPCOD
     .               LNODST,PROPST,RMAT1T,
     .               WORK1T(IFORCT(38)),                    !SHAPE
     .               STRANT,STRSGT,STRA0T,STRS0T,
     .               WORK1T(IFORCT( 1)),WORK1T(IFORCT( 2)),
     .               WORK1T(IFORCT( 3)),WORK1T(IFORCT( 4)),
     .               WORK1T(IFORCT( 5)),WORK1T(IFORCT( 6)),
     .               WORK1T(IFORCT( 7)),WORK1T(IFORCT( 8)),
     .               WORK1T(IFORCT( 9)),WORK1T(IFORCT(10)),
     .               WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .               WORK1T(IFORCT(27)),                         ! whape
     .                                  WORK1T(IFORCT(34)),      ! advem
     .               WORK1T(IFORCT(45)),WORK1T(IFORCT(48))) !teini,fpchl
      RETURN
C
C**** EVALUATE INTERNAL TRANSIENT RESISTING HEATS
C
   60 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATT=2
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(IFORDT(23)),WORK1T(IFORDT(24)),  ! cartd,dvolu
     .             WORK1T(IFORDT(25)),WORK1T(IFORDT(26)),  ! gpcod,shape
     .             WORK1T(IFORDT(27)),WORK1T(IFORDT(28)),  ! deriv,posgp
     .             WORK1T(IFORDT(29)),WORK1T(IFORDT(30)),  ! weigp,xjacm
     .             WORK1T(IFORDT(17)),WORK1T(IFORDT(31)),  ! elco1,emas1
     .                                WORK1T(IFORDT(32)),  ! eldi1
     .             WORK1T(IFORDT(35)),WORK1T(IFORDT(36)),  ! dvoli,elcoi
     .             WORK1T(IFORDT(37)),                     ! dispt
     .             0,ITASKT)
C
C**** EVALUATES WEIGHTING FUNCTIONS
C
      CALL PETR05T(WORK1T(IFORDT(17)),                     ! ELCOD
     .             PROPST,
     .             WORK1T(IFORDT(32)),                     ! ELDIS
     .             WORK1T(IFORDT(26)),WORK1T(IFORDT(23)),  ! SHAPE,CARTD
     .             WORK1T(IFORDT( 9)),WORK1T(IFORDT(10)),  ! posgp,weigp
     .             WORK1T(IFORDT(14)),WORK1T(IFORDT( 7)),  ! gpcod,shapd
     .             WORK1T(IFORDT(11)),WORK1T(IFORDT(12)),  ! derid,xjacm
     .             WORK1T(IFORDT(13)),                     ! cardd
     .             WORK1T(IFORDT( 2)),                     ! velcm
     .             WORK1T(IFORDT(18)),WORK1T(IFORDT(19)),  ! whape,hache
     .             WORK1T(IFORDT(11)),                     ! werid=derid
     .             WORK1T(IFORDT(21)),WORK1T(IFORDT(22)),  ! centr,advem
     .             WORK1T(IFORDT(33)),WORK1T(IFORDT(34)),  ! teini,fpchl
     .             EHISTT)
C
C**** EVALUATES INTERNAL TRANSIENT HEAT DUE TO SPECIFIC HEAT
C
      CALL FRDY05T(PROPST,WORK1T(IFORDT(2)),WORK1T(IFORDT(1)),
     .             WORK1T(IFORDT(24)),                     ! DVOLU
     .             WORK1T(IFORDT(26)),                     ! SHAPE
     .             EHISTT,
     .             WORK1T(IFORDT(32)),                     ! ELDIS
     .             WORK1T(IFORDT( 3)),WORK1T(IFORDT( 6)),
     .                                WORK1T(IFORDT(18)),  ! whape
     .             WORK1T(IFORDT(33)),WORK1T(IFORDT(34)),  ! teini,fpchl
     .             WORK1T(IFORDT(35)))                     ! dvoli
C
C**** EVALUATE PHASE CHANGE INTERNAL DYNAMIC RESISTING HEATS
C
      CALL CERO00T(WORK1T(IFORDT(16)),WORK1T(IFORDT(15)),NEVABT)
      CALL IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,NISDTT,
     .             WORK1T(IFORDT(32)),                     ! ELDIS
     .             PROPST,WORK1T(IFORDT(2)), 2,TEMTLT,
     .             TEMDLT)
C
      DO INPCCT=1,2
       IF(INPCCT.EQ.1) THEN
        MULPHT=MULPTT
        NPHCHT=NPHCTT
        NISOTT=NISTTT
        TEMPLT=TEMTLT
       ELSE
        MULPHT=MULPDT
        NPHCHT=NPHCDT
        NISOTT=NISDTT
        TEMPLT=TEMDLT
       ENDIF
C
       IF(MULPHT.EQ.1) THEN
        CALL FRUF05T(PROPST,WORK1T(IFORDT(2)),WORK1T(IFORDT(1)),
     .               WORK1T(IFORDT(24)),                    ! DVOLU
     .               WORK1T(IFORDT(26)),                    ! SHAPE
     .               EHISTT,
     .               WORK1T(IFORDT(32)),                    ! ELDIS
     .               WORK1T(IFORDT( 3)),
     .                                  WORK1T(IFORDT( 6)),
     .               WORK1T(IFORDT(15)),WORK1T(IFORDT(16)),
     .                                  WORK1T(IFORDT(18)), ! whape
     .               INPCCT,
     .              WORK1T(IFORDT(33)),WORK1T(IFORDT(34)),  !teini,fpchl
     .              WORK1T(IFORDT(35)))                     !dvoli
       ELSE
        IF(NISOTT.EQ.1) 
     .   CALL FRTF05T(PROPST,WORK1T(IFORDT(2)),WORK1T(IFORDT(1)),
     .                WORK1T(IFORDT(24)),                        ! DVOLU
     .                WORK1T(IFORDT(26)),                        ! SHAPE
     .                EHISTT,
     .                WORK1T(IFORDT(32)),                        ! ELDIS
     .                WORK1T(IFORDT( 3)),
     .                WORK1T(IFORDT( 6)),WORK1T(IFORDT( 7)),
     .                WORK1T(IFORDT( 8)),
     .                WORK1T(IFORDT(25)),                        ! GPCOD
     .                NPHCHT,
     .                WORK1T(IFORDT( 9)),
     .                WORK1T(IFORDT(17)),                        ! ELCOD
     .                WORK1T(IFORDT(10)),
     .                WORK1T(IFORDT(11)),WORK1T(IFORDT(12)),
     .                WORK1T(IFORDT(13)),WORK1T(IFORDT(14)),
     .                WORK1T(IFORDT(15)),WORK1T(IFORDT(16)),
     .               WORK1T(IFORDT(19)),WORK1T(IFORDT(20)), !hache,whade
     .               WORK1T(IFORDT(21)),WORK1T(IFORDT(22)), !centr,advem
     .               WORK1T(IFORDT(33)),WORK1T(IFORDT(34)), !teini,fpchl
     .               WORK1T(IFORDT(35)),                    !dvoli
     .                LNODST,INPCCT,TEMPLT)
        IF(NISOTT.EQ.2) 
     .   CALL FRMF05T(PROPST,WORK1T(IFORDT(2)),WORK1T(IFORDT(1)),
     .                WORK1T(IFORDT(24)),                        ! DVOLU
     .                WORK1T(IFORDT(26)),                        ! SHAPE
     .                EHISTT,
     .                WORK1T(IFORDT(32)),                        ! ELDIS
     .                WORK1T(IFORDT(3)),
     .                WORK1T(IFORDT(6)),
     .                WORK1T(IFORDT(7)),WORK1T(IFORDT(8)),
     .                WORK1T(IFORDT(25)),                        ! GPCOD
     .                WORK1T(IFORDT(9)),
     .                WORK1T(IFORDT(17)),                        ! ELCOD
     .                WORK1T(IFORDT(10)),
     .                WORK1T(IFORDT(11)),WORK1T(IFORDT(12)),
     .                WORK1T(IFORDT(13)),WORK1T(IFORDT(14)),
     .                WORK1T(IFORDT(15)),WORK1T(IFORDT(16)),
     .               WORK1T(IFORDT(19)),WORK1T(IFORDT(20)), !hache,whade
     .               WORK1T(IFORDT(21)),WORK1T(IFORDT(22)), !centr,advem
     .               WORK1T(IFORDT(33)),WORK1T(IFORDT(34)), !teini,fpchl
     .               WORK1T(IFORDT(35)),                    !dvoli
     .                LNODST,INPCCT)
       ENDIF                   ! mulpht.eq.1
      ENDDO                    ! inpcct=1,2
C
      CALL SUMMATT(WORK1T(IFORDT(1)),WORK1T(IFORDT(16)),
     .             WORK1T(IFORDT(15)),NEVABT)
C
C**** COMPUTES NODAL PHASE-CHANGE FUNCTION AT TIME t+dt
C     (useful to postprocess & to compute in the mechanical problem)
C
      CALL FPCN05T(WORK1T(IFORDT(32)),                     ! ELDIS
     .             WORK1T(IFORDT(2)),                      ! VELCM
     .             PROPST,                                 ! not used
     .             WORK1T(IFORDT(34)),     1,     0)       ! FPCHL
      RETURN
C
C**** EVALUATE EQUIVALENT VOLUME HEAT (CONSTANT !!!, i.e. AS EXT. HEAT)
C
   70 CONTINUE
      IF(ILDGRT.EQ.1) RETURN
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(ILOSUT(21)),WORK1T(ILOSUT(22)),  ! cartd,dvolu
     .             WORK1T(ILOSUT(23)),WORK1T(ILOSUT(24)),  ! gpcod,shape
     .             WORK1T(ILOSUT(25)),WORK1T(ILOSUT(26)),  ! deriv,posgp
     .             WORK1T(ILOSUT(27)),WORK1T(ILOSUT(28)),  ! weigp,xjacm
     .             WORK1T(ILOSUT(17)),WORK1T(ILOSUT(29)),  ! elco1,emas1
     .                                WORK1T(ILOSUT(31)),  ! eldi1
     .             WORK1T(ILOSUT(35)),WORK1T(ILOSUT(36)),  ! dvoli,elcoi
     .             WORK1T(ILOSUT(37)),                     ! dispt
     .             0,ITASKT)
C
C**** EVALUATES WEIGHTING FUNCTIONS
C
      CALL PETR05T(WORK1T(ILOSUT(17)),                     ! ELCOD
     .             PROPST,
     .             WORK1T(ILOSUT(31)),       ! TEINI=ELDIS for this case
     .             WORK1T(ILOSUT(24)),WORK1T(ILOSUT(21)),  ! SHAPE,CARTD
     .             WORK1T(ILOSUT(26)),WORK1T(ILOSUT(27)),  ! posgp,weigp
     .             WORK1T(ILOSUT( 4)),WORK1T(ILOSUT(10)),  ! gpcod,shapd
     .             WORK1T(ILOSUT( 2)),WORK1T(ILOSUT(13)),  ! derid,xjacm
     .             WORK1T(ILOSUT(14)),                     ! cardd
     .             WORK1T(ILOSUT(16)),                     ! velcm
     .             WORK1T(ILOSUT(18)),WORK1T(ILOSUT(19)),  ! whape,hache
     .             WORK1T(ILOSUT(25)),                     ! weriv=derid
     .             WORK1T(ILOSUT(20)),WORK1T(ILOSUT(32)),  ! centr,advem
     .             WORK1T(ILOSUT(31)),WORK1T(ILOSUT(34)),  ! teini,fpchl
     .             EHISTT)
C
C**** EVALUATE HEAT
C
      CALL LDGR05T(WORK1T(ILOSUT(22)),                    ! DVOLU
     .             PROPST,
     .             WORK1T(ILOSUT(24)),                    ! SHAPE
     .             WORK1T(ILOSUT( 1)),
     .             WORK1T(ILOSUT(31)),      ! TEINI=ELDIS for this case
     .             WORK1T(ILOSUT(16)),                    ! VELCM
     .             EHISTT,ILDGRT,
     .             WORK1T(ILOSUT(18)),WORK1T(ILOSUT(31)), ! whape,teini
     .             WORK1T(ILOSUT(34)),WORK1T(ILOSUT(35)), ! fpchl,dvoli
     .             WORK1T(ILOSUT(32)),WORK1T(ILOSUT(21)), ! advem,cartd
     .             WORK1T(IFORCT(13)))                    ! XJACM
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
      CALL CHEK01T(NDIMET,IELEMT,NNODLT,NRULET,NGAULT,     5,LUREST)
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
      CALL OUTG05T(EHISTT,STRANT,STRSGT)
      RETURN
C
C**** OUTPUT NODAL (SMOOTHED) HEAT FLUXES
C
  130 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      IF(ISMO1T.EQ.1) THEN
       NRUAUX=NRULET
       NRULET=5
      ENDIF
      CALL SMONOD(NDIMET,NNODLT,NNODST,NQUTRT)
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(IGSMOT(11)),WORK1T(IGSMOT(12)),  ! cartd,dvolu
     .             WORK1T(IGSMOT(13)),WORK1T(IGSMOT(14)),  ! gpcod,shape
     .             WORK1T(IGSMOT(15)),WORK1T(IGSMOT(16)),  ! deriv,posgp
     .             WORK1T(IGSMOT(17)),WORK1T(IGSMOT(18)),  ! weigp,xjacm
     .             WORK1T(IGSMOT(19)),WORK1T(IGSMOT(20)),  ! elco1,emas1
     .                                WORK1T(IGSMOT(21)),  ! eldi1
     .             WORK1T(IGSMOT(28)),WORK1T(IGSMOT(29)),  ! dvoli,elcoi
     .             WORK1T(IGSMOT(30)),                     ! dispt
     .             0,ITASKT)
C
      CALL SMOG05T(WORK1T(IGSMOT(12)),                     ! DVOLU
     .             WORK1T(IGSMOT(20)),                     ! EMASS
     .             WORK1T(IGSMOT(14)),                     ! SHAPE
     .             STRSGT,PROPST,
     .             WORK1T(IGSMOT( 5)),WORK1T(IGSMOT(4)))   ! STREB,STREA
      IF(ISMO1T.EQ.1) NRULET=NRUAUX
      RETURN
C
C**** OUTPUT NODAL (SMOOTHED) TEMP.GRADIENTS
C
  140 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      IF(ISMO1T.EQ.1) THEN
       NRUAUX=NRULET
       NRULET=5
      ENDIF
      CALL SMONOD(NDIMET,NNODLT,NNODST,NQUTRT)
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(IGSMOT(11)),WORK1T(IGSMOT(12)),  ! cartd,dvolu
     .             WORK1T(IGSMOT(13)),WORK1T(IGSMOT(14)),  ! gpcod,shape
     .             WORK1T(IGSMOT(15)),WORK1T(IGSMOT(16)),  ! deriv,posgp
     .             WORK1T(IGSMOT(17)),WORK1T(IGSMOT(18)),  ! weigp,xjacm
     .             WORK1T(IGSMOT(19)),WORK1T(IGSMOT(20)),  ! elco1,emas1
     .                                WORK1T(IGSMOT(21)),  ! eldi1
     .             WORK1T(IGSMOT(28)),WORK1T(IGSMOT(29)),  ! dvoli,elcoi
     .             WORK1T(IGSMOT(30)),                     ! dispt
     .             0,ITASKT)
C
      CALL SMOF05T(WORK1T(IGSMOT(12)),                     ! DVOLU
     .             WORK1T(IGSMOT(20)),                     ! EMASS
     .             WORK1T(IGSMOT(14)),                     ! SHAPE
     .             STRANT,PROPST,
     .             WORK1T(IGSMOT( 7)),WORK1T(IGSMOT(6)))   ! SFISB,SFISA
      IF(ISMO1T.EQ.1) NRULET=NRUAUX
      RETURN
C
C**** OUTPUT NODAL (SMOOTHED) INTERNAL VARIABLES
C
  150 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      IF(ISMO1T.EQ.1) THEN
       NRUAUX=NRULET
       NRULET=5
      ENDIF
      CALL SMONOD(NDIMET,NNODLT,NNODST,NQUTRT)
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(IGSMOT(11)),WORK1T(IGSMOT(12)),  ! cartd,dvolu
     .             WORK1T(IGSMOT(13)),WORK1T(IGSMOT(14)),  ! gpcod,shape
     .             WORK1T(IGSMOT(15)),WORK1T(IGSMOT(16)),  ! deriv,posgp
     .             WORK1T(IGSMOT(17)),WORK1T(IGSMOT(18)),  ! weigp,xjacm
     .             WORK1T(IGSMOT(19)),WORK1T(IGSMOT(20)),  ! elco1,emas1
     .                                WORK1T(IGSMOT(21)),  ! eldi1
     .             WORK1T(IGSMOT(28)),WORK1T(IGSMOT(29)),  ! dvoli,elcoi
     .             WORK1T(IGSMOT(30)),                     ! dispt
     .             0,ITASKT)
C
      CALL SMOI05T(WORK1T(IGSMOT(12)),                     ! DVOLU
     .             WORK1T(IGSMOT(13)),                     ! gpcod
     .             WORK1T(IGSMOT(20)),                     ! EMASS
     .             WORK1T(IGSMOT(14)),                     ! SHAPE
     .             EHISTT,PROPST,
     .             WORK1T(IGSMOT(10)),WORK1T(IGSMOT(9)))   ! SFIPB,SFIPA
      IF(ISMO1T.EQ.1) NRULET=NRUAUX
      RETURN
C
C**** OUTPOST
C
  160 CONTINUE
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(ISTART( 3)),WORK1T(ISTART( 4)),  ! cartd,dvolu
     .             WORK1T(ISTART( 5)),WORK1T(ISTART( 6)),  ! gpcod,shape
     .             WORK1T(ISTART( 7)),WORK1T(ISTART( 8)),  ! deriv,posgp
     .             WORK1T(ISTART( 9)),WORK1T(ISTART(10)),  ! weigp,xjacm
     .             WORK1T(ISTART(11)),WORK1T(ISTART(12)),  ! elco1,emas1
     .                                WORK1T(ISTART(13)),  ! eldi1
     .             WORK1T(ISTART(15)),WORK1T(ISTART(16)),  ! dvoli,elcoi
     .             WORK1T(ISTART(17)),                     ! dispt
     .             0,ITASKT)
      RETURN
C
C**** OUTPUT NODAL (SMOOTHED) POROSITY CRITERIA
C
  170 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
C     Note: NRULET is always 5 because the porosity criteria are defined
C           and stored at nodes (not at Gauss points).
C
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      NRUAUX=NRULET
      NRULET=5
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(IGSMOT(11)),WORK1T(IGSMOT(12)),  ! cartd,dvolu
     .             WORK1T(IGSMOT(13)),WORK1T(IGSMOT(14)),  ! gpcod,shape
     .             WORK1T(IGSMOT(15)),WORK1T(IGSMOT(16)),  ! deriv,posgp
     .             WORK1T(IGSMOT(17)),WORK1T(IGSMOT(18)),  ! weigp,xjacm
     .             WORK1T(IGSMOT(19)),WORK1T(IGSMOT(20)),  ! elco1,emas1
     .                                WORK1T(IGSMOT(21)),  ! eldi1
     .             WORK1T(IGSMOT(28)),WORK1T(IGSMOT(29)),  ! dvoli,elcoi
     .             WORK1T(IGSMOT(30)),                     ! dispt
     .             0,ITASKT)
C
      CALL PORO05T(WORK1T(IGSMOT(11)),WORK1T(IGSMOT(12)),  ! CARTD,DVOLU
     .             WORK1T(IGSMOT(20)),                     ! EMASS
     .             WORK1T(IGSMOT(14)),                     ! SHAPE
     .             PROPST,
     .             WORK1T(IGSMOT(24)),WORK1T(IGSMOT(23)),  ! SPIPB,SPIPA
     .             WORK1T(IGSMOT(25)),                     ! SPIPC
     .             WORK1T(IGSMOT(21)),WORK1T(IGSMOT(26)),  ! ELDI1,FPCHL
     .             WORK1T(IGSMOT(27)))                     ! VELCM
      NRULET=NRUAUX
      RETURN
C
C**** COMPUTES L*f_pc & ITS RATE (ALTERNATIVE OPTION FOR MICRO
C     ADVECTIVE EFFECTS, THAT IS, IEGFPC=1 & NMEMO3=1; see forcint.f)
C
  190 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      NEGATT=1
      IF(IMICR.EQ.1) NEGATT=2
      IF(ISMO1T.EQ.1) THEN
       NRUAUX=NRULET
       NRULET=5
      ENDIF
      CALL SMONOD(NDIMET,NNODLT,NNODST,NQUTRT)
      CALL SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(IFORCT(35)),WORK1T(IFORCT(36)),  ! cartd,dvolu
     .             WORK1T(IFORCT(37)),WORK1T(IFORCT(38)),  ! gpcod,shape
     .             WORK1T(IFORCT(39)),WORK1T(IFORCT(40)),  ! deriv,posgp
     .             WORK1T(IFORCT(41)),WORK1T(IFORCT(42)),  ! weigp,xjacm
     .             WORK1T(IFORCT(16)),WORK1T(IFORCT(43)),  ! elco1,emas1
     .                                WORK1T(IFORCT(44)),  ! eldi1
     .             WORK1T(IFORCT(49)),WORK1T(IFORCT(50)),  ! dvoli,elcoi
     .             WORK1T(IFORCT(13)),                     ! dispt
     .             0,ITASKT)
C
      CALL SMOX05T(WORK1T(IFORCT(36)),                     ! DVOLU
     .             WORK1T(IFORCT(37)),                     ! gpcod
     .             WORK1T(IFORCT(43)),                     ! EMASS
     .             WORK1T(IFORCT(38)),                     ! SHAPE
     .             EHISTT,PROPST,
     .             WORK1T(IFORCT(54)),WORK1T(IFORCT(53)))  ! SFIPB,SFIPA
      IF(ISMO1T.EQ.1) NRULET=NRUAUX
      RETURN
C
      END
