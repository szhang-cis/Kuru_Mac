      SUBROUTINE ELM101T(LNODST,PROELT,PROPST,WORK1T,
     .                   ELCODT,CARTDT,DVOLUT,GPCODT,SHAPET,EPMTXT,
     .                   RMAT1T,
     .                   EMASST,                                 !ELDAT
     .                   STRA0T,STRS0T,TEMPCT,                   !ELPRE
     .                   ELDIST,EHISTT,STRANT,STRSGT,            !ELVAR
     .                   CSTIFT,ESTIFT,WSTIFT,HSTIFT,PSTIFT,
     .                   QSTIFT,                                 !ELMAT
     .                   HTLODT,
     .                                                      ITASKT)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS COMPUTATIONS AT ELEMENT LEVEL FOR
C
C     ELEMENT NO. 101 : 1D/2D ELEMENT IN 2D/3D SPACE
C
C     FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  3   Evaluate Capacity Matrix, Save it in Jacobian
C                    Matrix
C       ITASK =  4   Evaluate Internal Resisting Heats
C       ITASK =  7   Evaluate equivalent Volume Heats
C       ITASK =  9   Check correctness of integration rule
C       ITASK = 16   Evaluate variables for outpost.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/NEGATITE/NEGATT
C
      DIMENSION LNODST(*), PROELT(*), PROPST(*), ELCODT(*),
     .          CARTDT(*), DVOLUT(*), GPCODT(*), SHAPET(*), 
     .          ELDIST(*), ESTIFT(*), WSTIFT(*), QSTIFT(*), HTLODT(*),
     .          WORK1T(*), STRSGT(*)
C
      DIMENSION EHISTT(NHISTT,*)
C
C**** BEGIN
C
      NNODLT=INT(PROELT( 2)) ! Select element properties
      NRULET=INT(PROELT( 3))
      NGAULT=INT(PROELT( 4))
      NTYPET=INT(PROELT( 6))
      THICKT=    PROELT( 7)
C
      NDIMLT=2
      IF(NTYPET.EQ.5.OR.NTYPET.EQ.3) NDIMLT=1
C
C**** Direct to appropriate processor
C
      GO TO ( 10,  1, 30, 40,  1,  1, 70,  1, 90,  1,
     .         1,  1,  1,  1,  1,160,  1,  1,  1), ITASKT
    1 RETURN ! Nothing
C
C**** Perform some needed initial computations
C
   10 CONTINUE
      NEGATT=2
      CALL SEI101T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(ISETMT( 8)),                     ! dvolu
     .             WORK1T(ISETMT( 9)),WORK1T(ISETMT(10)),  ! gpcod,shape
     .             WORK1T(ISETMT( 1)),WORK1T(ISETMT( 2)),  ! deriv,posgp
     .             WORK1T(ISETMT( 3)),WORK1T(ISETMT( 4)),  ! weigp,xjacm
     .             WORK1T(ISETMT( 6)),WORK1T(ISETMT(12)),  ! elco1,eldi1
     .             WORK1T(ISETMT(16)),WORK1T(ISETMT(17)),  ! dvoli,elcoi
     .             WORK1T(ISETMT(18)),                     ! dispt
     .             ITASKT)
      RETURN
C
C**** Evaluate Jacobian Matrix
C
   30 CONTINUE
C
C**** Perform some needed initial computations
C
      NEGATT=2
      CALL SEI101T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(ISTIFT(27)),                     ! dvolu
     .             WORK1T(ISTIFT(28)),WORK1T(ISTIFT(29)),  ! gpcod,shape
     .             WORK1T(ISTIFT(30)),WORK1T(ISTIFT(31)),  ! deriv,posgp
     .             WORK1T(ISTIFT(32)),WORK1T(ISTIFT(33)),  ! weigp,xjacm
     .             WORK1T(ISTIFT(18)),WORK1T(ISTIFT(39)),  ! elco1,eldi1
     .             WORK1T(ISTIFT(46)),WORK1T(ISTIFT(47)),  ! dvoli,elcoi
     .             WORK1T(ISTIFT(44)),                     ! dispt
     .             ITASKT)
C
      CALL ASIL00(WORK1T(ISTIFT(37)),WORK1T(ISTIFT(38)),   ! ESTIF,WSTIF
     .            NKOVAT,KDYNAT,
     .            WORK1T(ISTIFT(50)),                      ! FMATX
     .            NEVABT,     1)
C
      CALL MAS101T(WORK1T(ISTIFT(27)),                     ! DVOLU
     .             PROPST,
     .             WORK1T(ISTIFT(29)),                     ! SHAPE
     .             WORK1T(ISTIFT(37)),                     ! ESTIF
     .             EHISTT,
     .             WORK1T(ISTIFT(39)),                     ! ELDIS
     .             WORK1T(ISTIFT(13)),                     ! velcm
     .             WORK1T(ISTIFT(42)),WORK1T(ISTIFT(45)),  ! BOUCH,fpchl
     .             WORK1T(ISTIFT(48)),WORK1T(ISTIFT(28)))  ! EMATX,gpcod
C
C**** ASSEMBLY PROCESS (NMEMO6=1)
C
      CALL ASAL00(ESTIFT,WSTIFT,
     .            WORK1T(ISTIFT(37)),                      ! ESTIF
     .            WORK1T(ISTIFT(38)),                      ! WSTIF
     .            NKOVAT,KDYNAT,NMEMO6)
C
      CALL ASEL05T(CSTIFT,
     .             WORK1T(ISTIFT(37)),                     ! ESTIF
     .             WORK1T(ISTIFT(38)),                     ! WSTIF
     .             WORK1T(ISTIFT(41)))                     ! CSTI1
C
      RETURN
C
C**** Evaluate internal heats
C
   40 CONTINUE
C
C**** Perform some needed initial computations
C
      NEGATT=2
      CALL SEI101T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(IFORCT(36)),                     ! dvol1
     .             WORK1T(IFORCT(37)),WORK1T(IFORCT(38)),  ! gpco1,shap1
     .             WORK1T(IFORCT(39)),WORK1T(IFORCT(40)),  ! deriv,posgp
     .             WORK1T(IFORCT(41)),WORK1T(IFORCT(42)),  ! weigp,xjacm
     .             WORK1T(IFORCT(16)),WORK1T(IFORCT(44)),  ! elco1,eldi1
     .             WORK1T(IFORCT(49)),WORK1T(IFORCT(50)),  ! dvoli,elcoi
     .             WORK1T(IFORCT(13)),                     ! dispt
     .             ITASKT)
C
      CALL FRI101T(PROPST,
     .             WORK1T(IFORCT(44)),                     ! ELDIS
     .             WORK1T(IFORCT( 1)),
     .             WORK1T(IFORCT(36)),WORK1T(IFORCT(38)),  ! DVOLU,SHAPE
     .             EHISTT,
     .             WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .             WORK1T(IFORCT(46)),WORK1T(IFORCT(48)),  ! BOUCH,fpchl
     .             WORK1T(IFORCT(37)))                     ! gpco1
C
C**** Evaluate equivalent Volume Heats (not constant; see ldv101t.f)
C
      IF(ILDV2T.EQ.1)
     . CALL LDV101T(WORK1T(IFORCT(44)),                  ! ELDIS
     .              WORK1T(IFORCT(36)),                  ! DVOLU
     .              PROPST,
     .              WORK1T(IFORCT(38)),                  ! SHAPE
     .              WORK1T(IFORCT(42)),                  ! XJACM
     .              WORK1T(IFORCT( 1)),                  ! FORCE
     .              EHISTT,HTLODT,
     .              WORK1T(IFORCT(46)),WORK1T(IFORCT(48)),  !BOUCH,fpchl
     .              WORK1T(IFORCT(16)),WORK1T(IFORCT(37)),  !ELCO1,GPCO1
     .              WORK1T(IFORCT(50)),GPCODT)              !ELCOI
      RETURN
C
C**** Evaluate equivalent Volume Heats (constant but time dependent
C                                       heat; see ldv101t.f)
C
   70 CONTINUE
      IF(ILDV1T.EQ.0) RETURN
C
C**** Perform some needed initial computations
C
      NEGATT=2
      CALL SEI101T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(ILOSUT(22)),                     ! dvolu
     .             WORK1T(ILOSUT(23)),WORK1T(ILOSUT(24)),  ! gpcod,shape
     .             WORK1T(ILOSUT(25)),WORK1T(ILOSUT(26)),  ! deriv,posgp
     .             WORK1T(ILOSUT(27)),WORK1T(ILOSUT(28)),  ! weigp,xjacm
     .             WORK1T(ILOSUT(17)),WORK1T(ILOSUT(30)),  ! elco1,eldi1
     .             WORK1T(ILOSUT(35)),WORK1T(ILOSUT(36)),  ! dvoli,elcoi
     .             WORK1T(ILOSUT(37)),                     ! dispt
     .             ITASKT)
C
      CALL LDV101T(WORK1T(ILOSUT(31)),       ! TEINI=ELDIS for this case
     .             WORK1T(ILOSUT(22)),                     ! DVOLU
     .             PROPST,
     .             WORK1T(ILOSUT(24)),                     ! SHAPE
     .             WORK1T(ILOSUT(28)),                     ! XJACM
     .             WORK1T(ILOSUT( 1)),                     ! FORCE
     .             EHISTT,DUMMYT,
     .             WORK1T(ILOSUT(33)),WORK1T(ILOSUT(34)),  ! BOUCH,fpchl
     .             WORK1T(ILOSUT(17)),WORK1T(ILOSUT(23)),  ! ELCO1,GPCO1
     .             WORK1T(ILOSUT(36)),GPCODT)              ! ELCOI
      RETURN
C
C**** Check correctness of integration rule
C
   90 CONTINUE
      CALL CHEK01T(NDIMLT,IELEMT,NNODLT,NRULET,NGAULT,   101,LUREST)
      RETURN
C
C**** OUTPOST
C
  160 CONTINUE
      NEGATT=2
      CALL SEI101T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,
     .             WORK1T(ISTART( 4)),                     ! dvolu
     .             WORK1T(ISTART( 5)),WORK1T(ISTART( 6)),  ! gpcod,shape
     .             WORK1T(ISTART( 7)),WORK1T(ISTART( 8)),  ! deriv,posgp
     .             WORK1T(ISTART( 9)),WORK1T(ISTART(10)),  ! weigp,xjacm
     .             WORK1T(ISTART(11)),WORK1T(ISTART(13)),  ! elco1,eldi1
     .             WORK1T(ISTART(15)),WORK1T(ISTART(16)),  ! dvoli,elcoi
     .             WORK1T(ISTART(17)),                     ! dispt
     .             ITASKT)
      RETURN
C
      END
