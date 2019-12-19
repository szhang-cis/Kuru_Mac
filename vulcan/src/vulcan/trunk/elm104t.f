      SUBROUTINE ELM104T(LNODST,PROELT,PROPST,WORK1T,
     .                   ELCODT,CARTDT,DVOLUT,GPCODT,SHAPET,EPMTXT,
     .                   RMAT1T,
     .                   EMASST,VNORIT,                           !ELDAT
     .                   STRA0T,STRS0T,TEMPCT,                    !ELPRE
     .                   ELDIST,EHISTT,STRANT,STRSGT,             !ELVAR
     .                   CSTIFT,ESTIFT,WSTIFT,HSTIFT,PSTIFT,
     .                   QSTIFT,                                  !ELMAT
     .                                                      ITASKT)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS COMPUTATIONS AT ELEMENT LEVEL FOR
C
C     ELEMENT NO. 104 : 1D/2D ELEMENT IN 2D/3D SPACE (GAP ELEMENT)
C
C     FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  3   Evaluate Damping Matrix, Save it in Stiffness
C                    Matrix
C       ITASK =  4   Evaluate Internal Resisting Forces
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
      DIMENSION LNODST(*), PROELT(*), PROPST(*), ELCODT(*), CARTDT(*),
     .          DVOLUT(*), GPCODT(*), SHAPET(*), ELDIST(*), ESTIFT(*),
     .          WSTIFT(*), QSTIFT(*), WORK1T(*), ACBOUT(3), STRSGT(*),
     .          EHISTT(NHISTT,*)
      DIMENSION VNORIT(*)
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
      NOCOLT=INT(PROELT(11))
C
      NNODNT=NNODLT
      IF(NOCOLT.EQ.1) NNODNT=INT(PROELT(10))
C
C**** Direct to appropriate processor
C
      GO TO ( 10,  1, 30, 40,  1,  1,  1,  1, 90,  1,
     .         1,  1,  1,  1,  1,160,  1,  1,  1), ITASKT
    1 RETURN ! Nothing
C
C**** Perform some needed initial computations
C
   10 CONTINUE
      NEGATT=2
      CALL SEI104T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,VNORIT,
     .             WORK1T(ISETMT( 8)),                     ! dvolu
     .             WORK1T(ISETMT( 9)),WORK1T(ISETMT(10)),  ! gpcod,shape
     .             WORK1T(ISETMT( 1)),WORK1T(ISETMT( 2)),  ! deriv,posgp
     .             WORK1T(ISETMT( 3)),WORK1T(ISETMT( 4)),  ! weigp,xjacm
     .             WORK1T(ISETMT( 6)),WORK1T(ISETMT(12)),  ! elco1,eldi1
     .             WORK1T(ISETMT(15)),                     ! vnorl
     .             WORK1T(ISETMT(16)),WORK1T(ISETMT(17)),  ! dvoli,elcoi
     .             WORK1T(ISETMT(18)),                     ! dispt
     .             STRSGT,ITASKT)
      RETURN
C
C**** Evaluate Stiffness Matrix
C
   30 CONTINUE
C
C**** Perform some needed initial computations
C
      NEGATT=2
      CALL SEI104T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,VNORIT,
     .             WORK1T(ISTIFT(27)),                     ! dvolu
     .             WORK1T(ISTIFT(28)),WORK1T(ISTIFT(29)),  ! gpcod,shape
     .             WORK1T(ISTIFT(30)),WORK1T(ISTIFT(31)),  ! deriv,posgp
     .             WORK1T(ISTIFT(32)),WORK1T(ISTIFT(33)),  ! weigp,xjacm
     .             WORK1T(ISTIFT(18)),WORK1T(ISTIFT(39)),  ! elco1,eldi1
     .             WORK1T(ISTIFT(43)),                     ! vnorl
     .             WORK1T(ISTIFT(46)),WORK1T(ISTIFT(47)),  ! dvoli,elcoi
     .             WORK1T(ISTIFT(44)),                     ! dispt
     .             STRSGT,ITASKT)
C
      CALL ASIL00(WORK1T(ISTIFT(37)),WORK1T(ISTIFT(38)),   ! ESTIF,WSTIF
     .            NKOVAT,KDYNAT,
     .            WORK1T(ISTIFT(50)),                      ! FMATX
     .            NEVABT,     1)
C
      CALL MAS104T(WORK1T(ISTIFT(27)),                     ! DVOLU
     .             PROPST,
     .             WORK1T(ISTIFT(29)),                     ! SHAPE
     .             WORK1T(ISTIFT(37)),                     ! ESTIF
     .             EHISTT,
     .             WORK1T(ISTIFT(39)),                     ! ELDIS
     .             WORK1T(ISTIFT(13)),                     ! velcm
     .             WORK1T(ISTIFT(35)),                     ! preasl
     .             WORK1T(ISTIFT(36)),                     ! tgapsl
     .             WORK1T(ISTIFT(42)),WORK1T(ISTIFT(43)),  ! BOUCH,vnorl
     .             WORK1T(ISTIFT(44)),WORK1T(ISTIFT(45)),  ! DISPT,fpchl
     .             LNODST,
     .             WORK1T(ISTIFT(48)),WORK1T(ISTIFT(49)))  ! EMATX,FMATX
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
C**** Evaluate internal forces
C
   40 CONTINUE
C
C**** Perform some needed initial computations
C
      NEGATT=2
      CALL SEI104T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,VNORIT,
     .             WORK1T(IFORCT(36)),                     ! dvolu
     .             WORK1T(IFORCT(37)),WORK1T(IFORCT(38)),  ! gpcod,shape
     .             WORK1T(IFORCT(39)),WORK1T(IFORCT(40)),  ! deriv,posgp
     .             WORK1T(IFORCT(41)),WORK1T(IFORCT(42)),  ! weigp,xjacm
     .             WORK1T(IFORCT(16)),WORK1T(IFORCT(44)),  ! elco1,eldi1
     .             WORK1T(IFORCT(47)),                     ! vnorl
     .             WORK1T(IFORCT(49)),WORK1T(IFORCT(50)),  ! dvoli,elcoi
     .             WORK1T(IFORCT(13)),                     ! dispt
     .             STRSGT,ITASKT)
C
      CALL FRI104T(PROPST,
     .             WORK1T(IFORCT(44)),                     ! ELDIS
     .             WORK1T(IFORCT( 1)),
     .             WORK1T(IFORCT(36)),WORK1T(IFORCT(38)),  ! DVOLU,SHAPE
     .             EHISTT,
     .             WORK1T(IFORCT(13)),                     ! DISPT
     .             WORK1T(IFORCT(16)),                     ! ELCOD
     .             WORK1T(IFORCT(14)),WORK1T(IFORCT(15)),
     .             WORK1T(IFORCT(11)),WORK1T(IFORCT(12)),
     .             WORK1T(IFORCT(46)),WORK1T(IFORCT(47)),  ! BOUCH,vnorl
     .             WORK1T(IFORCT(48)),LNODST)              ! fpchl
      RETURN
C
C**** Check correctness of integration rule
C
   90 CONTINUE
      NNOBOT=NNODLT/2
      IF(NOCOLT.EQ.1) NNOBOT=NNODLT                ! non-coincident-mesh
      CALL CHEK01T(NDIMLT,IELEMT,NNOBOT,NRULET,NGAULT,   104,LUREST)
      RETURN
C
C**** OUTPOST
C
  160 CONTINUE
      NEGATT=2
      CALL SEI104T(DVOLUT,ELCODT,GPCODT,LNODST,PROPST,
     .             SHAPET,THICKT,ELDIST,VNORIT,
     .             WORK1T(ISTART( 4)),                     ! dvolu
     .             WORK1T(ISTART( 5)),WORK1T(ISTART( 6)),  ! gpcod,shape
     .             WORK1T(ISTART( 7)),WORK1T(ISTART( 8)),  ! deriv,posgp
     .             WORK1T(ISTART( 9)),WORK1T(ISTART(10)),  ! weigp,xjacm
     .             WORK1T(ISTART(11)),WORK1T(ISTART(13)),  ! elco1,eldi1
     .             WORK1T(ISTART(14)),                     ! vnorl
     .             WORK1T(ILOSUT(35)),WORK1T(ILOSUT(36)),  ! dvoli,elcoi
     .             WORK1T(ILOSUT(37)),                     ! dispt
     .             STRSGT,ITASKT)
      RETURN
C
      END
