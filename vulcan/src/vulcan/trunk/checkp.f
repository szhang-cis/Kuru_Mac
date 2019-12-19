      SUBROUTINE CHECKP(MATNO, PROEL, PROPS,
     .                  MATNOT,PROELT,PROPST)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS MECHANICAL & THERMAL PROPERTIES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION MATNO(NELEM),          PROEL(NPREL,NGRUP),
     .          PROPS(NPROP,NMATS)
      DIMENSION MATNOT(NELEMT),        PROELT(NPRELT,NGRUPT),
     .          PROPST(NPROPT,NMATST)
C
C**** OPTIONS TO CHECK (ICHEK):
C
C     1: EQUAL ELEMENT & MATERIAL NUMBERING
C     2: EQUAL MATERIAL NUMBERING
C
      ICHEK=2                     ! better as input
      GO TO (1,2) ICHEK
C
C**** OPTION 1
C
    1 CONTINUE
C
C**** LOOP OVER ELEMENTS
C
      DO 1000 IELEM=1,NELEMC      ! solid elements
      LGRUP=MATNO(IELEM)          ! set number
      LMATS=INT(PROEL(1,LGRUP))   ! material number
      LGRUPT=MATNOT(IELEM)
      LMATST=INT(PROELT(1,LGRUPT))
C
      IF(LGRUP.NE.LGRUPT)
     . CALL RUNEND('ERROR: LGRUP NE LGRUPT')
      IF(LMATS.NE.LMATST)
     . CALL RUNMEN('WARNING: LMATS NE LMATST')
C
C**** OBTAINS MATERIAL PROPERTIES
C
      CALL IDENPR(PROPS(1,LMATS))
      CALL IDENPRT(PROPST(1,LMATST))
C
C**** CHECKS NUMBER OF PHASE CHANGES
C
      IF(NPLAT.NE.NPLATM)
     . CALL RUNEND('ERROR: NPLAT NE NPLATM     ')
C
 1000 CONTINUE
      GO TO 100
C
C**** OPTION 2
C
    2 CONTINUE
C
C**** LOOP OVER MATERIALS
C
      NGRUX=NGRUP
      IF(NGRUPT.LT.NGRUP) NGRUX=NGRUPT
      DO 2000 IGRUP=1,NGRUX       ! sets
      LMATS=INT(PROEL(1,IGRUP))   ! material number
      LMATST=INT(PROELT(1,IGRUP))
C
      IF(LMATS.NE.LMATST)
     . CALL RUNMEN('WARNING: LMATS NE LMATST')   ! it should be RUNEND ?
C
      LTYPE=INT(PROEL(5,IGRUP))
      LTYPET=INT(PROELT(5,IGRUP))
C
C**** SOLIDS ELEMENTS
C
      IF(LTYPE.EQ.30.AND.LTYPET.EQ.5) THEN
C
C**** OBTAINS MATERIAL PROPERTIES
C
       CALL IDENPR(PROPS(1,LMATS))
       CALL IDENPRT(PROPST(1,LMATST))
C
C**** CHECKS NUMBER OF PHASE CHANGES
C
       IF(NPLAT.NE.NPLATM)
     .  CALL RUNEND('ERROR: NPLAT NE NPLATM     ')
C
C**** CHECKS FREE ENERGY MODELS
C
       IF(IFREN.NE.IFRENT) THEN
        IFRENT=IFREN
        CALL RUNMENT('WARNING: FREE ENERGY MODEL GIVEN BY MECH. INPUT')
       ENDIF
      ENDIF
C
C**** CONTACT ELEMENTS
C
      IF(LTYPE.EQ.4.OR.LTYPE.EQ.32.AND.LTYPET.EQ.104) THEN
C
C**** OBTAINS MATERIAL PROPERTIES
C
       CALL IDENPHT(PROPST(1,LMATST),    1)
C
C**** TRANSFERS NORMAL STIFFNESS FROM PROPS TO PROPST
C
       IF(NRADP.GT.0) THEN
        INRIG=1+
     .        1+2*NRADH1+1+2*NRADH2+1+2*NRADU+1+2*NRADP
        IF(IFILL.EQ.1)                  ! to be revised for filling
     .   INRIG=INRIG+
     .        1+2*NRADH1FI+1+2*NRADH2FI+1+2*NRADUFI+1+2*NRADPFI
        PROPST(INRIG+1,LMATST)=PROPS(2,LMATS)          ! RIGINT
        IF(LTYPE.EQ.4) PROPST(INRIG+2,LMATST)=1.0      ! IR432
        IF(LTYPE.EQ.32) PROPST(INRIG+2,LMATST)=2.0
       ENDIF
C
       NOCOL=INT(PROEL(22,IGRUP))
       NOCOLT=INT(PROELT(11,IGRUP))
       IF(NOCOL.EQ.1) THEN              ! NOCOLT is also 1; see couinc.f
C
C**** TRANSFERS FRICTIONAL HEATING INDEX & PARTITION COEFFICIENT FROM
C     PROPST TO PROPS
C     (it is not necessary in this case to obtain mechanical material
C      properties)
C
        PROPS(101,LMATS)=PROPST(INRIG+3,LMATST)        ! IGAFR
        PROPS(102,LMATS)=0.0D0
        IF(INT(PROPS(101,LMATS)).EQ.1)
     .   PROPS(102,LMATS)=PROPST(INRIG+4,LMATST)       ! AGAFP
       ENDIF
      ENDIF
C
      NTYPE=INT(PROEL(6,IGRUP))
      NTYPET=INT(PROELT(6,IGRUP))
C
      IF(NTYPE.NE.NTYPET) THEN
       INERRO=0
       IF((NTYPE.EQ.1.AND.NTYPET.EQ.2).OR.
     .    (NTYPE.EQ.2.AND.NTYPET.EQ.1)) INERRO=1
C
       IF(INERRO.EQ.0) CALL RUNEND('CHECKP: WRONG NTYPE OR NTYPET')
      ENDIF
C
 2000 CONTINUE
      GO TO 100
C
  100 CONTINUE
C
C**** OTHER CHECK
C
      IF(NFURES.NE.NFURESM) THEN
       CALL RUNMEN('WARNING: NFURES =/ NFURESM')
       CALL RUNMENT('WARNING: NFURES =/ NFURESM')
      ENDIF
      INIAUX=0
      INIAUXT=0
      IF(INITI.EQ.2) INIAUX=1
      IF(INITIT.EQ.2) INIAUXT=1
      IF(INIAUX.NE.INIAUXT)
     . CALL RUNEND('ERROR: MECH. OR THER. PREVIOUS RESULTS ARE LACKING')
C
C**** CHECKS NUMBER OF NODES PER ELEMENT (THERMAL & MECHANICAL DATA)
C
C
C     Notes:
C
C     NPOINC must be the same in the thermal & mechanical data files
C
C     NNODET corresponds to linear or n_dime-linear elements and NNODE
C     is transformed to n_dime-quadratic elements
C
C     Toblerones should be implemented
C
C
C**** LOOP OVER ELEMENTS
C
      DO 1001 IELEM=1,NELEMC      ! solid elements
      LGRUP=MATNO(IELEM)          ! set number
      LMATS=INT(PROEL(1,LGRUP))   ! material number
      NNODL=INT(PROEL(2,LGRUP))   ! element node number
      NRULE=INT(PROEL(3,LGRUP))   ! integration rule
      NGAUL=INT(PROEL(4,LGRUP))   ! Gauss point number
      LGRUPT=MATNOT(IELEM)
      LMATST=INT(PROELT(1,LGRUPT))
      NNODLT=INT(PROELT(2,LGRUPT))
      NRULET=INT(PROELT(3,LGRUPT))
      NGAULT=INT(PROELT(4,LGRUPT))
C
      IF(NNODL.NE.NNODLT) THEN
       INTERC=1
       INERRO=0
C
C**** 1D
C
       IF(NDIME.EQ.1) THEN
        IF(NNODL.EQ.3.AND.NNODLT.EQ.2) INERRO=1
       ENDIF
C
C**** 2D
C
       IF(NDIME.EQ.2) THEN
        IF(NNODL.EQ.6.AND.NNODLT.EQ.3) INERRO=1   ! triangles
        IF(NNODL.EQ.7.AND.NNODLT.EQ.3) INERRO=1
        IF(NNODL.EQ.4.AND.NNODLT.EQ.3) THEN
         IF(NRULE.EQ.3.OR.NRULE.EQ.6.OR.NRULE.EQ.7) INERRO=1
        ENDIF
        IF(NNODL.EQ.8.AND.NNODLT.EQ.4) INERRO=1   ! quadrilaterals
        IF(NNODL.EQ.9.AND.NNODLT.EQ.4) INERRO=1
       ENDIF
C
C**** 3D
C
       IF(NDIME.EQ.3) THEN
        IF(NNODL.EQ.10.AND.NNODLT.EQ.4) INERRO=1  ! tetrahedra
        IF(NNODL.EQ.15.AND.NNODLT.EQ.6) INERRO=0  ! toblerones
        IF(NNODL.EQ.20.AND.NNODLT.EQ.8) INERRO=1  ! bricks
        IF(NNODL.EQ.27.AND.NNODLT.EQ.8) INERRO=1
       ENDIF
C
       IF(INERRO.EQ.0) CALL RUNEND('CHECKP: WRONG NNODE OR NNODET')
C
      ENDIF
C
      IF(NGAUL.NE.NGAULT) THEN
       IF(NITERC.EQ.1.OR.NITERC.EQ.4)
     .  CALL RUNEND('CHECKP: WRONG NUMBER OF GAUSS POINTS FOR NITERC=1')
      ENDIF
C
 1001 CONTINUE
C
      RETURN
      END
