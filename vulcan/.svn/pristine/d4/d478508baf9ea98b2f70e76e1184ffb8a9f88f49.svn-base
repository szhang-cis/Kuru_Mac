      SUBROUTINE INMICRO(PROPST,ALPHAM)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES SOME MICROSTRUCTURAL PARAMETERS
C
C***********************************************************************
C
C             IPCFO=0 TOTAL PHASE-CHANGE FORMULATION
C                     f_pc(T) is bijective
C
C                     IPCMO=1 LINEAR
C                     IPCMO=2 IDEM 1 WITH f_s
C                     IPCMO=3 PARABOLIC
C                     IPCMO=4 CUBIC
C                     IPCMO=5 SCHEIL'S EQUATION
C                     IPCMO=6 LEVER'S EQUATION
C                     IPCMO=7 ....
C
C             IPCFO=1 RATE PHASE-CHANGE FORMULATION
C                     f_pc(\alpha_m_k,T) where:
C                     \alpha_m_k=microscopical variables
C                     => the computation of f^._pc is necessary
C
C                     IPCMO=1 MODEL ...
C                     IPCMO=2 MODEL ...
C
C***********************************************************************
C
C     Index of variables
C
C     BASMM = Density
C     BASCC = Capacity coefficient
C     BASKK = Conductivity
C     SOUR2T= L*Phase-change function (total)
C     DSOURT= L*Phase-change function rate (total)
C
C     ALPHAM= array of microstructural (microscopical) variables
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION PROPST(*), ALPHAM(NHISTM)
C
      INUPC=0
      INUPM=0                      ! number of microstructural ph-ch
      DO IPLAT=1,NPLAT
       IPCFO=INT(VPLAT(IPLAT,4))
       IPCMO=INT(VPLAT(IPLAT,5))
C
       IF(IPCFO.EQ.0) GO TO 101
C
C**** RATE PHASE-CHANGE FORMULATION (IPCFO=1) > MICROSCOPICAL MODELS
C
       INUPM=INUPM+1
       GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14), IPCMO
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 1
C
    1  CALL INMICRO1(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 2
C
    2  CALL INMICRO2(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 3
C
    3  CALL INMICRO3(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 4
C
    4  CALL INMICRO4(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 5
C
    5  CALL INMICRO5(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 6
C
    6  CALL INMICRO6(ALPHAM,INUPC,IPLAT,INUPM)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 7
C
    7  CALL INMICRO7(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 8
C
    8  CALL INMICRO8(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 9
C
    9  CALL INMICRO9(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 10
C
   10  CALL INMICRO10(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 11
C
   11  CALL INMICRO11(ALPHAM,INUPC,IPLAT,INUPM)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 12
C
   12  CALL INMICRO12(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 13
C
   13  CALL RUNENDT('ERROR: MODEL 13 NOT IMPLEMENTED YET')
       GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 14
C
   14  CALL INMICRO14(ALPHAM,INUPC,IPLAT)
       GO TO 100
C
  100  CONTINUE
C
  101  CONTINUE
C
      ENDDO         ! iplat=1,nplat
C
      RETURN
      END
