      SUBROUTINE DATPANT(ARRAYT,IFLAGT,ITASKT)
C***********************************************************************
C
C**** THIS ROUTINE INTERFACES VULCAN-DATA WITH PREVIOUS ANALYSIS-FILE
C     (*.pan)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_omt.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION ARRAYT(*)
C
C**** BEGIN
C
      NUNITT=LUPANT                               ! .pan
C
      GO TO (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     .        11                                      ),IFLAGT
C
    1 NRECDT=NRECGT+(IELEMT-1)*IDATPT(1,2)                      ! ELDAT
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,NDATAT*2,LENRCT,NRECDT)
      RETURN
C
    2 NRECDT=NRECCT+IDATCT(1)+(IELEMT-1)*(IDATPT(2,2)+
     .       IDATPT(3,2))                                       ! ELPRE
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,NPREVT*2,LENRCT,NRECDT)
      RETURN
C
    3 NRECDT=NRECCT+IDATCT(1)+(IELEMT-1)*(IDATPT(2,2)+
     .       IDATPT(3,2))+IDATPT(2,2)                           ! ELVAR
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,NSTATT*2,LENRCT,NRECDT)
      RETURN
C
    4 NRECDT=NRECCT+IDATCT(2)                                   ! DISTO
      IIII3=1+KDYNAT+NMEMO8
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,IIII3*NTOTVT*2,LENRCT,NRECDT)
      RETURN
C
    5 NRECDT=NRECCT+IDATCT(3)                                   ! DISPR
      IIII3=1+KDYNAT+NMEMO8  ! 1+KDYNAT+NMEMO8>1+NMEMO10 (see rssetct.f)
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,IIII3*NTOTVT*2,LENRCT,NRECDT)
      RETURN
C
    6 NRECDT=NRECCT+IDATCT(4)                                   ! HEADS
      KPORETC=0
      IF(KPORET.NE.0) KPORETC=1
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,KPORETC*NPOINT*4*2,LENRCT,
     .             NRECDT)
      RETURN
C
    7 NRECDT=NRECCT+IDATCT(5)                                   ! REFOR
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,NTOTVT*2,LENRCT,NRECDT)
      RETURN
C
    8 NRECDT=NRECCT+IDATCT(6)                                   ! TLOAD
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,NTOTVT*2*2,LENRCT,NRECDT)
      RETURN
C
    9 NRECDT=NRECCT+IDATCT(7)                                   ! HTLOD
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,NHLODT*NSUBFT*NFUNCT*2,
     .             LENRCT,NRECDT)
      RETURN
C
   10 NRECDT=NRECCT+IDATCT(8)                                   ! IFFIX
      KPORETA=1
      IF(KPORET-1.GT.0) KPORETA=2
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,NTOTVT*KPORETA*(1+NACTIT),
     .             LENRCT,NRECDT)
      RETURN
C
   11 NRECDT=NRECCT+IDATCT(9)                                   ! RLOAD
      CALL RESTDKT(NUNITT,ITASKT,ARRAYT,NTOTVT*2,LENRCT,NRECDT)
      RETURN
C
      END
