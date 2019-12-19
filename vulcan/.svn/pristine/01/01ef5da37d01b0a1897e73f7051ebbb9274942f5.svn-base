      SUBROUTINE DATRST(ARRAY,IFLAG,ITASK)
C***********************************************************************
C
C**** THIS ROUTINE INTERFACES VULCAN-DATA WITH RESTART-FILE
C     (*.rst or *.fan)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_om.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION ARRAY(*)
C
C**** BEGIN
C
      NUNIT=LURST                                      ! .rst
      IF(NFURESM.EQ.1) NUNIT=LUFAN                     ! .fan
C
      GO TO (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     .        11, 12                                  ),IFLAG
C
    1 NRECD=NRECG+(IELEM-1)*IDATP(1,2)                          ! ELDAT
      CALL RESTDK(NUNIT,ITASK,ARRAY,NDATA*2,LENRC,NRECD)
      RETURN
C
    2 NRECD=NRECC+IDATC(1)+(IELEM-1)*(IDATP(2,2)+IDATP(3,2))    ! ELPRE
      CALL RESTDK(NUNIT,ITASK,ARRAY,NPREV*2,LENRC,NRECD)
      RETURN
C
    3 NRECD=NRECC+IDATC(1)+(IELEM-1)*(IDATP(2,2)+IDATP(3,2))
     .                    +IDATP(2,2)                           ! ELVAR
      CALL RESTDK(NUNIT,ITASK,ARRAY,NSTAT*2,LENRC,NRECD)
      RETURN
C
    4 NRECD=NRECC+IDATC(2)                                      ! DISTO
      CALL RESTDK(NUNIT,ITASK,ARRAY,NDISO*NTOTV*2,LENRC,NRECD)
      RETURN
C
    5 NRECD=NRECC+IDATC(3)                                      ! DISPR
      CALL RESTDK(NUNIT,ITASK,ARRAY,NDISR*NTOTV*2,LENRC,NRECD)
      RETURN
C
    6 NRECD=NRECC+IDATC(4)                                      ! HEADS
      KPOREC=0
      IF(KPORE.NE.0) KPOREC=1
      CALL RESTDK(NUNIT,ITASK,ARRAY,KPOREC*NPOIN*4*2,LENRC,NRECD)
      RETURN
C
    7 NRECD=NRECC+IDATC(5)                                      ! REFOR
      CALL RESTDK(NUNIT,ITASK,ARRAY,NTOTV*2,LENRC,NRECD)
      RETURN
C
    8 NRECD=NRECC+IDATC(6)                                      ! TLOAD
      CALL RESTDK(NUNIT,ITASK,ARRAY,NTOTV*2*2,LENRC,NRECD)
      RETURN
C
    9 NRECD=NRECC+IDATC(7)                                      ! HTLOD
      CALL RESTDK(NUNIT,ITASK,ARRAY,NHLOD*NSUBF*NFUNC*2,LENRC,NRECD)
      RETURN
C
   10 NRECD=NRECC+IDATC(8)                                      ! IFFIX
      KPOREA=1
      IF(KPORE-1.GT.0) KPOREA=2
      CALL RESTDK(NUNIT,ITASK,ARRAY,NTOTV*KPOREA*(1+NACTI),LENRC,NRECD)
      RETURN
C
   11 NRECD=NRECC+IDATC(9)                                      ! RLOAD
      CALL RESTDK(NUNIT,ITASK,ARRAY,NTOTV*2,LENRC,NRECD)
      RETURN
C
   12 NRECD=NRECC+IDATC(10)                                     ! LACTI
      CALL RESTDK(NUNIT,ITASK,ARRAY,NELEM,LENRC,NRECD)
      RETURN
C
      END
