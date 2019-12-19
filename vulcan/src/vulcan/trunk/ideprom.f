      SUBROUTINE IDEPROM(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR MICROSTRUCTURAL PROBLEMS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'   ! thermal-microstructural
      INCLUDE 'nuee_om.f'   ! mechanical-microstructural
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION PROPS(*)
C
      I4=NPROPM+2                                   ! NPROPM+1 => KPOROT
      NPLATX=INT(PROPS(I4))
      IF(NPLATX.EQ.0) RETURN
C
      IF(NPLATM.NE.NPLATX)
     . CALL RUNEND('ERROR: NPLAT (mechanical) NE NPLAT (micro)')
C
      INUPC=0
      INUPM=0
C
C**** CORRECTS THE ASSIGMENT OF VPLATM GIVEN IN idenpr.f (ideipt.f, etc)
C
      IX=I4
      DO IPLAT=1,NPLATM
       I4A=IX+1
       IPCFOM=INT(PROPS(I4A))
       IF(IPCFOM.NE.0) THEN
C
        INUPM=INUPM+1
        VPLATM(IPLAT,6)=NNUPCM+INUPM
C
        I4B=I4A+1
        IPCMOM=INT(PROPS(I4B))
C
        GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14) IPCMOM
C
    1   CALL IDEPROM1(PROPS,IPLAT,IX)
        GO TO 100
C
    2   CALL IDEPROM2(PROPS,IPLAT,IX)
        GO TO 100
C
    3   CALL IDEPROM3(PROPS,IPLAT,IX)
        GO TO 100
C
    4   call runend('error in ideprom: ideprom4 not implemented yet')
        GO TO 100
C
    5   call runend('error in ideprom: ideprom5 not implemented yet')
        GO TO 100
C
    6   CALL IDEPROM6(PROPS,IPLAT,IX,INUPM)
        GO TO 100
C
    7   CALL IDEPROM7(PROPS,IPLAT,IX)
        GO TO 100
C
    8   CALL IDEPROM8(PROPS,IPLAT,IX)
        GO TO 100
C
    9   call runend('error in ideprom: ideprom9 not implemented yet')
        GO TO 100
C
   10   call runend('error in ideprom: ideprom10 not implemented yet')
        GO TO 100
C
   11   call runend('error in ideprom: ideprom11 not implemented yet')
        GO TO 100
C
   12   call runend('error in ideprom: ideprom12 not implemented yet')
        GO TO 100
C
   13   call runend('error in ideprom: ideprom13 not implemented yet')
        GO TO 100
C
   14   CALL IDEPROM14(PROPS,IPLAT,IX)
        GO TO 100
C
  100   CONTINUE
C
       ELSE
        IX=IX+1        ! 1=ipcfo
C
        INUPC=INUPC+1
        VPLATM(IPLAT,6)=INUPC
C
       ENDIF           ! ipcfom.ne.0
      ENDDO            ! do iplat=1,nplat
C
      RETURN
      END
