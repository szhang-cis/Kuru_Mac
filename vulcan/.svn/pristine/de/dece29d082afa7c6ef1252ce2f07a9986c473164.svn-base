      SUBROUTINE IDEPROS(PROPST)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR MICROSTRUCTURAL PROBLEMS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'                      ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION PROPST(*)
C
      I4=NPROPM+2                ! NPROPM+1 => KPOROT
      NPLAX=INT(PROPST(I4))
      IF(NPLAX.EQ.0) RETURN
C
      IF(NPLAT.NE.NPLAX)
     . CALL RUNENDT('ERROR: NPLAT (thermal) NE NPLAT (micro)')
C
      IX=I4
      INUPM=0                      ! number of microstructural ph-ch
      DO IPLAT=1,NPLAT
       I4A=IX+1
       IPCFO=INT(PROPST(I4A))
       IF(IPCFO.NE.0) THEN
        VPLAT(IPLAT,4)=PROPST(I4A)
        I4B=I4A+1
        IPCMO=INT(PROPST(I4B))
        VPLAT(IPLAT,5)=PROPST(I4B)
C
        INUPM=INUPM+1
        GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14), IPCMO
C
    1   CALL IDEPROS1(PROPST,IPLAT,IX)
        GO TO 100
C
    2   CALL IDEPROS2(PROPST,IPLAT,IX)
        GO TO 100
C
    3   CALL IDEPROS3(PROPST,IPLAT,IX)
        GO TO 100
C
    4   CALL IDEPROS4(PROPST,IPLAT,IX)
        GO TO 100
C
    5   CALL IDEPROS5(PROPST,IPLAT,IX)
        GO TO 100
C
    6   CALL IDEPROS6(PROPST,IPLAT,IX,INUPM)
        GO TO 100
C
    7   CALL IDEPROS7(PROPST,IPLAT,IX)
        GO TO 100
C
    8   CALL IDEPROS8(PROPST,IPLAT,IX)
        GO TO 100
C
    9   CALL IDEPROS9(PROPST,IPLAT,IX)
        GO TO 100
C
   10   CALL IDEPROS10(PROPST,IPLAT,IX)
        GO TO 100
C
   11   CALL IDEPROS11(PROPST,IPLAT,IX,INUPM)
        GO TO 100
C
   12   CALL IDEPROS12(PROPST,IPLAT,IX)
        GO TO 100
C
   13   CALL RUNENDT('ERROR: MODEL 13 NOT IMPLEMENTED YET')
        GO TO 100
C
   14   CALL IDEPROS14(PROPST,IPLAT,IX)
        GO TO 100
C
  100   CONTINUE
C
       ELSE
        IX=IX+1        ! 1=ipcfo
       ENDIF           ! ipcfo.ne.0
      ENDDO            ! do iplat=1,nplat
C
C**** PHASE-CHANGE ASSOCIATED TO THE POROSITY CRITERIA
C
      IF(NPOROT.GT.0) THEN
       IF(KPOROT.GT.0) THEN
        KPOROT=INT(PROPST(NPROPM+1))
       ENDIF
      ENDIF
C
      RETURN
      END
