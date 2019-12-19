      SUBROUTINE INPPROS(ITAPES,PROPSS,INDEXS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MICROSTRUCTURAL PROPERTIES
C
C     Notes:
C           NMATSS: number of materials
C           NGRUPS: number of sets
C           In general, NGRUPS ge NMATSS
C
C           INDEXS=0 => material properties read at GENERAL_DATA level
C           INDEXS=1 => material properties read at INTERVAL_DATA level
C
C           INDEXS=1 not implemented yet !!! (January 1997)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES (thermal & microstructural problems)
C
      INCLUDE 'nued_om.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      DIMENSION PROPSS(NPROPS,*)
C
      PARAMETER(MCOMDS=5)
      PARAMETER(MCOMES=2)
      CHARACTER*5 COMMDS(MCOMDS)
      CHARACTER*5 COMMES(MCOMES)
C
      DATA COMMDS/'STISO','STORT','STANI','BOUND','AIRGA'/
      DATA COMMES/'     ','FILLI'/
C
      IF(IEVFI.EQ.0) THEN
C
C**** INITIALIZATION
C
C     NSKIPS=number of properties from thermal to micro
C
       NSKIPS=1                                     ! KPOROT
       IF(INDEXS.EQ.0) THEN
        DO IPROPS=NSKIPS+1,NPROPS
         DO IMATSS=1,NMATSS
          PROPSS(IPROPS,IMATSS)=0.0D+00
         ENDDO
        ENDDO
       ENDIF                                        ! indexs.eq.0
C
C**** READ MICROSTRUCTURAL PROPERTIES
C
       ITAPES=LUDATS
       DO 11 IMATST=1,NMATSS
        NPRINT=0
        CALL LISTENS('INPPROS',NPRINT,ITAPES)
C
C**** LOOK FOR 'MATERIAL_DATA' CARD
C
        IF(WORDSS(1).NE.'MATER') GO TO 1110
C
        IF(PARAMS(1).EQ.0.0)
     .   CALL RUNENDS('INPPROS: WRONG MATERIAL NUMBER')
C
        NUMATS=INT(PARAMS(1))
C
        CALL PROPMIC(ITAPES,PROPSS(1,NUMATS))
C
   11  CONTINUE
C
C**** LOOK FOR 'END_MICROSTRUCTURAL_PROPERTIES' CARD
C
       NPRINT=0
       JTAPES=LUDATS
       CALL LISTENS('INPPROS',NPRINT,JTAPES)
 1110  IF(WORDSS(1).NE.'END_M')
     .  CALL RUNENDS('INPPROS: END_MICROSTRUCTURAL CARD NOT FOUND')
C
      ELSE          ! ievfi=1
C
c      IERR1S=1
c      IF(ITAPES.NE.LUDATS.AND.ITAPES.NE.LUMATS) GOTO 1000
c      IF(INDEXS.EQ.0) THEN
c       IF(ITAPES.EQ.LUMATS)                        ! external .mat file
c    .   OPEN(UNIT=LUMATS,FILE=CC1S,STATUS='OLD',ERR=1000)
c      ENDIF                                        ! indext.eq.0
c      IERR1S=0
C
c1000  IF(IERR1S.NE.0) THEN
c       IF(IERR1S.EQ.1) THEN
c        WRITE(LURESS,901)
c        CALL RUNENDS('ERROR IN OPENING FILES')
c       ENDIF
c      ENDIF
C
       WRITE(LURESS,900)
C
C**** INITIALIZATION
C
       IF(INDEXS.EQ.0) THEN
        DO IPROPS=1,NPROPS
         DO IMATSS=1,NMATSS
          PROPSS(IPROPS,IMATSS)=0.0D+00
         ENDDO
        ENDDO
       ENDIF                                        ! indexs.eq.0
C
C**** READ MATERIALS PROPERTIES
C
       DO 10 IMATSS=1,NMATSS
C
       NPRINS=0
       CALL LISTENS('INPPROS',NPRINS,ITAPES)
C
C**** LOOK FOR 'MATERIAL_DATA' CARD
C
       IF(INDEXS.EQ.0) THEN
        IF(WORDSS(1).NE.'MATER')
     .   CALL RUNENDS('INPPROS: WRONG MATERIAL_DATA CARD')
       ELSE
        IF(WORDSS(1).NE.'MATER') GO TO 2000
       ENDIF
C
       IF(PARAMS(1).EQ.0.0)
     .  CALL RUNENDS('INPPROS: WRONG MATERIAL NUMBER')
C
       NUMATS=INT(PARAMS(1))
C
C**** CONTROLS REPEATED MATERIALS
C
       IF(INDEXS.EQ.0) THEN
        DO IPROPS=1,NPROPS
         IF(PROPSS(IPROPS,NUMATS).NE.0.0)
     .    CALL RUNENDS('ERROR IN MATERIAL NUMBER')
        ENDDO
       ENDIF
C
       DO ICOMDS=1,MCOMDS
        IF(WORDSS(2).EQ.COMMDS(ICOMDS)) GOTO 100
       ENDDO
       CALL RUNENDS('INPPROS: WRONG MATERIAL CARD')
C
  100  GOTO(1,2,3,4,5), ICOMDS
C
    1  CONTINUE
       DO ICOMES=1,MCOMES
        IF(WORDSS(3).EQ.COMMES(ICOMES)) GOTO 110
       ENDDO
       CALL RUNENDS('INPPROS: WRONG MATERIAL CARD')
C
  110  GOTO(111,112), ICOMES
C
  111  CALL PROP05S(ITAPES,PROPSS(1,NUMATS))
       GOTO 10
C
  112  call runends('prop05fs not implemented')
c 112  CALL PROP05FS(ITAPES,PROPSS(1,NUMATS))
       GOTO 10
C
    2  CONTINUE
       DO ICOMES=1,MCOMES
        IF(WORDSS(3).EQ.COMMES(ICOMES)) GO TO 210
       ENDDO
       CALL RUNENDS('INPPROS: WRONG MATERIAL CARD')
C
  210  GOTO(211,212), ICOMES
C
  211  call runends('prop05so not implemented')
c 211  CALL PROP05SO(ITAPES,PROPSS(1,NUMATS))
       GOTO 10
C
  212  call runends('pro05fos not implemented')
c 212  CALL PRO05FOS(ITAPES,PROPSS(1,NUMATS))
       GOTO 10
C
    3  CONTINUE
       CALL RUNENDS('ERROR: FULLY ANISOTROPIC MATERIAL NOT IMPLEMENTED')
c      note: change VCONDS(500,2,9) in auxl_oms.f
       GOTO 10
C
    4  CONTINUE
       DO ICOMES=1,MCOMES
        IF(WORDSS(3).EQ.COMMES(ICOMES)) GOTO 410
       ENDDO
       CALL RUNENDS('INPPROS: WRONG MATERIAL CARD')
C
  410  GOTO(411,412), ICOMES
C
  411  call runends('pro101s not implemented')
c 411  CALL PRO101S(ITAPES,PROPSS(1,NUMATS))
       GOTO 10
C
  412  call runends('pro101fs not implemented')
c 412  CALL PRO101FS(ITAPES,PROPSS(1,NUMATS))
       GOTO 10
C
    5  call runends('pro104s not implemented')
c   5  CALL PRO104S(ITAPES,PROPSS(1,NUMATS))
       GOTO 10
C
   10  CONTINUE
C
C**** LOOK FOR 'END_PROPERTIES' CARD
C
       NPRINS=0
       JTAPES=LUDATS
       CALL LISTENS('INPPROS',NPRINS,JTAPES)
 2000  IF(WORDSS(1).NE.'END_P')
     .  CALL RUNENDS('INPPROS: END_PROPERTIES CARD NOT FOUND')
C
      ENDIF         ! ievfi.eq.0
C
      RETURN
  900 FORMAT(//,6X,18HELEMENT PROPERTIES,/,6X,18('-'))
  901 FORMAT(' ERROR IN OPENING PROPERTIES INPUT FILE 242 (52 LINUX)')
c 910 FORMAT(10I5)
c 920 FORMAT(//,
c    . 10X,40H MATERIAL NUMBER :                      ,10X,I5,/,
c    . 10X,16H ---------------                               ,/,
c    . 10X,40H MATERIAL MODEL  :                      ,10X,I5,/)
      END
