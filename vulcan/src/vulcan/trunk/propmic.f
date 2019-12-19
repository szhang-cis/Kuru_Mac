      SUBROUTINE PROPMIC(ITAPET,PROPSS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MICROSTRUCTURAL MATERIAL PROPERTIES FOR
C     SOLID MODEL (ELEMENT NO. 5)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
C
C**** MICROSCOPICAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      DIMENSION PROPSS(*)
C
C**** LOOKS FOR 'PHASE_CHANGE' CARD
C
C     Note: the condition NLINET.LE.5 can be changed if necessary;
C           see dimensions of VPLAT in auxl_omt.f & IPLAC in nued_om.f
C
      NPRINT=0
      CALL LISTENS('PROPMIC',NPRINT,ITAPET)
      IF(WORDSS(1).NE.'PHASE')
     . CALL RUNENDS('PROPMIC: WRONG MATERIAL_DATA CARD')
      NLINET=INT(PARAMS(1))
      IF(NLINET.EQ.0.OR.NLINET.GT.5)
     . CALL RUNENDS('PROPMIC: 1) WRONG NUMBER OF PHASE_CHANGES')
C
      KPOROX=0
      KPOROS=INT(PROPSS(1))
C
      NPOI1T=2                            ! 1=KPOROS
      NPOI2T=NPOI1T
      PROPSS(NPOI2T)=PARAMS(1)            ! nlinet
C
      IF(IMICR.GT.0) THEN
       NLINEA=NFPCH/2
       IF(NLINET.GT.NLINEA)
     .  CALL RUNENDS('PROPMIC: 2) WRONG NUMBER OF PHASE_CHANGES')
       CALL RUNMENS('WARNING: CHECK THE CORRESPONDANCE BETWEEN THERMAL &
     . MICROSTRUCTURAL PHASE-CHANGES    ')
      ENDIF
C
      DO ILINE=1,NLINET
       DO JLINE=1,NLINET
        DO JPLAC=1,NPLAC
         IPLAC(ILINE,JLINE,JPLAC)=0   ! initialization of coupling array
        ENDDO
       ENDDO
      ENDDO
C
      INUPC=0
      INUPM=0
C
      NPRINT=1
      DO LINEAT=1,NLINET
       WRITE(LURESS,810) LINEAT
       CALL LISTENS('PROPMIC',NPRINT,ITAPET)
C
       ICONTR1=0
       IF(WORDSS(1).EQ.'MACRO') THEN
        ICONTR1=1
        NPOI2T=NPOI2T+1
        PROPSS(NPOI2T)=0.0D0             ! IPCFO
        INUPC=INUPC+1
        WRITE(LURESS,811)
C
C**** PHASE-CHANGE ASSOCIATED TO THE POROSITY CRITERIA
C
        IF(NPOROS.GT.0) THEN
         IF(KPOROS.GT.0) THEN
          IF((INUPM+INUPC).EQ.KPOROS) THEN
           INDEXP=0
           KPOROX=INUPC
          ENDIF
         ENDIF
        ENDIF
       ENDIF                   ! 'macro'
C
       IF(WORDSS(1).EQ.'MICRO') THEN
        ICONTR1=1
        IF(WORDSS(2).EQ.'MODEL') THEN
         IPCMO=INT(PARAMS(1))
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=1.0D0  ! IPCFO
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1)
         INUPM=INUPM+1
         WRITE(LURESS,812)
C
         ICONTR2=0
         IF(IPCMO.EQ.1) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC1(ITAPET,PROPSS,NPOI2T,INUPM)
         ENDIF
C
         IF(IPCMO.EQ.2) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC2(ITAPET,PROPSS,NPOI2T,INUPM)
         ENDIF
C
         IF(IPCMO.EQ.3) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC3(ITAPET,PROPSS,NPOI2T,INUPM)
         ENDIF
C
         IF(IPCMO.EQ.4) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC4(ITAPET,PROPSS,NPOI2T,INUPM)
         ENDIF
C
         IF(IPCMO.EQ.5) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC5(ITAPET,PROPSS,NPOI2T,INUPM)
         ENDIF
C
         IF(IPCMO.EQ.6) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC6(ITAPET,PROPSS,NPOI2T,INUPM,NLINET)
         ENDIF
C
         IF(IPCMO.EQ.7) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC7(ITAPET,PROPSS,NPOI2T,INUPM)
         ENDIF
C
         IF(IPCMO.EQ.8) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          INUPM=INUPM-1                  ! no phase-change in this model
          CALL PROPMIC8(ITAPET,PROPSS,NPOI2T,INUPM)
         ENDIF
C
         IF(IPCMO.EQ.9) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC9(ITAPET,PROPSS,NPOI2T,INUPM)
         ENDIF
C
         IF(IPCMO.EQ.10) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC10(ITAPET,PROPSS,NPOI2T,INUPM,NLINET)
         ENDIF
C
         IF(IPCMO.EQ.11) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC11(ITAPET,PROPSS,NPOI2T,INUPM,NLINET)
         ENDIF
C
         IF(IPCMO.EQ.12) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC12(ITAPET,PROPSS,NPOI2T,INUPM,NLINET)
         ENDIF
C
         IF(IPCMO.EQ.13) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL RUNENDS('ERROR: MODEL 13 NOT IMPLEMENTED YET')
         ENDIF
C
         IF(IPCMO.EQ.14) THEN
          WRITE(LURESS,813) IPCMO
          ICONTR2=1
          CALL PROPMIC14(ITAPET,PROPSS,NPOI2T,INUPM,NLINET)
         ENDIF
C
         IF(ICONTR2.EQ.0)
     .    CALL RUNENDS('PROPMIC: WRONG MODEL NUMBER')
C
        ELSE
         CALL RUNENDS('PROPMIC: WRONG MODEL DATA')
        ENDIF
C
C**** PHASE-CHANGE ASSOCIATED TO THE POROSITY CRITERIA
C
        IF(NPOROS.GT.0) THEN
         IF(KPOROS.GT.0) THEN
          IF((INUPM+INUPC).EQ.KPOROS) THEN
           INDEXP=1
           KPOROX=INUPM
          ENDIF
         ENDIF
        ENDIF
       ENDIF                   ! 'micro'
C
       IF(ICONTR1.EQ.0)
     .  CALL RUNENDS('PROPMIC: WRONG MACRO-MICRO DATA')
C
      ENDDO                    ! lineat=1,nlinet
C
C**** CHECKS NUMBER OF PHASE-CHANGES
C
C     a) macroscopic phase-changes
C
      NLINEA=NFPCH/2
      IF(INUPC.GT.NNUPCS) NNUPCS=INUPC
      IF(NNUPCS.GT.NLINEA)
     . CALL RUNENDS('PROPMIC: 3) WRONG NUMBER OF PHASE_CHANGES')
C
C     b) microscopic phase-changes
C
      IF(INUPM.GT.NNUPM) NNUPM=INUPM
      IF(NNUPCS.GT.NLINEA)
     . CALL RUNENDS('PROPMIC: 4) WRONG NUMBER OF PHASE_CHANGES')
C
C     c) total phase-changes
C
      NNUPTS=NNUPCS+NNUPM
      IF(NNUPTS.GT.NLINEA)
     . CALL RUNENDS('PROPMIC: 5) WRONG NUMBER OF PHASE_CHANGES')
C
C**** PHASE-CHANGE ASSOCIATED TO THE POROSITY CRITERIA
C
      IF(NPOROS.GT.0) THEN
       IF(KPOROS.GT.0) THEN
        IF(INDEXP.EQ.1) THEN
         KPOROX=NNUPCS+KPOROX
        ENDIF
        PROPSS(1)=KPOROX
       ENDIF
      ENDIF
C
C**** DEALS WITH COUPLING BETWEEN MICROSTRUCTURAL MODELS
C
      CALL MICCOUP
C
C**** LOOKS FOR 'END_MATERIAL' CARD
C
      NPRINT=0
      CALL LISTENS('PROPMIC',NPRINT,ITAPET)
      IF(WORDSS(1).NE.'END_M')
     . CALL RUNENDS('PROPMIC: END_MATERIAL CARD NOT FOUND   ')
      WRITE(LURESS,809)
C
C**** CONTROLS DIMENSION OF PROPSS
C
      IF(NPOI2T.GT.NPROPM) THEN
       WRITE(LURESS,900) NPROPM,NPOI2T
       CALL RUNENDS('PROPMIC: TOO MANY PROPERTIES TO READ')
      ENDIF
C
      RETURN
c 804 FORMAT(E15.6,10X,E15.6)
  809 FORMAT(/)
  810 FORMAT(/,3X,'PHASE-CHANGE NUMBER  =',I5,/)
  811 FORMAT(/,3X,'MACROSCOPICAL PHASE-CHANGE',/)
  812 FORMAT(/,3X,'MICROSCOPICAL PHASE-CHANGE',/)
  813 FORMAT(3X,'MODEL  =',I5,/)
  900 FORMAT(//,'TOO MANY MATERIAL PROPERTIES TO READ:',/,
     .      20X,'NPROP =',I5,5X,'NUMBE =',I5,/)
      END
