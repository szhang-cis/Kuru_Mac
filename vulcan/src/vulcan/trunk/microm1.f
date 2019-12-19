      SUBROUTINE MICROM1(ROTET,DEROT,LIQUI,DELLS,
     .                   TEMPG,DTEMG,SHAPE,FPCHL,
     .                   ROTES,DEROS,LIQUS,DELSS,
     .                   YOUNG,
     .                   IPLAT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE MECHANICAL PROPERTIES AS A FUNCTION OF
C     THE MICROSTRUCTURAL VARIABLES FOR MODEL 1 (IPCMOM=1)
C
C     Note: this routine is not used in this model (IPCFOM & IPCMOM
C           are not redefined in ideprom1.f)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION SHAPE(*), FPCHL(NFPCH,*)
      DIMENSION FPAUX(27)
C
C**** RECOVERS THE PHASE-CHANGE FUNCTION
C
      ILSPC=INT(VPLATM(IPLAT,1))
      ISSPC=INT(VPLATM(IPLAT,2))
      DELEX=VPLATM(IPLAT,3)
      INDFC=INT(VPLATM(IPLAT,6))      ! index for fpchl
C
      IF(INTERC.EQ.1) THEN
       DO INODL=1,NNODL
        FPAUX(INODL)=FPCHL(INDFC,INODL)
       ENDDO
       CALL SMOMID(FPAUX,NDIME,NNODL,NQUTR)
       DO INODL=1,NNODL
        FPCHL(INDFC,INODL)=FPAUX(INODL)
       ENDDO
C
       DO INODL=1,NNODL
        FPAUX(INODL)=FPCHL(INDFC+NNUPTM,INODL)
       ENDDO
       CALL SMOMID(FPAUX,NDIME,NNODL,NQUTR)
       DO INODL=1,NNODL
        FPCHL(INDFC+NNUPTM,INODL)=FPAUX(INODL)
       ENDDO
      ENDIF
C
C**** CALCULATES THE PHASE-CHANGE FUNCTION AT GAUSS POINT
C
      IF(ILSPC.GT.0) THEN
       ROTET=0.0D+00
       DEROT=0.0D+00
       LIQUI=1
       DO INODE=1,NNODL
        ROTET=ROTET+SHAPE(INODE)*FPCHL(INDFC,INODE)
        DEROT=DEROT+SHAPE(INODE)*FPCHL(INDFC+NNUPTM,INODE)
       ENDDO
       ROTET=1.0D+00-ROTET
       IF(ROTET.GT.0.0) LIQUI=0
       DELLS=DELEX
      ENDIF
C
      IF(ISSPC.GT.0) THEN
       ROTES=0.0D+00
       DEROS=0.0D+00
       LIQUS=1
       DO INODE=1,NNODL
        ROTES=ROTES+SHAPE(INODE)*FPCHL(INDFC,INODE)
        DEROS=DEROS+SHAPE(INODE)*FPCHL(INDFC+NNUPTM,INODE)
       ENDDO
       ROTES=1.0D+00-ROTES
       IF(ROTES.GT.0.0) LIQUS=0
       DELSS=DELEX
      ENDIF
C
C**** THERMAL-MECHANICAL-MICROSCOPICAL PROPERTIES (see idepros1.f)
C
      TEINF=VPLATM(IPLAT, 8)          ! Ts
      TESUP=VPLATM(IPLAT, 9)          ! Tl
      HENER=VPLATM(IPLAT,10)          ! L
      IMECHMIC=INT(VPLATM(IPLAT,11))  ! mechanical-microstructural index
C
      RETURN
      END
