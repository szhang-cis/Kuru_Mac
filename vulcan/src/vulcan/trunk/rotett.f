      SUBROUTINE ROTETT(ROTET,DEROT,LIQUI,DELLS,
     .                  TEMPG,DTEMG,SHAPE,FPCHL,
     .                  ROTES,DEROS,LIQUS,DELSS,
     .                  STRRB,STRRC,
     .                  YOUNG,CCERO,CCEROM,CCEROP)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE SOLIDIFICATION VARIABLE
C
C
C     IPCFOM=0 TOTAL PHASE-CHANGE FORMULATION
C                     f_pc(T) is bijective
C
C                     IPCMOM=1 LINEAR
C                     IPCMOM=2 IDEM 1 WITH f_s
C                     IPCMOM=3 PARABOLIC
C                     IPCMOM=4 CUBIC
C                     IPCMOM=5 SCHEIL'S EQUATION
C                     IPCMOM=6 LEVER'S EQUATION
C                     IPCMOM=7 ....
C
C     IPCFOM=1 RATE PHASE-CHANGE FORMULATION
C                     f_pc(\alpha_m_k,T) where:
C                     \alpha_m_k=microscopical variables
C                     => the computation of f^._pc is necessary
C
C                     IPCMOM=1 MODEL ...
C                     IPCMOM=2 MODEL ...
C
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
      IF(NPLATM.EQ.0) RETURN
C
      STRRB=0.0D0
      STRRC=0.0D0
      DO IPLAT=1,NPLATM
       ITEMRO=INT(VPLATM(IPLAT,7))
       IF(ITEMRO.EQ.0) THEN
C
        IPCFOM=INT(VPLATM(IPLAT,4))
        IPCMOM=INT(VPLATM(IPLAT,5))
C
        IF(IPCFOM.EQ.0) THEN
C
C**** TOTAL PHASE-CHANGE FORMULATION (IPCFOM=0)
C
         ILSPC=INT(VPLATM(IPLAT,1))
         ISSPC=INT(VPLATM(IPLAT,2))
         DELEX=VPLATM(IPLAT,3)
         INDFC=INT(VPLATM(IPLAT,6))      ! index for fpchl
         ISEPC=INT(VPLATM(IPLAT,11))     ! sense of phase-change
C
C**** RECOVERS THE PHASE-CHANGE FUNCTION
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
         IF(IPCMOM.EQ.0) THEN                       ! constant expansion
C
          IF(ILSPC.GT.0) THEN
           ROTET=0.0D+00
           DEROT=0.0D+00
           DO INODE=1,NNODL
            ROTET=ROTET+SHAPE(INODE)*FPCHL(INDFC,INODE)
            DEROT=DEROT+SHAPE(INODE)*FPCHL(INDFC+NNUPTM,INODE)
           ENDDO
           IF(ISEPC.EQ.0) THEN                      ! cooling
            LIQUI=1
            ROTET=1.0D+00-ROTET
            IF(ROTET.GT.0.0D0) LIQUI=0
           ELSE                                     ! heating
            LIQUI=0
            IF(ROTET.EQ.1.0D0) LIQUI=1
           ENDIF
           IF(ROTET.LT.0.0D0) ROTET=0.0D0           ! controls
           IF(ROTET.GT.1.0D0) ROTET=1.0D0
           DELLS=DELEX
           STRRB=STRRB+1.0D+00/3.00D+00*ROTET*DELLS
           STRRC=STRRC+1.0D+00/3.00D+00*(ROTET-DEROT)*DELLS
          ENDIF
C
          IF(ISSPC.GT.0) THEN
           ROTES=0.0D+00
           DEROS=0.0D+00
           DO INODE=1,NNODL
            ROTES=ROTES+SHAPE(INODE)*FPCHL(INDFC,INODE)
            DEROS=DEROS+SHAPE(INODE)*FPCHL(INDFC+NNUPTM,INODE)
           ENDDO
           IF(ISEPC.EQ.0) THEN                      ! cooling
            LIQUS=1
            ROTES=1.0D+00-ROTES
            IF(ROTES.GT.0.0D0) LIQUS=0
           ELSE                                     ! heating
            LIQUS=0
            IF(ROTES.EQ.1.0D0) LIQUS=1
           ENDIF
           IF(ROTES.LT.0.0D0) ROTES=0.0D0           ! controls
           IF(ROTES.GT.1.0D0) ROTES=1.0D0
           DELSS=DELEX
           STRRB=STRRB+1.0D+00/3.00D+00*ROTES*DELSS
           STRRC=STRRC+1.0D+00/3.00D+00*(ROTES-DEROS)*DELSS
          ENDIF
C
         ELSE
          call runend('rotett: error with ipcfom=0 & ipcmom=1')
c         CALL MICPHASM
         ENDIF       ! ipcmom.eq.0
C
        ELSE         ! ipcfom.eq.0
C
C**** RATE PHASE-CHANGE FORMULATION (IPCFOM=1) > MICROSCOPICAL MODELS
C
         GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14), IPCMOM
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 1
C
    1    CALL MICROM1(ROTET,DEROT,LIQUI,DELLS,
     .                TEMPG,DTEMG,SHAPE,FPCHL,
     .                ROTES,DEROS,LIQUS,DELSS,
     .                YOUNG,
     .                IPLAT)
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 2
C
    2    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=2')
c        CALL MICROM2
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 3
C
    3    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=3')
c        CALL MICROM3
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 4
C
    4    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=4')
c        CALL MICROM4
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 5
C
    5    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=5')
c        CALL MICROM5
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 6
C
    6    CALL MICROM6(ROTET,DEROT,LIQUI,DELLS,
     .                TEMPG,DTEMG,SHAPE,FPCHL,
     .                ROTES,DEROS,LIQUS,DELSS,
     .                YOUNG,CCERO,CCEROM,CCEROP,STRRB,STRRC,
     .                IPLAT)
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 7
C
    7    CALL MICROM7(ROTET,DEROT,LIQUI,DELLS,
     .                TEMPG,DTEMG,SHAPE,FPCHL,
     .                ROTES,DEROS,LIQUS,DELSS,
     .                YOUNG,CCERO,CCEROM,CCEROP,
     .                IPLAT)
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 8
C
    8    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=8')
c        CALL MICROM8
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 9
C
    9    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=9')
c        CALL MICROM9
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 10
C
   10    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=10')
c        CALL MICROM10
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 11
C
   11    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=11')
c        CALL MICROM11
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 12
C
   12    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=12')
c        CALL MICROM12
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 13
C
   13    CONTINUE
         call runend('rotett: error with ipcfom=1 & ipcmom=13')
c        CALL MICROM13
         GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 14
C
   14    CALL MICROM14(ROTET,DEROT,LIQUI,DELLS,
     .                 TEMPG,DTEMG,SHAPE,FPCHL,
     .                 ROTES,DEROS,LIQUS,DELSS,
     .                 YOUNG,CCERO,CCEROM,CCEROP,STRRB,STRRC,
     .                 IPLAT)
         GO TO 100
C
  100    CONTINUE
C
        ENDIF        ! ipcfom.eq.0
C
       ELSE          ! itemro=1
C
C**** COMPUTES "ROTET" ACCORDING TO THE INPUT OF MECHANICAL DATA FILE
C
        TEMPS=VPLATM(IPLAT,8)
        TEMPL=VPLATM(IPLAT,9)
        CALL ROTEMR(ROTET,DEROT,LIQUI,TEMPG,DTEMG)
C
C**** ALWAYS TOTAL PHASE-CHANGE FORMULATION (IPCFOM=0) WITH IPCMOM=0
C
        ILSPC=INT(VPLATM(IPLAT,1))
        ISSPC=INT(VPLATM(IPLAT,2))
        DELEX=VPLATM(IPLAT,3)
C
        IF(ILSPC.GT.0) THEN
         DELLS=DELEX
        ENDIF
C
        IF(ISSPC.GT.0) THEN
         ROTES=ROTET
         DEROS=DEROT
         LIQUS=LIQUI
         DELSS=DELEX
        ENDIF
C
       ENDIF        ! itemro.eq.0
      ENDDO         ! iplat=1,nplatm
C
      RETURN
#endif
      END
