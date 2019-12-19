      SUBROUTINE MICROM14(ROTET,DEROT,LIQUI,DELLS,
     .                   TEMPG,DTEMG,SHAPE,FPCHL,
     .                   ROTES,DEROS,LIQUS,DELSS,
     .                   YOUNG,CCERO,CCEROM,CCEROP,STRRB,STRRC,
     .                   IPLAT)
C***********************************************************************
C
C     THIS ROUTINE EVALUATES THE RATE PHASE CHANGE FOR MICROS 14
C
C     Note: the mechanical properties as a function of the
C           microstructural are evaluated in micros14.f and transferred
C           to the mechanical problem as NNUPO variables
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

C--   coupling variables
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
      
C--   mechanical variables
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      
      DIMENSION SHAPE(*), FPCHL(NFPCH,*)
      DIMENSION FPAUX(27)

c--   transfers VPLATM to mechanicals parameters
      IPCFOM=INT(VPLATM(IPLAT,4))
      IPCMOM=INT(VPLATM(IPLAT,5))
      INDFC=INT(VPLATM(IPLAT,6)) !index for fpchl
      nmm=INT(VPLATM(IPLAT,8))  !mechanical model versio

cccccccccccccccccccccccccccccccccccccccccccccc
c--   mechanical model 1 (constant expansion)
cccccccccccccccccccccccccccccccccccccccccccccc
      if (nmm.eq.1) then

c--   transfers VPLATM to mechanicals parameters
         svc=VPLATM(IPLAT,9)
         
C--   phase change function
c     recovers
         IF(INTERC.EQ.1) THEN
            DO INODL=1,NNODL
               FPAUX(INODL)=FPCHL(INDFC,INODL)
            ENDDO
            CALL SMOMID(FPAUX,NDIME,NNODL,NQUTR)
            DO INODL=1,NNODL
               FPCHL(INDFC,INODL)=FPAUX(INODL)
            ENDDO
c     
            DO INODL=1,NNODL
               FPAUX(INODL)=FPCHL(INDFC+NNUPTM,INODL)
            ENDDO
            CALL SMOMID(FPAUX,NDIME,NNODL,NQUTR)
            DO INODL=1,NNODL
               FPCHL(INDFC+NNUPTM,INODL)=FPAUX(INODL)
            ENDDO
         ENDIF
         
c     calculate at gaus point
         ROTES=0.0D+00
cc         DEROS=0.0D+00
         DO INODE=1,NNODL
            ROTES=ROTES+SHAPE(INODE)*FPCHL(INDFC,INODE)
cc            DEROS=DEROS+SHAPE(INODE)*
cc     .           FPCHL(INDFC+NNUPTM,INODE)
         ENDDO
cc         LIQUS=0                ! heating
cc         IF(ROTES.EQ.1.0D0) LIQUS=1
cc         IF(ROTES.LT.0.0D0) ROTES=0.0D0 ! controls
cc         IF(ROTES.GT.1.0D0) ROTES=1.0D0
cc         DELSS=DELEX

c--   solve model 1
         STRRB=STRRB+1.0D+00/3.00D+00*ROTES*svc
cc         STRRC=STRRC+1.0D+00/3.00D+00*(ROTES-DEROS)*svc
      ENDIF
      
      RETURN
      END
      
