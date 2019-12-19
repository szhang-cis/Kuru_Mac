      SUBROUTINE MICROM6(ROTET,DEROT,LIQUI,DELLS,
     .                   TEMPG,DTEMG,SHAPE,FPCHL,
     .                   ROTES,DEROS,LIQUS,DELSS,
     .                   YOUNG,CCERO,CCEROM,CCEROP,STRRB,STRRC,
     .                   IPLAT)
C***********************************************************************
C
C     THIS ROUTINE EVALUATES THE RATE PHASE CHANGE FOR MICROS 6
C     VERSION 3
C
C     Note: the mechanical properties as a function of the
C           microstructural are evaluated in micros7.f and transferred
C           to the mechanical problem as NNUPO variables
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

C--   COUPLING VARIABLES
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'

C--   MECHANICAL VARIABLES
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION SHAPE(*), FPCHL(NFPCH,*)
      DIMENSION FPAUX(27)
      
      IPCFOM=INT(VPLATM(IPLAT,4))
      IPCMOM=INT(VPLATM(IPLAT,5))
      INDFC=INT(VPLATM(IPLAT,6)) !index for fpchl
      IVERSI=INT(VPLATM(IPLAT,8))

c*****************************************************************
c     
c     VERSION 3 (Austempered Ductile Iron model)
c     
c*****************************************************************
      if(IVERSI.eq.3) then
         nmm=INT(VPLATM(IPLAT,9)) !mechanical model version

cccccccccccccccccccccccccccccccccccccccccccccc
c--   mechanical model 1 (bhadeshia's model)
cccccccccccccccccccccccccccccccccccccccccccccc
         if (nmm.eq.1) then
            
c--   of vplatm to mechanical parameters
            TBF=VPLATM(IPLAT,10)
            alphaf=VPLATM(IPLAT,11)
            alphaa=VPLATM(IPLAT,12)
            
C--   variable from microstructural model
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
            
            vfBF=0.D0           !volumetric fraction of bainite ferrite
                                !(RECOVERS THE PHASE-CHANGE FUNCTION)
            DO INODE=1,NNODL
               vfBF=vfBF+SHAPE(INODE)*FPCHL(INDFC,INODE)
            ENDDO
            
c--   solve model 1
            if (vfBF.gt.0.D0) then
               
c     other varibles from microstructural model (RECOVERS THE ADDITIONAL
c     MICROSTRUCTURAL VARIABLES)
               IF(INTERC.EQ.1) THEN
                  DO INDFC=1,NNUPO
                     DO INODL=1,NNODL
                        FPAUX(INODL)=FPCHL(2*NNUPTM+INDFC,INODL)
                     ENDDO
                     CALL SMOMID(FPAUX,NDIME,NNODL,NQUTR)
                     DO INODL=1,NNODL
                        FPCHL(2*NNUPTM+INDFC,INODL)=FPAUX(INODL)
                     ENDDO
                  ENDDO
               ENDIF
               
               vfA=0.D0         !volumetric fraction of austenite
               DO INODE=1,NNODL
                  vfA=vfA+SHAPE(INODE)*FPCHL(2*NNUPTM+1,INODE)
               ENDDO
               
               cAo=0.D0         !initial carbon concentration of austenite
               DO INODE=1,NNODL
                  cAo=cAo+SHAPE(INODE)*FPCHL(2*NNUPTM+2,INODE)
               ENDDO
               
               cAR=0.D0         !carbon concentration of residual austenite
               DO INODE=1,NNODL
                  cAR=cAR+SHAPE(INODE)*FPCHL(2*NNUPTM+3,INODE) !uncomment this line
               ENDDO
               
c     control of temperature. It is necessary because the phase change strain
c     doesn't change when phase change finishes. It finishes when temp is less
c     than TBF
               TGAUST=TEMPG
               if (TEMPG.lt.TBF) TGAUST=TBF

c     lattice parameters at ambient temperature [Angstrom]
               aFo=2.873D0      !ferrite (bhadeshia)
c               aFo=2.865D0      !ferrite
               aAo=3.555D0+0.044D0*cAo !austenite
c               aAo=3.578D0+0.033D0*cA !austenite
               aeAo=3.555D0+0.044D0*cAR !residual austenite
c               aeAo=3.578D0+0.033D0*cAR !residual austenite
         
c     lattice parameters at TGAUST temperature
               aF=aFo*(1.D0+alphaf*(TGAUST-20.D0)) !ferrite
               aA=aAo*(1.D0+alphaa*(TGAUST-20.D0)) !austenite
               aeA=aeAo*(1.D0+alphaa*(TGAUST-20.D0)) !residual austenite

c     phase change strain 
               beta1=(2.D0*aF**3.D0-aeA**3.D0)/(1.D0+
     .              2.D0*vfA*aF**3.D0/(vfBF*aeA**3.D0))
               STRRB=STRRB+(1.D0-(aA**3.D0)/(beta1+aeA**3.D0))/3.D0
            endif
         endif
      endif
      
      RETURN
      END
