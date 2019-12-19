      SUBROUTINE PROPMIC1(ITAPET,PROPSS,NPOI2T,INUPM)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MICROSTRUCTURAL MATERIAL PROPERTIES FOR
C     MODEL 1
C                                                      
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
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
C**** THIS IS A TEST MODEL
C
      NPRINT=0
      CALL LISTENS('PROPMIC1',NPRINT,ITAPET)
      NPOI2T=NPOI2T+1
      PROPSS(NPOI2T)=PARAMS(1)       ! Ts
      NPOI2T=NPOI2T+1
      PROPSS(NPOI2T)=PARAMS(2)       ! Tl
      NPOI2T=NPOI2T+1
      PROPSS(NPOI2T)=PARAMS(3)       ! L
      NPOI2T=NPOI2T+1
      PROPSS(NPOI2T)=PARAMS(4)       ! IMECHMIC=0 for this model
C
      IF(PARAMS(1).GT.PARAMS(2))
     . CALL RUNENDS('PROPMIC1: Tl < Ts              ')
      IF(PARAMS(3).LT.0.0)
     . CALL RUNENDS('PROPMIC1: LATENT HEAT < 0.0    ')
C
      WRITE(LURESS,810) PROPSS(NPOI2T-2)
      WRITE(LURESS,811) PROPSS(NPOI2T-1)
      WRITE(LURESS,812) PROPSS(NPOI2T)
      WRITE(LURESS,809)
C
C**** OUTPUT (see pointes.f)
C
      KPLA1S=1
      IPLLLS(INUPM)=1
C
      RETURN
  809 FORMAT(/)
  810 FORMAT(3X,'SOLIDUS TEMPERATURE  =',E15.6)
  811 FORMAT(3X,'LIQUIDUS TEMPERATURE =',E15.6)
  812 FORMAT(3X,'SPECIFIC LATENT HEAT =',E15.6)
      END
