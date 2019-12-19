      SUBROUTINE ASELMTS(MATNOS,PROELS,PROPSS,CSTIFS,ESTIFS,WSTIFS,
     .                   PSTIFS,QSTIFS,HSTIFS)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES ELEMENTAL "EFFECTIVE" MATRIX   
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION MATNOS(*), PROELS(*), PROPSS(*)
      DIMENSION CSTIFS(*), ESTIFS(*), WSTIFS(*),
     .          PSTIFS(*), QSTIFS(*), HSTIFS(*) 
C
      CALL CPUTIMS(TIME1S)
C
      CALL ASELM3S(MATNOS,PROELS,PROPSS,CSTIFS,ESTIFS,WSTIFS)
C
      CALL CPUTIMS(TIME2S)
      CPUASS=CPUASS+(TIME2S-TIME1S)
C
      RETURN
      END
