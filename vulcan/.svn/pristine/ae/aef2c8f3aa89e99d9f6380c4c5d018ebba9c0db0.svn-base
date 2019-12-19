      SUBROUTINE ASELMTT(MATNOT,PROELT,PROPST,CSTIFT,ESTIFT,WSTIFT,
     .                   PSTIFT,QSTIFT,HSTIFT)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES ELEMENTAL "EFFECTIVE" MATRIX   
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION MATNOT(*), PROELT(*), PROPST(*)
      DIMENSION CSTIFT(*), ESTIFT(*), WSTIFT(*),
     .          PSTIFT(*), QSTIFT(*), HSTIFT(*) 
C
      CALL CPUTIMT(TIME1T)
C
      CALL ASELM3T(MATNOT,PROELT,PROPST,CSTIFT,ESTIFT,WSTIFT)
C
      CALL CPUTIMT(TIME2T)
      CPUAST=CPUAST+(TIME2T-TIME1T)
C
      RETURN
      END
