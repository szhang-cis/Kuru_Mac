      SUBROUTINE DECAPA(ERE0,ERE1,SUMA,SUMB,SUMC,ALENG,HCAPA,CAPAD,
     .                  DECAP,PROPS,GFM,GCM,YIELD,ENELI,DELMU,RETEN,
     .                  SIKMA)
C***********************************************************************
C
C***  ESTA RUTINA CALCULA LA VARIABLE INTERNA DE DANO POR DEGRADACION:
C               CAPAD
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION PROPS(*)
C
      GSUBF=PROPS(62)
      GSUBC=PROPS(67)
C
C***  CALCULO DE LA ENERERGIA ESPECIFICA A DISIPAR
C
      DCAPA=0.0D00
      CONS0=0.0D00
      CONS1=0.0D00
      GF=GSUBF/ALENG
      GC=GSUBC/ALENG
      cons0=ere0*dabs(yield/reten)/(gf*suma)
      cons1=ere1*dabs(yield)/(gc*suma)
c      IF(YIELD.EQ.0.0)THEN
c       GFM=GF
c       GCM=GC
c      ELSE
c       GFM=GF*RETEN*(SUMB/DABS(YIELD))
c       GCM=GC*(SUMC/DABS(YIELD))
c      ENDIF
c      IF(GFM.GT.0.000001)CONS0=ERE0/GFM
c      IF(GCM.GT.0.000001)CONS1=ERE1/GCM
C
C***  CALCULO DE LA FUNCION DE ESTADO DE CAPA
C
      HCAPA=CONS0+CONS1
      DECAP=HCAPA*ENELI*DELMU
      CAPA1=CAPAD+DECAP
      IF(CAPA1.GE.1.0D00)CAPA1=1.0
      CAPAD=CAPA1
C
      RETURN
      END
