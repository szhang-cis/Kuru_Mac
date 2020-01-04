      SUBROUTINE RSSETCS
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES  THE RECORD ADDRES FOR DIFFERENT ARRAYS
C     STORED IN THE CONVERGED RESULTS AREA
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
C
      GO TO (101,102,102,101,101,101,101,101), NMACHIS
C
  101 CONTINUE
      LENWDS=LENRCS/4
      GO TO 1000
C
  102 CONTINUE
      LENWDS=LENRCS
      GO TO 1000
C
 1000 CONTINUE
C
C**** EVALUATE FIRST RECORD FOR CONVERGED RESULTS AREA
C
      NRECGS=4
C
C............................................LPNTN(NPOIN)
      NMRECS=((NPOINS-1)/LENWDS+1)*KRENUS
      NRECGS=NRECGS+NMRECS
C............................................LNODS(NNODE,NELEM)
      NMRECS=(NNODES*NELEMS-1)/LENWDS+1
      NRECGS=NRECGS+NMRECS
C............................................MATNO(NELEM)
      NMRECS=(NELEMS-1)/LENWDS+1
      NRECGS=NRECGS+NMRECS
C............................................PROEL(NPREL,NGRUP)
      NMRECS=(NPRELS*NGRUPS*2-1)/LENWDS+1
      NRECGS=NRECGS+NMRECS
C............................................PROPS(NPROP,NMATS)
      NMRECS=(NPROPS*NMATSS*2-1)/LENWDS+1
      NRECGS=NRECGS+NMRECS
C
      NRECCS=NRECGS+IDATPS(1,2)*NELEMS
C
C**** EVALUATE IDATCS
C
C.....[  < ELPRE(NPREV) > + < ELVAR(NSTAT) >  ]  x NELEM
C
      IDATCS(1)=6  ! heading (1 record ) + commons 'inte_oms.f' (1r) +
                   ! common /ploter/ (4r)
C
C.....DISTOS(NTOTVS,1:3)
C
      IDATCS(2)=IDATCS(1)+(IDATPS(2,2)+IDATPS(3,2))*NELEMS
C
C.....DISPRS(NTOTVS,1:3)
C
      IIII3S=(1+KDYNAS+NMEMO8S)
      NMRECS=((NTOTVS*IIII3S)*2-1)/LENWDS+1
      IDATCS(3)=IDATCS(2)+NMRECS
C
C.....HEADSS(NPOINS,1:4)
C
      NMRECS=((NTOTVS*IIII3S)*2-1)/LENWDS+1
      IDATCS(4)=IDATCS(3)+NMRECS
C
C.....REFORS(NTOTVS)
C
      KPORESC=0
      IF(KPORES.NE.0) KPORESC=1
      NMRECS=((NPOINS*4)*2-1)/LENWDS+1
      IDATCS(5)=IDATCS(4)+NMRECS*KPORESC
C
C.....TLOADS(NTOTVS,1:2)
C
      NMRECS=(NTOTVS*2-1)/LENWDS+1
      IDATCS(6)=IDATCS(5)+NMRECS
C
C.....HTLODS(NHLODS,NSUBFS,NFUNCS)
C
      NMRECS=((NTOTVS*2)*2-1)/LENWDS+1
      IDATCS(7)=IDATCS(6)+NMRECS
C
C.....IFFIXS(NTOTVS,1:2)
C
      NMRECS=(NHLODS*NSUBFS*NFUNCS*2-1)/LENWDS+1
      IDATCS(8)=IDATCS(7)+NMRECS
C
C.....RLOADS(NTOTVS)
C
      KPORESA=1
      IF(KPORES-1.GT.0) KPORESA=2          ! smoothing of water pressure
      NMRECS=(NTOTVS*KPORESA-1)/LENWDS+1
      IDATCS(9)=IDATCS(8)+NMRECS
C
C.....LAST RECORD
C
      NMRECS=(NTOTVS*2-1)/LENWDS+1
      IDATCS(10)=IDATCS(9)+NMRECS
C
      NLENCS=IDATCS(10)
C
      RETURN
      END