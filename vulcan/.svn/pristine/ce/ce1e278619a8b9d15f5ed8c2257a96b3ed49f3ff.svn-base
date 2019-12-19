      SUBROUTINE RSSETC
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
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
C
      GO TO (101,102,102,101,101,101,101,101), NMACHIM
C
  101 CONTINUE
      LENWD=LENRC/4
      GO TO 1000
C
  102 CONTINUE
      LENWD=LENRC
      GO TO 1000
C
 1000 CONTINUE
C
C**** EVALUATE FIRST RECORD FOR CONVERGED RESULTS AREA
C
      NRECG=4
C
C............................................LPNTN(NPOIN)
      NMREC=((NPOIN-1)/LENWD+1)*KRENU
      NRECG=NRECG+NMREC
C............................................LNODS(NNODE,NELEM)
      NMREC=(NNODE*NELEM-1)/LENWD+1
      NRECG=NRECG+NMREC
C............................................MATNO(NELEM)
      NMREC=(NELEM-1)/LENWD+1
      NRECG=NRECG+NMREC
C............................................PROEL(NPREL,NGRUP)
      NMREC=(NPREL*NGRUP*2-1)/LENWD+1
      NRECG=NRECG+NMREC
C............................................PROPS(NPROP,NMATS)
      NMREC=(NPROP*NMATS*2-1)/LENWD+1
      NRECG=NRECG+NMREC
C
      NRECC=NRECG+IDATP(1,2)*NELEM
C
C**** EVALUATE IDATC
C
C.....[  < ELPRE(NPREV) > + < ELVAR(NSTAT) >  ]  x NELEM
C
      IDATC(1)=6   ! heading (1 record ) + commons 'inte_om.f' (1r) +
                   ! common /ploter/ (4r)
C
C.....DISTO(NTOTV,1:3)
C
      IDATC(2)=IDATC(1)+(IDATP(2,2)+IDATP(3,2))*NELEM
C
C.....DISPR(NTOTV,1:3)
C
      NMREC=((NTOTV*NDISO)*2-1)/LENWD+1
      IDATC(3)=IDATC(2)+NMREC
C
C.....HEADS(NPOIN,1:4)
C
      NMREC=((NTOTV*NDISR)*2-1)/LENWD+1
      IDATC(4)=IDATC(3)+NMREC
C
C.....REFOR(NTOTV)
C
      KPOREC=0
      IF(KPORE.NE.0) KPOREC=1
      NMREC=((NPOIN*4)*2-1)/LENWD+1
      IDATC(5)=IDATC(4)+NMREC*KPOREC
C
C.....TLOAD(NTOTV,1:2)
C
      NMREC=(NTOTV*2-1)/LENWD+1
      IDATC(6)=IDATC(5)+NMREC
C
C.....HTLOD(NHLOD,NSUBF,NFUNC)
C
      NMREC=((NTOTV*2)*2-1)/LENWD+1
      IDATC(7)=IDATC(6)+NMREC
C
C.....IFFIX(NTOTV,1:4)
C
      NMREC=(NHLOD*NSUBF*NFUNC*2-1)/LENWD+1
      IDATC(8)=IDATC(7)+NMREC
C
C.....RLOAD(NTOTV)
C
      KPOREA=1
      IF(KPORE-1.GT.0) KPOREA=2            ! smoothing of water pressure
      NMREC=(NTOTV*KPOREA*(1+NACTI)-1)/LENWD+1
      IDATC(9)=IDATC(8)+NMREC
C
C.....LACTI(NELEM)
C
      NMREC=(NTOTV*2-1)/LENWD+1
      IDATC(10)=IDATC(9)+NMREC
C
C.....LAST RECORD
C
      NMREC=(NELEM-1)/LENWD+1
      IDATC(11)=IDATC(10)+NMREC*NACTI
C
      NLENC=IDATC(11)
C
      RETURN
      END
