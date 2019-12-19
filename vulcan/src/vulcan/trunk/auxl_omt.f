C=============================================================== AUXLIAR
C**** INDWORTA
      INTEGER*4       IAQUAT,IASEMT,IFIXYT,IFORCT,
     .                IFORDT,IFORWT,IGSMOT,IIPCWT,
     .                ILOSUT,IQUAST,ISETMT,ISOLVT,
     .                ISTART,ISTIFT
C
      COMMON/INDWORTA/IAQUAT(10),IASEMT(10),IFIXYT(20),IFORCT(60),
     .                IFORDT(50),IFORWT(10),IGSMOT(50),IIPCWT(10),
     .                ILOSUT(50),IQUAST(10),ISETMT(30),ISOLVT(20),
     .                ISTART(20),ISTIFT(50)
C--------------------------------------------------------------- AUXLIAR
C**** MEMORYTA
      INTEGER*4       IADPRT,LPRINT,IADW1T,LWOR1T,IADSOT,LSOLVT,IRELET,
     .                LBYTST,IADDBT,LDABAT,IREL1T,LBYTBT,NBLIMT,LSTIFT
C
      COMMON/MEMORYTA/IADPRT,LPRINT,IADW1T,LWOR1T,IADSOT,LSOLVT,IRELET,
     .                LBYTST,IADDBT,LDABAT,IREL1T,LBYTBT,NBLIMT,LSTIFT
C--------------------------------------------------------------- AUXLIAR
C**** NEWINTTA
      INTEGER*4       NEWBOT,NEWLOT,NEWSTT,NEWFUT,NEWADT,NEWACT,
     .                IOFIXT,IOLOAT,              IOADVT,IOACTT,IOINIT,
     .                IOSTRT,ISUVOACT
C
      COMMON/NEWINTTA/NEWBOT,NEWLOT,NEWSTT,NEWFUT,NEWADT,NEWACT,
     .                IOFIXT,IOLOAT,              IOADVT,IOACTT,IOINIT,
     .                IOSTRT,ISUVOACT
C--------------------------------------------------------------- AUXLIAR
C**** PROCESTA
      INTEGER*4       IELEMT,LMATST,NDIMLT,NNODLT,NRULET,
     .                NGAULT,NTYPET,NCRITT,NKOSTT,NSTRET,NSTRST,NNODST,
     .                NQUTRT,LGRUPT
C
      COMMON/PROCESTA/IELEMT,LMATST,NDIMLT,NNODLT,NRULET,
     .                NGAULT,NTYPET,NCRITT,NKOSTT,NSTRET,NSTRST,NNODST,
     .                NQUTRT,LGRUPT
C--------------------------------------------------------------- AUXLIAR
C**** RESIDUTA
      REAL*8          AALPHT,GZEROT,GCURNT,AACOET,BBCOET
C
      COMMON/RESIDUTA/AALPHT,GZEROT,GCURNT,AACOET,BBCOET
C--------------------------------------------------------------- AUXLIAR
C**** RSTARTTA
C**** RSTARTTB
      INTEGER*4       INITIT,IRESTT,ISAVET,ISKIPT,KTIMET,KSTEPT,KTSTET,
     .                NWPOST,KPPCGT,KPPCNT
      REAL*8          TIMSTT,TLIMTT
C
      COMMON/RSTARTTA/INITIT,IRESTT,ISAVET,ISKIPT,KTIMET,KSTEPT,KTSTET,
     .                NWPOST,KPPCGT,KPPCNT
      COMMON/RSTARTTB/TIMSTT,TLIMTT
C--------------------------------------------------------------- AUXLIAR
C**** RUNTIMTA
      REAL*8          CPUINT,CPUDAT,CPUSTT,CPUSFT,CPUAST,CPUSOT,CPURET,
     .                CPURST,CPUOUT,CPURNT
C
      COMMON/RUNTIMTA/CPUINT,CPUDAT,CPUSTT,CPUSFT,CPUAST,CPUSOT,CPURET,
     .                CPURST,CPUOUT,CPURNT
C--------------------------------------------------------------- AUXLIAR
C**** SMIPTIT
      INTEGER*4  MNU4TX,MNU4TMX
      PARAMETER (MNU4TX=400)                  ! change if necessary
      PARAMETER (MNU4TMX=11*MNU4TX+46)        ! estimated; see setdatd.f
C
      INTEGER*4      IPLAST,IPLAOT,IPLANT,IPLAMT,
     .               IPLUAT,IPLUOT,IPLLLT,IPLXXT
      INTEGER*4      NNUINT,NNUNOT
      INTEGER*4      KPLA1T,KPLA2T,KPLA3T,KPLA4T,KPLA5T,
     .               KPLA6T,KPLA7T,KPLA8T,KPLA9T,KPLA10T
      INTEGER*4      NNUPC,NNUPT
C
      COMMON/SMIPT1T/IPLAST(MNU4TMX),IPLAOT(MNU4TMX),IPLANT(MNU4TMX),
     .               IPLAMT(MNU4TMX),IPLUAT(MNU4TMX),IPLUOT(MNU4TMX),
     .               IPLLLT(MNU4TMX),IPLXXT(MNU4TMX)
C
      COMMON/SMIPT2T/NNUINT,NNUNOT
      COMMON/SMIPT3T/KPLA1T,KPLA2T,KPLA3T,KPLA4T,KPLA5T,
     .               KPLA6T,KPLA7T,KPLA8T,KPLA9T,KPLA10T
      COMMON/SMIPMICT/NNUPC,NNUPT
C
      REAL*8          APLUOT
C
      COMMON/SMIPT1TA/APLUOT(MNU4TMX)
C--------------------------------------------------------------- AUXLIAR
C**** PROPERT1
      INTEGER*4       NDENS,NCAPA,NCOND,NPLAT,NSUBD,IFRENT,IFREKT,
     .                LINPC,ICACOT
C
      COMMON/PROPERT1/NDENS,NCAPA,NCOND,NPLAT,NSUBD(3,2),IFRENT,IFREKT,
     .                LINPC,ICACOT
C
C**** PROPERT2
      INTEGER*4  MPLAT,MVPLAT
      PARAMETER (MPLAT=5, MVPLAT=84)          ! change if necessary
C
      REAL*8          VDENS,VCAPA,VCOND,
     .                VPLAT,VPCFU,
     .                VPCAU,
     .                VGRO1,VGRO2,VGRO3,
     .                VGRO1E
C
      COMMON/PROPERT2/VDENS(20,2),VCAPA(500,2),VCOND(500,2,3),
     .                VPLAT(MPLAT,MVPLAT),VPCFU(MPLAT,200,3),
     .                VPCAU(MPLAT,200),
     .                VGRO1(50,2),VGRO2(100,2),VGRO3(100,2),
     .                VGRO1E(50,2)
C
C**** PROPERT3
      INTEGER*4       NRADH1,NRADH2,NRADU,NRADP,IR432,
     .                IGABO3,NRADH3,IGABO4
C
      COMMON/PROPERT3/NRADH1,NRADH2,NRADU,NRADP,IR432,
     .                IGABO3,NRADH3,IGABO4
C
C**** PROPERT4
      REAL*8          VRADH1,VRADH2,VRADU,VRADP,
     .                RIGINT,
     .                VRADH3,TENVI3,AGABO4
C
      COMMON/PROPERT4/VRADH1(22,2),VRADH2(20,2),VRADU(20,2),VRADP(20,2),
     .                RIGINT,
     .                VRADH3(20,2),TENVI3,AGABO4
C
C**** PROPERT5
      INTEGER*4       NSOTRT,ISOTRT,NDIMETO,NFILL,IFILL
C
      COMMON/PROPERT5/NSOTRT,ISOTRT,NDIMETO,NFILL,IFILL
C
C**** PROPERT6
      INTEGER*4       NDENSFI,NCAPAFI,NCONDFI,NPLATFI,NSUBDFI,
     .                IFRENFI,IFREKFI,LINPCFI
C
      COMMON/PROPERT6/NDENSFI,NCAPAFI,NCONDFI,NPLATFI,NSUBDFI(3,2),
     .                IFRENFI,IFREKFI,LINPCFI
C
C**** PROPERT7
      REAL*8          VDENSFI,VCAPAFI,VCONDFI,
     .                VPLATFI,VPCFUFI,
     .                VPCAUFI
C
      COMMON/PROPERT7/VDENSFI(20,2),VCAPAFI(500,2),VCONDFI(500,2,3),
     .                VPLATFI(MPLAT,50),VPCFUFI(MPLAT,200,3),
     .                VPCAUFI(MPLAT,200)
C
C**** PROPERT8
      INTEGER*4       NRADH1FI,NRADH2FI,NRADUFI,NRADPFI
C
      COMMON/PROPERT8/NRADH1FI,NRADH2FI,NRADUFI,NRADPFI
C
C**** PROPERT9
      REAL*8          VRADH1FI,VRADH2FI,VRADUFI,VRADPFI
C
      COMMON/PROPERT9/VRADH1FI(20,2),VRADH2FI(20,2),VRADUFI(20,2),
     .                VRADPFI(20,2)
C
C**** PROPERT10
      INTEGER*4        KPOROT,KPOROTFI
C
      COMMON/PROPERT10/KPOROT,KPOROTFI
C
C**** PROPERT11
      INTEGER*4        NNUM4T,NNUM4TM,INNUM4T,IWEUT,IMNMT
C
      COMMON/PROPERT11/NNUM4T,NNUM4TM,INNUM4T,IWEUT,IMNMT
C
C**** PROPERT12
      REAL*8           DNU, RNU, RNUZ1,
     .                 ANU, FNU, DPERC, PNUC, DELTAG, DELTAF, RAFGGN
C
      COMMON/PROPERT12/DNU(MNU4TX,3), RNU(MNU4TX),    RNUZ1(MNU4TX),
     .                 ANU(MNU4TX),   FNU(MNU4TX),    DPERC(MNU4TX,1),
     .                 PNUC(MNU4TX),  DELTAG(MNU4TX), DELTAF(MNU4TX),
     .                 RAFGGN(MNU4TX)
C--------------------------------------------------------------- AUXLIAR
C**** SMOOTHING
      INTEGER*4       ISMO1T,ISMO2T
      COMMON/SMOOTHBT/ISMO1T,ISMO2T
C--------------------------------------------------------------- AUXLIAR
C**** COSOLTRU
      INTEGER*4      NSOL1,NSOL2,NSOL3
      COMMON/COSOLIA/NSOL1,NSOL2,NSOL3
C
      REAL*8         TSOL1,TSOL2,TSOL3
      COMMON/COSOLIB/TSOL1,TSOL2,TSOL3
C
      INTEGER*4      ITRUC,ITRUCV
      COMMON/TRUCHOA/ITRUC,ITRUCV
C
      REAL*8         TRUCH,TRUCHB
      COMMON/TRUCHOB/TRUCH,TRUCHB
C
      INTEGER*4     NCETAT
      COMMON/TRUCHI/NCETAT
C=============================================================== AUXLIAR
