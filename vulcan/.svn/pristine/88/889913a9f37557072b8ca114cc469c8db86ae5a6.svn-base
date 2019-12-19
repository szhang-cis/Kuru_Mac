C=============================================================== PROBLEM
C**** DATBASTA
C
      INTEGER*4       IDATPT,IDATCT,
     .                LENRCT,NLENCT,NLENPT,NRECCT,NRECGT,NRECPT,NWORPT
C
      COMMON/DATBASTA/IDATPT(12,5),IDATCT(10),
     .                LENRCT,NLENCT,NLENPT,NRECCT,NRECGT,NRECPT,NWORPT
C--------------------------------------------------------------- PROBLEM
C**** DIMEMNTA
C
      INTEGER*4       NDIMET,NELEMT,NFUNCT,NGRUPT,NHLODT,NHISTT,NMATST,
     .                NPOINT,NPRELT,NPROPT,NSTR1T,NTOTVT,NSUBFT,NTOTGT,
     .                NTOTVM,NHIST1,ILDV1T,ILDV2T
C
      COMMON/DIMENNTA/NDIMET,NELEMT,NFUNCT,NGRUPT,NHLODT,NHISTT,NMATST,
     .                NPOINT,NPRELT,NPROPT,NSTR1T,NTOTVT,NSUBFT,NTOTGT,
     .                NTOTVM,NHIST1,ILDV1T,ILDV2T
C--------------------------------------------------------------- PROBLEM
C**** ELDATATA
C
      INTEGER*4       IDATAT,IPREVT,ISTATT,IMATXT
C
      COMMON/ELDATATA/IDATAT(11),IPREVT(3),ISTATT(4),IMATXT(7)
C--------------------------------------------------------------- PROBLEM
C***ELEMNTTA
C
      INTEGER*4       NDOFCT,NDOFNT,NEVABT,NEVACT,NGAUST,NKONDT,NKOVAT,
     .                NMOVAT,
     .                NNODET,NDATAT,NPREVT,NSTATT,NMATXT,
     .                NDOFCM
C
      COMMON/ELEMNTTA/NDOFCT,NDOFNT,NEVABT,NEVACT,NGAUST,NKONDT,NKOVAT,
     .                NMOVAT,
     .                NNODET,NDATAT,NPREVT,NSTATT,NMATXT,
     .                NDOFCM
C--------------------------------------------------------------- PROBLEM
C***HOURGLTA
C***HOURGLTB
C
      INTEGER*4       NHOURT,KELAST
      REAL*8          HPARAT
C
      COMMON/HOURGLTA/NHOURT,KELAST
      COMMON/HOURGLTB/HPARAT
C--------------------------------------------------------------- PROBLEM
C**** PROBLMTA
C
      INTEGER*4       KDYNAT,KPORET,KPOSTT,KPROBT,KSGAUT,KSMUST,KTEMPT,
     .                LARGET
C
      COMMON/PROBLMTA/KDYNAT,KPORET,KPOSTT,KPROBT,KSGAUT,KSMUST,KTEMPT,
     .                LARGET
C--------------------------------------------------------------- PROBLEM
C***SOLVERTA
C***SOLVERTB
C
      INTEGER*4       KRENUT,KSOLVT,KSYMMT,NWIDTT,MITCGT,NBUFAT,NPRIRT,
     .                MITGMT,MKRYLT,IPGMRT
      REAL*8          TOLCGT,TOLC1T,TOLGMT
C
      COMMON/SOLVERTA/KRENUT,KSOLVT,KSYMMT,NWIDTT,MITCGT,NBUFAT,NPRIRT,
     .                MITGMT,MKRYLT,IPGMRT
      COMMON/SOLVERTB/TOLCGT,TOLC1T,TOLGMT
C--------------------------------------------------------------- PROBLEM
C**** TITLESTA
C
      CHARACTER*8     TITLET,SUBTIT
C
      COMMON/TITLESTA/TITLET(8),SUBTIT(8)
C--------------------------------------------------------------- PROBLEM
C**** CONVECTION
C
      INTEGER*4     ICONVT,IGALET,IUPWI,IPERT,ISUPW,IEGFPC,ICUBIC,
     .              IUPWIC,IPERTC,ISUPWC,
     .              IUPWIG,IPERTG,ISUPWG,IGALFA
      COMMON/CONVE1/ICONVT,IGALET,IUPWI,IPERT,ISUPW,IEGFPC,ICUBIC,
     .              IUPWIC,IPERTC,ISUPWC,
     .              IUPWIG,IPERTG,ISUPWG,IGALFA
C
      REAL*8        EXTUP,GALFA,RAMON,EXTUPC,EXTUPG
      COMMON/CONVE2/EXTUP,GALFA,RAMON,EXTUPC,EXTUPG
C--------------------------------------------------------------- PROBLEM
C**** POROSITY
C
      INTEGER*4     NPREAT,NPOROT,NPRE1T
      COMMON/POROS1/NPREAT,NPOROT,NPRE1T
C--------------------------------------------------------------- PROBLEM
C**** ACTIVE ELEMENTS
C
      INTEGER*4      NACTIT
      COMMON/ACTIVT1/NACTIT
C--------------------------------------------------------------- PROBLEM
C**** CONTACT - NON-COINCIDENT MESH
C
      INTEGER*4       NOCOIT,NNODNT,LARGCT,NOCOLT
      COMMON/CONTNCMT/NOCOIT,NNODNT,LARGCT,NOCOLT
C--------------------------------------------------------------- PROBLEM
C**** NAMETA
C
C    To add names, see: frodest.f, frontst.f & skydest.f
C
      CHARACTER*80  CAT,CBT,CCT,CDT,CET,CFT,CGT,CHT,CIT,CJT,CKT,CLT,
     .              CMT,CNT,COT,CPT,CQT,CRT,CST,CTT,
     .              CA1T,CB1T,CC1T,CD1T,CE1T,CF1T,CG1T,CH1T,CI1T,
     .              C1T,C2T,C3T,C4T,C5T,C6T,C7T,C8T,C9T,C10T
C
      COMMON/NAMETA/CAT,CBT,CCT,CDT,CET,CFT,CGT,CHT,CIT,CJT,CKT,CLT,
     .              CMT,CNT,COT,CPT,CQT,CRT,CST,CTT,
     .              CA1T,CB1T,CC1T,CD1T,CE1T,CF1T,CG1T,CH1T,CI1T,
     .              C1T,C2T,C3T,C4T,C5T,C6T,C7T,C8T,C9T,C10T
C
      INTEGER*4     INUS1T,ILUS1T
      COMMON/USEFUL/INUS1T,ILUS1T
C
C**** LOGUNT
C
C     To add units, see: chek01t.f, froappt.f, frodest.f, froelit.f,
C                        frontst.f, frotapt.f,
C                        pcgitet.f, renum0t.f,
C                        skydest.f, skygrat.f, skyitet.f,
C                        solvert.f, updatet.f
C
      INTEGER*4     LUDTST,LUSOLT,LUFROT,LUFRHT,LUDATT,LUPRIT,LUREST,
     .              LUSO2T,LUFR2T,LUPOST,LURSTT,LUBFGT,LUPIPT,LUPANT,
     .              LUGEOT,LUSETT,LUMATT,LUINIT,LULOAT,LUFIXT,LUADVT,
     .              LUACTT,LUFANT,LUSTRT,LUINFT,
     .              LUCU1T,LUCU2T,LUCU3T,LUCU4T,LUCU5T,LUCU6T,LUCU7T,
     .              LUCU8T,LUCU9T,LUC10T
C
      COMMON/LOGUNT/LUDTST,LUSOLT,LUFROT,LUFRHT,LUDATT,LUPRIT,LUREST,
     .              LUSO2T,LUFR2T,LUPOST,LURSTT,LUBFGT,LUPIPT,LUPANT,
     .              LUGEOT,LUSETT,LUMATT,LUINIT,LULOAT,LUFIXT,LUADVT,
     .              LUACTT,LUFANT,LUSTRT,LUINFT,
     .              LUCU1T,LUCU2T,LUCU3T,LUCU4T,LUCU5T,LUCU6T,LUCU7T,
     .              LUCU8T,LUCU9T,LUC10T
C
C=============================================================== PROBLEM
