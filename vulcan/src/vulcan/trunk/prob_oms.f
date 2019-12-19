C=============================================================== PROBLEM
C**** DATBASSA
C
      INTEGER*4       IDATPS,IDATCS,
     .                LENRCS,NLENCS,NLENPS,NRECCS,NRECGS,NRECPS,NWORPS
C
      COMMON/DATBASSA/IDATPS(12,5),IDATCS(10),
     .                LENRCS,NLENCS,NLENPS,NRECCS,NRECGS,NRECPS,NWORPS
C--------------------------------------------------------------- PROBLEM
C**** DIMEMNSA
C
      INTEGER*4       NDIMES,NELEMS,NFUNCS,NGRUPS,NHLODS,NHISTS,NMATSS,
     .                NPOINS,NPRELS,NPROPS,NSTR1S,NTOTVS,NSUBFS,NTOTGS,
     .                NTOTVMS
C
      COMMON/DIMENNSA/NDIMES,NELEMS,NFUNCS,NGRUPS,NHLODS,NHISTS,NMATSS,
     .                NPOINS,NPRELS,NPROPS,NSTR1S,NTOTVS,NSUBFS,NTOTGS,
     .                NTOTVMS
C--------------------------------------------------------------- PROBLEM
C**** EVOLUTIO
C
      INTEGER*4       IEVFI
C
      COMMON/EVOLUTIO/IEVFI
C--------------------------------------------------------------- PROBLEM
C**** ELDATASA
C
      INTEGER*4       IDATAS,IPREVS,ISTATS,IMATXS
C
      COMMON/ELDATASA/IDATAS(11),IPREVS(3),ISTATS(4),IMATXS(7)
C--------------------------------------------------------------- PROBLEM
C***ELEMNTSA
C
      INTEGER*4       NDOFCS,NDOFNS,NEVABS,NEVACS,NGAUSS,NKONDS,NKOVAS,
     .                NMOVAS,
     .                NNODES,NDATAS,NPREVS,NSTATS,NMATXS,
     .                NDOFCMS
C
      COMMON/ELEMNTSA/NDOFCS,NDOFNS,NEVABS,NEVACS,NGAUSS,NKONDS,NKOVAS,
     .                NMOVAS,
     .                NNODES,NDATAS,NPREVS,NSTATS,NMATXS,
     .                NDOFCMS
C--------------------------------------------------------------- PROBLEM
C**** HOURGLSA
C**** HOURGLSB
C
      INTEGER*4       NHOURS,KELASS
      REAL*8          HPARAS
C
      COMMON/HOURGLSA/NHOURS,KELASS
      COMMON/HOURGLSB/HPARAS
C--------------------------------------------------------------- PROBLEM
C**** PROBLMSA
C
      INTEGER*4       KDYNAS,KPORES,KPOSTS,KPROBS,KSGAUS,KSMUSS,KTEMPS,
     .                LARGES
C
      COMMON/PROBLMSA/KDYNAS,KPORES,KPOSTS,KPROBS,KSGAUS,KSMUSS,KTEMPS,
     .                LARGES
C--------------------------------------------------------------- PROBLEM
C**** SOLVERSA
C**** SOLVERSB
C
      INTEGER*4       KRENUS,KSOLVS,KSYMMS,NWIDTS,MITCGS,NBUFAS,NPRIRS
      REAL*8          TOLCGS,TOLC1S
C
      COMMON/SOLVERSA/KRENUS,KSOLVS,KSYMMS,NWIDTS,MITCGS,NBUFAS,NPRIRS
      COMMON/SOLVERSB/TOLCGS,TOLC1S
C--------------------------------------------------------------- PROBLEM
C**** TITLESSA
C
      CHARACTER*8     TITLES,SUBTIS
C
      COMMON/TITLESSA/TITLES(8),SUBTIS(8)
C--------------------------------------------------------------- PROBLEM
C**** CONVECTIONS
C
      INTEGER*4      ICONVS,IGALES,IUPWIS,IPERTS,ISUPWS,ICUBICS,IGALFAS
      COMMON/CONVE1S/ICONVS,IGALES,IUPWIS,IPERTS,ISUPWS,ICUBICS,IGALFAS
C
      REAL*8         EXTUPS
      COMMON/CONVE2S/EXTUPS
C--------------------------------------------------------------- PROBLEM
C**** POROSITYS
C
      INTEGER*4      NPREAS,NPOROS,NPRE1S
      COMMON/POROS1S/NPREAS,NPOROS,NPRE1S
C--------------------------------------------------------------- PROBLEM
C**** ACTIVE ELEMENTS
C
      INTEGER*4      NACTIS
      COMMON/ACTIVS1/NACTIS
C--------------------------------------------------------------- PROBLEM
C**** NAMESA
C
C    To add names, see: frodess.f, frontss.f & skydess.f
C
      CHARACTER*80  CAS,CBS,CCS,CDS,CES,CFS,CGS,CHS,CIS,CJS,CKS,CLS,
     .              CMS,CNS,COS,CPS,CQS,CRS,CSS,CTS,
     .              CA1S,CB1S,CC1S,CD1S,CE1S,CF1S,CG1S,
     .              C1S,C2S,C3S,C4S,C5S,C6S,C7S,C8S,C9S,C10S
 
      COMMON/NAMESA/CAS,CBS,CCS,CDS,CES,CFS,CGS,CHS,CIS,CJS,CKS,CLS,
     .              CMS,CNS,COS,CPS,CQS,CRS,CSS,CTS,
     .              CA1S,CB1S,CC1S,CD1S,CE1S,CF1S,CG1S,
     .              C1S,C2S,C3S,C4S,C5S,C6S,C7S,C8S,C9S,C10S
c
c     CHARACTER*80  CBS,CDS,CES,CFS,CGS,CJS,
c    .              C1S
C
c     COMMON/NAMESA/CBS,CDS,CES,CFS,CGS,CJS,
c    .              C1S
C
      INTEGER*4      INUS1S,ILUS1S
      COMMON/USEFULS/INUS1S,ILUS1S
C
C**** LOGUNS
C
C     To add units, see: chek01s.f, froapps.f, frodess.f, froelis.f,
C                        frontss.f, frotaps.f, inpceks.f,
C                        pcgites.f, renum0s.f,
C                        skydess.f, sktgras.f, skyites.f,
C                        solvers.f, updates.f
C
c     INTEGER*4     LUDTSS,LUSOLS,LUFROS,LUFRHS,LUDATS,LUPRIS,LURESS,
c    .              LUSO2S,LUFR2S,LUPOSS,LURSTS,LUBFGS,LUPIPS,LUPANS,
c    .              LUGEOS,LUSETS,LUMATS,LUINIS,LULOAS,LUFIXS,LUADVS,
c    .              LUACTS,LUFANS,
c    .              LUCU1S,LUCU2S,LUCU3S,LUCU4S,LUCU5S,LUCU6S,LUCU7S,
c    .              LUCU8S,LUCU9S,LUC10S
C
c     COMMON/LOGUNT/LUDTSS,LUSOLS,LUFROS,LUFRHS,LUDATS,LUPRIS,LURESS,
c    .              LUSO2S,LUFR2S,LUPOSS,LURSTS,LUBFGS,LUPIPS,LUPANS,
c    .              LUGEOS,LUSETS,LUMATS,LUINIS,LULOAS,LUFIXS,LUADVS,
c    .              LUACTS,LUFANS,
c    .              LUCU1S,LUCU2S,LUCU3S,LUCU4S,LUCU5S,LUCU6S,LUCU7S,
c    .              LUCU8S,LUCU9S,LUC10S
c
      INTEGER*4     LUSOLS,LUFRHS,LUDATS,LUPRIS,LURESS,LUPOSS,
     .              LUCU1S,LUINFS
C
      COMMON/LOGUNS/LUSOLS,LUFRHS,LUDATS,LUPRIS,LURESS,LUPOSS,
     .              LUCU1S,LUINFS
C
C=============================================================== PROBLEM
