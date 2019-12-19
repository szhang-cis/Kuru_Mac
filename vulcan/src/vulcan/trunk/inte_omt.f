C=============================================================== INTERVL
C**** CONSTNTA
      REAL*8          DITERT,DTIMET,GRAVYT,GVECTT,STIFIT,
     .                TALFAT,TBETAT,TOLERT,XTIMET,WLUMPT
C
      COMMON/CONSTNTA/DITERT,DTIMET,GRAVYT,GVECTT(3),STIFIT,
     .                TALFAT,TBETAT,TOLERT,XTIMET,WLUMPT
C--------------------------------------------------------------- INTERVL
C**** CURRENTA
C**** CURRENTB
      INTEGER*4       ITIMET,ISTEPT,IITERT,KDTIMT,KRESLT,KSTIFT,KUNLDT,
     .                NCHEKT
      REAL*8          ARCLNT,FACTOT,FACPRT,PITERT,PVALUT,STICUT,
     .                STIINT,TFACTT,TTIMET
C
      COMMON/CURRENTA/ITIMET,ISTEPT,IITERT,KDTIMT,KRESLT,KSTIFT,KUNLDT,
     .                NCHEKT
      COMMON/CURRENTB/ARCLNT,FACTOT,FACPRT,PITERT,PVALUT,STICUT,
     .                STIINT,TFACTT,TTIMET
C--------------------------------------------------------------- INTERVL
C**** INTERVTA
      INTEGER*4       IMPOST,KALGOT,KARCLT,KINTET,KOPTIT,KCONVT,KSAVET,
     .                LACCET,LAUTOT,LINEST,MITERT,NALGOT,NBACKT,NCDIST,
     .                NSTEPT,NOUTPT
C
      COMMON/INTERVTA/IMPOST,KALGOT,KARCLT,KINTET,KOPTIT,KCONVT,KSAVET,
     .                LACCET,LAUTOT,LINEST,MITERT,NALGOT,NBACKT,NCDIST,
     .                NSTEPT,NOUTPT(50)
C=============================================================== INTERVL
