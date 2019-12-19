C=============================================================== INP_OUT
C**** LISTENTA:
C**** LISTENTB:
C**** LISTENTD:
        INTEGER*4  MAXWPT
        PARAMETER  (MAXWPT=50)
C
        INTEGER*4       NNWORT,NNPART
        REAL*8          PARAMT(MAXWPT),DPARAT(MAXWPT)
        CHARACTER*5     WORDST(MAXWPT),DWORDT(MAXWPT)
C
        COMMON/LISTENTA/NNWORT,NNPART
        COMMON/LISTENTB/PARAMT,DPARAT
        COMMON/LISTENTD/WORDST,DWORDT
C--------------------------------------------------------------- INP_OUT
C**** PLOTERTA:
C**** PLOTERTB:
        INTEGER*4    MMCURT,MSPLOT        ! Max curves & Size of MPLOTT
        PARAMETER    (MMCURT=40,MSPLOT=240) ! MSPLOT=MMCURT*6
C
        INTEGER*4     NCOLDT,NCURVT,NPONTT(MMCURT),MPLOTT(MMCURT,2,3)
        CHARACTER*8   FORMAT
C
        INTEGER*4     IPRCOT
C
        COMMON/PLOTERTA/NCOLDT,NCURVT,NPONTT,MPLOTT
        COMMON/PLOTERTB/FORMAT
        COMMON/PLOTERTC/IPRCOT
C--------------------------------------------------------------- INP_OUT
C**** PRIOUTTA:
        INTEGER*4       KFEMVT,KPRI0T,KPRI1T,KPRI2T,KPRI3T,KPRI4T,
     .                  KPRI5T,KPRI6T,KPRI7T,KPRI8T,KPRI9T,
     .                  ICALPO
C
        COMMON/PRIOUTTA/KFEMVT,KPRI0T,KPRI1T,KPRI2T,KPRI3T,KPRI4T,
     .                  KPRI5T,KPRI6T,KPRI7T,KPRI8T,KPRI9T,
     .                  ICALPO
C=============================================================== INP_OUT
