      SUBROUTINE FROTAPS(KFRON,LOCEL,NACVA,NDEST,NEVAB,NFRON)
C***********************************************************************
C
C**** THIS ROUTINE READS FROM TAPE THE DESTINATION OF THE CURRENT
C     ELEMENTAL EQUATION INTO FRONT
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
c     COMMON/LOGUNT/LUDTST,LUSOLT,LUFROT,LUFRHT,LUDATT,LUPRIT,LUREST,
c    .              LUSO2T,LUFR2T,LUPOST,LURSTT,LUBFGT,LUPIPT,LUPANT,
c    .              LUGEOT,LUSETT,LUMATT,LUINIT,LULOAT,LUFIXT,LUADVT,
c    .              LUACTS,LUFANS,
c    .              LUCU1T,LUCU2T,LUCU3T,LUCU4T,LUCU5T,LUCU6T,LUCU7T,
c    .              LUCU8T,LUCU9T,LUC10T
      COMMON/LOGUNS/LUSOLS,LUFRHS,LUDATS,LUPRIS,LURESS,LUPOSS,
     .              LUCU1S
C
      DIMENSION LOCEL(NEVAB), NACVA(NFRON), NDEST(NEVAB)
C
      READ(LUFRHS) LOCEL,NDEST,KFRON
      READ(LUFRHS) (NACVA(IACVA), IACVA=1,KFRON)
C
      DO IACVA=KFRON+1,NFRON
       NACVA(IACVA)=0
      ENDDO
C
      RETURN
      END