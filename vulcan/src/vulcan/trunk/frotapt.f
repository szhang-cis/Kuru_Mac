      SUBROUTINE FROTAPT(KFRON,LOCEL,NACVA,NDEST,NEVAB,NFRON)
C***********************************************************************
C
C**** THIS ROUTINE READS FROM TAPE THE DESTINATION OF THE CURRENT
C     ELEMENTAL EQUATION INTO FRONT
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/LOGUNT/LUDTST,LUSOLT,LUFROT,LUFRHT,LUDATT,LUPRIT,LUREST,
     .              LUSO2T,LUFR2T,LUPOST,LURSTT,LUBFGT,LUPIPT,LUPANT,
     .              LUGEOT,LUSETT,LUMATT,LUINIT,LULOAT,LUFIXT,LUADVT,
     .              LUACTT,LUFANT,LUSTRT,
     .              LUCU1T,LUCU2T,LUCU3T,LUCU4T,LUCU5T,LUCU6T,LUCU7T,
     .              LUCU8T,LUCU9T,LUC10T
C
      DIMENSION LOCEL(NEVAB), NACVA(NFRON), NDEST(NEVAB)
C
      READ(LUFRHT) LOCEL,NDEST,KFRON
      READ(LUFRHT) (NACVA(IACVA), IACVA=1,KFRON)
C
      DO IACVA=KFRON+1,NFRON
       NACVA(IACVA)=0
      ENDDO
C
      RETURN
      END
