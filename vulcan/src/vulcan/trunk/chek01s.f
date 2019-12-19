      SUBROUTINE CHEK01S(NDIME,IELEM,KNODE,KRULE,KGAUS,ITYPE)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS THE INTEGRATION RULE SELECTED FOR ELEMENT NO.1
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/ERRORSS/NERORS
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
      IF(NDIME.EQ.1) THEN
       IF(KRULE.EQ.1.OR.KRULE.EQ.2) THEN
        if((itype.eq.101.or.itype.eq.104).and.knode.le.2.and.
     .      kgaus.eq.1) return          ! possibility to use subinteg.
       IF((KNODE.EQ.2.OR.KNODE.EQ.3).AND.
     .    (KGAUS.EQ.2.OR.KGAUS.EQ.3))                 RETURN
       ENDIF
C
C**** USER'S DEFINITION
C
       IF(KRULE.EQ.5) THEN
        IF(KGAUS.NE.1) RETURN               ! for KGAUS=1, use KRULE=1-4
       ENDIF
C
C**** SPECIAL INTEGRATION RULE FOR NEGATIVE TEMPERATURE CONTROL
C
       IF(KRULE.EQ.11) RETURN
       IF(KRULE.EQ.13) RETURN
C
      ELSE IF(NDIME.EQ.2) THEN
       IF(KRULE.GT.0.AND.KRULE.LE.3) THEN
        IF((KNODE.EQ.4.OR.KNODE.EQ.8.OR.KNODE.EQ.9).AND.
     .     (KRULE.LE.2.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.4.OR.KGAUS.EQ.9))) RETURN
        IF((KNODE.EQ.3.OR.KNODE.EQ.6).AND.
     .     (KRULE.EQ.3.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.3.OR.KGAUS.EQ.4.
     .                 OR.KGAUS.EQ.6.OR.KGAUS.EQ.7))) RETURN
       ENDIF
C
C**** USER'S DEFINITION
C
       IF(KRULE.EQ.5) THEN
        IF(KGAUS.NE.1) RETURN               ! for KGAUS=1, use KRULE=1-4
       ENDIF
C
       IF(KRULE.EQ.6) THEN                  ! other option
        IF(KGAUS.EQ.3)                                RETURN
       ENDIF
       IF(KRULE.EQ.7) THEN                  ! other option
        IF(KGAUS.EQ.3)                                RETURN
       ENDIF
C
C**** SPECIAL INTEGRATION RULE FOR NEGATIVE TEMPERATURE CONTROL
C
       IF(KRULE.EQ.11) RETURN
       IF(KRULE.EQ.13) RETURN
      ELSE                                       ! ndime=3
       IF(KRULE.GT.0.AND.KRULE.LE.4) THEN
        IF((KNODE.EQ.8.OR.KNODE.EQ.20.OR.KNODE.EQ.27).AND.
     .     (KRULE.LE.2.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.8.OR.KGAUS.EQ.27)))  RETURN
        IF((KNODE.EQ.8.OR.KNODE.EQ.20.OR.KNODE.EQ.27).AND.
     .     (KRULE.EQ.4.AND.
     .     (KGAUS.NE.6.OR.KGAUS.NE.14.OR.KGAUS.NE.15))) RETURN
        IF((KNODE.EQ.4.OR.KNODE.EQ.10).AND.
     .     (KRULE.EQ.3.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.4.OR.KGAUS.EQ.5)))   RETURN
       ENDIF
C
C**** USER'S DEFINITION
C
       IF(KRULE.EQ.5) THEN
        IF(KGAUS.NE.1) RETURN               ! for KGAUS=1, use KRULE=1-4
       ENDIF
C
       IF(KRULE.EQ.8) THEN                  ! other option
        IF(KGAUS.EQ.2)                                  RETURN
       ENDIF
C
C**** SPECIAL INTEGRATION RULE FOR NEGATIVE TEMPERATURE CONTROL 
C
       IF(KRULE.EQ.11) RETURN
       IF(KRULE.EQ.13) RETURN
C
      ENDIF
C
      WRITE(LURESS,900) IELEM,KNODE,KRULE,KGAUS
      NERORS=NERORS+1
      RETURN
C
  900 FORMAT(
     .' INCORRECT INTEGRATION FOR ELEMENT NO. ',I5,/,
     .10X,'NNODE =',I5,2X,'NRULE =',I5,2X,'NGAUS =',I5,//)
      END
