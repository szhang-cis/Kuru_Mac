      SUBROUTINE CHEK01T(NDIME,IELEM,KNODE,KRULE,KGAUS,ITYPE,LURES)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS THE INTEGRATION RULE SELECTED FOR ELEMENT NO.1
C
C     Notes:
C
C     For KGAUS=1, use KRULE=1-4
C
C     KRULE=11,13 are be used to obtain negative temperature control in
C     presence of large temperature gradients
C
C     restricitions over KRULE=5,11,13 need to be implemented !! Ag/99
C     (also in chek01t.f & chek01s.f)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/ERRORST/NERORT
C
C**** CHECKS INTEGRATION RULE
C
      IF(NDIME.EQ.1) THEN
       IF(KRULE.EQ.1.OR.KRULE.EQ.2) THEN
        if((itype.eq.101.or.itype.eq.104).and.knode.le.2.and.
     .      kgaus.eq.1) return          ! possibility to use subinteg.
       IF((KNODE.EQ.2.OR.KNODE.EQ.3).AND.
     .    (KGAUS.EQ.2.OR.KGAUS.EQ.3))                 RETURN
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                RETURN
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                RETURN
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                RETURN
       ENDIF
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
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                RETURN
       ENDIF
       IF(KRULE.EQ.6) THEN                  ! other option
        IF(KGAUS.EQ.3)                                RETURN
       ENDIF
       IF(KRULE.EQ.7) THEN                  ! other option
        IF(KGAUS.EQ.3)                                RETURN
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                RETURN
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                RETURN
       ENDIF
C
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
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                  RETURN
       ENDIF
       IF(KRULE.EQ.8) THEN                  ! other option
        IF(KGAUS.EQ.2)                                  RETURN
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                  RETURN
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                  RETURN
       ENDIF
C
      ENDIF
C
      WRITE(LURES,900) IELEM,KNODE,KRULE,KGAUS
      NERORT=NERORT+1
      RETURN
C
  900 FORMAT(
     .' INCORRECT INTEGRATION FOR ELEMENT NO. ',I5,/,
     .10X,'NNODE =',I5,2X,'NRULE =',I5,2X,'NGAUS =',I5,//)
      END
