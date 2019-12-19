      SUBROUTINE CHEK01(NDIME,IELEM,KNODE,KRULE,KGAUS,NQUTR)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS THE INTEGRATION RULE SELECTED FOR ELEMENT NO.1
C
C     Notes:
C
C     For KGAUS=1, use KRULE=1-4
C
C     KRULE=11,13 are be used to obtain lumped mass matrix
C
C     restricitions over KRULE=5,11,13 need to be implemented !! Ag/99
C     (also in chek01t.f & chek01s.f)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/ERRORS/NEROR,NEROB,NEROT,NEROI
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      COMMON/LOGUN/LUDTS,LUSOL,LUFRO,LUFRH,LUDAT,LUPRI,LURES,
     .             LUSO2,LUFR2,LUPOS,LURST,LUBFG,LUPIP,LUPAN,LUINF,
     .             LUGEO,LUSET,LUMAT,LUINI,LULOA,LUFIX,LUIN1,
     .             LUTUN,LUCON,LUACT,LUFAN,
     .             LUCU1,LUCU2,LUCU3,LUCU4,LUCU5,LUCU6,LUCU7,
     .             LUCU8,LUCU9,LUC10
C
C**** CHECKS INTEGRATION RULE
C
      IRETUR=0
C
      IF(NDIME.EQ.1) THEN
       IF(KRULE.EQ.1.OR.KRULE.EQ.2) THEN
        IF((KNODE.EQ.2.OR.KNODE.EQ.3).AND.
     .     (KGAUS.EQ.2.OR.KGAUS.EQ.3))                IRETUR=1
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUR=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUR=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUR=1
       ENDIF
C
      ELSE IF(NDIME.EQ.2) THEN
       IF(KRULE.GT.0.AND.KRULE.LE.3) THEN
        IF((KNODE.EQ.4.OR.KNODE.EQ.8.OR.KNODE.EQ.9).AND.
     .     (KRULE.LE.2.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.4.OR.KGAUS.EQ.9))) IRETUR=1
        IF((KNODE.EQ.4).AND.(NQUTR.EQ.2).AND.
     .     (KRULE.EQ.3.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.3.OR.KGAUS.EQ.4.OR.
     .                    KGAUS.EQ.6.OR.KGAUS.EQ.7.OR.
     .                    KGAUS.EQ.13)))              IRETUR=1
        IF((KNODE.EQ.5).AND.
     .     (KRULE.EQ.1.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.4.OR.KGAUS.EQ.9))) IRETUR=1
        IF((KNODE.EQ.3.OR.KNODE.EQ.6).AND.
     .     (KRULE.EQ.3.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.3.OR.KGAUS.EQ.4.OR.
     .                    KGAUS.EQ.6.OR.KGAUS.EQ.7.OR.
     .                    KGAUS.EQ.13)))              IRETUR=1
        IF((KNODE.EQ.3).AND.(NQUTR.EQ.2).AND.
     .     (KRULE.LE.2.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.4.OR.KGAUS.EQ.9))) IRETUR=1
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUR=1
       ENDIF
       IF(KRULE.EQ.6) THEN                  ! other option
        IF(KGAUS.EQ.3)                                IRETUR=1
       ENDIF
       IF(KRULE.EQ.7) THEN                  ! other option
        IF(KGAUS.EQ.3)                                IRETUR=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUR=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUR=1
       ENDIF
C
      ELSE                                  ! ndime=3
       IF(KRULE.GT.0.AND.KRULE.LE.4) THEN
        IF((KNODE.EQ.8.OR.KNODE.EQ.20.OR.KNODE.EQ.27).AND.
     .     (KRULE.LE.2.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.8.OR.KGAUS.EQ.27)))  IRETUR=1
        IF((KNODE.EQ.8.OR.KNODE.EQ.20.OR.KNODE.EQ.27).AND.
     .     (KRULE.EQ.4.AND.
     .     (KGAUS.NE.6.OR.KGAUS.NE.14.OR.KGAUS.NE.15))) IRETUR=1
        IF((KNODE.EQ.4.OR.KNODE.EQ.10).AND.
     .     (KRULE.EQ.3.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.4.OR.KGAUS.EQ.5)))   IRETUR=1
        IF((KNODE.EQ.4).AND.(NQUTR.EQ.2).AND.
     .     (KRULE.LE.2.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.8.OR.KGAUS.EQ.27)))  IRETUR=1
        IF((KNODE.EQ.5).AND.
     .     (KRULE.EQ.3.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.4.OR.KGAUS.EQ.5)))   IRETUR=1
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUR=1
       ENDIF
       IF(KRULE.EQ.8) THEN                  ! other option
        IF(KGAUS.EQ.2)                                IRETUR=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUR=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUR=1
       ENDIF
C
      ENDIF
C
      IF(IRETUR.EQ.0) THEN
       WRITE(LURES,900) IELEM,KNODE,KRULE,KGAUS
       NEROR=NEROR+1
      ENDIF
C
C**** CHECKS COMPATIBILITY OF KPARB AND INTEGRATION RULE
C
      KPARB=INT(PPARB)
      IF(KPARB.LT.0) KPARB=-KPARB
C
      IF(KPARB.NE.0) THEN
C
      IRETUB=0
C
      IF(NDIME.EQ.1) THEN
       IF(KRULE.EQ.1.OR.KRULE.EQ.2) THEN
        IF(KNODE.EQ.2)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUB=1
       ENDIF
C
      ELSE IF(NDIME.EQ.2) THEN
       IF(KRULE.GT.0.AND.KRULE.LE.3) THEN
        IF(KNODE.EQ.4)                                IRETUB=1
        IF(KNODE.EQ.5)                                IRETUB=1
        IF((KNODE.EQ.8.OR.KNODE.EQ.9).AND.
     .     (KRULE.LE.2.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.4))) THEN
          IF(KPARB.GE.1.AND.KPARB.LE.3)               IRETUB=1
        ENDIF
        IF(KNODE.EQ.3)                                IRETUB=1
        IF(KNODE.EQ.6.AND.
     .    (KRULE.EQ.3.AND.
     .    (KGAUS.EQ.1.OR.KGAUS.EQ.3.OR.KGAUS.EQ.4.
     .                OR.KGAUS.EQ.6.OR.KGAUS.EQ.7.
     .                OR.KGAUS.EQ.13)))               IRETUB=1
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.6) THEN                  ! other option
        IF(KGAUS.EQ.3)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.7) THEN                  ! other option
        IF(KGAUS.EQ.3)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUB=1
       ENDIF
C
      ELSE
       IF(KRULE.GT.0.AND.KRULE.LE.4) THEN
        IF(KNODE.EQ.8)                                IRETUB=1
        IF((KNODE.EQ.20.OR.KNODE.EQ.27).AND.
     .     (KRULE.LE.2.AND.
     .     (KGAUS.EQ.1.OR.KGAUS.EQ.8))) THEN
         IF(KPARB.GE.1.AND.KPARB.LE.3)                IRETUB=1
        ENDIF
        IF((KNODE.EQ.8.OR.KNODE.EQ.20.OR.KNODE.EQ.27).AND.
     .     (KRULE.EQ.4.AND.
     .     (KGAUS.NE.6.OR.KGAUS.NE.14.OR.KGAUS.NE.15))) IRETUB=1
        IF(KNODE.EQ.4)                                  IRETUB=1
        IF(KNODE.EQ.5)                                  IRETUB=1
        IF(KNODE.EQ.10.AND.
     .    (KRULE.EQ.3.AND.
     .    (KGAUS.EQ.1.OR.KGAUS.EQ.4.OR.KGAUS.EQ.5)))    IRETUB=1
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.8) THEN                  ! other option
        IF(KGAUS.EQ.2)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUB=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUB=1
       ENDIF
C
      ENDIF
C
      IF(IRETUB.EQ.0) THEN
       WRITE(LURES,901) IELEM,KNODE,KRULE,KGAUS,KPARB
       NEROB=NEROB+1
      ENDIF
C
      ENDIF                                 ! kparb.ne.0
C
C**** CHECKS THE CONVENIENT THE USE OF INDEX KPART
C
      KPART=INT(PPART)
      KPARI=INT(PPARI)
C
      IRETUT=0
C
      IF(NDIME.EQ.1) THEN
       IF(KRULE.EQ.1.OR.KRULE.EQ.2) THEN
        IF(KNODE.EQ.2) THEN
         IF(KPART.EQ.0)                               IRETUT=1
         IF(KPART.EQ.1.AND.KPARI.EQ.1)                IRETUT=1
        ENDIF
        IF(KNODE.EQ.3) THEN                               
         IF(KPART.EQ.1)                               IRETUT=1
        ENDIF
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUT=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUT=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUT=1
       ENDIF
C
      ELSE IF(NDIME.EQ.2) THEN
       IF(KRULE.GT.0.AND.KRULE.LE.3) THEN
        IF(KNODE.EQ.4) THEN
         IF(KPART.EQ.0)                               IRETUT=1
         IF(KPART.EQ.1.AND.KPARI.EQ.1)                IRETUT=1
        ENDIF
        IF(KNODE.EQ.8.OR.KNODE.EQ.9) THEN
         IF(KPART.EQ.1)                               IRETUT=1
        ENDIF
        IF(KNODE.EQ.3) THEN
         IF(KPART.EQ.0)                               IRETUT=1
         IF(KPART.EQ.1)                               IRETUT=1
        ENDIF
        IF(KNODE.EQ.6.AND.
     .    (KRULE.EQ.3.AND.
     .    (KGAUS.EQ.1.OR.KGAUS.EQ.3.OR.KGAUS.EQ.4.
     .                OR.KGAUS.EQ.6.OR.KGAUS.EQ.7.
     .                OR.KGAUS.EQ.13)))               IRETUT=1
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUT=1
       ENDIF
       IF(KRULE.EQ.6) THEN                  ! other option
        IF(KGAUS.EQ.3)                                IRETUT=1
       ENDIF       
       IF(KRULE.EQ.7) THEN                  ! other option
        IF(KGAUS.EQ.3)                                IRETUT=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUT=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUT=1
       ENDIF
C
      ELSE
       IF(KRULE.GT.0.AND.KRULE.LE.4) THEN
        IF(KNODE.EQ.8) THEN
         IF(KPART.EQ.0)                               IRETUT=1
         IF(KPART.EQ.1.AND.KPARI.EQ.1)                IRETUT=1
        ENDIF
        IF(KNODE.EQ.20.OR.KNODE.EQ.27) THEN
         IF(KPART.EQ.1)                               IRETUT=1
        ENDIF
        IF((KNODE.EQ.8.OR.KNODE.EQ.20.OR.KNODE.EQ.27).AND.
     .     (KRULE.EQ.4.AND.
     .     (KGAUS.NE.6.OR.KGAUS.NE.14.OR.KGAUS.NE.15))) IRETUT=1
        IF(KNODE.EQ.4) THEN
         IF(KPART.EQ.0)                                 IRETUT=1
         IF(KPART.EQ.1)                                 IRETUT=1
        ENDIF
        IF(KNODE.EQ.5) THEN
         IF(KPART.EQ.0)                                 IRETUT=1
         IF(KPART.EQ.1)                                 IRETUT=1
        ENDIF
        IF(KNODE.EQ.10) THEN
         IF(KPART.EQ.1)                                 IRETUT=1
        ENDIF
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUT=1
       ENDIF
       IF(KRULE.EQ.8) THEN                  ! other option
        IF(KGAUS.EQ.2)                                IRETUT=1
       ENDIF       
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUT=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUT=1
       ENDIF
C
      ENDIF
C
      IF(IRETUT.EQ.0) THEN
       WRITE(LURES,902) IELEM,KNODE,KRULE,KGAUS,KPART
       NEROT=NEROT+1
      ENDIF
C
C**** CHECKS THE CONVENIENT THE USE OF INDEX KPARI
C
      IRETUI=0
C
      IF(NDIME.EQ.1) THEN
       IF(KRULE.EQ.1.OR.KRULE.EQ.2) THEN
        IF(KNODE.EQ.2) THEN
         IF(KPARI.EQ.1)                               IRETUI=1
        ENDIF
        IF(KNODE.EQ.3) THEN
         IF(KPARI.EQ.0)                               IRETUI=1
        ENDIF
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUI=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUI=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUI=1
       ENDIF
C
      ELSE IF(NDIME.EQ.2) THEN
       IF(KRULE.GT.0.AND.KRULE.LE.3) THEN
        IF(KNODE.EQ.4) THEN
         IF(KPARI.EQ.1)                               IRETUI=1
        ENDIF
        IF(KNODE.EQ.8.OR.KNODE.EQ.9) THEN
         IF(KPARI.EQ.0)                               IRETUI=1
        ENDIF
        IF(KNODE.EQ.3) THEN
         IF(KPARI.EQ.0)                               IRETUI=1
         IF(KPARI.EQ.1)                               IRETUI=1
        ENDIF
        IF(KNODE.EQ.6.AND.
     .    (KRULE.EQ.3.AND.
     .    (KGAUS.EQ.1.OR.KGAUS.EQ.3.OR.KGAUS.EQ.4.
     .                OR.KGAUS.EQ.6.OR.KGAUS.EQ.7.
     .                OR.KGAUS.EQ.13)))               IRETUI=1
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUI=1
       ENDIF
       IF(KRULE.EQ.6) THEN                  ! other option
        IF(KGAUS.EQ.3)                                IRETUI=1
       ENDIF
       IF(KRULE.EQ.7) THEN                  ! other option
        IF(KGAUS.EQ.3)                                IRETUI=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUI=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUI=1
       ENDIF
C
      ELSE
       IF(KRULE.GT.0.AND.KRULE.LE.4) THEN
        IF(KNODE.EQ.8) THEN
         IF(KPARI.EQ.1)                               IRETUI=1
        ENDIF
        IF(KNODE.EQ.20.OR.KNODE.EQ.27) THEN
         IF(KPARI.EQ.0)                               IRETUI=1
        ENDIF
        IF((KNODE.EQ.8.OR.KNODE.EQ.20.OR.KNODE.EQ.27).AND.
     .     (KRULE.EQ.4.AND.
     .     (KGAUS.NE.6.OR.KGAUS.NE.14.OR.KGAUS.NE.15))) IRETUI=1
        IF(KNODE.EQ.4) THEN
         IF(KPARI.EQ.0)                                 IRETUI=1
         IF(KPARI.EQ.1)                                 IRETUI=1
        ENDIF
        IF(KNODE.EQ.5) THEN
         IF(KPARI.EQ.0)                                 IRETUI=1
         IF(KPARI.EQ.1)                                 IRETUI=1
        ENDIF
        IF(KNODE.EQ.10) THEN
         IF(KPARI.EQ.0)                                 IRETUI=1
        ENDIF
       ENDIF
       IF(KRULE.EQ.5) THEN                  ! user's definition
        IF(KGAUS.NE.1)                                IRETUI=1
       ENDIF
       IF(KRULE.EQ.8) THEN                  ! other option
        IF(KGAUS.EQ.2)                                IRETUI=1
       ENDIF
       IF(KRULE.EQ.11) THEN                 ! special integration
        IF(KGAUS.NE.1)                                IRETUI=1
       ENDIF
       IF(KRULE.EQ.13) THEN
        IF(KGAUS.NE.1)                                IRETUI=1
       ENDIF
C
      ENDIF
C
      IF(IRETUI.EQ.0) THEN
       WRITE(LURES,903) IELEM,KNODE,KRULE,KGAUS,KPARI
       NEROI=NEROI+1
      ENDIF
C
      RETURN
C
  900 FORMAT(
     .' INCORRECT INTEGRATION FOR ELEMENT NO. ',I7,/,
     .10X,'NNODE =',I5,2X,'NRULE =',I5,2X,'NGAUS =',I5,//)
  901 FORMAT(
     .' INCORRECT INDEX KPARB FOR ELEMENT NO. ',I7,/,
     .10X,'NNODE =',I5,2X,'NRULE =',I5,2X,'NGAUS =',I5,'KPARB =',I5,//)
  902 FORMAT(
     .' INCONVENIENT INDEX KPART FOR ELEMENT NO. ',I7,/,
     .10X,'NNODE =',I5,2X,'NRULE =',I5,2X,'NGAUS =',I5,'KPART =',I5,//)
  903 FORMAT(
     .' INCONVENIENT INDEX KPARI FOR ELEMENT NO. ',I7,/,
     .10X,'NNODE =',I5,2X,'NRULE =',I5,2X,'NGAUS =',I5,'KPARI =',I5,//)
      END
