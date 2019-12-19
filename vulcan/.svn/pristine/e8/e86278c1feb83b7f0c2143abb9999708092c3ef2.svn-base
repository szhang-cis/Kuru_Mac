      SUBROUTINE CPOLYT(TEMPG,PROPS,
     .                  CPOL1,CPOL2,CPOL3,EPOL1,EPOL2,EPOL3,
     .                  CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .                  CPOL9,CPOL10,
     .                  CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .                  CEMIN,CEMAX,CEETA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE POLYMER CONSTANTS COEFFICIENT
C     (TEMPERATURE-DEPENDENT)
C     
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*)
      DIMENSION CPOL9(*), CPOL10(*)
C
      IPEP2=INT(PROPS(2))
      IF(IPEP2.NE.40) RETURN
C
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN
#ifndef restricted
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58)
     .  CPOL1=VPOL1(1,1)
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58)
     .  CPOL2=VPOL2(1,1)
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58)
     .  CPOL3=VPOL3(1,1)
       IF(IFREN.EQ.53) THEN
        EPOL1=VPOE1(1,1)
        EPOL2=VPOE2(1,1)
        EPOL3=VPOE3(1,1)
       ENDIF
       IF(IFREN.EQ.56.OR.IFREN.EQ.57) THEN
        CPOL4=VPOL4(1,1)
       ENDIF
       IF(IFREN.EQ.56) THEN
        CPOL5=VPOL5(1,1)
        CPOL6=VPOL6(1,1)
        CPOL7=VPOL7(1,1)
        CPOL8=VPOL8(1,1)
       ENDIF
#endif
       CPENI=VPENI(1,1)
       RETURN
      ENDIF                ! ipep3.eq.1
C
#ifndef restricted
      IPEP4=INT(PROPS(4))
      IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58) THEN
        IF(TEMPG.LE.VPOL1(1,2)) CPOL1=VPOL1(1,1)
        DO IPOL1=2,NPOL1
         I1=IPOL1-1
         I2=IPOL1
         IF(TEMPG.GE.VPOL1(I1,2).AND.TEMPG.LE.VPOL1(I2,2)) 
     .    CPOL1=(VPOL1(I2,1)-VPOL1(I1,1))/(VPOL1(I2,2)-VPOL1(I1,2))*
     .          (TEMPG-VPOL1(I1,2))+VPOL1(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOL1(NPOL1,2)) CPOL1=VPOL1(NPOL1,1)
C
        DPOL1=0.0D0        ! temperature derivative of C1
       ENDIF               ! ifren.eq.51....
C
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58) THEN
        IF(TEMPG.LE.VPOL2(1,2)) CPOL2=VPOL2(1,1)
        DO IPOL2=2,NPOL2
         I1=IPOL2-1
         I2=IPOL2
         IF(TEMPG.GE.VPOL2(I1,2).AND.TEMPG.LE.VPOL2(I2,2))
     .    CPOL2=(VPOL2(I2,1)-VPOL2(I1,1))/(VPOL2(I2,2)-VPOL2(I1,2))*
     .          (TEMPG-VPOL2(I1,2))+VPOL2(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOL2(NPOL2,2)) CPOL2=VPOL2(NPOL2,1)
C
        DPOL2=0.0D0        ! temperature derivative of C2
       ENDIF               ! ifren.eq.51....
C
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58) THEN
        IF(TEMPG.LE.VPOL3(1,2)) CPOL3=VPOL3(1,1)
        DO IPOL3=2,NPOL3
         I1=IPOL3-1
         I2=IPOL3
         IF(TEMPG.GE.VPOL3(I1,2).AND.TEMPG.LE.VPOL3(I2,2))
     .    CPOL3=(VPOL3(I2,1)-VPOL3(I1,1))/(VPOL3(I2,2)-VPOL3(I1,2))*
     .          (TEMPG-VPOL3(I1,2))+VPOL3(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOL3(NPOL3,2)) CPOL3=VPOL3(NPOL3,1)
C
        DPOL3=0.0D0        ! temperature derivative of C3
       ENDIF               ! ifren.eq.51....
C
       IF(IFREN.EQ.53) THEN
        IF(TEMPG.LE.VPOE1(1,2)) EPOL1=VPOE1(1,1)
        DO IPOE1=2,NPOE1
         I1=IPOE1-1
         I2=IPOE1
         IF(TEMPG.GE.VPOE1(I1,2).AND.TEMPG.LE.VPOE1(I2,2)) 
     .    EPOL1=(VPOE1(I2,1)-VPOE1(I1,1))/(VPOE1(I2,2)-VPOE1(I1,2))*
     .          (TEMPG-VPOE1(I1,2))+VPOE1(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOE1(NPOE1,2)) EPOL1=VPOE1(NPOE1,1)
C
        DEPOL1=0.0D0       ! temperature derivative of E1
C
        IF(TEMPG.LE.VPOE2(1,2)) EPOL2=VPOE2(1,1)
        DO IPOE2=2,NPOE2
         I1=IPOE2-1
         I2=IPOE2
         IF(TEMPG.GE.VPOE2(I1,2).AND.TEMPG.LE.VPOE2(I2,2))
     .    EPOL2=(VPOE2(I2,1)-VPOE2(I1,1))/(VPOE2(I2,2)-VPOE2(I1,2))*
     .          (TEMPG-VPOE2(I1,2))+VPOE2(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOE2(NPOE2,2)) EPOL2=VPOE2(NPOE2,1)
C
        DEPOL2=0.0D0       ! temperature derivative of E2
C
        IF(TEMPG.LE.VPOE3(1,2)) EPOL3=VPOE3(1,1)
        DO IPOE3=2,NPOE3
         I1=IPOE3-1
         I2=IPOE3
         IF(TEMPG.GE.VPOE3(I1,2).AND.TEMPG.LE.VPOE3(I2,2))
     .    EPOL3=(VPOE3(I2,1)-VPOE3(I1,1))/(VPOE3(I2,2)-VPOE3(I1,2))*
     .          (TEMPG-VPOE3(I1,2))+VPOE3(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOE3(NPOE3,2)) EPOL3=VPOE3(NPOE3,1)
C
        DEPOL3=0.0D0       ! temperature derivative of E3
       ENDIF               ! ifren.eq.53
C
       IF(IFREN.EQ.56.OR.IFREN.EQ.57) THEN
        IF(TEMPG.LE.VPOL4(1,2)) CPOL4=VPOL4(1,1)
        DO IPOL4=2,NPOL4
         I1=IPOL4-1
         I2=IPOL4
         IF(TEMPG.GE.VPOL4(I1,2).AND.TEMPG.LE.VPOL4(I2,2))
     .    CPOL4=(VPOL4(I2,1)-VPOL4(I1,1))/(VPOL4(I2,2)-VPOL4(I1,2))*
     .          (TEMPG-VPOL4(I1,2))+VPOL4(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOL4(NPOL4,2)) CPOL4=VPOL4(NPOL4,1)
C
        DPOL4=0.0D0        ! temperature derivative of C4
       ENDIF               ! ifren.eq.56.or.ifren.eq.57
C
       IF(IFREN.EQ.56) THEN
        IF(TEMPG.LE.VPOL5(1,2)) CPOL5=VPOL5(1,1)
        DO IPOL5=2,NPOL5
         I1=IPOL5-1
         I2=IPOL5
         IF(TEMPG.GE.VPOL5(I1,2).AND.TEMPG.LE.VPOL5(I2,2))
     .    CPOL5=(VPOL5(I2,1)-VPOL5(I1,1))/(VPOL5(I2,2)-VPOL5(I1,2))*
     .          (TEMPG-VPOL5(I1,2))+VPOL5(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOL5(NPOL5,2)) CPOL5=VPOL5(NPOL5,1)
C
        DPOL5=0.0D0        ! temperature derivative of C5
C
        IF(TEMPG.LE.VPOL6(1,2)) CPOL6=VPOL6(1,1)
        DO IPOL6=2,NPOL6
         I1=IPOL6-1
         I2=IPOL6
         IF(TEMPG.GE.VPOL6(I1,2).AND.TEMPG.LE.VPOL6(I2,2))
     .    CPOL6=(VPOL6(I2,1)-VPOL6(I1,1))/(VPOL6(I2,2)-VPOL6(I1,2))*
     .          (TEMPG-VPOL6(I1,2))+VPOL6(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOL6(NPOL6,2)) CPOL6=VPOL6(NPOL6,1)
C
        DPOL6=0.0D0        ! temperature derivative of C6
C
        IF(TEMPG.LE.VPOL7(1,2)) CPOL7=VPOL7(1,1)
        DO IPOL7=2,NPOL7
         I1=IPOL7-1
         I2=IPOL7
         IF(TEMPG.GE.VPOL7(I1,2).AND.TEMPG.LE.VPOL7(I2,2))
     .    CPOL7=(VPOL7(I2,1)-VPOL7(I1,1))/(VPOL7(I2,2)-VPOL7(I1,2))*
     .          (TEMPG-VPOL7(I1,2))+VPOL7(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOL7(NPOL7,2)) CPOL7=VPOL7(NPOL7,1)
C
        DPOL7=0.0D0        ! temperature derivative of C7
C
        IF(TEMPG.LE.VPOL8(1,2)) CPOL8=VPOL8(1,1)
        DO IPOL8=2,NPOL8
         I1=IPOL8-1
         I2=IPOL8
         IF(TEMPG.GE.VPOL8(I1,2).AND.TEMPG.LE.VPOL8(I2,2))
     .    CPOL8=(VPOL8(I2,1)-VPOL8(I1,1))/(VPOL8(I2,2)-VPOL8(I1,2))*
     .          (TEMPG-VPOL8(I1,2))+VPOL8(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOL8(NPOL8,2)) CPOL8=VPOL8(NPOL8,1)
C
        DPOL8=0.0D0        ! temperature derivative of C8
       ENDIF               ! ifren.eq.56
C
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58) THEN
        IF(TEMPG.LE.VPENI(1,2)) CPENI=VPENI(1,1)
        DO IPENI=2,NPENI
         I1=IPENI-1
         I2=IPENI
         IF(TEMPG.GE.VPENI(I1,2).AND.TEMPG.LE.VPENI(I2,2))
     .    CPENI=(VPENI(I2,1)-VPENI(I1,1))/(VPENI(I2,2)-VPENI(I1,2))*
     .          (TEMPG-VPENI(I1,2))+VPENI(I1,1)
        ENDDO
        IF(TEMPG.GE.VPENI(NPENI,2)) CPENI=VPENI(NPENI,1)
C
        DPENI=0.0D0        ! temperature derivative of PENIN
       ENDIF               ! ifren.eq.51....
C
       IF(NCHAI.GT.0) THEN
        DO IN=1,NCHAI
         IF(TEMPG.LE.VPOL9(IN,1,2)) CPOL9(IN)=VPOL9(IN,1,1)
         DO IPOL9=2,NPOL9
          I1=IPOL9-1
          I2=IPOL9
          IF(TEMPG.GE.VPOL9(IN,I1,2).AND.TEMPG.LE.VPOL9(IN,I2,2))
     .     CPOL9(IN)=(VPOL9(IN,I2,1)-VPOL9(IN,I1,1))/
     .               (VPOL9(IN,I2,2)-VPOL9(IN,I1,2))*
     .               (TEMPG-VPOL9(IN,I1,2))+VPOL9(IN,I1,1)
         ENDDO
         IF(TEMPG.GE.VPOL9(IN,NPOL9,2)) CPOL9(IN)=VPOL9(IN,NPOL9,1)
C
         DPOL9=0.0D0        ! temperature derivative of tau parameter
C
         IF(TEMPG.LE.VPOL10(IN,1,2)) CPOL10(IN)=VPOL10(IN,1,1)
         DO IPOL10=2,NPOL10
          I1=IPOL10-1
          I2=IPOL10
          IF(TEMPG.GE.VPOL10(IN,I1,2).AND.TEMPG.LE.VPOL10(IN,I2,2))
     .     CPOL10(IN)=(VPOL10(IN,I2,1)-VPOL10(IN,I1,1))/
     .                (VPOL10(IN,I2,2)-VPOL10(IN,I1,2))*
     .                (TEMPG-VPOL10(IN,I1,2))+VPOL10(IN,I1,1)
         ENDDO
         IF(TEMPG.GE.VPOL10(IN,NPOL10,2)) CPOL10(IN)=VPOL10(IN,NPOL10,1)
C
         DPOL10=0.0D0      ! temperature derivative of beta parameter
        ENDDO              ! in=1,nchai
       ENDIF               ! nchai gt 0
C
       IF(NFIBN.GT.0) THEN
        IF(NFIBM.EQ.1.OR.NFIBM.EQ.2) THEN
         IF(TEMPG.LE.VFIB1(1,2)) CFIB1=VFIB1(1,1)
         DO IFIB1=2,NFIB1
          I1=IFIB1-1
          I2=IFIB1
          IF(TEMPG.GE.VFIB1(I1,2).AND.TEMPG.LE.VFIB1(I2,2))
     .     CFIB1=(VFIB1(I2,1)-VFIB1(I1,1))/(VFIB1(I2,2)-VFIB1(I1,2))*
     .           (TEMPG-VFIB1(I1,2))+VFIB1(I1,1)
         ENDDO
         IF(TEMPG.GE.VFIB1(NFIB1,2)) CFIB1=VFIB1(NFIB1,1)
C
         DFIB1=0.0D0       ! temperature derivative of C1
C
         IF(TEMPG.LE.VFIB2(1,2)) CFIB2=VFIB2(1,1)
         DO IFIB2=2,NFIB2
          I1=IFIB2-1
          I2=IFIB2
          IF(TEMPG.GE.VFIB2(I1,2).AND.TEMPG.LE.VFIB2(I2,2))
     .     CFIB2=(VFIB2(I2,1)-VFIB2(I1,1))/(VFIB2(I2,2)-VFIB2(I1,2))*
     .           (TEMPG-VFIB2(I1,2))+VFIB2(I1,1)
         ENDDO
         IF(TEMPG.GE.VFIB2(NFIB2,2)) CFIB2=VFIB2(NFIB2,1)
C
         DFIB2=0.0D0       ! temperature derivative of C2
C
         IF(TEMPG.LE.VFIB3(1,2)) CFIB3=VFIB3(1,1)
         DO IFIB3=2,NFIB3
          I1=IFIB3-1
          I2=IFIB3
          IF(TEMPG.GE.VFIB3(I1,2).AND.TEMPG.LE.VFIB3(I2,2))
     .     CFIB3=(VFIB3(I2,1)-VFIB3(I1,1))/(VFIB3(I2,2)-VFIB3(I1,2))*
     .           (TEMPG-VFIB3(I1,2))+VFIB3(I1,1)
         ENDDO
         IF(TEMPG.GE.VFIB3(NFIB3,2)) CFIB3=VFIB3(NFIB3,1)
C
         DFIB3=0.0D0       ! temperature derivative of C3
C
         IF(TEMPG.LE.VFIB4(1,2)) CFIB4=VFIB4(1,1)
         DO IFIB4=2,NFIB4
          I1=IFIB4-1
          I2=IFIB4
          IF(TEMPG.GE.VFIB4(I1,2).AND.TEMPG.LE.VFIB4(I2,2))
     .     CFIB4=(VFIB4(I2,1)-VFIB4(I1,1))/(VFIB4(I2,2)-VFIB4(I1,2))*
     .           (TEMPG-VFIB4(I1,2))+VFIB4(I1,1)
         ENDDO
         IF(TEMPG.GE.VFIB4(NFIB4,2)) CFIB4=VFIB4(NFIB4,1)
C
         DFIB4=0.0D0       ! temperature derivative of C4
        ENDIF              ! nfibm eq 1,2
C
        IF(NFIBM.EQ.2) THEN
         IF(TEMPG.LE.VFIB5(1,2)) CFIB5=VFIB5(1,1)
         DO IFIB5=2,NFIB5
          I1=IFIB5-1
          I2=IFIB5
          IF(TEMPG.GE.VFIB5(I1,2).AND.TEMPG.LE.VFIB5(I2,2))
     .     CFIB5=(VFIB5(I2,1)-VFIB5(I1,1))/(VFIB5(I2,2)-VFIB5(I1,2))*
     .           (TEMPG-VFIB5(I1,2))+VFIB5(I1,1)
         ENDDO
         IF(TEMPG.GE.VFIB5(NFIB5,2)) CFIB5=VFIB5(NFIB5,1)
C
         DFIB5=0.0D0       ! temperature derivative of C5
        ENDIF              ! nfibm eq 2
       ENDIF               ! nfibn gt 0
C
       IF(IDAMG.EQ.1) THEN
        IF(TEMPG.LE.VEMIN(1,2)) CEMIN=VEMIN(1,1)
        DO IEMIN=2,NEMIN
         I1=IEMIN-1
         I2=IEMIN
         IF(TEMPG.GE.VEMIN(I1,2).AND.TEMPG.LE.VEMIN(I2,2))
     .    CEMIN=(VEMIN(I2,1)-VEMIN(I1,1))/(VEMIN(I2,2)-VEMIN(I1,2))*
     .          (TEMPG-VEMIN(I1,2))+VEMIN(I1,1)
        ENDDO
        IF(TEMPG.GE.VEMIN(NEMIN,2)) CEMIN=VEMIN(NEMIN,1)
C
        DEMIN=0.0D0        ! temperature derivative of E_MIN
C
        IF(TEMPG.LE.VEMAX(1,2)) CEMAX=VEMAX(1,1)
        DO IEMAX=2,NEMAX
         I1=IEMAX-1
         I2=IEMAX
         IF(TEMPG.GE.VEMAX(I1,2).AND.TEMPG.LE.VEMAX(I2,2))
     .    CEMAX=(VEMAX(I2,1)-VEMAX(I1,1))/(VEMAX(I2,2)-VEMAX(I1,2))*
     .          (TEMPG-VEMAX(I1,2))+VEMAX(I1,1)
        ENDDO
        IF(TEMPG.GE.VEMAX(NEMAX,2)) CEMAX=VEMAX(NEMAX,1)
C
        DEMAX=0.0D0        ! temperature derivative of E_MAX
C
        IF(TEMPG.LE.VEETA(1,2)) CEETA=VEETA(1,1)
        DO IEETA=2,NEETA
         I1=IEETA-1
         I2=IEETA
         IF(TEMPG.GE.VEETA(I1,2).AND.TEMPG.LE.VEETA(I2,2))
     .    CEETA=(VEETA(I2,1)-VEETA(I1,1))/(VEETA(I2,2)-VEETA(I1,2))*
     .          (TEMPG-VEETA(I1,2))+VEETA(I1,1)
        ENDDO
        IF(TEMPG.GE.VEETA(NEETA,2)) CEETA=VEETA(NEETA,1)
C
        DEETA=0.0D0        ! temperature derivative of ETA
       ENDIF               ! idamg eq 1
      ENDIF                ! ipep4.eq.0.or.ipep4.eq.2
C
#endif
      RETURN
      END
