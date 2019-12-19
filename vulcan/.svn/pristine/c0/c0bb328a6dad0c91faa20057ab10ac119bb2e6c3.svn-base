      SUBROUTINE OUTG30(EHIST,STRAN,STRSG,DVOLU,KPARI)
C***********************************************************************
C
C**** THIS ROUTINE PRINT THE GAUSSIAN VARIABLES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION EHIST(NHIST,*), STRAN(NSTR1,*), 
     .          STRSG(NSTR1,*)
      DIMENSION STRSP(3)
      DIMENSION DVOLU(*),       STSMO(6)
C
      KELGS=0
C
      WRITE(LURES,900) IELEM
C
      IF(NTYPE.EQ.1.OR.NTYPE.EQ.2) THEN
       IF(KPRI2.EQ.1) WRITE(LURES,960)
       IF(KPRI3.EQ.1) WRITE(LURES,962)
      ELSE IF(NTYPE.EQ.3) THEN
       IF(KPRI2.EQ.1) WRITE(LURES,965)
       IF(KPRI3.EQ.1) WRITE(LURES,966)
      ELSE IF(NTYPE.EQ.4) THEN
       IF(KPRI2.EQ.1) WRITE(LURES,970)
       IF(KPRI3.EQ.1) WRITE(LURES,971)
      ELSE IF(NTYPE.EQ.5) THEN
       IF(KPRI2.EQ.1) WRITE(LURES,972)
      ENDIF
C
C**** CALCULATE SMOOTHED STRESSES FOR THE SOLIDIFICATION PROBLEM
C
      GO TO(1,2,3) NDIME
C
    1 CONTINUE
C
C**** 1D CASE (not implemented yet)
C
      CALL RUNEND('ERROR IN OUTG30; NDIME = 1         ')
      GO TO 4
C
    2 CONTINUE
C
C**** 2D CASE FOR 3-NODED, 4-NODED & 8-NODED ELEMENTS
C
      ICHEQ=0
      IF(NNODL.EQ.3.AND.NGAUL.EQ.1) THEN
       ICHEQ=1
       VOLIN=1.0/DVOLU(1)
      ENDIF
      IF(NNODL.EQ.4.AND.NGAUL.EQ.4) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4))
       if(kpari.eq.2) volin=0.25D0
      ENDIF
      IF(NNODL.EQ.8.AND.NGAUL.EQ.4) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4))
      ENDIF
      IF(NNODL.EQ.8.AND.NGAUL.EQ.9) THEN                      ! pichuleo
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4)+
     .            DVOLU(5)+DVOLU(6)+DVOLU(7)+DVOLU(8)+DVOLU(9))
       ICHEQ=1
      ENDIF
      IF(ICHEQ.EQ.0)
     .  CALL RUNEND('ERROR IN OUTG30; NDIME = 2         ')
      GO TO 4
C
    3 CONTINUE
C
C**** 3D CASE FOR 4-NODED & 8-NODED ELEMENTS
C
      ICHEQ=0
      IF(NNODL.EQ.4.AND.NGAUL.EQ.1) THEN
       ICHEQ=1
       VOLIN=1.0/DVOLU(1)
      ENDIF
      IF(NNODL.EQ.8.AND.NGAUL.EQ.8) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4)+
     .            DVOLU(5)+DVOLU(6)+DVOLU(7)+DVOLU(8))
      ENDIF
      IF(ICHEQ.EQ.0)
     . CALL RUNEND('ERROR IN OUTG30; NDIME = 3         ')
C
    4 CONTINUE
C
      DO ISTR1=1,NSTR1
       STSMO(ISTR1)=0.0D0
       DO IGAUS=1,NGAUL
        if(kpari.eq.1)
     .  STSMO(ISTR1)=STSMO(ISTR1)+STRSG(ISTR1,IGAUS)*DVOLU(IGAUS)
        if(kpari.eq.2)
     .  STSMO(ISTR1)=STSMO(ISTR1)+STRSG(ISTR1,IGAUS)
       ENDDO
       STSMO(ISTR1)=STSMO(ISTR1)*VOLIN
      ENDDO
C
C**** LOOP OVER GAUSS POINTS
C
      DO 30 IGAUS=1,NGAUL
      KELGS=KELGS+1
C
C**** CALCULATE THE PRINCIPAL STRESSES
C
      IF(KPRI3.EQ.1) THEN
       IF(NTYPE.LE.3) CALL PRIDIR(NTYPE,STRSG(1,IGAUS),STRSP)
       IF(NTYPE.EQ.4) CALL PRINSI(STRSG(1,IGAUS),STRSP)
      ENDIF
C
      IF(NTYPE.NE.5) THEN
       EPSIV=STRAN(1,IGAUS)+STRAN(2,IGAUS)+STRAN(4,IGAUS)
      ELSE
       EPSIV=STRAN(1,IGAUS)
      ENDIF
C
C**** PRINT STRESSES
C
      IF(NTYPE.NE.4) THEN
       IF(KPRI2.EQ.1) WRITE(LURES,905) KELGS,
     .                                (STSMO(ISTR1),ISTR1=1,NSTR1),
     .                                 EPSIV    ! STSMO instead of STRSG
       IF(KPRI3.EQ.1.AND.NTYPE.NE.5)
     .  WRITE(LURES,905) KELGS,(STRSP(ISTRE),ISTRE=1,3)
      ENDIF
C
      IF(NTYPE.EQ.4) THEN
       IF(KPRI2.EQ.1)
     .  WRITE(LURES,905) KELGS,(STRSG(ISTR1,IGAUS),ISTR1=1,NSTR1),
     .                   EPSIV     ! it should be STSMO instead of STRSG
       IF(KPRI3.EQ.1)
     .  WRITE(LURES,905) KELGS,(STRSP(ISTR1),ISTR1=1,3)
      ENDIF
   30 CONTINUE
C
      IF(NTYPE.EQ.1.OR.NTYPE.EQ.2) THEN
       IF(KPRI4.EQ.1) WRITE(LURES,994)
      ELSE IF(NTYPE.EQ.3) THEN
       IF(KPRI4.EQ.1) WRITE(LURES,994)
      ELSE IF(NTYPE.EQ.4) THEN
       IF(KPRI4.EQ.1) WRITE(LURES,994)
      ENDIF
C
C**** PRINT STATE VARIABLES
C
C     Note: although IPLAO should be used, it is not used in order to
C           print separately the different variables (see pointe.f)
C
      IF(KPRI4.EQ.1) THEN
       DO 40 IGAUS=1,NGAUL
       WRITE(LURES,981) IGAUS,           ! plastic strain tensor (STRAP)
     .                        (EHIST(INDEX,IGAUS),INDEX=1,NSTR1)
       INDEX=IPLAS(3)                    ! effective plastic deformation
       WRITE(LURES,981) IGAUS,EHIST(INDEX,IGAUS)
C
       INDEY=INDEX
       IF(KPLA1.EQ.1) THEN               ! isotropic hardening (EBASE)
        WRITE(LURES,981) IGAUS,EHIST(INDEY+1,IGAUS)
        INDEY=INDEY+1
       ENDIF
       IF(KPLA2.EQ.1) THEN               ! kinematic hardening (EBASE)
        WRITE(LURES,981) IGAUS,(EHIST(INDEY+INDEX,IGAUS),INDEX=1,NSTR1)
        INDEY=INDEY+NSTR1
       ENDIF
       IF(KPLA4.EQ.1) THEN               ! porosity+tot. hard. (EBASE)

       ENDIF
       IF(KPLA7.EQ.1) THEN               ! dot e^p (EBASE)

       ENDIF
       IF(KPLA8.EQ.1) THEN               ! recrystall. (S & L) (EBASE)
        WRITE(LURES,981) IGAUS,(EHIST(INDEY+INDEX,IGAUS),INDEX=1,2)
        INDEY=INDEY+2
       ENDIF
       IF(KPLA9.EQ.1) THEN               ! old equiv. stress (EBASE)

       ENDIF
       IF(KPLA10.EQ.1) THEN              ! viscous stress (EBASE)

       ENDIF
       IF(KPLA11.EQ.1) THEN              ! not used now!

       ENDIF
C
       IF(KPLA3.EQ.1) THEN               ! damage (DBASE)
        WRITE(LURES,981) IGAUS,(EHIST(INDEX,IGAUS),
     .                                        INDEX=IPLAS(4),IPLAS(5)-1)
       ENDIF
       IF(KPLA5.EQ.1) THEN               ! shrinkage (SHRIN)

       ENDIF
       IF(KPLA6.EQ.1) THEN               ! fatigue (FATIG)

       ENDIF
       IF(KPLA12.EQ.1) THEN              ! dual phase steel (VDUAL)
        WRITE(LURES,981) IGAUS,(EHIST(INDEX,IGAUS),
     .                       INDEX=IPLAS(8)+2*NSTR1+2*NKOST,
     .                             IPLAS(8)+2*NSTR1+2*NKOST+NSTR1-1)
        WRITE(LURES,981) IGAUS,(EHIST(INDEX,IGAUS),
     .                       INDEX=IPLAS(8)+3*NSTR1+2*NKOST,
     .                             IPLAS(8)+3*NSTR1+2*NKOST+NSTR1-1)
        WRITE(LURES,981) IGAUS,(EHIST(INDEX,IGAUS),
     .                       INDEX=IPLAS(8)+4*NSTR1+2*NKOST,
     .                             IPLAS(8)+4*NSTR1+2*NKOST+1-1)
        WRITE(LURES,981) IGAUS,(EHIST(INDEX,IGAUS),
     .                       INDEX=IPLAS(8)+4*NSTR1+2*NKOST+1,
     .                             IPLAS(8)+4*NSTR1+2*NKOST+2-1)
        WRITE(LURES,981) IGAUS,(EHIST(INDEX,IGAUS),
     .                       INDEX=IPLAS(8)+4*NSTR1+2*NKOST+2,
     .                             IPLAS(8)+4*NSTR1+2*NKOST+3-1)
        WRITE(LURES,981) IGAUS,(EHIST(INDEX,IGAUS),
     .                       INDEX=IPLAS(8)+4*NSTR1+2*NKOST+3,
     .                             IPLAS(8)+4*NSTR1+2*NKOST+4-1)
       ENDIF
   40  CONTINUE
      ENDIF
C
      IF(NCRIT.EQ.5) THEN
       DO 50 IGAUS=1,NGAUL
       IF(NTYPE.NE.4)
     .  WRITE(LURES,991) IGAUS, (EHIST(INDE1,IGAUS),INDE1=16,17),
     .                      (EHIST(INDE2,IGAUS),INDE2=13,14),
     .                      (EHIST(INDE3,IGAUS),INDE3= 1, 4)
       IF(NTYPE.EQ.4)
     .  WRITE(LURES,992) IGAUS, (EHIST(INDE1,IGAUS),INDE1=16,18),
     .                      (EHIST(INDE2,IGAUS),INDE2=13,15),
     .                      (EHIST(INDE3,IGAUS),INDE3= 1, 9)
   50  CONTINUE
      ENDIF
C
      IF(NCRIT.EQ.21.OR.NCRIT.EQ.22.OR.NCRIT.EQ.23) THEN
       DO 60 IGAUS=1,NGAUL
       WRITE(LURES,993) IGAUS, (EHIST(INDE1,IGAUS),INDE1=10,12),
     .                      EHIST(13,IGAUS),EHIST(17,IGAUS),
     .                      EHIST(19,IGAUS)
   60  CONTINUE
      ENDIF
C
      RETURN
C
  900 FORMAT(5X,13HELEMENT NO. =,I5)
  905 FORMAT(5X,I5,8E15.6)
  960 FORMAT(5X,5H G.P.,4X,10H XX-STRESS,5X,10H YY-STRESS,
     .       5X,10H XY-STRESS,5X,10H ZZ-STRESS,3X,12H DVOL.STRAIN)
  962 FORMAT(5X,5H G.P.,5X,10H  MAX P.S.,5X,10H  MIN P.S.,
     .          10X,5HANGLE)
  965 FORMAT(5X,5H G.P.,5X,10H RR-STRESS,5X,10H ZZ-STRESS,
     .       5X,10H RZ-STRESS,5X,10H TT-STRESS,3X,12H DVOL.STRAIN)
  966 FORMAT(5X,5H G.P.,5X,10H  MAX P.S.,5X,10H  MIN P.S.,
     .       10X,5HANGLE)
  970 FORMAT(5X,5H G.P.,4X,10H XX-STRESS,5X,10H YY-STRESS,
     .       5X,10H XY-STRESS,5X,10H ZZ-STRESS,5X,10H XZ-STRESS,
     .       5X,10H YZ-STRESS,3X,12H DVOL.STRAIN)
  971 FORMAT(5X,5H G.P.,3X,10H     P1   ,5X,10H     P2   ,
     .       5X,10H     P3   ,4X,12H DVOL.STRAIN)
  972 FORMAT(5X,5H G.P.,4X,10H XX-STRESS,3X,12H DVOL.STRAIN)
  980 FORMAT(5X,5H G.P.,4X,11H Q-CO. CEN., 4X,11H P-CO. CEN.,
     .                  4X,11H MJ. S-AXIS, 4X,11H EQ.ME.PR. ,
     .                  4X,11H  SPEC.VOL., 2X,13H PLA.VOL.STR.)
  981 FORMAT(5X,I5,6E15.6)
  991 FORMAT(I5,2E15.6,5X,2E15.6,//,2(40X,2E15.6,/))
  992 FORMAT(I5,3E15.6,5X,3E15.6,//,3(55X,3E15.6,/))
  993 FORMAT(I5,6E15.6)
  994 FORMAT(5X,' G.P.    TOTAL PLASTIC DEFORMATION TENSOR',/,
     .      14X,'EFFECTIVE PLASTIC DEFORMATION',/,
     .      14X,'ISOTROPIC HARDENING FUNCTION (IF CONSIDERED)',/,
     .      14X,'KINEMATIC HARDENING TENSOR (IF CONSIDERED)',/,
     .      14X,'POROSITY AND TOTAL HARDENING (IF CONSIDERED)',/,
     .      14X,'EFFECTIVE PLASTIC DEFORMATION RATE (IF CONSIDERED)',/,
     .      14X,'RECRYSTALLIZATION VARIABLES (IF CONSIDERED)',/,
     .      14X,'VISCOUS VARIABLES (IF CONSIDERED)',/,
     .      14X,'DAMAGE VARIABLE (IF CONSIDERED)',/,
     .      14X,'SHRINKAGE (IF CONSIDERED)',/,
     .      14X,'FATIGUE VARIABLES (IF CONSIDERED)',/,
     .      14X,'DUAL PHASE STEEL VARIABLES (IF CONSIDERED)')
      END
