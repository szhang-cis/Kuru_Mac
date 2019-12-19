      SUBROUTINE MICPHAS(TGAUST,TGAUAT,ILAHET,
     .                   IPLAT,HENER,TSOE1,TSOE2,
     .                   TEINF,TESUP,DELEE,IPCMO,IAUXX)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION FOR MORE GENERAL
C     FORMS OF LATENT HEAT RELEASE
C
C***********************************************************************
C
C     TOTAL PHASE-CHANGE FORMULATION
C
C                     f_pc(T) is bijective
C
C                     IPCMO=1 PIECE WISE LINEAR
C                     IPCMO=2 IDEM 1 WITH f_s
C                     IPCMO=3 PARABOLIC
C                     IPCMO=4 CUBIC
C                     IPCMO=5 SCHEIL'S EQUATION
C                     IPCMO=6 LEVER EQUATION
C                     IPCMO=7 ....
C
C
C***********************************************************************
C
C     Index of variables:
C
C     TGAUST=current temperature at Gauss point
C     TGAUAT=last converged temperature at Gauss point
C     PROPST=array containing the thermal & microstructural material
C            properties
C     NPLAT =number of phase-changes
C     HENER =specific latent heat
C     TSOE1 =array contaning the phase-change function values for each
C            phase-change for the current step
C     TSOE2 =array contaning the phase-change function values for each
C            phase-change for the last converged step
C     TEINF ="solidus" temperature
C     TESUP ="liquidus" temperature
C     DELEE =temperature range "liquidus-solidus"
C
C     IAUXX = index for standard (=0) or filling (=1) material
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION TSOE1(5), TSOE2(5)
C
      GO TO (1,2) IAUXX+1
C
C**** DEALS WITH STANDARD MATERIAL
C
    1 HENE1=HENER
      HENE2=HENER
C
C**** COMPUTES f_pc
C
      IF(ILAHET.EQ.1) THEN
       GO TO (10,10,20,30,40,50), IPCMO
C
   10  CONTINUE                               ! piece-wise linear
       NPCFU=INT(VPLAT(IPLAT,6))
       IF(TGAUST.LE.VPCFU(IPLAT,1,2)) TSOE1(IPLAT)=VPCFU(IPLAT,1,1)
       DO IPCFU=2,NPCFU
        I1=IPCFU-1
        I2=IPCFU
        IF(TGAUST.GT.VPCFU(IPLAT,I1,2).AND.TGAUST.LE.VPCFU(IPLAT,I2,2))
     .   TSOE1(IPLAT)=(VPCFU(IPLAT,I2,1)-VPCFU(IPLAT,I1,1))/
     .                (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .                (TGAUST-VPCFU(IPLAT,I1,2))+VPCFU(IPLAT,I1,1)
       ENDDO
       IF(TGAUST.GT.VPCFU(IPLAT,NPCFU,2))
     .  TSOE1(IPLAT)=VPCFU(IPLAT,NPCFU,1)
C
       IF(TGAUAT.LE.VPCFU(IPLAT,1,2)) TSOE2(IPLAT)=VPCFU(IPLAT,1,1)
       DO IPCFU=2,NPCFU
        I1=IPCFU-1
        I2=IPCFU
        IF(TGAUAT.GT.VPCFU(IPLAT,I1,2).AND.TGAUAT.LE.VPCFU(IPLAT,I2,2))
     .   TSOE2(IPLAT)=(VPCFU(IPLAT,I2,1)-VPCFU(IPLAT,I1,1))/
     .                (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .                (TGAUAT-VPCFU(IPLAT,I1,2))+VPCFU(IPLAT,I1,1)
       ENDDO
       IF(TGAUAT.GT.VPCFU(IPLAT,NPCFU,2))
     .  TSOE2(IPLAT)=VPCFU(IPLAT,NPCFU,1)
C
       IF(IPCMO.EQ.2) THEN         ! piece-wise linear f_s (f_pc=1-f_s)
        TSOE1(IPLAT)=1.0D0-TSOE1(IPLAT)
        TSOE2(IPLAT)=1.0D0-TSOE2(IPLAT)
       ENDIF
C
       NLATE=INT(VPLAT(IPLAT,7))
       IF(NLATE.EQ.1) THEN                 ! temp.-dep. L
        NTYPC=INT(VPLAT(IPLAT,8))
        IF(NTYPC.EQ.0) THEN                ! secant L
         IF(TGAUST.LE.VPCFU(IPLAT,1,2)) HENE1=VPCFU(IPLAT,1,3)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUST.GT.VPCFU(IPLAT,I1,2).AND.
     .       TGAUST.LE.VPCFU(IPLAT,I2,2))
     .     HENE1=(VPCFU(IPLAT,I2,3)-VPCFU(IPLAT,I1,3))/
     .           (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .           (TGAUST-VPCFU(IPLAT,I1,2))+VPCFU(IPLAT,I1,3)
         ENDDO
         IF(TGAUST.GT.VPCFU(IPLAT,NPCFU,2)) HENE1=VPCFU(IPLAT,NPCFU,3)
C
         IF(TGAUAT.LE.VPCFU(IPLAT,1,2)) HENE2=VPCFU(IPLAT,1,3)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUAT.GT.VPCFU(IPLAT,I1,2).AND.
     .       TGAUAT.LE.VPCFU(IPLAT,I2,2))
     .     HENE2=(VPCFU(IPLAT,I2,3)-VPCFU(IPLAT,I1,3))/
     .           (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .           (TGAUAT-VPCFU(IPLAT,I1,2))+VPCFU(IPLAT,I1,3)
         ENDDO
         IF(TGAUAT.GT.VPCFU(IPLAT,NPCFU,2)) HENE2=VPCFU(IPLAT,NPCFU,3)
        ELSE                               ! tangent L
C
         IF(TGAUST.LE.VPCFU(IPLAT,1,2)) HENE1=VPCAU(IPLAT,1)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUST.GT.VPCFU(IPLAT,I1,2).AND.
     .       TGAUST.LE.VPCFU(IPLAT,I2,2))
     .     HENE1=(VPCAU(IPLAT,I2)-VPCAU(IPLAT,I1))/
     .           (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .           (TGAUST-VPCFU(IPLAT,I1,2))+VPCAU(IPLAT,I1)
         ENDDO
         IF(TGAUST.GT.VPCFU(IPLAT,NPCFU,2)) HENE1=VPCAU(IPLAT,NPCFU)
C
         IF(TGAUAT.LE.VPCFU(IPLAT,1,2)) HENE2=VPCAU(IPLAT,1)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUAT.GT.VPCFU(IPLAT,I1,2).AND.
     .       TGAUAT.LE.VPCFU(IPLAT,I2,2))
     .     HENE2=(VPCAU(IPLAT,I2)-VPCAU(IPLAT,I1))/
     .           (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .           (TGAUAT-VPCFU(IPLAT,I1,2))+VPCAU(IPLAT,I1)
         ENDDO
         IF(TGAUAT.GT.VPCFU(IPLAT,NPCFU,2)) HENE2=VPCAU(IPLAT,NPCFU)
        ENDIF             ! ntypc.eq.0
       ENDIF              ! nlate.eq.1
       GO TO 1000
C
   20  CONTINUE                               ! parabolic
       IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=((TGAUST-TEINF)/DELEE)**2.0D0
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
       IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=((TGAUAT-TEINF)/DELEE)**2.0D0
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
       GO TO 1000
C
   30  CONTINUE                               ! cubic
       CALL RUNENDT('MICPHAS=ERROR: IPCMO EQ 3')
       GO TO 1000
C
   40  CONTINUE                               ! Scheil's equation
       TEMEL=VPLAT(IPLAT,6)                   ! Melting temperature
       TEUTE=VPLAT(IPLAT,7)                   ! Eutectic temperature
       PARTI=(TEMEL-TESUP)/(TEMEL-TEINF)      ! partition coefficient
       PARTX=1.0D0/(PARTI-1.0D0)              ! Exponent
       IF(TGAUST.LE.TEUTE) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEUTE.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=((TGAUST-TEMEL)/(TESUP-TEMEL))**PARTX
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
       IF(TGAUAT.LE.TEUTE) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEUTE.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=((TGAUAT-TEMEL)/(TESUP-TEMEL))**PARTX
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
       GO TO 1000
C
   50  CONTINUE                               ! Lever equation
       TEMEL=VPLAT(IPLAT,6)                   ! Melting temperature
       PARTI=(TEMEL-TESUP)/(TEMEL-TEINF)      ! partition coefficient
       PARTX=1.0D0/(1.0D0-PARTI)              ! Factor
c      IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0 ! first form
c      IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
c    .  TSOE1(IPLAT)=1.0-((TGAUST-TESUP)/(TGAUST-TEMEL))*PARTX
c      IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
       PARTX=((TESUP-TEMEL)/(TGAUST-TEMEL))   ! second (equivalent) form
       IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=((TGAUST-TEINF)/DELEE)*PARTX
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
       PARTX=((TESUP-TEMEL)/(TGAUAT-TEMEL))
       IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=((TGAUAT-TEINF)/DELEE)*PARTX
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
       GO TO 1000
C
 1000  CONTINUE
C
       TSOE1(IPLAT)=TSOE1(IPLAT)*HENE1
       TSOE2(IPLAT)=TSOE2(IPLAT)*HENE2
C
      ENDIF                   ! ilahet.eq.1
C
C**** COMPUTES df_pc/dT
C
      IF(ILAHET.EQ.2) THEN
       GO TO (11,11,21,31,41,51), IPCMO
C
   11  CONTINUE                               ! piece-wise linear
       NPCFU=INT(VPLAT(IPLAT,6))
       IF(TGAUST.LE.VPCFU(IPLAT,1,2)) TSOE1(IPLAT)=0.0D0
       DO IPCFU=2,NPCFU
        I1=IPCFU-1
        I2=IPCFU
        IF(TGAUST.GT.VPCFU(IPLAT,I1,2).AND.TGAUST.LE.VPCFU(IPLAT,I2,2))
     .   TSOE1(IPLAT)=(VPCFU(IPLAT,I2,1)-VPCFU(IPLAT,I1,1))/
     .                (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))
       ENDDO
       IF(TGAUST.GT.VPCFU(IPLAT,NPCFU,2))
     .  TSOE1(IPLAT)=0.0D0
C
       IF(TGAUAT.LE.VPCFU(IPLAT,1,2)) TSOE2(IPLAT)=0.0D0
       DO IPCFU=2,NPCFU
        I1=IPCFU-1
        I2=IPCFU
        IF(TGAUAT.GT.VPCFU(IPLAT,I1,2).AND.TGAUAT.LE.VPCFU(IPLAT,I2,2))
     .   TSOE2(IPLAT)=(VPCFU(IPLAT,I2,1)-VPCFU(IPLAT,I1,1))/
     .                (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))
       ENDDO
       IF(TGAUAT.GT.VPCFU(IPLAT,NPCFU,2))
     .  TSOE2(IPLAT)=0.0D0
C
       IF(IPCMO.EQ.2) THEN   ! piece-wise linear f_s (df_pc/dT=-df_s/dT)
        TSOE1(IPLAT)=-TSOE1(IPLAT)
        TSOE2(IPLAT)=-TSOE2(IPLAT)
       ENDIF
C
       NLATE=INT(VPLAT(IPLAT,7))
       IF(NLATE.EQ.1) THEN                 ! temp.-dep. L
        NTYPC=INT(VPLAT(IPLAT,8))
        IF(NTYPC.EQ.0) THEN                ! secant L
         IF(TGAUST.LE.VPCFU(IPLAT,1,2)) HENE1=VPCFU(IPLAT,1,3)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUST.GT.VPCFU(IPLAT,I1,2).AND.
     .       TGAUST.LE.VPCFU(IPLAT,I2,2))
     .     HENE1=(VPCFU(IPLAT,I2,3)-VPCFU(IPLAT,I1,3))/
     .           (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .           (TGAUST-VPCFU(IPLAT,I1,2))+VPCFU(IPLAT,I1,3)
         ENDDO
         IF(TGAUST.GT.VPCFU(IPLAT,NPCFU,2)) HENE1=VPCFU(IPLAT,NPCFU,3)
C
         IF(TGAUAT.LE.VPCFU(IPLAT,1,2)) HENE2=VPCFU(IPLAT,1,3)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUAT.GT.VPCFU(IPLAT,I1,2).AND.
     .       TGAUAT.LE.VPCFU(IPLAT,I2,2))
     .     HENE2=(VPCFU(IPLAT,I2,3)-VPCFU(IPLAT,I1,3))/
     .           (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .           (TGAUAT-VPCFU(IPLAT,I1,2))+VPCFU(IPLAT,I1,3)
         ENDDO
         IF(TGAUAT.GT.VPCFU(IPLAT,NPCFU,2)) HENE2=VPCFU(IPLAT,NPCFU,3)
        ELSE                               ! tangent L
C
c        IF(IPCMO.EQ.2) THEN      ! see ideprot.f
c         DO IPCFU=1,NPCFU
c          VPCFU(IPLAT,IPCFU,1)=1.0D+00-VPCFU(IPLAT,IPCFU,1)
c         ENDDO
c        ENDIF
C
         IF(TGAUST.LE.VPCFU(IPLAT,1,2)) HENE1=VPCAU(IPLAT,1)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUST.GT.VPCFU(IPLAT,I1,2).AND.
     .       TGAUST.LE.VPCFU(IPLAT,I2,2))
     .     HENE1=(VPCAU(IPLAT,I2)-VPCAU(IPLAT,I1))/
     .           (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .           (TGAUST-VPCFU(IPLAT,I1,2))+VPCAU(IPLAT,I1)
         ENDDO
         IF(TGAUST.GT.VPCFU(IPLAT,NPCFU,2)) HENE1=VPCAU(IPLAT,NPCFU)
C
         IF(TGAUAT.LE.VPCFU(IPLAT,1,2)) HENE2=VPCAU(IPLAT,1)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUAT.GT.VPCFU(IPLAT,I1,2).AND.
     .       TGAUAT.LE.VPCFU(IPLAT,I2,2))
     .     HENE2=(VPCAU(IPLAT,I2)-VPCAU(IPLAT,I1))/
     .           (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .           (TGAUAT-VPCFU(IPLAT,I1,2))+VPCAU(IPLAT,I1)
         ENDDO
         IF(TGAUAT.GT.VPCFU(IPLAT,NPCFU,2)) HENE2=VPCAU(IPLAT,NPCFU)
        ENDIF             ! ntypc.eq.0
       ENDIF              ! nlate.eq.1
       GO TO 1001
C
   21  CONTINUE                               ! parabolic
       IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=2.0D0*(TGAUST-TEINF)/(DELEE**2.0D0)
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0D0
C
       IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=2.0D0*(TGAUAT-TEINF)/(DELEE**2.0D0)
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0D0
       GO TO 1001
C
   31  CONTINUE                               ! cubic
       CALL RUNENDT('MICPHAS=ERROR: IPCMO EQ 3')
       GO TO 1001
C
   41  CONTINUE                               ! Scheil's equation
       TEMEL=VPLAT(IPLAT,6)                   ! Melting temperature
       TEUTE=VPLAT(IPLAT,7)                   ! Eutectic temperature
       PARTI=(TEMEL-TESUP)/(TEMEL-TEINF)      ! partition coefficient
       PARTX=1.0D0/(PARTI-1.0D0)              ! Exponent
       IF(TGAUST.LE.TEUTE) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEUTE.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=PARTX*
     .       ((TGAUST-TEMEL)/(TESUP-TEMEL))**(PARTX-1.0D0)/(TESUP-TEMEL)
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0D0
C
       IF(TGAUAT.LE.TEUTE) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEUTE.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=PARTX*
     .       ((TGAUAT-TEMEL)/(TESUP-TEMEL))**(PARTX-1.0D0)/(TESUP-TEMEL)
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0D0
       GO TO 1001
C
   51  CONTINUE                               ! Lever equation
       TEMEL=VPLAT(IPLAT,6)
       PARTX=((TESUP-TEMEL)/(TGAUST-TEMEL))
       IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=((TEMEL-TESUP)*(TEMEL-TEINF))/
     .                                   (((TGAUST-TEMEL)**2.0D0)*DELEE)
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0D0
C
       IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=((TEMEL-TESUP)*(TEMEL-TEINF))/
     .                                   (((TGAUAT-TEMEL)**2.0D0)*DELEE)
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0D0
       GO TO 1001
C
 1001  CONTINUE
C
       TSOE1(IPLAT)=TSOE1(IPLAT)*HENE1
       TSOE2(IPLAT)=TSOE2(IPLAT)*HENE2
C
      ENDIF                   ! ilahet.eq.2
C
      RETURN
C
C**** DEALS WITH FILLING MATERIAL
C
    2 HENE1=HENER
      HENE2=HENER
C
C**** COMPUTES f_pc
C
      IF(ILAHET.EQ.1) THEN
       GO TO (100,100,200,300,400,500), IPCMO
C
  100  CONTINUE                               ! piece-wise linear
       NPCFU=INT(VPLATFI(IPLAT,6))
       IF(TGAUST.LE.VPCFUFI(IPLAT,1,2))
     .  TSOE1(IPLAT)=VPCFUFI(IPLAT,1,1)
       DO IPCFU=2,NPCFU
        I1=IPCFU-1
        I2=IPCFU
        IF(TGAUST.GT.VPCFUFI(IPLAT,I1,2).AND.
     .     TGAUST.LE.VPCFUFI(IPLAT,I2,2))
     .   TSOE1(IPLAT)=(VPCFUFI(IPLAT,I2,1)-VPCFU(IPLAT,I1,1))/
     .                (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .                (TGAUST-VPCFUFI(IPLAT,I1,2))+VPCFUFI(IPLAT,I1,1)
       ENDDO
       IF(TGAUST.GT.VPCFUFI(IPLAT,NPCFU,2))
     .  TSOE1(IPLAT)=VPCFUFI(IPLAT,NPCFU,1)
C
       IF(TGAUAT.LE.VPCFUFI(IPLAT,1,2))
     .  TSOE2(IPLAT)=VPCFUFI(IPLAT,1,1)
       DO IPCFU=2,NPCFU
        I1=IPCFU-1
        I2=IPCFU
        IF(TGAUAT.GT.VPCFUFI(IPLAT,I1,2).AND.
     .     TGAUAT.LE.VPCFUFI(IPLAT,I2,2))
     .   TSOE2(IPLAT)=(VPCFUFI(IPLAT,I2,1)-VPCFUFI(IPLAT,I1,1))/
     .                (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .                (TGAUAT-VPCFUFI(IPLAT,I1,2))+VPCFUFI(IPLAT,I1,1)
       ENDDO
       IF(TGAUAT.GT.VPCFUFI(IPLAT,NPCFU,2))
     .  TSOE2(IPLAT)=VPCFUFI(IPLAT,NPCFU,1)
C
       IF(IPCMO.EQ.2) THEN         ! piece-wise linear f_s (f_pc=1-f_s)
        TSOE1(IPLAT)=1.0D0-TSOE1(IPLAT)
        TSOE2(IPLAT)=1.0D0-TSOE2(IPLAT)
       ENDIF
C
       NLATE=INT(VPLATFI(IPLAT,7))
       IF(NLATE.EQ.1) THEN                 ! temp.-dep. L
        NTYPC=INT(VPLATFI(IPLAT,8))
        IF(NTYPC.EQ.0) THEN                ! secant L
         IF(TGAUST.LE.VPCFUFI(IPLAT,1,2)) HENE1=VPCFUFI(IPLAT,1,3)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUST.GT.VPCFUFI(IPLAT,I1,2).AND.
     .       TGAUST.LE.VPCFUFI(IPLAT,I2,2))
     .     HENE1=(VPCFUFI(IPLAT,I2,3)-VPCFUFI(IPLAT,I1,3))/
     .           (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .           (TGAUST-VPCFUFI(IPLAT,I1,2))+VPCFUFI(IPLAT,I1,3)
         ENDDO
         IF(TGAUST.GT.VPCFUFI(IPLAT,NPCFU,2))
     .    HENE1=VPCFUFI(IPLAT,NPCFU,3)
C
         IF(TGAUAT.LE.VPCFUFI(IPLAT,1,2)) HENE2=VPCFUFI(IPLAT,1,3)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUAT.GT.VPCFUFI(IPLAT,I1,2).AND.
     .       TGAUAT.LE.VPCFUFI(IPLAT,I2,2))
     .     HENE2=(VPCFUFI(IPLAT,I2,3)-VPCFUFI(IPLAT,I1,3))/
     .           (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .           (TGAUAT-VPCFUFI(IPLAT,I1,2))+VPCFUFI(IPLAT,I1,3)
         ENDDO
         IF(TGAUAT.GT.VPCFUFI(IPLAT,NPCFU,2))
     .    HENE2=VPCFUFI(IPLAT,NPCFU,3)
        ELSE                               ! tangent L
C
         IF(TGAUST.LE.VPCFUFI(IPLAT,1,2)) HENE1=VPCAUFI(IPLAT,1)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUST.GT.VPCFUFI(IPLAT,I1,2).AND.
     .       TGAUST.LE.VPCFUFI(IPLAT,I2,2))
     .     HENE1=(VPCAUFI(IPLAT,I2)-VPCAUFI(IPLAT,I1))/
     .           (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .           (TGAUST-VPCFUFI(IPLAT,I1,2))+VPCAUFI(IPLAT,I1)
         ENDDO
         IF(TGAUST.GT.VPCFUFI(IPLAT,NPCFU,2)) HENE1=VPCAUFI(IPLAT,NPCFU)
C
         IF(TGAUAT.LE.VPCFUFI(IPLAT,1,2)) HENE2=VPCAUFI(IPLAT,1)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUAT.GT.VPCFUFI(IPLAT,I1,2).AND.
     .       TGAUAT.LE.VPCFUFI(IPLAT,I2,2))
     .     HENE2=(VPCAUFI(IPLAT,I2)-VPCAUFI(IPLAT,I1))/
     .           (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .           (TGAUAT-VPCFUFI(IPLAT,I1,2))+VPCAUFI(IPLAT,I1)
         ENDDO
         IF(TGAUAT.GT.VPCFUFI(IPLAT,NPCFU,2)) HENE2=VPCAUFI(IPLAT,NPCFU)
        ENDIF             ! ntypc.eq.0
       ENDIF              ! nlate.eq.1
       GO TO 1002
C
  200  CONTINUE                               ! parabolic
       IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=((TGAUST-TEINF)/DELEE)**2.0D0
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
       IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=((TGAUAT-TEINF)/DELEE)**2.0D0
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
       GO TO 1002
C
  300  CONTINUE                               ! cubic
       CALL RUNENDT('MICPHAS=ERROR: IPCMO EQ 3')
       GO TO 1002
C
  400  CONTINUE                               ! Scheil's equation
       TEMEL=VPLATFI(IPLAT,6)                 ! Melting temperature
       TEUTE=VPLATFI(IPLAT,7)                 ! Eutectic temperature
       PARTI=(TEMEL-TESUP)/(TEMEL-TEINF)      ! partition coefficient
       PARTX=1.0D0/(PARTI-1.0D0)              ! Exponent
       IF(TGAUST.LE.TEUTE) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEUTE.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=((TGAUST-TEMEL)/(TESUP-TEMEL))**PARTX
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
       IF(TGAUAT.LE.TEUTE) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEUTE.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=((TGAUAT-TEMEL)/(TESUP-TEMEL))**PARTX
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
       GO TO 1002
C
  500  CONTINUE                               ! Lever equation
       TEMEL=VPLATFI(IPLAT,6)
       PARTX=((TESUP-TEMEL)/(TGAUST-TEMEL))
       IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=((TGAUST-TEINF)/DELEE)*PARTX
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
       PARTX=((TESUP-TEMEL)/(TGAUAT-TEMEL))
       IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=((TGAUAT-TEINF)/DELEE)*PARTX
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
       GO TO 1002
C
 1002  CONTINUE
C
       TSOE1(IPLAT)=TSOE1(IPLAT)*HENE1
       TSOE2(IPLAT)=TSOE2(IPLAT)*HENE2
C
      ENDIF                   ! ilahet.eq.1
C
C**** COMPUTES df_pc/dT
C
      IF(ILAHET.EQ.2) THEN
       GO TO (110,110,210,310,410,510), IPCMO
C
  110  CONTINUE                               ! piece-wise linear
       NPCFU=INT(VPLATFI(IPLAT,6))
       IF(TGAUST.LE.VPCFUFI(IPLAT,1,2)) TSOE1(IPLAT)=0.0D0
       DO IPCFU=2,NPCFU
        I1=IPCFU-1
        I2=IPCFU
        IF(TGAUST.GT.VPCFUFI(IPLAT,I1,2).AND.
     .     TGAUST.LE.VPCFUFI(IPLAT,I2,2))
     .   TSOE1(IPLAT)=(VPCFUFI(IPLAT,I2,1)-VPCFUFI(IPLAT,I1,1))/
     .                (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))
       ENDDO
       IF(TGAUST.GT.VPCFUFI(IPLAT,NPCFU,2))
     .  TSOE1(IPLAT)=0.0D0
C
       IF(TGAUAT.LE.VPCFUFI(IPLAT,1,2)) TSOE2(IPLAT)=0.0D0
       DO IPCFU=2,NPCFU
        I1=IPCFU-1
        I2=IPCFU
        IF(TGAUAT.GT.VPCFUFI(IPLAT,I1,2).AND.
     .     TGAUAT.LE.VPCFUFI(IPLAT,I2,2))
     .   TSOE2(IPLAT)=(VPCFUFI(IPLAT,I2,1)-VPCFUFI(IPLAT,I1,1))/
     .                (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))
       ENDDO
       IF(TGAUAT.GT.VPCFUFI(IPLAT,NPCFU,2))
     .  TSOE2(IPLAT)=0.0D0
C
       IF(IPCMO.EQ.2) THEN   ! piece-wise linear f_s (df_pc/dT=-df_s/dT)
        TSOE1(IPLAT)=-TSOE1(IPLAT)
        TSOE2(IPLAT)=-TSOE2(IPLAT)
       ENDIF
C
       NLATE=INT(VPLATFI(IPLAT,7))
       IF(NLATE.EQ.1) THEN                 ! temp.-dep. L
        NTYPC=INT(VPLATFI(IPLAT,8))
        IF(NTYPC.EQ.0) THEN                ! secant L
         IF(TGAUST.LE.VPCFUFI(IPLAT,1,2)) HENE1=VPCFUFI(IPLAT,1,3)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUST.GT.VPCFUFI(IPLAT,I1,2).AND.
     .       TGAUST.LE.VPCFUFI(IPLAT,I2,2))
     .     HENE1=(VPCFUFI(IPLAT,I2,3)-VPCFUFI(IPLAT,I1,3))/
     .           (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .           (TGAUST-VPCFUFI(IPLAT,I1,2))+VPCFUFI(IPLAT,I1,3)
         ENDDO
         IF(TGAUST.GT.VPCFUFI(IPLAT,NPCFU,2))
     .    HENE1=VPCFUFI(IPLAT,NPCFU,3)
C
         IF(TGAUAT.LE.VPCFUFI(IPLAT,1,2)) HENE2=VPCFUFI(IPLAT,1,3)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUAT.GT.VPCFUFI(IPLAT,I1,2).AND.
     .       TGAUAT.LE.VPCFUFI(IPLAT,I2,2))
     .     HENE2=(VPCFUFI(IPLAT,I2,3)-VPCFUFI(IPLAT,I1,3))/
     .           (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .           (TGAUAT-VPCFUFI(IPLAT,I1,2))+VPCFUFI(IPLAT,I1,3)
         ENDDO
         IF(TGAUAT.GT.VPCFUFI(IPLAT,NPCFU,2))
     .    HENE2=VPCFUFI(IPLAT,NPCFU,3)
        ELSE                               ! tangent L
C
c        IF(IPCMO.EQ.2) THEN          ! see ideprot.f
c         DO IPCFU=1,NPCFU
c          VPCFUFI(IPLAT,IPCFU,1)=1.0D+00-VPCFUFI(IPLAT,IPCFU,1)
c         ENDDO
c        ENDIF
C
         IF(TGAUST.LE.VPCFUFI(IPLAT,1,2)) HENE1=VPCAUFI(IPLAT,1)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUST.GT.VPCFUFI(IPLAT,I1,2).AND.
     .       TGAUST.LE.VPCFUFI(IPLAT,I2,2))
     .     HENE1=(VPCAUFI(IPLAT,I2)-VPCAUFI(IPLAT,I1))/
     .           (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .           (TGAUST-VPCFUFI(IPLAT,I1,2))+VPCAUFI(IPLAT,I1)
         ENDDO
         IF(TGAUST.GT.VPCFUFI(IPLAT,NPCFU,2)) HENE1=VPCAUFI(IPLAT,NPCFU)
C
         IF(TGAUAT.LE.VPCFUFI(IPLAT,1,2)) HENE2=VPCAUFI(IPLAT,1)
         DO IPCFU=2,NPCFU
          I1=IPCFU-1
          I2=IPCFU
          IF(TGAUAT.GT.VPCFUFI(IPLAT,I1,2).AND.
     .       TGAUAT.LE.VPCFUFI(IPLAT,I2,2))
     .     HENE2=(VPCAUFI(IPLAT,I2)-VPCAUFI(IPLAT,I1))/
     .           (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .           (TGAUAT-VPCFUFI(IPLAT,I1,2))+VPCAUFI(IPLAT,I1)
         ENDDO
         IF(TGAUAT.GT.VPCFUFI(IPLAT,NPCFU,2)) HENE2=VPCAUFI(IPLAT,NPCFU)
        ENDIF             ! ntypc.eq.0
       ENDIF              ! nlate.eq.1
       GO TO 1003
C
  210  CONTINUE                               ! parabolic
       IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=2.0D0*(TGAUST-TEINF)/(DELEE**2.0D0)
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0D0
C
       IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=2.0D0*(TGAUAT-TEINF)/(DELEE**2.0D0)
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0D0
       GO TO 1003
C
  310  CONTINUE                               ! cubic
       CALL RUNENDT('MICPHAS=ERROR: IPCMO EQ 3')
       GO TO 1003
C
  410  CONTINUE                               ! Scheil's equation
       TEMEL=VPLATFI(IPLAT,6)                 ! Melting temperature
       TEUTE=VPLATFI(IPLAT,7)                 ! Eutectic temperature
       PARTI=(TEMEL-TESUP)/(TEMEL-TEINF)      ! partition coefficient
       PARTX=1.0D0/(PARTI-1.0D0)              ! Exponent
       IF(TGAUST.LE.TEUTE) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEUTE.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=PARTX*
     .       ((TGAUST-TEMEL)/(TESUP-TEMEL))**(PARTX-1.0D0)/(TESUP-TEMEL)
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0D0
C
       IF(TGAUAT.LE.TEUTE) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEUTE.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=PARTX*
     .       ((TGAUAT-TEMEL)/(TESUP-TEMEL))**(PARTX-1.0D0)/(TESUP-TEMEL)
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0D0
       GO TO 1003
C
  510  CONTINUE                               ! Lever equation
       TEMEL=VPLAT(IPLAT,6)
       PARTX=((TESUP-TEMEL)/(TGAUST-TEMEL))
       IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
       IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .  TSOE1(IPLAT)=((TEMEL-TESUP)*(TEMEL-TEINF))/
     .                                   (((TGAUST-TEMEL)**2.0D0)*DELEE)
       IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0D0
C
       IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
       IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .  TSOE2(IPLAT)=((TEMEL-TESUP)*(TEMEL-TEINF))/
     .                                   (((TGAUAT-TEMEL)**2.0D0)*DELEE)
       IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0D0
       GO TO 1003
C
 1003  CONTINUE
C
       TSOE1(IPLAT)=TSOE1(IPLAT)*HENE1
       TSOE2(IPLAT)=TSOE2(IPLAT)*HENE2
C
      ENDIF                   ! ilahet.eq.2
C
      RETURN
C
      END
