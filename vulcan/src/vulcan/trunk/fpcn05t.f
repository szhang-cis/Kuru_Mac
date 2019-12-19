      SUBROUTINE FPCN05T(TEINIT,VELCMT,PROPST,FPCHLT,IN1FPC,IN2FPC)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES:
C
C     IN2FPC=0 => THE PHASE-CHANGE FUNCTION (IN1FPC=0,1) AND ITS RATE
C                 (IN1FPC=1) AT NODAL POINTS
C                 (USEFUL FOR MECHANICAL PROBLEM)
C
C     IN2FPC=1 => THE PHASE-CHANGE FUNCTION TIMES THE LATENT HEAT
C                 (USEFUL FOR ALTERNATIVE OPTION FOR MICRO ADV. EFFECTS)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
      INCLUDE 'nuef_om.f'   ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION TEINIT(NDOFCT,*), VELCMT(*), PROPST(*), FPCHLT(NFPCH,*)
C
      DIMENSION TSOE1(5), TSOE2(5)
      DIMENSION TSOE1FI(5), TSOE2FI(5)
C
      IF(IN2FPC.EQ.0) THEN
C
C**** INITIALIZATION
C
       DO IFPCH=1,NNUPC       ! do not remove micro, additional & pseudo
        DO INODLT=1,NNODLT
         FPCHLT(IFPCH,INODLT)=0.0D0
         FPCHLT(IFPCH+NNUPT,INODLT)=0.0D0
        ENDDO
       ENDDO
C
       IF(NPLAT.EQ.0) GOTO 1002
C
C**** TEMPERATURE AT NODE
C
       DO INODLT=1,NNODLT
        TGAUST=TEINIT(1,INODLT)
        IF(IN1FPC.EQ.1) THEN
         DTEMPT=VELCMT(INODLT)
         TGAUAT=TGAUST-DTEMPT*DTIMET
        ENDIF
C
        INUPC=0
        DO IPLAT=1,NPLAT
         TSOE1(IPLAT)=0.0D0
         TSOE2(IPLAT)=0.0D0
         TEINF=VPLAT(IPLAT,1)
         TESUP=VPLAT(IPLAT,2)
         IPCFO=INT(VPLAT(IPLAT,4))
         IPCMO=INT(VPLAT(IPLAT,5))
         DELEE=TESUP-TEINF
C
C**** COMPUTES f_pc
C
         IF(IPCFO.EQ.1) GO TO 1001     ! skip microscopical phase-change
C
         IF(IPCMO.EQ.-1) GO TO 1001    ! skip no real phase-change
C
         INUPC=INUPC+1
         IF(IPCMO.EQ.0.OR.IPCMO.EQ.-2) THEN                ! f_pc linear
C
          IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .     TSOE1(IPLAT)=1.0D0/DELEE*(TGAUST-TEINF)
          IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
          IF(ICONVT.EQ.1) THEN
           IF(ICUBIC.EQ.0) THEN   ! fpc cubical (coherent with ILAHET=2)
            IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
            IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
             TEMRO=TEINF
             TEMPL=TESUP
             TEMPG=TGAUST
             AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                  3.0D0*TEMRO*TEMRO*(TEMPL-TEMRO)-
     .                     1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
             ALP11=-1.0D0/AAMAY
             ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
             ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
             ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-
     .                                                       ALP33*TEMPL
             ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+
     .                                                 ALP33*TEMPG+ALP44
             DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
             TSOE1(IPLAT)=1.0D0-ROTET
            ENDIF
            IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
           ENDIF                              ! icubic.eq.0
          ENDIF                               ! iconvt.eq.1
C
          IF(IN1FPC.EQ.1) THEN
           IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .      TSOE2(IPLAT)=1.0D0/DELEE*(TGAUAT-TEINF)
           IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
C
           IF(ICONVT.EQ.1) THEN
            IF(ICUBIC.EQ.0) THEN  ! fpc cubical (coherent with ILAHET=2)
             IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
             IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
              TEMRO=TEINF
              TEMPL=TESUP
              TEMPG=TGAUAT
              AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                  3.0D0*TEMRO*TEMRO*(TEMPL-TEMRO)-
     .                     1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
              ALP11=-1.0D0/AAMAY
              ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
              ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
              ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-
     .                                                       ALP33*TEMPL
              ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+
     .                                                 ALP33*TEMPG+ALP44
              DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
              TSOE2(IPLAT)=1.0D0-ROTET
             ENDIF
             IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
            ENDIF                             ! icubic.eq.0
           ENDIF                              ! iconvt.eq.1
          ENDIF                               ! in1fpc.eq.1
C
         ELSE                                    ! f_pc general
          GO TO (10,10,20,30,40,50), IPCMO
C
   10     CONTINUE                               ! piece-wise linear
          NPCFU=INT(VPLAT(IPLAT,6))
          IF(TGAUST.LE.VPCFU(IPLAT,1,2)) TSOE1(IPLAT)=VPCFU(IPLAT,1,1)
          DO IPCFU=2,NPCFU
           I1=IPCFU-1
           I2=IPCFU
           IF(TGAUST.GT.VPCFU(IPLAT,I1,2).AND.
     .        TGAUST.LE.VPCFU(IPLAT,I2,2))
     .      TSOE1(IPLAT)=(VPCFU(IPLAT,I2,1)-VPCFU(IPLAT,I1,1))/
     .                   (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .                   (TGAUST-VPCFU(IPLAT,I1,2))+VPCFU(IPLAT,I1,1)
          ENDDO
          IF(TGAUST.GT.VPCFU(IPLAT,NPCFU,2))
     .     TSOE1(IPLAT)=VPCFU(IPLAT,NPCFU,1)
C
          IF(IN1FPC.EQ.1) THEN
           IF(TGAUAT.LE.VPCFU(IPLAT,1,2)) TSOE2(IPLAT)=VPCFU(IPLAT,1,1)
           DO IPCFU=2,NPCFU
            I1=IPCFU-1
            I2=IPCFU
            IF(TGAUAT.GT.VPCFU(IPLAT,I1,2).AND.
     .         TGAUAT.LE.VPCFU(IPLAT,I2,2))
     .       TSOE2(IPLAT)=(VPCFU(IPLAT,I2,1)-VPCFU(IPLAT,I1,1))/
     .                    (VPCFU(IPLAT,I2,2)-VPCFU(IPLAT,I1,2))*
     .                    (TGAUAT-VPCFU(IPLAT,I1,2))+VPCFU(IPLAT,I1,1)
           ENDDO
           IF(TGAUAT.GT.VPCFU(IPLAT,NPCFU,2))
     .      TSOE2(IPLAT)=VPCFU(IPLAT,NPCFU,1)
          ENDIF         ! in1fpc.eq.1
C
          IF(IPCMO.EQ.2) THEN       ! piece-wise linear f_s (f_pc=1-f_s)
           TSOE1(IPLAT)=1.0D0-TSOE1(IPLAT)
           IF(IN1FPC.EQ.1) TSOE2(IPLAT)=1.0D0-TSOE2(IPLAT)
          ENDIF
          GO TO 1000
C
   20     CONTINUE                               ! parabolic
          IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .     TSOE1(IPLAT)=((TGAUST-TEINF)/DELEE)**2.0D0
          IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
          IF(IN1FPC.EQ.1) THEN
           IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .      TSOE2(IPLAT)=((TGAUAT-TEINF)/DELEE)**2.0D0
           IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
          ENDIF        ! in1fpc.eq.1
          GO TO 1000
C
   30     CONTINUE                               ! cubic
          CALL RUNENDT('FPCN05T=ERROR: IPCMO EQ 3')
          GO TO 1000
C
   40     CONTINUE                               ! Scheil's equation
          PARTC=VPLAT(IPLAT,6)
          PARTX=-1.0D0/(PARTC-1.0D0)
          IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .     TSOE1(IPLAT)=((TGAUST-TEINF)/DELEE)**PARTX
          IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
          IF(IN1FPC.EQ.1) THEN
           IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .      TSOE2(IPLAT)=((TGAUAT-TEINF)/DELEE)**PARTX
           IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
          ENDIF        ! in1fpc.eq.1
          GO TO 1000
C
   50     CONTINUE                               ! Lever equation
          TEMEL=VPLAT(IPLAT,6)
          PARTX=((TESUP-TEMEL)/(TGAUST-TEMEL))
          IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .     TSOE1(IPLAT)=((TGAUST-TEINF)/DELEE)*PARTX
          IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
          IF(IN1FPC.EQ.1) THEN
           PARTX=((TESUP-TEMEL)/(TGAUAT-TEMEL))
           IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .      TSOE2(IPLAT)=((TGAUAT-TEINF)/DELEE)*PARTX
           IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
          ENDIF                 ! in1fpc.eq.1
          GO TO 1000
C
 1000     CONTINUE
         ENDIF                      ! ipcmo.eq.0
C
C**** COMPUTES PSEUDO
C
         PSEUDO=1.0D0               ! metal
         IF(IFILL.EQ.1) THEN
          IF(IMICR.EQ.0) THEN
           IPSEU=2*NNUPT+1
          ELSE
           IPSEU=2*NNUPT+NNUPO+1
          ENDIF
          PSEUDO=FPCHLT(IPSEU,INODLT)
         ENDIF
         IF(IN1FPC.EQ.0) PSEUDO=1.0D0   ! assumption in the fpc_0 output
C
C**** COMPUTES FPCHLT AND ITS RATE
C
         FPCHLT(INUPC,INODLT)=PSEUDO*TSOE1(IPLAT)
         IF(IN1FPC.EQ.1)
     .    FPCHLT(INUPC+NNUPT,INODLT)=PSEUDO*(TSOE1(IPLAT)-TSOE2(IPLAT))/
     .                               DTIMET
C
 1001    CONTINUE
        ENDDO                       ! iplat=1,nplat
       ENDDO                        ! inodlt=1,nnodlt
C
 1002  CONTINUE
C
C**** DEALS WITH FILLING MATERIAL
C
       IF(IFILL.EQ.0) RETURN
C
       IF(NPLATFI.EQ.0) RETURN
C
C**** TEMPERATURE AT NODE
C
       DO INODLT=1,NNODLT
        TGAUST=TEINIT(1,INODLT)
        IF(IN1FPC.EQ.1) THEN
         DTEMPT=VELCMT(INODLT)
         TGAUAT=TGAUST-DTEMPT*DTIMET
        ENDIF
C
        INUPC=0
        DO IPLAT=1,NPLATFI
         TSOE1FI(IPLAT)=0.0D0
         TSOE2FI(IPLAT)=0.0D0
         TEINF=VPLATFI(IPLAT,1)
         TESUP=VPLATFI(IPLAT,2)
         IPCFO=INT(VPLATFI(IPLAT,4))
         IPCMO=INT(VPLATFI(IPLAT,5))
         DELEE=TESUP-TEINF
C
C**** COMPUTES f_pc
C
         IF(IPCFO.EQ.1) GO TO 2001     ! skip microscopical phase-change
C
         IF(IPCMO.EQ.-1) GO TO 2001    ! skip no real phase-change
C
         INUPC=INUPC+1
         IF(IPCMO.EQ.0.OR.IPCMO.EQ.-2) THEN                ! f_pc linear
C
          IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0D0
          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .     TSOE1FI(IPLAT)=1.0D0/DELEE*(TGAUST-TEINF)
          IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=1.0D0
C
          IF(ICONVT.EQ.1) THEN
           IF(ICUBIC.EQ.0) THEN   ! fpc cubical (coherent with ILAHET=2)
            IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0D0
            IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
             TEMRO=TEINF
             TEMPL=TESUP
             TEMPG=TGAUST
             AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                  3.0D0*TEMRO*TEMRO*(TEMPL-TEMRO)-
     .                     1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
             ALP11=-1.0D0/AAMAY
             ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
             ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
             ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-
     .                                                       ALP33*TEMPL
             ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+
     .                                                 ALP33*TEMPG+ALP44
             DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
             TSOE1FI(IPLAT)=1.0D0-ROTET
            ENDIF
            IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=1.0D0
           ENDIF                              ! icubic.eq.0
          ENDIF                               ! iconvt.eq.1
C
          IF(IN1FPC.EQ.1) THEN
           IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .      TSOE2FI(IPLAT)=1.0D0/DELEE*(TGAUAT-TEINF)
           IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=1.0D0
C
           IF(ICONVT.EQ.1) THEN
            IF(ICUBIC.EQ.0) THEN  ! fpc cubical (coherent with ILAHET=2)
             IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0D0
             IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
              TEMRO=TEINF
              TEMPL=TESUP
              TEMPG=TGAUAT
              AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                  3.0D0*TEMRO*TEMRO*(TEMPL-TEMRO)-
     .                     1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
              ALP11=-1.0D0/AAMAY
              ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
              ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
              ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-
     .                                                       ALP33*TEMPL
              ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+
     .              ALP33*TEMPG+ALP44
              DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
              TSOE2FI(IPLAT)=1.0D0-ROTET
             ENDIF
             IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=1.0D0
            ENDIF                             ! icubic.eq.0
           ENDIF                              ! iconvt.eq.1
          ENDIF                               ! in1fpc.eq.1
C
         ELSE                                    ! f_pc general
          GO TO (110,110,120,130,140,150), IPCMO
C
  110     CONTINUE                               ! piece-wise linear
          NPCFU=INT(VPLATFI(IPLAT,6))
          IF(TGAUST.LE.VPCFUFI(IPLAT,1,2)) TSOE1FI(IPLAT)=
     .                                     VPCFUFI(IPLAT,1,1)
          DO IPCFU=2,NPCFU
           I1=IPCFU-1
           I2=IPCFU
           IF(TGAUST.GT.VPCFUFI(IPLAT,I1,2).AND.
     .        TGAUST.LE.VPCFUFI(IPLAT,I2,2))
     .      TSOE1FI(IPLAT)=(VPCFUFI(IPLAT,I2,1)-VPCFUFI(IPLAT,I1,1))/
     .                  (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .                  (TGAUST-VPCFUFI(IPLAT,I1,2))+VPCFUFI(IPLAT,I1,1)
          ENDDO
          IF(TGAUST.GT.VPCFUFI(IPLAT,NPCFU,2))
     .     TSOE1FI(IPLAT)=VPCFUFI(IPLAT,NPCFU,1)
C
          IF(IN1FPC.EQ.1) THEN
           IF(TGAUAT.LE.VPCFUFI(IPLAT,1,2)) TSOE2FI(IPLAT)=
     .                                      VPCFUFI(IPLAT,1,1)
           DO IPCFU=2,NPCFU
            I1=IPCFU-1
            I2=IPCFU
            IF(TGAUAT.GT.VPCFUFI(IPLAT,I1,2).AND.
     .         TGAUAT.LE.VPCFUFI(IPLAT,I2,2))
     .       TSOE2FI(IPLAT)=(VPCFUFI(IPLAT,I2,1)-VPCFUFI(IPLAT,I1,1))/
     .                  (VPCFUFI(IPLAT,I2,2)-VPCFUFI(IPLAT,I1,2))*
     .                  (TGAUAT-VPCFUFI(IPLAT,I1,2))+VPCFUFI(IPLAT,I1,1)
           ENDDO
           IF(TGAUAT.GT.VPCFUFI(IPLAT,NPCFU,2))
     .      TSOE2FI(IPLAT)=VPCFUFI(IPLAT,NPCFU,1)
          ENDIF         ! in1fpc.eq.1
C
          IF(IPCMO.EQ.2) THEN       ! piece-wise linear f_s (f_pc=1-f_s)
           TSOE1FI(IPLAT)=1.0D0-TSOE1FI(IPLAT)
           IF(IN1FPC.EQ.1) TSOE2FI(IPLAT)=1.0-TSOE2FI(IPLAT)
          ENDIF
          GO TO 2000
C
  120     CONTINUE                               ! parabolic
          IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0D0
          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .     TSOE1FI(IPLAT)=((TGAUST-TEINF)/DELEE)**2.0D0
          IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=1.0D0
C
          IF(IN1FPC.EQ.1) THEN
           IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .      TSOE2FI(IPLAT)=((TGAUAT-TEINF)/DELEE)**2.0D0
           IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=1.0D0
          ENDIF        ! in1fpc.eq.1
          GO TO 2000
C
  130     CONTINUE                               ! cubic
          CALL RUNENDT('FPCN05T=ERROR: IPCMO EQ 3')
          GO TO 2000
C
  140     CONTINUE                               ! Scheil's equation
          PARTC=VPLATFI(IPLAT,6)
          PARTX=-1.0D0/(PARTC-1.0D0)
          IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0D0
          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .     TSOE1FI(IPLAT)=((TGAUST-TEINF)/DELEE)**PARTX
          IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=1.0D0
C
          IF(IN1FPC.EQ.1) THEN
           IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .      TSOE2FI(IPLAT)=((TGAUAT-TEINF)/DELEE)**PARTX
           IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=1.0D0
          ENDIF        ! in1fpc.eq.1
          GO TO 2000
C
  150     CONTINUE                               ! Lever equation
          TEMEL=VPLATFI(IPLAT,6)
          PARTX=((TESUP-TEMEL)/(TGAUST-TEMEL))
          IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0D0
          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .     TSOE1FI(IPLAT)=((TGAUST-TEINF)/DELEE)*PARTX
          IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=1.0D0
C
          IF(IN1FPC.EQ.1) THEN
           PARTX=((TESUP-TEMEL)/(TGAUAT-TEMEL))
           IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .      TSOE2FI(IPLAT)=((TGAUAT-TEINF)/DELEE)*PARTX
           IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=1.0D0
          ENDIF                 ! in1fpc.eq.1
          GO TO 2000
C
 2000     CONTINUE
         ENDIF                      ! ipcmo.eq.0
C
C**** COMPUTES PSEUDO 
C
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO=FPCHLT(IPSEU,INODLT)
         IF(IN1FPC.EQ.0) PSEUDO=1.0D0   ! assumption in the fpc_0 output
C
C**** COMPUTES FPCHLT AND ITS RATE
C
         FPCHLT(INUPC,INODLT)=FPCHLT(INUPC,INODLT)+
     .                                     (1.0D0-PSEUDO)*TSOE1FI(IPLAT)
         IF(IN1FPC.EQ.1)
     .    FPCHLT(INUPC+NNUPT,INODLT)=FPCHLT(INUPC+NNUPT,INODLT)+
     .             (1.0D0-PSEUDO)*(TSOE1FI(IPLAT)-TSOE2FI(IPLAT))/DTIMET
C
 2001    CONTINUE
        ENDDO                       ! iplat=1,nplatfi
       ENDDO                        ! inodlt=1,nnodlt
      ENDIF                         ! in2fpc.eq.0
C
      IF(IN2FPC.EQ.1) THEN
C
C**** INITIALIZATION
C
       IFPCH=2*NNUPT+NNUPO+NFILL+IGALFA+1
       DO INODLT=1,NNODLT
        FPCHLT(IFPCH,INODLT)=0.0D0                  ! only remove L*f_pc
        IF(KDYNAT.EQ.1) FPCHLT(IFPCH+1,INODLT)=0.0D0        ! & its rate
       ENDDO
C
       IF(NPLAT.EQ.0) GOTO 3002
C
C**** TEMPERATURE AT NODE
C
       DO INODLT=1,NNODLT
        TGAUST=TEINIT(1,INODLT)          ! at t + alpha dt is not needed
        TGAUSX=TGAUST                    ! at time t + dt
        IF(KDYNAT.EQ.1) THEN
         DTEMPT=VELCMT(INODLT)
        ENDIF
        TGAUIT=0.0D0                     ! not used
C
C**** COMPUTES PSEUDO
C
        PSEUDO=1.0D0               ! metal
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO=FPCHLT(IPSEU,INODLT)
        ENDIF
C
        ILAHET=1                               ! note: nmemo3=0 (imicr0)
        ISINRT=2
C
        CALL CAPCOFT(BASMM, BASCC, PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
C
        FPCHLT(IFPCH,INODLT)=SOUR2T
        IF(KDYNAT.EQ.1) FPCHLT(IFPCH+1,INODLT)=DSOURT/DTIMET
C
       ENDDO                        ! inodlt=1,nnodlt
 3002  CONTINUE
      ENDIF                         ! in2fpc.eq.1
C
      RETURN
      END
