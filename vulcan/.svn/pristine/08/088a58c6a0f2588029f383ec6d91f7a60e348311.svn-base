      SUBROUTINE IDEPROS3(PROPST,IPLAT,IA1)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 3 (IPCMO=3) OF RATE PHASE-CHANGE FORMULATIONS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION PROPST(*)
C
      IA2=IA1+2         ! 2=ipcfo,ipcmo
C
      SICON=PROPST(IA2+1)
      SBCON=PROPST(IA2+2)
      SICONE=PROPST(IA2+3)
C
      I1=IA2+4
      NLICO=INT(PROPST(I1))
      DO ILICO=1,NLICO
       IA=I1-1+2*ILICO
       IB=IA+1
       VGRO2(ILICO,1)=PROPST(IA)
       VGRO2(ILICO,2)=PROPST(IB)
      ENDDO
C
      DO ILICO=1,NLICO                 ! orders according to composition
       VGRO2(NLICO+ILICO,1)=VGRO2(ILICO,1)
       VGRO2(NLICO+ILICO,2)=VGRO2(ILICO,2)
      ENDDO
      DO ILICO=1,NLICO
       VMIN=VGRO2(ILICO,1)
       DO JLICO=ILICO,NLICO
        IF(VGRO2(JLICO,1).LT.VMIN) THEN
         VMIN=VGRO2(JLICO,1)
         VAUX1=VGRO2(NLICO+ILICO,1)
         VAUX2=VGRO2(NLICO+ILICO,2)
         VGRO2(NLICO+ILICO,1)=VGRO2(NLICO+JLICO,1)
         VGRO2(NLICO+ILICO,2)=VGRO2(NLICO+JLICO,2)
         VGRO2(NLICO+JLICO,1)=VAUX1
         VGRO2(NLICO+JLICO,2)=VAUX2
        ENDIF
       ENDDO
      ENDDO
C
      I1=I1+2*NLICO+1
      PARCO=PROPST(I1)
C
C**** SOLIDUS LINE:
C
C     TWO OPTIONS: 1) solidus line read as input (like the liquidus)
C                  2) solidus line obtained from the liquidus line and
C                     the partition coefficient (adopted Sept/98)
C
      IOPTION=2                        ! better as input
      IF(IOPTION.EQ.2) THEN
      DO ILICO=1,NLICO                 ! computes solidus line
       VGRO3(ILICO,1)=VGRO2(ILICO,1)*PARCO
       VGRO3(ILICO,2)=VGRO2(ILICO,2)
      ENDDO
C
      DO ILICO=1,NLICO                 ! orders according to composition
       VGRO3(NLICO+ILICO,1)=VGRO3(ILICO,1)
       VGRO3(NLICO+ILICO,2)=VGRO3(ILICO,2)
      ENDDO
      DO ILICO=1,NLICO
       VMIN=VGRO3(ILICO,1)
       DO JLICO=ILICO,NLICO
        IF(VGRO3(JLICO,1).LT.VMIN) THEN
         VMIN=VGRO3(JLICO,1)
         VAUX1=VGRO3(NLICO+ILICO,1)
         VAUX2=VGRO3(NLICO+ILICO,2)
         VGRO3(NLICO+ILICO,1)=VGRO3(NLICO+JLICO,1)
         VGRO3(NLICO+ILICO,2)=VGRO3(NLICO+JLICO,2)
         VGRO3(NLICO+JLICO,1)=VAUX1
         VGRO3(NLICO+JLICO,2)=VAUX2
        ENDIF
       ENDDO
      ENDDO
      ENDIF                              ! ioption.eq.2
C
      IMODECOM=3+1+2*NLICO+1             ! total prop. of comp.+equil.
C
      I1=I1+1
      IMNUC=INT(PROPST(I1))              ! dendritic nucleation model
C
      IF(IMNUC.EQ.1.OR.IMNUC.EQ.2) THEN
       TEMAV=PROPST(I1+1)
       TEMDE=PROPST(I1+2)
       GRDEM=PROPST(I1+3)
       IMODENUC=4                        ! 1+3=total prop. of nuc. model
      ENDIF
C
      I1=I1+IMODENUC
      IMGRO=INT(PROPST(I1))              ! dendritic growth model
C
      IF(IMGRO.EQ.1.OR.IMGRO.EQ.2.OR.IMGRO.EQ.3.OR.IMGRO.EQ.4.OR.
     .   IMGRO.EQ.5) THEN
       I1=I1+1
       NRAUN=INT(PROPST(I1))
       DO IRAUN=1,NRAUN
        IA=I1-1+2*IRAUN
        IB=IA+1
        VGRO1(IRAUN,1)=PROPST(IA)
        VGRO1(IRAUN,2)=PROPST(IB)
       ENDDO
C
       I1=I1+2*NRAUN+1
       DIFLI=PROPST(I1)
       I1=I1+1
       HENER=PROPST(I1)
       IMODEGRO=2+2*NRAUN+2              ! total prop. of growth model
      ENDIF
C
      I1=I1+1
      IMNUCE=INT(PROPST(I1))             ! eutectic nucleation model
C
      IF(IMNUCE.EQ.1.OR.IMNUCE.EQ.2) THEN
       TEMAVE=PROPST(I1+1)
       TEMDEE=PROPST(I1+2)
       GRDEME=PROPST(I1+3)
       IMODENUCE=4                       ! 1+3=total prop. of nuc. model
      ENDIF
C
      I1=I1+IMODENUCE
      IMGROE=INT(PROPST(I1))             ! eutectic growth model
C
      IF(IMGROE.EQ.1.OR.IMGROE.EQ.2.OR.IMGROE.EQ.3.OR.
     .   IMGROE.EQ.4.OR.IMGROE.EQ.5.OR.IMGROE.EQ.6) THEN
       I1=I1+1
       NRAUNE=INT(PROPST(I1))
       DO IRAUN=1,NRAUNE
        IA=I1-1+2*IRAUN
        IB=IA+1
        VGRO1E(IRAUN,1)=PROPST(IA)
        VGRO1E(IRAUN,2)=PROPST(IB)
       ENDDO
C
       I1=I1+2*NRAUNE+1
       HENERE=PROPST(I1)
       IMODEGROE=2+2*NRAUNE+1            ! total prop. of growth model
      ENDIF
C
      I1=I1+1
      IFPCDT= INT(PROPST(I1))
      IAFLOJ= INT(PROPST(I1+1))
      NSUGUEM=INT(PROPST(I1+2))
      EXTGUE=     PROPST(I1+3)
      IMODEFPC=4                         ! total prop. of numer. strat.
C
      VPLAT(IPLAT,6)=SICON
      VPLAT(IPLAT,7)=SBCON
      VPLAT(IPLAT,8)=SICONE
      VPLAT(IPLAT,9)=FLOAT(NLICO)
      VPLAT(IPLAT,10)=PARCO
C
      VPLAT(IPLAT,11)=FLOAT(IMNUC)
      VPLAT(IPLAT,12)=FLOAT(IMGRO)
      VPLAT(IPLAT,13)=FLOAT(IMNUCE)
      VPLAT(IPLAT,14)=FLOAT(IMGROE)
C
      IF(IMNUC.EQ.1.OR.IMNUC.EQ.2) THEN
       VPLAT(IPLAT,14+1)=TEMAV
       VPLAT(IPLAT,14+2)=TEMDE
       VPLAT(IPLAT,14+3)=GRDEM
       IVN=3
      ENDIF
C
      IF(IMGRO.EQ.1.OR.IMGRO.EQ.2.OR.IMGRO.EQ.3.OR.IMGRO.EQ.4.OR.
     .   IMGRO.EQ.5) THEN
       VPLAT(IPLAT,14+IVN+1)=FLOAT(NRAUN)
       VPLAT(IPLAT,14+IVN+2)=DIFLI
       VPLAT(IPLAT,14+IVN+3)=HENER
       IVG=3
      ENDIF
C
      IF(IMNUCE.EQ.1.OR.IMNUCE.EQ.2) THEN
       VPLAT(IPLAT,14+IVN+IVG+1)=TEMAVE
       VPLAT(IPLAT,14+IVN+IVG+2)=TEMDEE
       VPLAT(IPLAT,14+IVN+IVG+3)=GRDEME
       IVNE=3
      ENDIF
C
      IF(IMGROE.EQ.1.OR.IMGROE.EQ.2.OR.IMGROE.EQ.3.OR.
     .   IMGROE.EQ.4.OR.IMGROE.EQ.5.OR.IMGROE.EQ.6) THEN
       VPLAT(IPLAT,14+IVN+IVG+IVNE+1)=FLOAT(NRAUNE)
       VPLAT(IPLAT,14+IVN+IVG+IVNE+2)=HENERE
       IVGE=2
      ENDIF
C
      VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+1)=FLOAT(IFPCDT)
      VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+2)=FLOAT(IAFLOJ)
      VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+3)=FLOAT(NSUGUEM)
      VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+4)=      EXTGUE
      IVFPC=4
C
      IF((14+IVN+IVG+IVNE+IVGE+IVFPC).GT.50)
     . CALL RUNENDT('ERROR: INDEX OF VPLAT GT 50 - CHANGE auxl_omt.f')
C
      IMODE=                      ! imode=total num. of prop. of model 3
     .      IMODECOM+IMODENUC+IMODEGRO+
     .      IMODENUCE+IMODEGROE+IMODEFPC
C
      IA1=IA2+IMODE
C
      RETURN
      END
