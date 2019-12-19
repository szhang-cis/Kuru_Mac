      SUBROUTINE IDEIPT(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THE ISOTROPIC, PLASTIC &
C     TEMPERATURE-DEPENDENT MODEL
C
C***********************************************************************
C
C     NYOUN
C     VYOUN(1)=Young's Modulus                at VYOUN(2)=temperature
C
C     NPOIS
C     VPOIS(1)=Poisson's Ratio                at VPOIS(2)=temperature
C
C     NALPH
C     VALPH(1)=Thermal Dilatation Coefficient at VALPH(2)=temperature
C
C     NALPU=indicates secant or tangent Dilatation Coefficient
C
C     NCCER
C     VCCER(1)=Hardening Function             at VCCER(2)=temperature
C
C
C     ISOTT=indicator of isotropic hardening model
C
C     NCCOE
C     VCCOE(1)=Hardening Coefficient          at VCCOE(2)=temperature
C
C     NCCOB
C     VCCOB(1)=Coefficient b                  at VCCOE(2)=temperature
C
C     NCCOQ
C     VCCOQ(1)=Coefficient Q                  at VCCOE(2)=temperature
C
C     IPORO=indicator of porosity model (Gurson)
C
C     NPCOE
C     VPCOE(1)=Porosity Coefficient           at VPCOE(2)=temperature
C
C     NPCO1
C     VPCO1(1)=Q1 Coefficient                 at VPCOE(2)=temperature
C
C     NPCO2
C     VPCO2(1)=Q2 Coefficient                 at VPCOE(2)=temperature
C
C     TEREF=Reference Temperature
C
C     TEMPS=Solidus Temperature input in mechanical data file
C     TEMPL=Liquidus Temperature input in mechanical data file
C
C     ISHRI=Shrinkage indicator
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*)
C
      NCRIT=INT(PROPS(36))
      NCRIP=INT(PROPS(52))
      NASNA=INT(PROPS(47))
C
C**** YOUNG MODULUS
C
      I1=61
      NYOUN=INT(PROPS(I1))
C
      DO IYOUN=1,NYOUN
       IA=I1-1+2*IYOUN
       IB=IA+1
       VYOUN(IYOUN,1)=PROPS(IA)
       VYOUN(IYOUN,2)=PROPS(IB)
      ENDDO
C
C**** POISSON RATIO
C
      I2=I1+2*NYOUN+1
      NPOIS=INT(PROPS(I2))
C
      DO IPOIS=1,NPOIS
       IA=I2-1+2*IPOIS
       IB=IA+1
       VPOIS(IPOIS,1)=PROPS(IA)
       VPOIS(IPOIS,2)=PROPS(IB)
      ENDDO
C
C**** REFERENCE TEMPERATURE
C
      TEREF=PROPS(11)        ! see propipt.f
C
C**** THERMAL DILATATION
C
      I3=I2+2*NPOIS+1
      NALPH=INT(PROPS(I3))
C
      I3=I3+1
      NALPU=INT(PROPS(I3))
      VALPU(1,2)=PROPS(I3)
C
      DO IALPH=1,NALPH
       IA=I3-1+2*IALPH
       IB=IA+1
       VALPH(IALPH,1)=PROPS(IA)
       VALPH(IALPH,2)=PROPS(IB)
      ENDDO
C
      IF(NALPU.EQ.1) THEN                ! tangent alpha
       VALPU(1,1)=VALPH(1,1)             ! computes alpha^s
       DO IALPU=2,NALPH
        VALPU(IALPU,1)=VALPH(IALPU,1)    ! alpha^s=alpha
        VAUXD=VALPH(IALPU,2)-TEREF
        IF(VAUXD.GT.0) THEN
         VAUXX=0.0D+00
         DO JALPU=2,IALPU
          I1=JALPU-1
          I2=JALPU
          VAUXX=VAUXX+(VALPH(I2,1)+VALPH(I1,1))*
     .                (VALPH(I2,2)-VALPH(I1,2))/
     .                (VALPH(IALPU,2)-TEREF)
         ENDDO
         VALPU(IALPU,1)=0.5*VAUXX+VALPH(1,1)*(VALPH(1,2)-TEREF)/
     .                                    (VALPH(IALPU,2)-TEREF)
        ENDIF
       ENDDO
      ENDIF                             ! nalpu.eq.1
C
C**** THERMAL HARDENING FUNCTION
C
      I4=I3+2*NALPH+1
      NCCER=INT(PROPS(I4))
C
      DO ICCER=1,NCCER
       IA=I4-1+2*ICCER
       IB=IA+1
       VCCER(ICCER,1)=PROPS(IA)
       VCCER(ICCER,2)=PROPS(IB)
      ENDDO
C
C**** ISOTROPIC HARDENING INDICATOR
C
      I5=I4+2*NCCER+1
      ISOTT=INT(PROPS(I5))
c     IF(ISOTT.EQ.0) THEN
       I6=I5
       NCCOE=0
       NCCOEX=NCCOE
c     ENDIF
C
C**** HARDENING COEFFICIENT
C
      IF(ISOTT.EQ.1.OR.ISOTT.EQ.2.OR.ISOTT.EQ.4.OR.ISOTT.EQ.5.OR.
     .   ISOTT.EQ.8.OR.ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.
     .   ISOTT.EQ.13) THEN
       I6=I5+1
       NCCOE=INT(PROPS(I6))
C
       DO ICCOE=1,NCCOE
        IA=I6-1+2*ICCOE
        IB=IA+1
        VCCOE(ICCOE,1)=PROPS(IA)
        VCCOE(ICCOE,2)=PROPS(IB)
       ENDDO
       NCCOEX=NCCOE
      ENDIF
C
C**** B COEFFICIENT
C
      IF(ISOTT.EQ.3.OR.ISOTT.EQ.6.OR.ISOTT.EQ.7.OR.ISOTT.EQ.8.OR.
     .   ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.ISOTT.EQ.13) THEN
c      I6=I5+1
       I6=I6+2*NCCOE+1
       NCCOB=INT(PROPS(I6))
C
       DO ICCOB=1,NCCOB
        IA=I6-1+2*ICCOB
        IB=IA+1
        VCCOB(ICCOB,1)=PROPS(IA)
        VCCOB(ICCOB,2)=PROPS(IB)
       ENDDO
C
C**** Q COEFFICIENT
C
       I6A=I6+2*NCCOB+1
       NCCOQ=INT(PROPS(I6A))
C
       DO ICCOQ=1,NCCOQ
        IA=I6A-1+2*ICCOQ
        IB=IA+1
        VCCOQ(ICCOQ,1)=PROPS(IA)
        VCCOQ(ICCOQ,2)=PROPS(IB)
       ENDDO
C
       I6=I6A
       NCCOEX=NCCOQ
      ENDIF
C
C**** KINEMATIC HARDENING INDICATOR
C
      I7=I6+2*NCCOEX+1
      IKINE=INT(PROPS(I7))
      IF(IKINE.EQ.0) THEN
       I8=I7
       NKCOE=0
       NKCOEX=NKCOE
      ENDIF
C
C**** HARDENING COEFFICIENT
C
      IF(IKINE.EQ.1.OR.IKINE.EQ.2) THEN
       I8=I7+1
       NKCOE=INT(PROPS(I8))
C
       DO IKCOE=1,NKCOE
        IA=I8-1+2*IKCOE
        IB=IA+1
        VKCOE(IKCOE,1)=PROPS(IA)
        VKCOE(IKCOE,2)=PROPS(IB)
       ENDDO
       NKCOEX=NKCOE
      ENDIF
C
C**** B COEFFICIENT
C
      IF(IKINE.EQ.3.OR.IKINE.EQ.4) THEN
       I8=I7+1
       NKCOB=INT(PROPS(I8))
C
       DO IKCOB=1,NKCOB
        IA=I8-1+2*IKCOB
        IB=IA+1
        VKCOB(IKCOB,1)=PROPS(IA)
        VKCOB(IKCOB,2)=PROPS(IB)
       ENDDO
C
C**** Q COEFFICIENT
C
       I8A=I8+2*NKCOB+1
       NKCOQ=INT(PROPS(I8A))
C
       DO IKCOQ=1,NKCOQ
        IA=I8A-1+2*IKCOQ
        IB=IA+1
        VKCOQ(IKCOQ,1)=PROPS(IA)
        VKCOQ(IKCOQ,2)=PROPS(IB)
       ENDDO
C
       I8=I8A
       NKCOEX=NKCOQ
      ENDIF
C
      IVIFL=0
      IRECR=0
C
C**** COEFFICIENT 1
C
      IF(NCRIT.EQ.39) THEN
       I9=I8+2*NKCOEX+1
       NCOE1=INT(PROPS(I9))
C
       DO ICOE1=1,NCOE1
        IA=I9-1+2*ICOE1
        IB=IA+1
        VCOE1(ICOE1,1)=PROPS(IA)
        VCOE1(ICOE1,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** COEFFICIENT 2
C
      IF(NCRIT.EQ.39) THEN
       I9=I9+2*NCOE1+1
       NCOE2=INT(PROPS(I9))
C
       DO ICOE2=1,NCOE2
        IA=I9-1+2*ICOE2
        IB=IA+1
        VCOE2(ICOE2,1)=PROPS(IA)
        VCOE2(ICOE2,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** COEFFICIENT 3
C
      IF(NCRIT.EQ.39) THEN
       I9=I9+2*NCOE2+1
       NCOE3=INT(PROPS(I9))
C
       DO ICOE3=1,NCOE3
        IA=I9-1+2*ICOE3
        IB=IA+1
        VCOE3(ICOE3,1)=PROPS(IA)
        VCOE3(ICOE3,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** COEFFICIENT 4
C
      IF(NCRIT.EQ.39) THEN
       I9=I9+2*NCOE3+1
       NCOE4=INT(PROPS(I9))
C
       DO ICOE4=1,NCOE4
        IA=I9-1+2*ICOE4
        IB=IA+1
        VCOE4(ICOE4,1)=PROPS(IA)
        VCOE4(ICOE4,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** COEFFICIENT 5
C
      IF(NCRIT.EQ.39) THEN
       I9=I9+2*NCOE4+1
       NCOE5=INT(PROPS(I9))
C
       DO ICOE5=1,NCOE5
        IA=I9-1+2*ICOE5
        IB=IA+1
        VCOE5(ICOE5,1)=PROPS(IA)
        VCOE5(ICOE5,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** COEFFICIENT 6
C
      IF(NCRIP.EQ.40) THEN
       I9=I9+2*NCOE5+1
       NCOE6=INT(PROPS(I9))
C
       DO ICOE6=1,NCOE6
        IA=I9-1+2*ICOE6
        IB=IA+1
        VCOE6(ICOE6,1)=PROPS(IA)
        VCOE6(ICOE6,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** COEFFICIENT 7
C
      IF(NCRIP.EQ.40) THEN
       I9=I9+2*NCOE6+1
       NCOE7=INT(PROPS(I9))
C
       DO ICOE7=1,NCOE7
        IA=I9-1+2*ICOE7
        IB=IA+1
        VCOE7(ICOE7,1)=PROPS(IA)
        VCOE7(ICOE7,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** COEFFICIENT 8
C
      IF(NCRIP.EQ.40) THEN
       I9=I9+2*NCOE7+1
       NCOE8=INT(PROPS(I9))
C
       DO ICOE8=1,NCOE8
        IA=I9-1+2*ICOE8
        IB=IA+1
        VCOE8(ICOE8,1)=PROPS(IA)
        VCOE8(ICOE8,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** LANKFORD COEFFICIENT
C
      IF(NCRIT.EQ.42) THEN
       I9=I8+2*NKCOEX+1
       NLANK=INT(PROPS(I9))
C
       DO ILANK=1,NLANK
        IA=I9-1+2*ILANK
        IB=IA+1
        VLANK(ILANK,1)=PROPS(IA)
        VLANK(ILANK,2)=PROPS(IB)
       ENDDO
C
       I8=I9
       NKCOEX=NLANK
      ENDIF
C
C**** POROSITY INDICATOR
C
      NPCOE=NKCOEX
      IPORO=0
      IF(NCRIT.EQ.41.OR.NCRIP.EQ.41) THEN
       I7=I8+2*NKCOEX+1
       IPORO=INT(PROPS(I7))
C
C**** POROSITY COEFFICIENT
C
       IF(IPORO.EQ.1.OR.IPORO.EQ.2.OR.IPORO.EQ.3.OR.IPORO.EQ.4) THEN
        I8=I7+1
        NPCOE=INT(PROPS(I8))
C
        DO IPCOE=1,NPCOE
         IA=I8-1+2*IPCOE
         IB=IA+1
         VPCOE(IPCOE,1)=PROPS(IA)
         VPCOE(IPCOE,2)=PROPS(IB)
        ENDDO
       ENDIF
C
       IF(IPORO.EQ.4) THEN
        I8=I8+2*NPCOE+1
        NPCO1=INT(PROPS(I8))
C
        DO IPCO1=1,NPCO1
         IA=I8-1+2*IPCO1
         IB=IA+1
         VPCO1(IPCO1,1)=PROPS(IA)
         VPCO1(IPCO1,2)=PROPS(IB)
        ENDDO
C
        I8=I8+2*NPCO1+1
        NPCO2=INT(PROPS(I8))
C
        DO IPCO2=1,NPCO2
         IA=I8-1+2*IPCO2
         IB=IA+1
         VPCO2(IPCO2,1)=PROPS(IA)
         VPCO2(IPCO2,2)=PROPS(IB)
        ENDDO
        NPCOE=NPCO2
       ENDIF
      ENDIF
C
C**** PHASE-CHANGES
C
      I9=I8+2*NPCOE+1
      NPLATM=INT(PROPS(I9))
C
      IX=I9
      IF(NPLATM.NE.0) THEN
       DO IPLAT=1,NPLATM
        I9A=IX+1
C
        ITEMRO=INT(PROPS(I9A))
        VPLATM(IPLAT,7)=PROPS(I9A)
        I9A=I9A+1
        TEMPS=PROPS(I9A)
        VPLATM(IPLAT,8)=PROPS(I9A)
        I9A=I9A+1
        TEMPL=PROPS(I9A)
        VPLATM(IPLAT,9)=PROPS(I9A)
        I9A=I9A+1
        ISHRI=INT(PROPS(I9A))
        VPLATM(IPLAT,10)=PROPS(I9A)
        I9A=I9A+1
        ISEPC=INT(PROPS(I9A))
        VPLATM(IPLAT,11)=PROPS(I9A)
        I9A=I9A+1
        NTEMRO=5
C
        IPCFOM=INT(PROPS(I9A))
        VPLATM(IPLAT,4)=PROPS(I9A)
        I9B=I9A+1
        IPCMOM=INT(PROPS(I9B))
        VPLATM(IPLAT,5)=PROPS(I9B)
C
        NPCPCM=2                                  ! 2=ipcfo,ipcmo
        IF(IPCFOM.EQ.0) THEN
C
         IF(IPCMOM.EQ.0) THEN
          IA=I9B+1
          IB=IA+1
          IC=IB+1
          VPLATM(IPLAT,1)=PROPS(IA)               ! ILSPC
          VPLATM(IPLAT,2)=PROPS(IB)               ! ISSPC
          VPLATM(IPLAT,3)=PROPS(IC)               ! EXPAN (constant)
          VPLATM(IPLAT,6)=IPLAT                   ! index for f_pc
          IX=IX+NTEMRO+NPCPCM+3                   ! 3=ILSPC,ISSPC,EXPAN
         ENDIF
C
         IF(IPCMOM.EQ.1.OR.IPCMOM.EQ.2) THEN
          call runend('error in ideipt')
          I9C=I9B+1
          NPCFUM=INT(PROPS(I9C))
          VPLATM(IPLAT,6)=PROPS(I9C)
          I9D=I9C+1
          NLATEM=INT(PROPS(I9D))
          VPLATM(IPLAT,10)=PROPS(I9D)             ! vplat(10) is busy!!
          IF(NLATEM.EQ.0) THEN                    ! constant expansion
           DO IPCFU=1,NPCFUM
            IA=I9D-1+2*IPCFU
            IB=IA+1
            EXPAN(IPLAT,IPCFU,1)=PROPS(IA)
            EXPAN(IPLAT,IPCFU,2)=PROPS(IB)
           ENDDO
           VPLATM(IPLAT,1)=EXPAN(IPLAT,1,2)       ! Ts
           VPLATM(IPLAT,2)=EXPAN(IPLAT,NPCFUM,2)  ! Tl
           I9E=IB+1
           VPLATM(IPLAT,3)=PROPS(I9E)
           IX=IX+NTEMRO+2*NPCFUM+NPCPCM+3         ! 3=npcfu,nlate,L
          ELSE                                    ! temp.-dep. L
           I9D=I9D+1
           NTYPCM=INT(PROPS(I9D))                 ! secant or tangent L
           VPLATM(IPLAT,11)=PROPS(I9D)
           DO IPCFU=1,NPCFUM
            IA=I9D-2+3*IPCFU
            IB=IA+1
            IC=IB+1
            EXPAN(IPLAT,IPCFU,1)=PROPS(IA)
            EXPAN(IPLAT,IPCFU,2)=PROPS(IB)
            EXPAN(IPLAT,IPCFU,3)=PROPS(IC)
           ENDDO
C
           IF(NTYPCM.EQ.1) THEN                   ! tangent L
            EXPAU(IPLAT,1)=EXPAN(IPLAT,1,3)       ! computes L^s
            DO IPCFU=2,NPCFUM
             VAUXX=0.0D+00
             DO JPCFU=2,IPCFU
              I1=JPCFU-1
              I2=JPCFU
              IF(IPCMOM.EQ.1) THEN
               VAUXX=VAUXX+(EXPAN(IPLAT,I2,3)+EXPAN(IPLAT,I1,3))*
     .               (EXPAN(IPLAT,I2,1)-EXPAN(IPLAT,I1,1))/
     .               EXPAN(IPLAT,IPCFU,1)
              ELSE                                ! considers f_pc
               VAUXX=VAUXX+(EXPAN(IPLAT,I2,3)+EXPAN(IPLAT,I1,3))*
     .               (-EXPAN(IPLAT,I2,1)+EXPAN(IPLAT,I1,1))/
     .               (1.0-EXPAN(IPLAT,IPCFU,1))
              ENDIF
             ENDDO
             EXPAU(IPLAT,IPCFU)=VAUXX-EXPAU(IPLAT,1)
            ENDDO
           ENDIF                                  ! ntypcm.eq.1
C
           VPLATM(IPLAT,1)=EXPAN(IPLAT,1,2)       ! Ts
           VPLATM(IPLAT,2)=EXPAN(IPLAT,NPCFUM,2)  ! Tl
           VPLATM(IPLAT,3)=EXPAN(IPLAT,NPCFUM,3)  ! last L
           IX=IX+NTEMRO+3*NPCFUM+NPCPCM+3         ! 3=npcfu,nlate,ntypc
          ENDIF
         ENDIF
C
         IF(IPCMOM.EQ.3) THEN
          call runend('error in ideipt')
          IA=I9B+1
          IB=IA+1
          IC=IB+1
          VPLATM(IPLAT,1)=PROPS(IA)
          VPLATM(IPLAT,2)=PROPS(IB)
          VPLATM(IPLAT,3)=PROPS(IC)
          IX=IX+NTEMRO+NPCPCM+3                   ! 3=Ts,Tl,L
         ENDIF
C
         IF(IPCMOM.EQ.4) THEN
          call runend('error in ideipt')
         ENDIF
C
         IF(IPCMOM.EQ.5) THEN
          call runend('error in ideipt')
         ENDIF
C
         IF(IPCMOM.EQ.6) THEN
          call runend('error in ideipt')
          IA=I9B+1
          IB=IA+1
          IC=IB+1
          ID=IC+1
          VPLATM(IPLAT,1)=PROPS(IA)               ! Ts
          VPLATM(IPLAT,2)=PROPS(IB)               ! Tl
          VPLATM(IPLAT,6)=PROPS(IC)               ! Tm
          VPLATM(IPLAT,3)=PROPS(ID)               ! L
          IX=IX+NTEMRO+NPCPCM+4                   ! 4=Ts,Tl,Tm,L
         ENDIF
C
         IF(IPCMOM.EQ.7) THEN
          call runend('error in ideipt')
         ENDIF
C
        ELSE            ! ipcfom=1
         CALL RUNEND('IDEIPT: ERROR WITH IPCFOM=1 >> NOT POSSIBLE')
        ENDIF           ! ipcfom.eq.0
       ENDDO            ! do iplat=1,nplatm
      ENDIF             ! nplat.ne.0
C
C**** DAMAGE INDICATOR
C
      I10=IX+1
      IDAMG=INT(PROPS(I10))
      IF(IDAMG.EQ.0) THEN
       I12=I10
       NFRAC=0
      ENDIF
C
C**** DAMAGE STRESS
C
      IF(IDAMG.EQ.1.OR.IDAMG.EQ.2) THEN
       I11=I10+1
       NDAMA=INT(PROPS(I11))
C
       DO IDAMA=1,NDAMA
        IA=I11-1+2*IDAMA
        IB=IA+1
        VDAMA(IDAMA,1)=PROPS(IA)
        VDAMA(IDAMA,2)=PROPS(IB)
       ENDDO
C
C**** FRACTURE ENERGY
C
       I12=I11+2*NDAMA+1
       NFRAC=INT(PROPS(I12))
C
       DO IFRAC=1,NFRAC
        IA=I12-1+2*IFRAC
        IB=IA+1
        VFRAC(IFRAC,1)=PROPS(IA)
        VFRAC(IFRAC,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** COEFFICIENT 9
C
      IF(IDAMG.EQ.3) THEN
       I11=I10+1
       NCOE9=INT(PROPS(I11))
C
       DO ICOE9=1,NCOE9
        IA=I11-1+2*ICOE9
        IB=IA+1
        VCOE9(ICOE9,1)=PROPS(IA)
        VCOE9(ICOE9,2)=PROPS(IB)
       ENDDO
      ENDIF
C
C**** GAMMA PARAMETERS
C
      IF(IDAMG.EQ.21) THEN
       I11=I10+1
       GAMMAP=PROPS(I11)
       I12=I11+1
       GAMMAM=PROPS(I12)
       NFRAC=0
      ENDIF
C
C**** FREE ENERGY MODEL
C
      I13=I12+2*NFRAC+1
      IFREN=INT(PROPS(I13))
C
      RETURN
      END
