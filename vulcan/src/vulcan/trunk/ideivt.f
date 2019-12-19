      SUBROUTINE IDEIVT(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THE ISOTROPIC,
C     VISCOPLASTIC & TEMPERATURE-DEPENDENT MODEL WITH DAMAGE
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
C     VCCER(1)=Thermal Hardening Function     at VCCER(2)=temperature
C
C     IVIFL=indicator of viscous law model
C
C     NVISC
C     VVISC(1)=Viscosity Function             at VVISC(2)=temperature
C
C     NZABA
C     VZABA(1)=Coefficient A                  at VZABA(2)=temperature
C
C     NZABB
C     VZABB(1)=Coefficient B                  at VZABB(2)=temperature
C
C     NZABC
C     VZABC(1)=Coefficient C                  at VZABC(2)=temperature
C
C     NEXPO
C     VEXPO(1)=Exponent Function              at VEXPO(2)=temperature
C
C     ISOTT=indicator of isotropic hardening model
C
C     NCCOE
C     VCCOE(1)=Hardening Coefficient          at VCCOE(2)=temperature
C
C     NCCOB
C     VCCOB(1)=Coefficient b                  at VCCOB(2)=temperature
C
C     NCCOQ
C     VCCOQ(1)=Coefficient Q                  at VCCOQ(2)=temperature
C
C     IPORO=indicator of porosity model (Gurson)
C
C     NPCOE
C     VPCOE(1)=Porosity Coefficient           at VPCOE(2)=temperature
C
C     NPCO1
C     VPCO1(1)=Q1 Coefficient                 at VPCO1(2)=temperature
C
C     NPCO2
C     VPCO2(1)=Q2 Coefficient                 at VPCO2(2)=temperature
C
C     TEREF=Reference Temperature
C
C     TEMPS=Solidus Temperature input in mechanical data file
C     TEMPL=Liquidus Temperature input in mechanical data file
C
C     ISHRI=Shrinkage indicator
C
C     IDAMG=indicator of damage model
C
C     NDAMA
C     VDAMA(1)=Damage stress                  at VDAMA(2)=temperature
C
C     NFRAC
C     VFRAC(1)=Fracture energy                at VFRAC(2)=temperature
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
C
C**** YOUNG MODULUS
C
      IX=61
      NYOUN=INT(PROPS(IX))
C
      DO IYOUN=1,NYOUN
       IA=IX-1+2*IYOUN
       IB=IA+1
       VYOUN(IYOUN,1)=PROPS(IA)
       VYOUN(IYOUN,2)=PROPS(IB)
      ENDDO
      IX=IX+2*NYOUN
C
C**** POISSON RATIO
C
      IX=IX+1
      NPOIS=INT(PROPS(IX))
C
      DO IPOIS=1,NPOIS
       IA=IX-1+2*IPOIS
       IB=IA+1
       VPOIS(IPOIS,1)=PROPS(IA)
       VPOIS(IPOIS,2)=PROPS(IB)
      ENDDO
      IX=IX+2*NPOIS
C
C**** REFERENCE TEMPERATURE
C
      TEREF=PROPS(11)        ! see propipt.f
C
C**** THERMAL DILATATION
C
      IX=IX+1
      NALPH=INT(PROPS(IX))
C
      IX=IX+1
      NALPU=INT(PROPS(IX))
      VALPU(1,2)=PROPS(IX)
C
      DO IALPH=1,NALPH
       IA=IX-1+2*IALPH
       IB=IA+1
       VALPH(IALPH,1)=PROPS(IA)
       VALPH(IALPH,2)=PROPS(IB)
      ENDDO
      IX=IX+2*NALPH
C
      IF(NALPU.EQ.1) THEN                ! tangent alpha
       VALPU(1,1)=VALPH(1,1)             ! computes alpha^s
       DO IALPU=2,NALPH
        VALPU(IALPU,1)=VALPH(IALPU,1)    ! alpha^s=alpha
        VAUXD=VALPH(IALPU,2)-TEREF
        IF(VAUXD.GT.0) THEN
         VAUXX=0.0D0
         DO JALPU=2,IALPU
          I1=JALPU-1
          I2=JALPU
          VAUXX=VAUXX+(VALPH(I2,1)+VALPH(I1,1))*
     .                (VALPH(I2,2)-VALPH(I1,2))/
     .                (VALPH(IALPU,2)-TEREF)
         ENDDO
         VALPU(IALPU,1)=0.5D0*VAUXX+VALPH(1,1)*(VALPH(1,2)-TEREF)/
     .                                     (VALPH(IALPU,2)-TEREF)
        ENDIF
       ENDDO
      ENDIF                             ! nalpu.eq.1
C
C**** THERMAL HARDENING FUNCTION
C
      IX=IX+1
      NCCER=INT(PROPS(IX))
C
      DO ICCER=1,NCCER
       IA=IX-1+2*ICCER
       IB=IA+1
       VCCER(ICCER,1)=PROPS(IA)
       VCCER(ICCER,2)=PROPS(IB)
      ENDDO
      IX=IX+2*NCCER
C
C**** FLOW PARAMETERS (VISCOUS LAW MODEL)
C
      IX=IX+1
      IVIFL=INT(PROPS(IX))
C
C**** VISCOSITY FUNCTION OR SMA COEFFICIENT
C
      IF(IVIFL.EQ.1.OR.IVIFL.EQ.2.OR.IVIFL.EQ.3.OR.IVIFL.EQ.6) THEN
       IX=IX+1
       NVISC=INT(PROPS(IX))
C
       DO IVISC=1,NVISC
        IA=IX-1+2*IVISC
        IB=IA+1
        VVISC(IVISC,1)=PROPS(IA)
        VVISC(IVISC,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NVISC
      ENDIF          ! ivifl.eq.1.or.ivifl.eq.2.or.ivifl.eq.3.or. ...
C
C**** ZABARAS' COEFFICIENTS (A, B & C) OR OR SMA COEFFICIENTS
C
      IF(IVIFL.EQ.4.OR.IVIFL.EQ.6) THEN
C
C**** COEFFICIENT A
C
       IX=IX+1
       NZABA=INT(PROPS(IX))
C
       DO IZABA=1,NZABA
        IA=IX-1+2*IZABA
        IB=IA+1
        VZABA(IZABA,1)=PROPS(IA)
        VZABA(IZABA,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NZABA
C
C**** COEFFICIENT B
C
       IX=IX+1
       NZABB=INT(PROPS(IX))
C
       DO IZABB=1,NZABB
        IA=IX-1+2*IZABB
        IB=IA+1
        VZABB(IZABB,1)=PROPS(IA)
        VZABB(IZABB,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NZABB
C
C**** COEFFICIENT C
C
       IX=IX+1
       NZABC=INT(PROPS(IX))
C
       DO IZABC=1,NZABC
        IA=IX-1+2*IZABC
        IB=IA+1
        VZABC(IZABC,1)=PROPS(IA)
        VZABC(IZABC,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NZABC
      ENDIF          ! ivifl.eq.4.or.ivifl.eq.6
C
C**** EXPONENT FUNCTION OR SMA COEFFICIENT
C
      IF((IVIFL.GE.1.AND.IVIFL.LE.4).OR.IVIFL.EQ.6) THEN
       IX=IX+1
       NEXPO=INT(PROPS(IX))
C
       DO IEXPO=1,NEXPO
        IA=IX-1+2*IEXPO
        IB=IA+1
        VEXPO(IEXPO,1)=PROPS(IA)
        VEXPO(IEXPO,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NEXPO
      ENDIF          ! (ivifl.ge.1.or.ivifl.le.4).or. ...
C
C**** BUSSO'S MODEL
C
      IF(IVIFL.EQ.5) THEN
       NBUS1=6
       DO I=1,NBUS1
        IA=IX+I
        RECRY(I)=PROPS(IA)
       ENDDO
       IX=IX+NBUS1
      ENDIF          ! ivifl.eq.5
C
C**** ISOTROPIC HARDENING INDICATOR
C
      IX=IX+1
      ISOTT=INT(PROPS(IX))
C
C**** HARDENING COEFFICIENT
C
      IF(ISOTT.EQ.1.OR.ISOTT.EQ.2.OR.ISOTT.EQ.4.OR.ISOTT.EQ.5.OR.
     .   ISOTT.EQ.8.OR.ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.
     .   ISOTT.EQ.13) THEN
       IX=IX+1
       NCCOE=INT(PROPS(IX))
C
       DO ICCOE=1,NCCOE
        IA=IX-1+2*ICCOE
        IB=IA+1
        VCCOE(ICCOE,1)=PROPS(IA)
        VCCOE(ICCOE,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCCOE
      ENDIF
C
C**** B COEFFICIENT
C
      IF(ISOTT.EQ.3.OR.ISOTT.EQ.6.OR.ISOTT.EQ.7.OR.ISOTT.EQ.8.OR.
     .   ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.ISOTT.EQ.13) THEN
       IX=IX+1
       NCCOB=INT(PROPS(IX))
C
       DO ICCOB=1,NCCOB
        IA=IX-1+2*ICCOB
        IB=IA+1
        VCCOB(ICCOB,1)=PROPS(IA)
        VCCOB(ICCOB,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCCOB
C
C**** Q COEFFICIENT
C
       IX=IX+1
       NCCOQ=INT(PROPS(IX))
C
       DO ICCOQ=1,NCCOQ
        IA=IX-1+2*ICCOQ
        IB=IA+1
        VCCOQ(ICCOQ,1)=PROPS(IA)
        VCCOQ(ICCOQ,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCCOQ
      ENDIF
C
C**** KINEMATIC HARDENING INDICATOR
C
      IX=IX+1
      IKINE=INT(PROPS(IX))
C
C**** HARDENING COEFFICIENT
C
      IF(IKINE.EQ.1.OR.IKINE.EQ.2) THEN
       IX=IX+1
       NKCOE=INT(PROPS(IX))
C
       DO IKCOE=1,NKCOE
        IA=IX-1+2*IKCOE
        IB=IA+1
        VKCOE(IKCOE,1)=PROPS(IA)
        VKCOE(IKCOE,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NKCOE
      ENDIF
C
C**** B COEFFICIENT
C
      IF(IKINE.EQ.3.OR.IKINE.EQ.4) THEN
       IX=IX+1
       NKCOB=INT(PROPS(IX))
C
       DO IKCOB=1,NKCOB
        IA=IX-1+2*IKCOB
        IB=IA+1
        VKCOB(IKCOB,1)=PROPS(IA)
        VKCOB(IKCOB,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NKCOB
C
C**** Q COEFFICIENT
C
       IX=IX+1
       NKCOQ=INT(PROPS(IX))
C
       DO IKCOQ=1,NKCOQ
        IA=IX-1+2*IKCOQ
        IB=IA+1
        VKCOQ(IKCOQ,1)=PROPS(IA)
        VKCOQ(IKCOQ,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NKCOQ
      ENDIF
C
C**** RECRYSTALLIZATION MODEL
C
      IX=IX+1
      IRECR=INT(PROPS(IX))
C
C**** INITIAL, DEFORMATION RESISTANCE AND GRAIN GROWTH COEFFICIENTS
C
      IF(IRECR.EQ.1) THEN
       NBUS2=14
       DO I=1,NBUS2
        IA=IX+I
        RECRY(NBUS1+I)=PROPS(IA)
       ENDDO
       IX=IX+NBUS2
      ENDIF          ! irecr.eq.1
C
C**** COEFFICIENT 1 (GREEN SAND)
C
      IF(NCRIT.EQ.39) THEN
       IX=IX+1
       NCOE1=INT(PROPS(IX))
C
       DO ICOE1=1,NCOE1
        IA=IX-1+2*ICOE1
        IB=IA+1
        VCOE1(ICOE1,1)=PROPS(IA)
        VCOE1(ICOE1,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCOE1
C
C**** COEFFICIENT 2 (GREEN SAND)
C
       IX=IX+1
       NCOE2=INT(PROPS(IX))
C
       DO ICOE2=1,NCOE2
        IA=IX-1+2*ICOE2
        IB=IA+1
        VCOE2(ICOE2,1)=PROPS(IA)
        VCOE2(ICOE2,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCOE2
C
C**** COEFFICIENT 3 (GREEN SAND)
C
       IX=IX+1
       NCOE3=INT(PROPS(IX))
C
       DO ICOE3=1,NCOE3
        IA=IX-1+2*ICOE3
        IB=IA+1
        VCOE3(ICOE3,1)=PROPS(IA)
        VCOE3(ICOE3,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCOE3
C
C**** COEFFICIENT 4 (GREEN SAND)
C
       IX=IX+1
       NCOE4=INT(PROPS(IX))
C
       DO ICOE4=1,NCOE4
        IA=IX-1+2*ICOE4
        IB=IA+1
        VCOE4(ICOE4,1)=PROPS(IA)
        VCOE4(ICOE4,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCOE4
C
C**** COEFFICIENT 5 (GREEN SAND)
C
       IX=IX+1
       NCOE5=INT(PROPS(IX))
C
       DO ICOE5=1,NCOE5
        IA=IX-1+2*ICOE5
        IB=IA+1
        VCOE5(ICOE5,1)=PROPS(IA)
        VCOE5(ICOE5,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCOE5
      ENDIF
C
C**** COEFFICIENT 6 (GREEN SAND)
C
      IF(NCRIP.EQ.40) THEN
       IX=IX+1
       NCOE6=INT(PROPS(IX))
C
       DO ICOE6=1,NCOE6
        IA=IX-1+2*ICOE6
        IB=IA+1
        VCOE6(ICOE6,1)=PROPS(IA)
        VCOE6(ICOE6,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCOE6
C
C**** COEFFICIENT 7 (GREEN SAND)
C
       IX=IX+1
       NCOE7=INT(PROPS(IX))
C
       DO ICOE7=1,NCOE7
        IA=IX-1+2*ICOE7
        IB=IA+1
        VCOE7(ICOE7,1)=PROPS(IA)
        VCOE7(ICOE7,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCOE7
C
C**** COEFFICIENT 8
C
       IX=IX+1
       NCOE8=INT(PROPS(IX))
C
       DO ICOE8=1,NCOE8
        IA=IX-1+2*ICOE8
        IB=IA+1
        VCOE8(ICOE8,1)=PROPS(IA)
        VCOE8(ICOE8,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCOE8
      ENDIF
C
C**** POROSITY INDICATOR
C
      IPORO=0
      IF(NCRIT.EQ.41.OR.NCRIP.EQ.41) THEN
       IX=IX+1
       IPORO=INT(PROPS(IX))
C
C**** POROSITY COEFFICIENT
C
       IF(IPORO.EQ.1.OR.IPORO.EQ.2.OR.IPORO.EQ.3.OR.IPORO.EQ.4) THEN
        IX=IX+1
        NPCOE=INT(PROPS(IX))
C
        DO IPCOE=1,NPCOE
         IA=IX-1+2*IPCOE
         IB=IA+1
         VPCOE(IPCOE,1)=PROPS(IA)
         VPCOE(IPCOE,2)=PROPS(IB)
        ENDDO
        IX=IX+2*NPCOE
       ENDIF
C
       IF(IPORO.EQ.4) THEN
        IX=IX+1
        NPCO1=INT(PROPS(IX))
C
        DO IPCO1=1,NPCO1
         IA=IX-1+2*IPCO1
         IB=IA+1
         VPCO1(IPCO1,1)=PROPS(IA)
         VPCO1(IPCO1,2)=PROPS(IB)
        ENDDO
        IX=IX+2*NPCO1
C
        IX=IX+1
        NPCO2=INT(PROPS(IX))
C
        DO IPCO2=1,NPCO2
         IA=IX-1+2*IPCO2
         IB=IA+1
         VPCO2(IPCO2,1)=PROPS(IA)
         VPCO2(IPCO2,2)=PROPS(IB)
        ENDDO
        IX=IX+2*NPCO2
       ENDIF
      ENDIF
C
C**** PHASE-CHANGES
C
      IX=IX+1
      NPLATM=INT(PROPS(IX))
C
      IF(NPLATM.NE.0) THEN
       DO IPLAT=1,NPLATM
        IX=IX+1
C
        ITEMRO=INT(PROPS(IX))
        VPLATM(IPLAT,7)=PROPS(IX)
        IX=IX+1
        TEMPS=PROPS(IX)
        VPLATM(IPLAT,8)=PROPS(IX)
        IX=IX+1
        TEMPL=PROPS(IX)
        VPLATM(IPLAT,9)=PROPS(IX)
        IX=IX+1
        ISHRI=INT(PROPS(IX))
        VPLATM(IPLAT,10)=PROPS(IX)
        IX=IX+1
        ISEPC=INT(PROPS(IX))
        VPLATM(IPLAT,11)=PROPS(IX)
        IX=IX+1
        NTEMRO=5
C
        IPCFOM=INT(PROPS(IX))
        VPLATM(IPLAT,4)=PROPS(IX)
        IX=IX+1
        IPCMOM=INT(PROPS(IX))
        VPLATM(IPLAT,5)=PROPS(IX)
C
        NPCPCM=2                              ! 2=ipcfo,ipcmo
        IF(IPCFOM.EQ.0) THEN
C
         IF(IPCMOM.EQ.0) THEN
          IA=IX+1
          IB=IA+1
          IC=IB+1
          VPLATM(IPLAT,1)=PROPS(IA)               ! ILSPC
          VPLATM(IPLAT,2)=PROPS(IB)               ! ISSPC
          VPLATM(IPLAT,3)=PROPS(IC)               ! EXPAN (constant)
          VPLATM(IPLAT,6)=IPLAT                   ! index for f_pc
          IX=IX+3                                 ! 3=ILSPC,ISSPC,EXPAN
         ENDIF
C
         IF(IPCMOM.EQ.1.OR.IPCMOM.EQ.2) THEN      ! to be revised
          IX=IX+1
          NPCFUM=INT(PROPS(IX))
          VPLATM(IPLAT,6)=PROPS(IX)
          IX=IX+1
          NLATEM=INT(PROPS(IX))
          VPLATM(IPLAT,7)=PROPS(IX)
          IF(NLATEM.EQ.0) THEN                    ! constant expansion
           DO IPCFU=1,NPCFUM
            IA=IX-1+2*IPCFU
            IB=IA+1
            EXPAN(IPLAT,IPCFU,1)=PROPS(IA)
            EXPAN(IPLAT,IPCFU,2)=PROPS(IB)
           ENDDO
           VPLATM(IPLAT,1)=EXPAN(IPLAT,1,2)       ! Ts
           VPLATM(IPLAT,2)=EXPAN(IPLAT,NPCFUM,2)  ! Tl
           IC=IB+1
           VPLATM(IPLAT,3)=PROPS(IC)
           IX=IX+NTEMRO+2*NPCFUM+NPCPCM+3         ! 3=npcfu,nlate,L
          ELSE                                    ! temp.-dep. L
           IX=IX+1
           NTYPCM=INT(PROPS(IX))                  ! secant or tangent L
           VPLATM(IPLAT,8)=PROPS(IX)
           DO IPCFU=1,NPCFUM
            IA=IX-2+3*IPCFU
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
             VAUXX=0.0D0
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
     .               (1.0D0-EXPAN(IPLAT,IPCFU,1))
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
         IF(IPCMOM.EQ.3) THEN                     ! to be revised
          IA=IX+1
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
         IF(IPCMOM.EQ.6) THEN                     ! to be revised
          IA=IX+1
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
         CALL RUNEND('IDEIVT: ERROR WITH IPCFOM=1 >> NOT POSSIBLE')
        ENDIF           ! ipcfom.eq.0
       ENDDO            ! do iplat=1,nplatm
      ENDIF             ! nplat.ne.0
C
C**** DAMAGE INDICATOR
C
      IX=IX+1
      IDAMG=INT(PROPS(IX))
C
C**** DAMAGE STRESS
C
      IF(IDAMG.EQ.1.OR.IDAMG.EQ.2) THEN
       IX=IX+1
       NDAMA=INT(PROPS(IX))
C
       DO IDAMA=1,NDAMA
        IA=IX-1+2*IDAMA
        IB=IA+1
        VDAMA(IDAMA,1)=PROPS(IA)
        VDAMA(IDAMA,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NDAMA
C
C**** FRACTURE ENERGY
C
       IX=IX+1
       NFRAC=INT(PROPS(IX))
C
       DO IFRAC=1,NFRAC
        IA=IX-1+2*IFRAC
        IB=IA+1
        VFRAC(IFRAC,1)=PROPS(IA)
        VFRAC(IFRAC,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NFRAC
      ENDIF
C
C**** COEFFICIENT 9
C
      IF(IDAMG.EQ.3) THEN
       IX=IX+1
       NCOE9=INT(PROPS(IX))
C
       DO ICOE9=1,NCOE9
        IA=IX-1+2*ICOE9
        IB=IA+1
        VCOE9(ICOE9,1)=PROPS(IA)
        VCOE9(ICOE9,2)=PROPS(IB)
       ENDDO
       IX=IX+2*NCOE9
      ENDIF
C
C**** GAMMA PARAMETERS
C
      IF(IDAMG.EQ.21) THEN
       IX=IX+1
       GAMMAP=PROPS(IX)
       IX=IX+1
       GAMMAM=PROPS(IX)
      ENDIF
C
C**** FREE ENERGY MODEL
C
      IX=IX+1
      IFREN=INT(PROPS(IX))
C
      RETURN
      END
