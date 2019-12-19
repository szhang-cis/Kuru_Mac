      SUBROUTINE IDEIVT1(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THE ISOTROPIC,
C     VISCOPLASTIC & TEMPERATURE-DEPENDENT MODEL FOR SG CAST IRON
C
C***********************************************************************
C
C     NYOUNF
C     VYOUNF(1)=Young's Modulus of ferrite      at VYOUNF(2)=temperature
C
C     NYOUNA
C     VYOUNA(1)=Young's Modulus of austenite    at VYOUNA(2)=temperature
C
C     NPOISF
C     VPOISF(1)=Poisson's Ratio of ferrite      at VPOISF(2)=temperature
C
C     NPOISA
C     VPOISA(1)=Poisson's Ratio of austenite    at VPOISA(2)=temperature
C
C     NALPHF
C     VALPHF(1)=Thermal Dilatation Coefficient  at VALPHF(2)=temperature
C
C     NALPUF=indicates secant or tangent Dilatation Coefficient
C
C     NALPHA
C     VALPHA(1)=Thermal Dilatation Coefficient  at VALPHA(2)=temperature
C
C     NALPUA=indicates secant or tangent Dilatation Coefficient
C
C     NCCER
C     VCCER(1)=Thermal Hardening Function     at VCCER(2)=temperature
C
C     IVIFL=indicator of viscous law model
C
C     NVISC
C     VVISC(1)=Viscosity Function             at VVISC(2)=temperature
C
C     NZABAF
C     VZABAF(1)=Coefficient A                 at VZABAF(2)=temperature
C
C     NZABBF
C     VZABBF(1)=Coefficient B                 at VZABBF(2)=temperature
C
C     NZABCF
C     VZABCF(1)=Coefficient C                 at VZABCF(2)=temperature
C
C     NZABAA
C     VZABAA(1)=Coefficient A                 at VZABAA(2)=temperature
C
C     NZABBA
C     VZABBA(1)=Coefficient B                 at VZABBA(2)=temperature
C
C     NZABCA
C     VZABCA(1)=Coefficient C                 at VZABCA(2)=temperature
C
C     NEXPOF
C     VEXPOF(1)=Exponent Function             at VEXPOF(2)=temperature
C
C     NEXPOA
C     VEXPOA(1)=Exponent Function             at VEXPOA(2)=temperature
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
C**** YOUNG MODULUS
C
      I1F=61
      NYOUNF=INT(PROPS(I1F))
C
      DO IYOUN=1,NYOUNF
       IA=I1F-1+2*IYOUN
       IB=IA+1
       VYOUNF(IYOUN,1)=PROPS(IA)
       VYOUNF(IYOUN,2)=PROPS(IB)
      ENDDO
C
      I1A=I1F+2*NYOUNF+1
      NYOUNA=INT(PROPS(I1A))
C
      DO IYOUN=1,NYOUNA
       IA=I1A-1+2*IYOUN
       IB=IA+1
       VYOUNA(IYOUN,1)=PROPS(IA)
       VYOUNA(IYOUN,2)=PROPS(IB)
      ENDDO
C
C**** POISSON RATIO
C
      I2F=I1A+2*NYOUNA+1
      NPOISF=INT(PROPS(I2F))
C
      DO IPOIS=1,NPOISF
       IA=I2F-1+2*IPOIS
       IB=IA+1
       VPOISF(IPOIS,1)=PROPS(IA)
       VPOISF(IPOIS,2)=PROPS(IB)
      ENDDO
C
      I2A=I2F+2*NPOISF+1
      NPOISA=INT(PROPS(I2A))
C
      DO IPOIS=1,NPOISA
       IA=I2A-1+2*IPOIS
       IB=IA+1
       VPOISA(IPOIS,1)=PROPS(IA)
       VPOISA(IPOIS,2)=PROPS(IB)
      ENDDO
C
C**** REFERENCE TEMPERATURE
C
      TEREF=PROPS(11)        ! see propivt1.f
C
C**** THERMAL DILATATION
C
      I3F=I2A+2*NPOISA+1
      NALPHF=INT(PROPS(I3F))
C
      I3F=I3F+1
      NALPUF=INT(PROPS(I3F))
      VALPUF(1,2)=PROPS(I3F)
C
      DO IALPH=1,NALPHF
       IA=I3F-1+2*IALPH
       IB=IA+1
       VALPHF(IALPH,1)=PROPS(IA)
       VALPHF(IALPH,2)=PROPS(IB)
      ENDDO
C
      IF(NALPUF.EQ.1) THEN                ! tangent alpha
       VALPUF(1,1)=VALPHF(1,1)            ! computes alpha^s
       DO IALPU=2,NALPHF
        VALPUF(IALPU,1)=VALPHF(IALPU,1)   ! alpha^s=alpha
        VAUXD=VALPHF(IALPU,2)-TEREF
        IF(VAUXD.GT.0) THEN
         VAUXX=0.0D+00
         DO JALPU=2,IALPU
          I1=JALPU-1
          I2=JALPU
          VAUXX=VAUXX+(VALPHF(I2,1)+VALPHF(I1,1))*
     .                (VALPHF(I2,2)-VALPHF(I1,2))/
     .                (VALPHF(IALPU,2)-TEREF)
         ENDDO
         VALPUF(IALPU,1)=0.5*VAUXX+VALPHF(1,1)*(VALPHF(1,2)-TEREF)/
     .                                      (VALPHF(IALPU,2)-TEREF)
        ENDIF
       ENDDO
      ENDIF                             ! nalpuf.eq.1
C
      I3A=I3F+2*NALPHF+1
      NALPHA=INT(PROPS(I3A))
C
      I3A=I3A+1
      NALPUA=INT(PROPS(I3A))
      VALPUA(1,2)=PROPS(I3A)
C
      DO IALPH=1,NALPHA
       IA=I3A-1+2*IALPH
       IB=IA+1
       VALPHA(IALPH,1)=PROPS(IA)
       VALPHA(IALPH,2)=PROPS(IB)
      ENDDO
C
      IF(NALPUA.EQ.1) THEN                ! tangent alpha
       VALPUA(1,1)=VALPHA(1,1)            ! computes alpha^s
       DO IALPU=2,NALPHA
        VALPUA(IALPU,1)=VALPHA(IALPU,1)   ! alpha^s=alpha
        VAUXD=VALPHA(IALPU,2)-TEREF
        IF(VAUXD.GT.0) THEN
         VAUXX=0.0D+00
         DO JALPU=2,IALPU
          I1=JALPU-1
          I2=JALPU
          VAUXX=VAUXX+(VALPHA(I2,1)+VALPHA(I1,1))*
     .                (VALPHA(I2,2)-VALPHA(I1,2))/
     .                (VALPHA(IALPU,2)-TEREF)
         ENDDO
         VALPUA(IALPU,1)=0.5*VAUXX+VALPHA(1,1)*(VALPHA(1,2)-TEREF)/
     .                                      (VALPHA(IALPU,2)-TEREF)
        ENDIF
       ENDDO
      ENDIF                             ! nalpua.eq.1
C
C**** THERMAL HARDENING FUNCTION
C
      I4F=I3A+2*NALPHA+1
      NCCERF=INT(PROPS(I4F))
C
      DO ICCER=1,NCCERF
       IA=I4F-1+2*ICCER
       IB=IA+1
       VCCERF(ICCER,1)=PROPS(IA)
       VCCERF(ICCER,2)=PROPS(IB)
      ENDDO
C
      I4A=I4F+2*NCCERF+1
      NCCERA=INT(PROPS(I4A))
C
      DO ICCER=1,NCCERA
       IA=I4A-1+2*ICCER
       IB=IA+1
       VCCERA(ICCER,1)=PROPS(IA)
       VCCERA(ICCER,2)=PROPS(IB)
      ENDDO
C
C**** FLOW PARAMETERS (VISCOUS LAW MODEL)
C
      I5=I4A+2*NCCERA+1
      IVIFL=INT(PROPS(I5))
C
C**** VISCOSITY FUNCTION
C
      IF(IVIFL.EQ.1.OR.IVIFL.EQ.2.OR.IVIFL.EQ.3) THEN
       I6F=I5+1
       NVISCF=INT(PROPS(I6F))
C
       DO IVISC=1,NVISCF
        IA=I6F-1+2*IVISC
        IB=IA+1
        VVISCF(IVISC,1)=PROPS(IA)
        VVISCF(IVISC,2)=PROPS(IB)
       ENDDO
C
       I6A=I6F+2*NVISCF+1
       NVISCA=INT(PROPS(I6A))
C
       DO IVISC=1,NVISCA
        IA=I6A-1+2*IVISC
        IB=IA+1
        VVISCA(IVISC,1)=PROPS(IA)
        VVISCA(IVISC,2)=PROPS(IB)
       ENDDO
C
      ENDIF          ! ivifl.eq.1.or.ivifl.eq.2.or.ivifl.eq.3
C
C**** ZABARAS' COEFFICIENTS (A, B & C)
C
      IF(IVIFL.EQ.4)
     . CALL RUNEND('IVIFL=4 not implemented - see ideivt.f')
C
C**** EXPONENT FUNCTION
C
      I7F=I6A+2*NVISCA+1
      NEXPOF=INT(PROPS(I7F))
C
      DO IEXPO=1,NEXPOF
       IA=I7F-1+2*IEXPO
       IB=IA+1
       VEXPOF(IEXPO,1)=PROPS(IA)
       VEXPOF(IEXPO,2)=PROPS(IB)
      ENDDO
C
      I7A=I7F+2*NEXPOF+1
      NEXPOA=INT(PROPS(I7A))
C
      DO IEXPO=1,NEXPOA
       IA=I7A-1+2*IEXPO
       IB=IA+1
       VEXPOA(IEXPO,1)=PROPS(IA)
       VEXPOA(IEXPO,2)=PROPS(IB)
      ENDDO
C
C**** ISOTROPIC HARDENING INDICATOR
C
      I8=I7A+2*NEXPOA+1
      ISOTT=INT(PROPS(I8))
c     IF(ISOTT.EQ.0) THEN
       I9=I8
       NCCOE=0
       NCCOEX=NCCOE
c     ENDIF
C
C**** HARDENING COEFFICIENT
C
      IF(ISOTT.EQ.1.OR.ISOTT.EQ.2.OR.ISOTT.EQ.4.OR.ISOTT.EQ.5.OR.
     .   ISOTT.EQ.8.OR.ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
       I9F=I8+1
       NCCOEF=INT(PROPS(I9F))
C
       DO ICCOE=1,NCCOEF
        IA=I9F-1+2*ICCOE
        IB=IA+1
        VCCOEF(ICCOE,1)=PROPS(IA)
        VCCOEF(ICCOE,2)=PROPS(IB)
       ENDDO
C
       I9A=I9F+2*NCCOEF+1
       NCCOEA=INT(PROPS(I9A))
C
       DO ICCOE=1,NCCOEA
        IA=I9A-1+2*ICCOE
        IB=IA+1
        VCCOEA(ICCOE,1)=PROPS(IA)
        VCCOEA(ICCOE,2)=PROPS(IB)
       ENDDO
C
       I9=I9A
       NCCOEX=NCCOEA
      ENDIF
C
C**** B COEFFICIENT
C
      IF(ISOTT.EQ.3.OR.ISOTT.EQ.6.OR.ISOTT.EQ.7.OR.ISOTT.EQ.8.OR.
     .   ISOTT.EQ.10.OR.ISOTT.EQ.11) THEN
c      I9F=I8+1
       I9F=I9+2*NCCOEX+1
       NCCOBF=INT(PROPS(I9F))
C
       DO ICCOB=1,NCCOBF
        IA=I9F-1+2*ICCOB
        IB=IA+1
        VCCOBF(ICCOB,1)=PROPS(IA)
        VCCOBF(ICCOB,2)=PROPS(IB)
       ENDDO
C
       I9A=I9F+2*NCCOBF+1
       NCCOBA=INT(PROPS(I9A))
C
       DO ICCOB=1,NCCOBA
        IA=I9A-1+2*ICCOB
        IB=IA+1
        VCCOBA(ICCOB,1)=PROPS(IA)
        VCCOBA(ICCOB,2)=PROPS(IB)
       ENDDO
C
C**** Q COEFFICIENT
C
       I9AF=I9A+2*NCCOBA+1
       NCCOQF=INT(PROPS(I9AF))
C
       DO ICCOQ=1,NCCOQF
        IA=I9AF-1+2*ICCOQ
        IB=IA+1
        VCCOQF(ICCOQ,1)=PROPS(IA)
        VCCOQF(ICCOQ,2)=PROPS(IB)
       ENDDO
C
       I9AA=I9AF+2*NCCOQF+1
       NCCOQA=INT(PROPS(I9AA))
C
       DO ICCOQ=1,NCCOQA
        IA=I9AA-1+2*ICCOQ
        IB=IA+1
        VCCOQA(ICCOQ,1)=PROPS(IA)
        VCCOQA(ICCOQ,2)=PROPS(IB)
       ENDDO
C
       I9=I9AA
       NCCOEX=NCCOQA
      ENDIF
C
C**** KINEMATIC HARDENING INDICATOR
C
      I7=I9+2*NCCOEX+1
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
       I8F=I7+1
       NKCOEF=INT(PROPS(I8F))
C
       DO IKCOE=1,NKCOEF
        IA=I8F-1+2*IKCOE
        IB=IA+1
        VKCOEF(IKCOE,1)=PROPS(IA)
        VKCOEF(IKCOE,2)=PROPS(IB)
       ENDDO
C
       I8A=I8F+2*NKCOEF+1
       NKCOEA=INT(PROPS(I8A))
C
       DO IKCOE=1,NKCOEA
        IA=I8A-1+2*IKCOE
        IB=IA+1
        VKCOEA(IKCOE,1)=PROPS(IA)
        VKCOEA(IKCOE,2)=PROPS(IB)
       ENDDO
C
       I8=I8A
       NKCOEX=NKCOEA
      ENDIF
C
C**** B COEFFICIENT
C
      IF(IKINE.EQ.3.OR.IKINE.EQ.4) THEN
       I8F=I7+1
       NKCOBF=INT(PROPS(I8F))
C
       DO IKCOB=1,NKCOBF
        IA=I8F-1+2*IKCOB
        IB=IA+1
        VKCOBF(IKCOB,1)=PROPS(IA)
        VKCOBF(IKCOB,2)=PROPS(IB)
       ENDDO
C
       I8A=I8F+2*NKCOBF+1
       NKCOBA=INT(PROPS(I8A))
C
       DO ICCOB=1,NKCOBA
        IA=I8A-1+2*ICCOB
        IB=IA+1
        VKCOBA(ICCOB,1)=PROPS(IA)
        VKCOBA(ICCOB,2)=PROPS(IB)
       ENDDO
C
C**** Q COEFFICIENT
C
       I8AF=I8A+2*NKCOBF+1
       NKCOQF=INT(PROPS(I8AF))
C
       DO IKCOQ=1,NKCOQF
        IA=I8AF-1+2*IKCOQ
        IB=IA+1
        VKCOQF(IKCOQ,1)=PROPS(IA)
        VKCOQF(IKCOQ,2)=PROPS(IB)
       ENDDO
C
       I8AA=I8AF+2*NKCOQF+1
       NKCOQA=INT(PROPS(I8AA))
C
       DO IKCOQ=1,NKCOQA
        IA=I8AA-1+2*IKCOQ
        IB=IA+1
        VKCOQA(IKCOQ,1)=PROPS(IA)
        VKCOQA(IKCOQ,2)=PROPS(IB)
       ENDDO
C
       I8=I8AA
       NKCOEX=NKCOQA
      ENDIF
C
      IVIFL=0
      IRECR=0
C
      IPORO=0                 ! no porosity
C
C**** PHASE-CHANGES
C
      I9=I8+2*NKCOEX+1
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
          VPLATM(IPLAT,10)=PROPS(I9D)             ! vplat(19) is busy!!
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
         CALL RUNEND('IDEIVT1: ERROR WITH IPCFOM=1 >> NOT POSSIBLE')
        ENDIF           ! ipcfom.eq.0
       ENDDO            ! do iplat=1,nplatm
      ENDIF             ! nplat.ne.0
C
C**** DAMAGE
C
      IDAMG=0           ! no damage
C
C**** FREE ENERGY MODEL
C
      I11=IX+1
      IFREN=INT(PROPS(I11))
C
      RETURN
      END
