      SUBROUTINE IDEOPT(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THE ORTHOTROPIC, PLASTIC &
C     TEMPERATURE-DEPENDENT MODEL
C
C***********************************************************************
C
C     NYOUN
C     VYOUN(1) =Young's Modulus, component xx   at VYOUN(2) =temperature
C     VYOUNF(1)=Young's modulus, component yy   at VYOUNF(2)=temperature
C     VYOUNA(1)=Young's modulus, component zz   at VYOUNA(2)=temperature
C
C     VPOIS(1) =Poisson's ratio, component xy   at VPOIS(2) =temperature
C     VPOISF(1)=Poisson's ratio, component xz   at VPOISF(2)=temperature
C     VPOISA(1)=Poisson's ratio, component yz   at VPOISA(2)=temperature
C
C     VSHEA(1) =Shear modulus, component xy     at VSHEA(2) =temperature
C     VSHEAF(1)=Shear modulus, component xz     at VSHEAF(2)=temperature
C     VSHEAA(1)=Shear modulus, component yz     at VSHEAA(2)=temperature
C
C     NALPU=indicates secant or tangent Dilatation Coefficient
C
C     NALPH
C     VALPH(1) =Thermal dilation, component xx  at VALPH(2) =temperature
C     VALPHF(1)=Thermal dilation, component xx  at VALPHF(2)=temperature
C     VALPHA(1)=Thermal dilation, component xx  at VALPHA(2)=temperature
C
C     NCCER
C     VCCER(1)=Hardening Function               at VCCER(2)=temperature
C
C
C     ISOTT=indicator of isotropic hardening model
C
C     NCCOE
C     VCCOE(1)=Hardening Coefficient            at VCCOE(2)=temperature
C
C     NCCOB
C     VCCOB(1)=Coefficient b                    at VCCOE(2)=temperature
C
C     NCCOQ
C     VCCOQ(1)=Coefficient Q                    at VCCOE(2)=temperature
C
C     IPORO=indicator of porosity model (Gurson)
C
C     NPCOE
C     VPCOE(1)=Porosity Coefficient             at VPCOE(2)=temperature
C
C     NPCO1
C     VPCO1(1)=Q1 Coefficient                   at VPCOE(2)=temperature
C
C     NPCO2
C     VPCO2(1)=Q2 Coefficient                   at VPCOE(2)=temperature
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
C**** YOUNG MODULI
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
      I1=I1+2*NYOUN+1
      NYOUNF=INT(PROPS(I1))
C
      DO IYOUN=1,NYOUNF
       IA=I1-1+2*IYOUN
       IB=IA+1
       VYOUNF(IYOUN,1)=PROPS(IA)
       VYOUNF(IYOUN,2)=PROPS(IB)
      ENDDO
C
      I1=I1+2*NYOUNF+1
      NYOUNA=INT(PROPS(I1))
C
      DO IYOUN=1,NYOUNA
       IA=I1-1+2*IYOUN
       IB=IA+1
       VYOUNA(IYOUN,1)=PROPS(IA)
       VYOUNA(IYOUN,2)=PROPS(IB)
      ENDDO
C
C**** POISSON RATII
C
      I2=I1+2*NYOUNA+1
      NPOIS=INT(PROPS(I2))
C
      DO IPOIS=1,NPOIS
       IA=I2-1+2*IPOIS
       IB=IA+1
       VPOIS(IPOIS,1)=PROPS(IA)
       VPOIS(IPOIS,2)=PROPS(IB)
      ENDDO
C
      I2=I2+2*NPOIS+1
      NPOISF=INT(PROPS(I2))
C
      DO IPOIS=1,NPOISF
       IA=I2-1+2*IPOIS
       IB=IA+1
       VPOISF(IPOIS,1)=PROPS(IA)
       VPOISF(IPOIS,2)=PROPS(IB)
      ENDDO
C
      I2=I2+2*NPOISF+1
      NPOISA=INT(PROPS(I2))
C
      DO IPOIS=1,NPOISA
       IA=I2-1+2*IPOIS
       IB=IA+1
       VPOISA(IPOIS,1)=PROPS(IA)
       VPOISA(IPOIS,2)=PROPS(IB)
      ENDDO
C
C**** SHEAR MODULI
C
      I2A=I2+2*NPOISA+1
      NSHEA=INT(PROPS(I2A))
C
      DO ISHEA=1,NSHEA
       IA=I2A-1+2*ISHEA
       IB=IA+1
       VSHEA(ISHEA,1)=PROPS(IA)
       VSHEA(ISHEA,2)=PROPS(IB)
      ENDDO
C
      NSHEAX=NSHEA
      IF(NDIME.EQ.3) THEN
       I2A=I2A+2*NSHEA+1
       NSHEAF=INT(PROPS(I2A))
C
       DO ISHEA=1,NSHEAF
        IA=I2A-1+2*ISHEA
        IB=IA+1
        VSHEAF(ISHEA,1)=PROPS(IA)
        VSHEAF(ISHEA,2)=PROPS(IB)
       ENDDO
C
       I2A=I2A+2*NSHEAF+1
       NSHEAA=INT(PROPS(I2A))
C
       DO ISHEA=1,NSHEAA
        IA=I2A-1+2*ISHEA
        IB=IA+1
        VSHEAA(ISHEA,1)=PROPS(IA)
        VSHEAA(ISHEA,2)=PROPS(IB)
       ENDDO
       NSHEAX=NSHEAA
      ENDIF
C
C**** REFERENCE TEMPERATURE
C
      TEREF=PROPS(11)        ! see propopt.f
C
C**** THERMAL DILATATION
C
      I3=I2A+2*NSHEAX+1
      NALPU=INT(PROPS(I3))
      VALPU(1,2)=PROPS(I3)
C
      I3=I3+1
      NALPH=INT(PROPS(I3))
C
      DO IALPH=1,NALPH
       IA=I3-1+2*IALPH
       IB=IA+1
       VALPH(IALPH,1)=PROPS(IA)
       VALPH(IALPH,2)=PROPS(IB)
      ENDDO
C
      I3=I3+2*NALPH+1
      NALPHF=INT(PROPS(I3))
C
      DO IALPH=1,NALPHF
       IA=I3-1+2*IALPH
       IB=IA+1
       VALPHF(IALPH,1)=PROPS(IA)
       VALPHF(IALPH,2)=PROPS(IB)
      ENDDO
C
      I3=I3+2*NALPHF+1
      NALPHA=INT(PROPS(I3))
C
      DO IALPH=1,NALPHA
       IA=I3-1+2*IALPH
       IB=IA+1
       VALPHA(IALPH,1)=PROPS(IA)
       VALPHA(IALPH,2)=PROPS(IB)
      ENDDO
C
      IF(NALPU.EQ.1) THEN                ! tangent alpha
c     not implemented yet (alpha^s for each component must be computed)
       call runend('ideopt: nalpu=1 not implemented')
c      VALPU(1,1)=VALPH(1,1)             ! computes alpha^s
c      DO IALPU=2,NALPH
c       VALPU(IALPU,1)=VALPH(IALPU,1)    ! alpha^s=alpha
c       VAUXD=VALPH(IALPU,2)-TEREF
c       IF(VAUXD.GT.0) THEN
c        VAUXX=0.0D+00
c        DO JALPU=2,IALPU
c         I1=JALPU-1
c         I2=JALPU
c         VAUXX=VAUXX+(VALPH(I2,1)+VALPH(I1,1))*
c    .                (VALPH(I2,2)-VALPH(I1,2))/
c    .                (VALPH(IALPU,2)-TEREF)
c        ENDDO
c        VALPU(IALPU,1)=0.5*VAUXX+VALPH(1,1)*(VALPH(1,2)-TEREF)/
c    .                                    (VALPH(IALPU,2)-TEREF)
c       ENDIF
c      ENDDO
      ENDIF                             ! nalpu.eq.1
C
C**** THERMAL HARDENING FUNCTION
C
      I4=I3+2*NALPHA+1
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
      I6=I5
      NCCOE=0
      NCCOEX=NCCOE
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
     .   ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.
     .   ISOTT.EQ.13) THEN
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
      I8=I7
      NKCOE=0
      NKCOEX=NKCOE
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
C**** PARAMETERS OF HILL 48 YIELD FUNCTION
C
      IF(NCRIT.EQ.42) THEN
       NC42T=INT(PROPS(45))
C
       IF(NC42T.EQ.1) THEN          ! GENERAL MODEL
C
C**** HILL COEFFICIENTS (F to N)
C
        I9=I8+2*NKCOEX+1
        NLAN1(1)=INT(PROPS(I9))     ! F
C
        DO ILANK=1,NLAN1(1)
         IA=I9-1+2*ILANK
         IB=IA+1
         VLAN1(ILANK,1)=PROPS(IA)
         VLAN1(ILANK,2)=PROPS(IB)
        ENDDO
C
        I9=I9+2*NLAN1(1)+1
        NLAN2(1)=INT(PROPS(I9))     ! G
C
        DO ILANK=1,NLAN2(1)
         IA=I9-1+2*ILANK
         IB=IA+1
         VLAN2(ILANK,1)=PROPS(IA)
         VLAN2(ILANK,2)=PROPS(IB)
        ENDDO
C
        I9=I9+2*NLAN2(1)+1
        NLAN3(1)=INT(PROPS(I9))     ! H
C
        DO ILANK=1,NLAN3(1)
         IA=I9-1+2*ILANK
         IB=IA+1
         VLAN3(ILANK,1)=PROPS(IA)
         VLAN3(ILANK,2)=PROPS(IB)
        ENDDO
C
        I9=I9+2*NLAN3(1)+1
        NLAN4(1)=INT(PROPS(I9))     ! N
C
        DO ILANK=1,NLAN4(1)
         IA=I9-1+2*ILANK
         IB=IA+1
         VLAN4(ILANK,1)=PROPS(IA)
         VLAN4(ILANK,2)=PROPS(IB)
        ENDDO
C
        I8=I9
        NKCOEX=NLAN4(1)
C
        IF(NTYPE.EQ.4) THEN
         I9=I9+2*NLAN4(1)+1
         NLAN5(1)=INT(PROPS(I9))    ! M
C
         DO ILANK=1,NLAN5(1)
          IA=I9-1+2*ILANK
          IB=IA+1
          VLAN5(ILANK,1)=PROPS(IA)
          VLAN5(ILANK,2)=PROPS(IB)
         ENDDO
C
         I9=I9+2*NLAN5(1)+1
         NLAN6(1)=INT(PROPS(I9))    ! L
C
         DO ILANK=1,NLAN6(1)
          IA=I9-1+2*ILANK
          IB=IA+1
          VLAN6(ILANK,1)=PROPS(IA)
          VLAN6(ILANK,2)=PROPS(IB)
         ENDDO
C
         I8=I9
         NKCOEX=NLAN6(1)
        ENDIF
       ENDIF
C
       IF(NC42T.EQ.2.OR.            ! PLANE STRESS & S-BASED MODEL
     .    NC42T.EQ.3) THEN          ! PLANE STRESS & R-BASED MODEL
C
C**** STRENGTH COEFFICIENTS (S_45, S_90 & S_BIAXIAL) OR
C     LANKFORD COEFFICIENTS (R0, R45 & R90)
C
        I9=I8+2*NKCOEX+1
        NLAN1(1)=INT(PROPS(I9))     ! S45 OR R0
C
        DO ILANK=1,NLAN1(1)
         IA=I9-1+2*ILANK
         IB=IA+1
         VLAN1(ILANK,1)=PROPS(IA)
         VLAN1(ILANK,2)=PROPS(IB)
        ENDDO
C
        I8=I9
        NKCOEX=NLAN1(1)
C
        IF(NTYPE.NE.3) THEN
         I9=I9+2*NLAN1(1)+1
         NLAN2(1)=INT(PROPS(I9))    ! S90 OR R45
C
         DO ILANK=1,NLAN2(1)
          IA=I9-1+2*ILANK
          IB=IA+1
          VLAN2(ILANK,1)=PROPS(IA)
          VLAN2(ILANK,2)=PROPS(IB)
         ENDDO
C
         I9=I9+2*NLAN2(1)+1
         NLAN3(1)=INT(PROPS(I9))    ! SB OR R90
C
         DO ILANK=1,NLAN3(1)
          IA=I9-1+2*ILANK
          IB=IA+1
          VLAN3(ILANK,1)=PROPS(IA)
          VLAN3(ILANK,2)=PROPS(IB)
         ENDDO
C
         I8=I9
         NKCOEX=NLAN3(1)
        ENDIF
       ENDIF              ! nc42t.eq.2,3
      ENDIF               ! ncrit.eq.42
C
C**** PARAMETERS OF HILL 48 FLOW POTENTIAL
C
      IF(NCRIP.EQ.42) THEN
       NC42P=INT(PROPS(46))
C
       IF(NC42P.EQ.1) THEN          ! GENERAL MODEL
        IF(NASNA.EQ.0) THEN         ! associate plasticity
         DO ILANK=1,NLAN1(1)        ! F
          VLAN1(ILANK,3)=VLAN1(ILANK,1)
          VLAN1(ILANK,4)=VLAN1(ILANK,2)
         ENDDO
         DO ILANK=1,NLAN2(1)        ! G
          VLAN2(ILANK,3)=VLAN2(ILANK,1)
          VLAN2(ILANK,4)=VLAN2(ILANK,2)
         ENDDO
         DO ILANK=1,NLAN3(1)        ! H
          VLAN3(ILANK,3)=VLAN3(ILANK,1)
          VLAN3(ILANK,4)=VLAN3(ILANK,2)
         ENDDO
         DO ILANK=1,NLAN4(1)        ! N
          VLAN4(ILANK,3)=VLAN4(ILANK,1)
          VLAN4(ILANK,4)=VLAN4(ILANK,2)
         ENDDO
         IF(NTYPE.EQ.4) THEN
          DO ILANK=1,NLAN4(1)       ! M
           VLAN5(ILANK,3)=VLAN5(ILANK,1)
           VLAN5(ILANK,4)=VLAN5(ILANK,2)
          ENDDO
          DO ILANK=1,NLAN6(1)       ! L
           VLAN6(ILANK,3)=VLAN6(ILANK,1)
           VLAN6(ILANK,4)=VLAN6(ILANK,2)
          ENDDO
         ENDIF
        ELSE                        ! non-associate plasticity
C
C**** HILL COEFFICIENTS (F to N)
C
         I9=I8+2*NKCOEX+1
         NLAN1(2)=INT(PROPS(I9))    ! F
C
         DO ILANK=1,NLAN1(2)
          IA=I9-1+2*ILANK
          IB=IA+1
          VLAN1(ILANK,3)=PROPS(IA)
          VLAN1(ILANK,4)=PROPS(IB)
         ENDDO
C
         I9=I9+2*NLAN1(2)+1
         NLAN2(2)=INT(PROPS(I9))    ! G
C
         DO ILANK=1,NLAN2(2)
          IA=I9-1+2*ILANK
          IB=IA+1
          VLAN2(ILANK,3)=PROPS(IA)
          VLAN2(ILANK,4)=PROPS(IB)
         ENDDO
C
         I9=I9+2*NLAN2(2)+1
         NLAN3(2)=INT(PROPS(I9))    ! H
C
         DO ILANK=1,NLAN3(2)
          IA=I9-1+2*ILANK
          IB=IA+1
          VLAN3(ILANK,3)=PROPS(IA)
          VLAN3(ILANK,4)=PROPS(IB)
         ENDDO
C
         I9=I9+2*NLAN3(2)+1
         NLAN4(2)=INT(PROPS(I9))    ! N
C
         DO ILANK=1,NLAN4(2)
          IA=I9-1+2*ILANK
          IB=IA+1
          VLAN4(ILANK,3)=PROPS(IA)
          VLAN4(ILANK,4)=PROPS(IB)
         ENDDO
C
         I8=I9
         NKCOEX=NLAN4(2)
C
         IF(NTYPE.EQ.4) THEN
          I9=I9+2*NLAN4(2)+1
          NLAN5(2)=INT(PROPS(I9))   ! M
C
          DO ILANK=1,NLAN5(2)
           IA=I9-1+2*ILANK
           IB=IA+1
           VLAN5(ILANK,3)=PROPS(IA)
           VLAN5(ILANK,4)=PROPS(IB)
          ENDDO
C
          I9=I9+2*NLAN5(2)+1
          NLAN6(2)=INT(PROPS(I9))   ! L
C
          DO ILANK=1,NLAN6(2)
           IA=I9-1+2*ILANK
           IB=IA+1
           VLAN6(ILANK,3)=PROPS(IA)
           VLAN6(ILANK,4)=PROPS(IB)
          ENDDO
C
          I8=I9
          NKCOEX=NLAN6(2)
         ENDIF
        ENDIF             ! nasna.eq.0
       ENDIF              ! nc42p.eq.1
C
       IF(NC42P.EQ.3) THEN          ! PLANE STRESS & R-BASED MODEL
C
C**** LANKFORD COEFFICIENTS (R0, R45 & R90)
C
        IF(NASNA.EQ.0) THEN         ! associate plasticity
         DO ILANK=1,NLAN1(1)        ! R0
          VLAN1(ILANK,3)=VLAN1(ILANK,1)
          VLAN1(ILANK,4)=VLAN1(ILANK,2)
         ENDDO
         DO ILANK=1,NLAN2(1)        ! R45
          VLAN2(ILANK,3)=VLAN2(ILANK,1)
          VLAN2(ILANK,4)=VLAN2(ILANK,2)
         ENDDO
         DO ILANK=1,NLAN3(1)        ! R90
          VLAN3(ILANK,3)=VLAN3(ILANK,1)
          VLAN3(ILANK,4)=VLAN3(ILANK,2)
         ENDDO
        ELSE                        ! non-associate plasticity
C
C**** LANKFORD COEFFICIENTS
C
         I9=I8+2*NKCOEX+1
         NLAN1(2)=INT(PROPS(I9))    ! R0
C
         DO ILANK=1,NLAN1(2)
          IA=I9-1+2*ILANK
          IB=IA+1
          VLAN1(ILANK,3)=PROPS(IA)
          VLAN1(ILANK,4)=PROPS(IB)
         ENDDO
C
         I8=I9
         NKCOEX=NLAN1(2)
C
         IF(NTYPE.NE.3) THEN
          I9=I9+2*NLAN1(2)+1
          NLAN2(2)=INT(PROPS(I9))   ! R45
C
          DO ILANK=1,NLAN2(2)
           IA=I9-1+2*ILANK
           IB=IA+1
           VLAN2(ILANK,3)=PROPS(IA)
           VLAN2(ILANK,4)=PROPS(IB)
          ENDDO
C
          I9=I9+2*NLAN2(2)+1
          NLAN3(2)=INT(PROPS(I9))   ! R90
C
          DO ILANK=1,NLAN3(2)
           IA=I9-1+2*ILANK
           IB=IA+1
           VLAN3(ILANK,3)=PROPS(IA)
           VLAN3(ILANK,4)=PROPS(IB)
          ENDDO
C
          I8=I9
          NKCOEX=NLAN3(2)
         ENDIF
        ENDIF             ! nasna.eq.0
       ENDIF              ! nc42p.eq.3
      ENDIF               ! ncrip.eq.42
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
        NTEMRO=4
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
