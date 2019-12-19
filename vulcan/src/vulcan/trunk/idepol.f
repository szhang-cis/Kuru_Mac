      SUBROUTINE IDEPOL(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THE POLYMER SOLID MODEL
C     (TEMPERATURE-DEPENDENT MODEL)
C
C***********************************************************************
C
C     NPOL1
C     VPOL1(1)=Constant 1                at VPOL1(2)=temperature
C
C     ....
C     ....                               ....
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
C**** FREE ENERGY MODEL
C
      IFREN=INT(PROPS(61))
      IFREM=INT(PROPS(58))
C
C**** CONSTANT C1
C
      I1=62
      IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .   IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58) THEN
       NPOL1=INT(PROPS(I1))
C
       DO IPOL1=1,NPOL1
        IA=I1-1+2*IPOL1
        IB=IA+1
        VPOL1(IPOL1,1)=PROPS(IA)
        VPOL1(IPOL1,2)=PROPS(IB)
       ENDDO
      ENDIF                  ! ifren.eq.51....
C
C**** CONSTANT C2
C
      IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .   IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58) THEN
       I2=I1+2*NPOL1+1
       NPOL2=INT(PROPS(I2))
C
       DO IPOL2=1,NPOL2
        IA=I2-1+2*IPOL2
        IB=IA+1
        VPOL2(IPOL2,1)=PROPS(IA)
        VPOL2(IPOL2,2)=PROPS(IB)
       ENDDO
       I3=I2
       NPOLX=NPOL2
      ENDIF                  ! ifren.eq.51....
C
C**** CONSTANT C3
C
      IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.
     .   IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58) THEN
       I3=I2+2*NPOL2+1
       NPOL3=INT(PROPS(I3))
C
       DO IPOL3=1,NPOL3
        IA=I3-1+2*IPOL3
        IB=IA+1
        VPOL3(IPOL3,1)=PROPS(IA)
        VPOL3(IPOL3,2)=PROPS(IB)
       ENDDO
       NPOLX=NPOL3
      ENDIF                  ! ifren.eq.51....
C
C**** EXPONENT 1
C
      IF(IFREN.EQ.53) THEN
       I1=I3+2*NPOL3+1
       NPOE1=INT(PROPS(I1))
C
       DO IPOE1=1,NPOE1
        IA=I1-1+2*IPOE1
        IB=IA+1
        VPOE1(IPOE1,1)=PROPS(IA)
        VPOE1(IPOE1,2)=PROPS(IB)
       ENDDO
C
C**** EXPONENT 2
C
       I2=I1+2*NPOE1+1
       NPOE2=INT(PROPS(I2))
C
       DO IPOE2=1,NPOE2
        IA=I2-1+2*IPOE2
        IB=IA+1
        VPOE2(IPOE2,1)=PROPS(IA)
        VPOE2(IPOE2,2)=PROPS(IB)
       ENDDO
C
C**** EXPONENT 3
C
       I3=I2+2*NPOE2+1
       NPOE3=INT(PROPS(I3))
C
       DO IPOE3=1,NPOE3
        IA=I3-1+2*IPOE3
        IB=IA+1
        VPOE3(IPOE3,1)=PROPS(IA)
        VPOE3(IPOE3,2)=PROPS(IB)
       ENDDO
       NPOLX=NPOE3
      ENDIF                  ! ifren.eq.53
C
C**** CONSTANT C4
C
      IF(IFREN.EQ.56.OR.IFREN.EQ.57) THEN
       I4=I3+2*NPOL3+1
       NPOL4=INT(PROPS(I4))
C
       DO IPOL4=1,NPOL4
        IA=I4-1+2*IPOL4
        IB=IA+1
        VPOL4(IPOL4,1)=PROPS(IA)
        VPOL4(IPOL4,2)=PROPS(IB)
       ENDDO
       NPOLX=NPOL4
       I3=I4
      ENDIF                  ! ifren.eq.56.or.ifren.eq.57
C
C**** CONSTANT C5
C
      IF(IFREN.EQ.56) THEN
       I5=I4+2*NPOL4+1
       NPOL5=INT(PROPS(I5))
C
       DO IPOL5=1,NPOL5
        IA=I5-1+2*IPOL5
        IB=IA+1
        VPOL5(IPOL5,1)=PROPS(IA)
        VPOL5(IPOL5,2)=PROPS(IB)
       ENDDO
       NPOLX=NPOL5
C
C**** CONSTANT C6
C
       I6=I5+2*NPOL5+1
       NPOL6=INT(PROPS(I6))
C
       DO IPOL6=1,NPOL6
        IA=I6-1+2*IPOL6
        IB=IA+1
        VPOL6(IPOL6,1)=PROPS(IA)
        VPOL6(IPOL6,2)=PROPS(IB)
       ENDDO
       NPOLX=NPOL6
C
C**** CONSTANT C7
C
       I7=I6+2*NPOL6+1
       NPOL7=INT(PROPS(I7))
C
       DO IPOL7=1,NPOL7
        IA=I7-1+2*IPOL7
        IB=IA+1
        VPOL7(IPOL7,1)=PROPS(IA)
        VPOL7(IPOL7,2)=PROPS(IB)
       ENDDO
       NPOLX=NPOL7
C
C**** CONSTANT C8
C
       I8=I7+2*NPOL7+1
       NPOL8=INT(PROPS(I8))
C
       DO IPOL8=1,NPOL8
        IA=I8-1+2*IPOL8
        IB=IA+1
        VPOL8(IPOL8,1)=PROPS(IA)
        VPOL8(IPOL8,2)=PROPS(IB)
       ENDDO
       NPOLX=NPOL8
       I3=I8
      ENDIF                  ! ifren.eq.56
C
C**** VOLUMETRIC TERM MODEL
C
      IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .   IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58) THEN
       I4=I3+2*NPOLX+1
       NMODI=INT(PROPS(I4))
C
C**** PENALTY FOR INCOMPRESSIBILITY
C
       I4=I4+1
       NPENI=INT(PROPS(I4))
C
       DO IPENI=1,NPENI
        IA=I4-1+2*IPENI
        IB=IA+1
        VPENI(IPENI,1)=PROPS(IA)
        VPENI(IPENI,2)=PROPS(IB)
       ENDDO
C
       IX=I4
       NX=NPENI
      ENDIF                  ! ifren.eq.51....
C
C**** REFERENCE TEMPERATURE
C
      TEREF=PROPS(11)        ! see proppol.f
C
C**** THERMAL DILATATION
C
      I5=IX+2*NX+1
      NALPH=INT(PROPS(I5))
C
      I5=I5+1
      NALPU=INT(PROPS(I5))
      VALPU(1,2)=PROPS(I5)
C
      DO IALPH=1,NALPH
       IA=I5-1+2*IALPH
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
         VALPU(IALPU,1)=0.5D0*VAUXX+VALPH(1,1)*(VALPH(1,2)-TEREF)/
     .                                     (VALPH(IALPU,2)-TEREF)
        ENDIF
       ENDDO
      ENDIF                             ! nalpu.eq.1
C
C**** THERMAL HARDENING FUNCTION
C
      I6=I5+2*NALPH+1
      NCCER=INT(PROPS(I6))
C
      DO ICCER=1,NCCER
       IA=I6-1+2*ICCER
       IB=IA+1
       VCCER(ICCER,1)=PROPS(IA)
       VCCER(ICCER,2)=PROPS(IB)
      ENDDO
C
C**** VISCOELASTIC PARAMETERS
C
      I7=I6+2*NCCER+1
      NCHAI=INT(PROPS(I7))
      NCCERX=0
      IF(NCHAI.GT.0) THEN
       DO IN=1,NCHAI
        I7A=I7+1+NCCERX
        NPOL9=INT(PROPS(I7A))
C
        DO IPOL9=1,NPOL9
         IA=I7A-1+2*IPOL9
         IB=IA+1
         VPOL9(IN,IPOL9,1)=PROPS(IA)
         VPOL9(IN,IPOL9,2)=PROPS(IB)
        ENDDO
C
        I7B=I7A+2*NPOL9+1
        NPOL10=INT(PROPS(I7B))
C      
        DO IPOL10=1,NPOL10
         IA=I7B-1+2*IPOL10
         IB=IA+1
         VPOL10(IN,IPOL10,1)=PROPS(IA)
         VPOL10(IN,IPOL10,2)=PROPS(IB)
        ENDDO
        NCCERX=NCCERX+2*NPOL9+1+2*NPOL10+1
       ENDDO
      ENDIF
C
C**** FIBER-REINFORCED PARAMETERS
C
      I8=I7+NCCERX+1
      NFIBN=INT(PROPS(I8))
      IF(NFIBN.GT.0) THEN
       I8A=I8+1
       NFIBM=INT(PROPS(I8A))
       IF(NFIBM.EQ.1.OR.NFIBM.EQ.2) THEN
C
C**** CONSTANT C1
C
        I8B=I8A+1
        NFIB1=INT(PROPS(I8B))
C
        DO IFIB1=1,NFIB1
         IA=I8B-1+2*IFIB1
         IB=IA+1
         VFIB1(IFIB1,1)=PROPS(IA)
         VFIB1(IFIB1,2)=PROPS(IB)
        ENDDO
C
C**** CONSTANT C2
C
        I8C=I8B+2*NFIB1+1
        NFIB2=INT(PROPS(I8C))
C
        DO IFIB2=1,NFIB2
         IA=I8C-1+2*IFIB2
         IB=IA+1
         VFIB2(IFIB2,1)=PROPS(IA)
         VFIB2(IFIB2,2)=PROPS(IB)
        ENDDO
C
C**** CONSTANT C3
C
        I8D=I8C+2*NFIB2+1
        NFIB3=INT(PROPS(I8D))
C
        DO IFIB3=1,NFIB3
         IA=I8D-1+2*IFIB3
         IB=IA+1
         VFIB3(IFIB3,1)=PROPS(IA)
         VFIB3(IFIB3,2)=PROPS(IB)
        ENDDO
C
C**** CONSTANT C4
C
        I8E=I8D+2*NFIB3+1
        NFIB4=INT(PROPS(I8E))
C
        DO IFIB4=1,NFIB4
         IA=I8E-1+2*IFIB4
         IB=IA+1
         VFIB4(IFIB4,1)=PROPS(IA)
         VFIB4(IFIB4,2)=PROPS(IB)
        ENDDO
        NCCERX=5+2*(NFIB1+NFIB2+NFIB3+NFIB4)
       ENDIF
C
       IF(NFIBM.EQ.2) THEN
C
C**** CONSTANT C5
C
        I8F=I8E+2*NFIB4+1
        NFIB5=INT(PROPS(I8F))
C
        DO IFIB5=1,NFIB5
         IA=I8F-1+2*IFIB5
         IB=IA+1
         VFIB5(IFIB5,1)=PROPS(IA)
         VFIB5(IFIB5,2)=PROPS(IB)
        ENDDO
        NCCERX=5+2*(NFIB1+NFIB2+NFIB3+NFIB4+NFIB5)
       ENDIF
C
       IF(NFIBM.EQ.3) THEN

       ENDIF
      ELSE
       NCCERX=0
      ENDIF
C
      I9=I8+NCCERX+1
      ISOTT=INT(PROPS(I9))
C
      I10=I9+1
      IKINE=INT(PROPS(I10))
C
      I11=I10+1
      NPLATM=INT(PROPS(I11))
C
      I12=I11+1
      ISHRI=INT(PROPS(I12))
C
      I13=I12+1
      IDAMG=INT(PROPS(I13))
      IF(IDAMG.EQ.1) THEN
C
C**** STRAIN ENERGY AT INITIAL DAMAGE (E_MIN)
C
       I14=I13+1
       NEMIN=INT(PROPS(I14))
C
       DO IEMIN=1,NEMIN
        IA=I14-1+2*IEMIN
        IB=IA+1
        VEMIN(IEMIN,1)=PROPS(IA)
        VEMIN(IEMIN,2)=PROPS(IB)
       ENDDO
C
C**** STRAIN ENERGY AT TOTAL DAMAGE (E_MAX)
C
       I15=I14+2*NEMIN+1
       NEMAX=INT(PROPS(I15))
C
       DO IEMAX=1,NEMAX
        IA=I15-1+2*IEMAX
        IB=IA+1
        VEMAX(IEMAX,1)=PROPS(IA)
        VEMAX(IEMAX,2)=PROPS(IB)
       ENDDO
C
C**** ETA PARAMETER
C
       I16=I15+2*NEMAX+1
       NEETA=INT(PROPS(I16))
C
       DO IEETA=1,NEETA
        IA=I16-1+2*IEETA
        IB=IA+1
        VEETA(IEETA,1)=PROPS(IA)
        VEETA(IEETA,2)=PROPS(IB)
       ENDDO
      ENDIF
C
      RETURN
      END
