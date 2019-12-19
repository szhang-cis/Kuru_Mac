      SUBROUTINE VISLAW(VISCO,FPARA,KVILA,CRIPS,TVELO,
     .                  GRAVE,TEMPE,PROAC,PSCON,NDIME,
     .                  FSLPG,IERROR,KINF3,KINF4)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE VISCOSITY AS A DEPENDENT VARIABLE
C
C
C     Notes:
C
C     VISCO: current viscosity
C
C     E:E=Eij*Eij=1/4*(dUi/dXj+dUj/dXi)*(dUi/dXj+dUj/dXi)
C             =1/4*[dUi/dXj*(dUi/dXj+dUj/dXi)+dUj/dXi*(dUi/dXj+dUj/dXi)]
C     Eij=Eji
C                =1/2*[dUi/dXj*(dUi/dXj+dUj/dXi)]
C     1/2*(E:E) = 1/4*EPSIJ
C
C     I2(E)=1/2(E:E)  --------------->   4*I2(E)=EPSIJ
C
C     When convergence problems appear, EPSIL at time t (instead of at
C     time t+dt) can be used (see updpro.f & nsmatx.f)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      REAL*8  TVELO(NDIME),    GRAVE(NDIME,NDIME),
     .        PROAC(    6),    FPARA(20)
C
      IERROR=0
C
C---------------------------------------------------- constant viscosity
      IF(KVILA.EQ.0) THEN
       VISCO=PROAC(2)
C
C**** TEMPERATURE & STRAIN RATE VISCOSITY MODELS
C
C------------------------- temperature & strain rate dependent viscosity
      ELSE IF(KVILA.EQ.1) THEN
       EPSIJ=0.0D0
       DO I=1,NDIME
        DO J=1,NDIME
         EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
        END DO
       END DO
       EPSIL=DSQRT(EPSIJ)
       ZETA =EPSIL*DEXP(FPARA(4)/FPARA(5)/(273.0D0+TEMPE))
       SIGMA=DLOG((ZETA/FPARA(2))**(1.0D0/FPARA(3))+
     .       DSQRT((ZETA/FPARA(2))**(2.0D0/FPARA(3))+1.0D0))/FPARA(1)
       IF(EPSIL.GT.1.0D-4) THEN
        VISCO=SIGMA/DSQRT(3.0D0)/EPSIL
        IF(VISCO.GT.PROAC(2)) VISCO=PROAC(2)
       ELSE
        VISCO=PROAC(2)
       END IF
C
C------------------------------------------------------------- power law
      ELSE IF(KVILA.EQ.2) THEN
       EPSIJ=0.0D0
       DO I=1,NDIME
        DO J=1,NDIME
         EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
        END DO
       END DO
C
       EPSIL=DSQRT(EPSIJ)                ! ORI=1
C      EPSIL=DSQRT(2.0D0*EPSIJ/3.0D0)    ! ORI=2
       IF(EPSIL.GT.1.0D-8) THEN
C       VISCO=FPARA(1)*DEXP(-FPARA(3)*TEMPE)/(EPSIL**FPARA(2))   ! ORI=1
        VISCO=FPARA(1)*DEXP(FPARA(3)/(TEMPE+273.0D0))/
     .                                         (EPSIL**FPARA(2)) ! ORI=2
        IF(VISCO.GT.PROAC(2)) VISCO=PROAC(2)
       ELSE
        VISCO=PROAC(2)
       END IF
C
C------------------------------- LDPE: law 1 Sastrohartono et al. (1995)
      ELSE IF(KVILA.EQ.3) THEN
       EPSIJ=0.0D0
       DO I=1,NDIME
        DO J=1,NDIME
         EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
        END DO
       END DO
       EPSIL=DSQRT(EPSIJ)
       IF(EPSIL.GT.1.0D-8) THEN
        VISCO=FPARA(1)*(EPSIL/FPARA(2))**(FPARA(3)-1.0D0)*
     .        DEXP(FPARA(4)*(FPARA(5)-TEMPE))
        IF(VISCO.GT.PROAC(2)) VISCO=PROAC(2)
       ELSE
        VISCO=PROAC(2)
       END IF
C
C------------------------------- LDPE: law 2 Sastrohartono et al. (1995)
      ELSE IF(KVILA.EQ.4) THEN
       EPSIJ=0.0D0
       DO I=1,NDIME
        DO J=1,NDIME
         EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
        END DO
       END DO
       EPSIL=DSQRT(EPSIJ)
       IF(EPSIL.GT.1.0D-8) THEN
        VISCO=FPARA(1)*DEXP(FPARA(2)/(TEMPE+273.0D0))/
     .        (1.0D0+FPARA(3)*(FPARA(1)*EXP(FPARA(2)/
     .        (TEMPE+273.0D0))*EPSIL)**(1-FPARA(4)))/10.0D0
        IF(VISCO.GT.PROAC(2)) VISCO=PROAC(2)
       ELSE
        VISCO=PROAC(2)
       END IF
C
C------------------------------- LDPE: law 3 Sastrohartono et al. (1995)
      ELSE IF(KVILA.EQ.5) THEN
       VISCO=FPARA(1)*
     .       DEXP(FPARA(4)*(FPARA(5)-TEMPE))
       IF(VISCO.GT.PROAC(2)) VISCO=PROAC(2)
C
C------------------------------ Perzyna's model for metal (Tesis Codina)
      ELSE IF(KVILA.EQ.6) THEN
       EPSIJ=0.0D0
       DO I=1,NDIME
        DO J=1,NDIME
         EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
        END DO
       END DO
       EPSIL=DSQRT(EPSIJ/3.0D0)
       IF(EPSIL.GT.1.0D-8) THEN
        DEFOR=DSQRT(3.0D0)*EPSIL*DEXP(FPARA(4)/
     .                                       (FPARA(5)*(TEMPE+273.0D0)))
        FLUEN=1.0D0/FPARA(1)*((DEFOR/FPARA(2))**(1.0D0/FPARA(3))+
     .                  DSQRT((DEFOR/FPARA(2))**(2.0D0/FPARA(3))+1.0D0))
        IF(FPARA(6).GT.10.0D0**5.0D0) THEN
         VISCO=FLUEN/(3.0D0*EPSIL)
        ELSE
         VISCO=(FLUEN+(EPSIL/FPARA(6))**(1.0D0/FPARA(7)))
     .                                                    /(3.0D0*EPSIL)
        END IF
        IF(VISCO.GT.PROAC(2)) VISCO=PROAC(2)
       ELSE
        VISCO=PROAC(2)
       END IF
C
C------------------------------- water viscosity: Ishikawa et al. (2000)
      ELSE IF(KVILA.EQ.7) THEN
C      DENSE=FPARA(5)/(1+FPARA(1)*TEMPE+FPARA(2)*TEMPE**2+
C    .       FPARA(3)*TEMPE**3+FPARA(4)*TEMPE**4)
       DENSE=999.8395D0/(1.0D0-0.67896452D-4*TEMPE+
     .       0.907294338D-5*TEMPE**2.0D0-
     .       0.964568125D-7*TEMPE**3.0D0+
     .       0.873702983D-9*TEMPE**4.0D0)
       VISCO=DENSE*(1.79023D-6-6.07496D-8*TEMPE+
     .              1.46296D-9*TEMPE**2.0D0-
     .              2.27242D-11*TEMPE**3.0D0+
     .              1.61111D-13*TEMPE**4.0D0)
C
C--------- water viscosity with constant density: Ishikawa et al. (2000)
      ELSE IF(KVILA.EQ.8) THEN         ! Water viscosity Ishikawa et al.
       VISCO=999.8395D0*(1.79023D-6-6.07496D-8*TEMPE+
     .                              1.46296D-9*TEMPETEMPE-
     .                             2.27242D-11*TEMPE*TEMPE*TEMPE+
     .                             1.61111D-13*TEMPE*TEMPE*TEMPE*TEMPE)
C
C-------------- water viscosity (external input): Ishikawa et al. (2000)
      ELSE IF(KVILA.EQ.9) THEN
       DENSE=FPARA(5)/(1.0D0+FPARA(1)*TEMPE+FPARA(2)*TEMPE**2.0D0+
     .       FPARA(3)*TEMPE**3.0D0+FPARA(4)*TEMPE**4.0D0)
       VISCO=DENSE*(FPARA(6)+FPARA(7)*TEMPE+FPARA(8)*TEMPE**2.0D0+
     .       FPARA(9)*TEMPE**3.0D0+FPARA(10)*TEMPE**4.0D0)
C
C- water viscosity with rho=cte (external input): Ishikawa et al. (2000)
      ELSE IF(KVILA.EQ.10) THEN
       VISCO=PROAC(1)*(FPARA(6)+FPARA(7)*TEMPE+FPARA(8)*TEMPE**2.0D0+
     .       FPARA(9)*TEMPE**3.0D0+FPARA(10)*TEMPE**4.0D0)
C
C-----------------------------------------------------------------------
      ELSE IF(KVILA.EQ.11) THEN    ! linear variation with T
       IF(PSCON.LE.CRIPS) THEN     ! identification index=1 for medium 1
        VISCO=FPARA(4)-
     .        (FPARA(4)-FPARA(3))/(FPARA(2)-FPARA(1))*(FPARA(2)-TEMPE)
        IF(TEMPE.LT.FPARA(1)) VISCO=FPARA(3)
        IF(TEMPE.GT.FPARA(2)) VISCO=FPARA(4)
       ELSE
        VISCO=PROAC(2)         ! to be generalized !!!!
       END IF
C
C-----------------------------------------------------------------------
      ELSE IF(KVILA.EQ.12) THEN    ! Non-Newtonian flow with front
       EPSIJ=0.0D0
       DO I=1,NDIME
        DO J=1,NDIME
         EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
        END DO
       END DO
       EPSIL=0.0D0
       IF(EPSIJ.GT.1.0D-8) EPSIL=DSQRT(0.5D0*EPSIJ)
       VISCO=PROAC(2)
       IF(PSCON.LE.CRIPS) THEN     ! air (PROAC(2) < FPARA(3))
        IF(FPARA(3).GT.0.0D0) THEN
         VISCO=FPARA(3)-FPARA(4)*EPSIL
         IF(VISCO.LT.PROAC(2)) VISCO=PROAC(2)
        END IF
       ELSE                        ! metal (PROAC(2) < FPARA(1))
        IF(FPARA(1).GT.0.0D0) THEN
         VISCO=FPARA(1)-FPARA(2)*EPSIL
         IF(VISCO.LT.PROAC(2)) VISCO=PROAC(2)
        END IF
       END IF
C
C-----------------------------------------------------------------------
      ELSE IF(KVILA.EQ.13) THEN    ! Bingham plastic model
       EPSIJ=0.0D0
       DO I=1,NDIME
        DO J=1,NDIME
         EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
        END DO
       END DO
       EPSIL=0.0D0
       IF(EPSIJ.GT.1.0D-8) EPSIL=DSQRT(0.5D0*EPSIJ)
       VISCO=PROAC(2)   ! FPARA(1): tau_0; FPARA(2): mu_inf, FPARA(3): n
       IF((FPARA(2)*EPSIL).LT.FPARA(1)) VISCO=FPARA(2)  ! constant
c      IF((FPARA(2)*EPSIL).LT.(FPARA(3)*FPARA(1)))      ! linear interp.
c    .  VISCO=-(FPARA(2)-PROAC(2))/(FPARA(3)*FPARA(1)/FPARA(2))*EPSIL+
c    .                                                          FPARA(2)
       IF((FPARA(2)*EPSIL).GE.FPARA(1))                 ! secant
     .  VISCO=(FPARA(1)+PROAC(2)*(EPSIL-FPARA(1)/FPARA(2)))/EPSIL
      END IF
C-----------------------------------------------------------------------
C
C**** VISCOSITY MODELS FOR PHASE-CHANGE PROBLEMS
C
C--------------------------------------------------------- mushy model 1
      IF(KVILA.EQ.101) THEN
c      FSTOL=1.0D-06
       FSTOL=PROAC(2)/FPARA(1)
       IF(FSLPG.GE.FSTOL) THEN
        VISCO=PROAC(2)/FSLPG
       ELSE
        VISCO=FPARA(1)
       END IF
C
C---------------------------------------------------- isothermic model 1
      ELSE IF(KVILA.EQ.102) THEN
       IF(FSLPG.GE.0.5D0) THEN
        VISCO=PROAC(2)
       ELSE
        VISCO=FPARA(1)
       END IF
C
C--------------------------------------------------------- mushy model 2
      ELSE IF(KVILA.EQ.103) THEN
       FSTOL=DSQRT(PROAC(2)/FPARA(1))
       IF(FSLPG.GE.FSTOL) THEN
        VISCO=PROAC(2)/(FSLPG*FSLPG)
       ELSE
        VISCO=FPARA(1)
       END IF
C
C---------------------------------------------------- isothermic model 2
      ELSE IF(KVILA.EQ.104) THEN
       FSTOL=PROAC(2)/FPARA(1)
       IF(FSLPG.GT.FSTOL) THEN
        VISCO=PROAC(2)
       ELSE
        VISCO=FPARA(1)
       END IF
C
C----------- isothermic model 3 (twin-roll process): Chang & Weng (1995)
      ELSE IF(KVILA.EQ.105) THEN
       EPSIJ=0.0D0
       DO I=1,NDIME
        DO J=1,NDIME
         EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
        END DO
       END DO
       EPSIL=0.0D0
       IF(EPSIJ.GT.1.0D-8) EPSIL=DSQRT(EPSIJ/3.0D0)
       IF(FSLPG.LT.0.5D0) THEN
        IF(EPSIJ.GT.1.0D-8) THEN
         CONS1=FPARA(1)*DEXP(FPARA(2)/TEMPE)
         CONS2=FPARA(3)*TEMPE**2.0D0+FPARA(4)*TEMPE+FPARA(5)
         VISCO=CONS1/3.0D0*EPSIL**(CONS2-1.0D0)
        ELSE
         IF(CONS2.GT.1.0D0) THEN
          VISCO=0.0D0
         ELSE
          VISCO=FPARA(6)
         ENDIF
        ENDIF
        IF(VISCO.GT.FPARA(6)) VISCO=FPARA(6)
       ELSE
        VISCO=PROAC(2)+DSQRT(2.0D0)*PROAC(1)*FPARA(7)*FPARA(7)*EPSIL
       END IF
C
C---------------- mushy model 3 (twin-roll process): Chang & Weng (1995)
      ELSE IF(KVILA.EQ.106) THEN                   ! to be revised !!!!!
       FSTOL=1.0D-06
       IF(FSLPG.LT.FSTOL) THEN
        EPSIJ=0.0D0
        DO I=1,NDIME
         DO J=1,NDIME
          EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
         END DO
        END DO
        EPSIL=DSQRT(EPSIJ/3.0D0)
        CONS1=FPARA(1)*DEXP(FPARA(2)/TEMPE)
        CONS2=FPARA(3)*TEMPE**2.0D0+FPARA(4)*TEMPE+FPARA(5)
        FLUEN=CONS1
        VISCO=FLUEN/3.0D0
        IF(EPSIL.GT.1.0E-8) THEN
         FLUEN=CONS1*EPSIL**CONS2
         VISCO=FLUEN/(3.0D0*EPSIL)
        ENDIF
        IF(VISCO.GT.FPARA(6)) VISCO=FPARA(6)
       ELSE IF(FLSPG.LE.(1.0D0-FSTOL)) THEN
        VISCO=PROAC(2)/FSLPG
       ELSE IF(FLSPG.GT.(1.0D0-FSTOL)) THEN
        VISCO=PROAC(2)
       END IF
C
C-------------------------- isothermic model 4 for tin: Raw & Lee (1991)
      ELSE IF(KVILA.EQ.107) THEN
       IF(FSLPG.GE.0.5D0) THEN
        VISCO=3.004D-3-5.53D-6*TEMPE+3.74D-9*TEMPE*TEMPE
       ELSE
        VISCO=FPARA(1)
       END IF
C
C-------- isothermic model 5 for water (rho=cte): Ishikawa et al. (2000)
      ELSE IF(KVILA.EQ.108) THEN
       IF(FSLPG.GE.0.5D0) THEN
        VISCO=999.8395D0*(1.79023D-6-6.07496D-8*TEMPE+
     .                               1.46296D-9*TEMPE*TEMPE-
     .                              2.27242D-11*TEMPE*TEMPE*TEMPE+
     .                              1.61111D-13*TEMPE*TEMPE*TEMPE*TEMPE)
       ELSE
        VISCO=FPARA(1)
       END IF
C
C-------------------- isothermic model 6: idem KVILA=108 with turbulence
      ELSE IF(KVILA.EQ.109) THEN
       IF(FSLPG.GE.0.5D0) THEN
        EPSIJ=0.0D0
        DO I=1,NDIME
         DO J=1,NDIME
          EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
         END DO
        END DO
        EPSIL=0.0D0
        IF(EPSIJ.GT.1.0D-8) EPSIL=DSQRT(0.5D0*EPSIJ)
        VISCO=999.8395D0*(1.79023D-6-6.07496D-8*TEMPE+
     .                               1.46296D-9*TEMPE*TEMPE-
     .                              2.27242D-11*TEMPE*TEMPE*TEMPE+
     .                              1.61111D-13*TEMPE*TEMPE*TEMPE*TEMPE)
        VISCO=VISCO+DSQRT(2.0D0)*999.8395D0*FPARA(2)*FPARA(2)*EPSIL
       ELSE
        VISCO=FPARA(1)
       END IF
C
C--------------- isothermic model 7 (the contrary of isothermic model 2)
      ELSE IF(KVILA.EQ.110) THEN
       FSTOL=PROAC(2)/FPARA(1)
       IF(FSLPG.GT.(1.0D0-FSTOL)) THEN
        VISCO=PROAC(2)
       ELSE
        VISCO=FPARA(1)
       END IF
C
C------------------ mushy model 3 (generalization of mushy models 1 & 2)
      ELSE IF(KVILA.EQ.111) THEN
       IF(FPARA(2).EQ.0.0D0) IERROR=1
       FSTOL=(PROAC(2)/FPARA(1))**(1.0D0/FPARA(2))
       IF(FSLPG.GT.FSTOL) THEN
        VISCO=PROAC(2)/FSLPG**FPARA(2)
       ELSE
        VISCO=FPARA(1)
       END IF
C
C-----------------------------------------------------------------------
      END IF
C
C**** ALGEBRAIC TURBULENCE VISCOSITY MODELS FOR FILLING PROBLEMS
C
C--------------------------------------- metal: turbulent - air: laminar
      IF(KVILA.EQ.200) THEN
       IF(PSCON.LE.CRIPS) THEN
        VISCO=PROAC(2)
       ELSE
        EPSIJ=0.0D0
        DO I=1,NDIME
         DO J=1,NDIME
          EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
         END DO
        END DO
        EPSIL=0.0D0
        IF(EPSIJ.GT.1.0D-8) EPSIL=DSQRT(0.5D0*EPSIJ)
        VISCO=PROAC(2)+DSQRT(2.0D0)*PROAC(1)*FPARA(1)*FPARA(1)*EPSIL
       END IF
      END IF
C
C------------------------------------------------ metal & air: turbulent
      IF(KVILA.EQ.201) THEN
       EPSIJ=0.0D0
       DO I=1,NDIME
        DO J=1,NDIME
         EPSIJ=EPSIJ+GRAVE(I,J)*(GRAVE(I,J)+GRAVE(J,I))
        END DO
       END DO
       EPSIL=0.0D0
       IF(EPSIJ.GT.1.0D-8) EPSIL=DSQRT(0.5D0*EPSIJ)
       VISCO=PROAC(2)
       JFRE1=2*(KINF3-1)+1         ! indexes for free parameters
       JFRE2=2*(KINF4-1)+1
       IF(PSCON.LE.CRIPS) THEN     ! default material ("air")
        IF(FPARA(JFRE2).GT.0.0D0)
     .   VISCO=PROAC(2)+           ! turbulent viscosity
     .             DSQRT(2.0D0)*PROAC(1)*FPARA(JFRE2)*FPARA(JFRE2)*EPSIL
        IF(FPARA(JFRE2+1).GT.0.0D0) THEN
         IF(VISCO.GT.FPARA(JFRE2+1)) VISCO=FPARA(JFRE2+1)      ! cut-off
        END IF
       ELSE                        ! inner material ("metal")
        IF(FPARA(JFRE1).GT.0.0D0)
     .   VISCO=PROAC(2)+           ! turbulent viscosity
     .             DSQRT(2.0D0)*PROAC(1)*FPARA(JFRE1)*FPARA(JFRE1)*EPSIL
        IF(FPARA(JFRE1+1).GT.0.0D0) THEN
         IF(VISCO.GT.FPARA(JFRE1+1)) VISCO=FPARA(JFRE1+1)      ! cut-off
        END IF
       END IF
      END IF
C
C-----------------------------------------------------------------------
C
      END
