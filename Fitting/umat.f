C======================================================================
      MODULE GLOBVAR
        integer, parameter :: MaxNodes = 2730
        integer, parameter :: MaxElem = 2000
        real, parameter :: thickness = 2.0e-3    ! wall thickness
        integer, parameter :: diff = 2499   !difference of element No. between two seccessive layer
        integer, parameter :: nl = 4             ! number of layers through the thickness
        integer, parameter :: NNL = 25           ! number of elements in longitudinal direction
        integer, parameter :: NT = 20            ! number of elements in circumferencial direction
  
        integer, dimension(MaxElem,9) :: conectivity  ! 
        real(8), dimension(MaxNodes,4) :: coordinate  !
        real(8), dimension(MaxElem) :: Xc, Yc, Zc     ! centeroid of elements
  
        real(8), dimension(3*MaxElem,3) :: Rel
        real(8), dimension(8*3*MaxElem,3) :: DFGin
        real(8), dimension(8*3*MaxElem,3) :: DFGnew
        real(8), dimension(MaxElem) :: volume
        real(8), dimension(nl,3) :: av_Fl
c Homeostatic stress 
        real(8), dimension(MaxElem,8) :: S11,S22,S33,S12,S13,S23,S_smc
        real(8), dimension(MaxElem,8,4) :: S_col
      END MODULE
C======================================================================
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
C
      USE GLOBVAR
  
      INCLUDE 'aba_param.inc'
  
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C definition of parameters	
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! Variables passed in for information
      integer :: NDI            ! number of direct stress components at this point
      integer :: NSHR           ! number of shear stress components
      integer :: NTENS          ! size of the stress or strain component array (ndi + nshr)
      integer :: NSTATV         ! number of state variables
      integer :: NPT            ! Integration point number (punto de gauss)
      integer :: NOEL           ! Element number
  
      real(8), dimension(NTENS) :: STRESS  ! Cauchy stress
      real(8), dimension(NSTATV) :: STATEV    ! state dependent variables
      real(8), dimension(NTENS,NTENS) :: DDSDDE   ! spatial stiffness tensor
      real(8), dimension(NTENS) :: DDSDDT  ! variation of stress wrt temperature
      real(8), dimension(NTENS) :: DRPLDE  ! variation of heat generation wrt strain
      real(8), dimension(NTENS) :: STRAN  ! total strain at the beginning of the increment
      real(8), dimension(NTENS) :: DSTRAN  ! strain increments
      real(8), dimension(1) :: PREDEF  ! predefined field variables at the start of the increment
      real(8), dimension(1) :: DPRED  ! predefined field increments
      real(8), dimension(NPROPS) :: PROPS  ! user-defined material constants
      real(8), dimension(3) :: COORDS  ! coordinates of this point
      real(8), dimension(3,3) :: DROT  ! rotation increment matrix
      real(8), dimension(3,3) :: DFGRD0   ! deformation gradient at the beginning of the increment
      real(8), dimension(3,3) :: DFGRD1   ! deformation gradient at the end of the increment
      real(8), dimension(2) :: TIME    ! Actual time
      real(8) :: DTIME
      character*8 CMNAME  ! user defined material name

! General assignments
      integer, parameter :: noff=4 !number of fiber family
      real(8), parameter :: zero=0.d0, half=0.5d0, one=1.d0, 
     1         two=2.d0, three=3.d0, four=4.d0, five=5.d0,
     2         six=6.d0, seven=7.d0, eight=8.d0, nine=9.d0
      real(8), parameter :: pi = 3.141592653589d0
! for initialization
      real(8) :: t_step,t_tot
! Some importatn variables
      integer :: reference,GR
      real(8), dimension(3,3) :: R33,identity        ! Identity matrix
      real(8), dimension(3,3) :: F,C,inv_F          ! Mixture deformation gradient
      real(8) :: DET

C Incompressible parameters
      real(8), dimension(ntens) :: STRESS_u  ! Cauchy stress of volumetric part
      real(8), dimension(ntens,ntens) :: DSDE_u
      real(8) :: kapa, U, EK, PR

C Elastic phenomena
      !Elastin parameters
      real(8), dimension(ntens,ntens) :: DSDE_e      ! spatial stiffness tensor of elastin 
      real(8), dimension(ntens) :: STRESS_e      ! Cauchy stress of elastin 
      real(8) :: C10,EG,EG23,DET_e,Ibar1
      real(8), dimension(3,3) :: G_e,F_e,Fbar_e,C_e
      real(8), dimension(3,3) :: Bbar_e,inv_C_e
      real(8), dimension(6) :: BBAR   !Left Cauchy-Green tensor of isochoric part of elastic
      !Collagen parameters
      real(8), dimension(noff) :: alpha
      real(8), dimension(ntens,ntens) :: DSDE_c,DSDE_ci
      real(8), dimension(ntens) :: STRESS_c,STRESS_ci
      real(8), dimension(3,3) :: F_ci,C_ci,m0m0
      real(8), dimension(3) :: m0,m0_ci,mi
      real(8) :: c2t,c3t,c2c,c3c,gc,c2,c3
      real(8) :: Ibar4_c,lambda_ci
      real(8) :: coef0,coef1,coef2
      !SMC parameters
      real(8), dimension(ntens,ntens) :: DSDE_smc
      real(8), dimension(ntens) :: STRESS_smc
      real(8), dimension(3,3) :: F_smc,C_smc
      real(8), dimension(3) :: m0_smc
      real(8) :: c2t_smc,c3t_smc,c2c_smc,c3c_smc,gsmc
      real(8) :: alpha_smc,Ibar4_smc,lambda_smc

C Densities and Fractions mass
      real(8) :: den_e,den_smc,den_tot
      real(8) :: rho_e,rho_smc
      real(8),dimension(noff) :: den_c,rho_c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Initialization of effects
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      t_step = TIME(1) + DTIME
      t_tot = TIME(2) + DTIME
C Growth and Remodeling
      GR = 2
C Densities
      den_e = 241.5d0    ! elastin density
      den_c(1) = 65.1d0      !alfa = 0
      den_c(2) = den_c(1)  !alfa = 90
      den_c(3) = 260.4d0     !alfa = 45
      den_c(4) = den_c(3)  !alfa = -45
      den_smc = 157.5d0  ! smooth muscle cells density
      den_tot = den_e + den_c(1) + den_c(2)
     1          + den_c(3) + den_c(4) + den_smc
      rho_e = den_e/den_tot  ! elastin mass fraction
      do i=1,noff
        rho_c(i) = den_c(i)/den_tot  ! collagen mass fraction
      enddo
      rho_smc = den_smc/den_tot
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C General assignments
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      identity = reshape((/ one, zero, zero, zero, 
     1 one, zero, zero, zero, one /), shape(identity)) ! identity matrix
C
      R33 = identity
C 
C Calculate deformation gradient
C
      F = DFGRD1
      call determinant(F,DET)
      C = MATMUL(transpose(F),F)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Calculate the stiffness and cauchy sttress of compressibility in mixture
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      C10 = PROPS(2)                !in my model mu=2c_e ()
      kapa = PROPS(3)

      kapa = rho_e*C10*1.0e4

      EK = two*kapa*(two*DET - one)
      PR = two*kapa*(DET - one)

      STRESS_u = zero
      DO K1 = 1,NDI
        STRESS_u(K1) = PR
      END DO 
  
      DSDE_u = zero
      DSDE_u(1,1) = EK
      DSDE_u(2,2) = EK
      DSDE_u(3,3) = EK
      DSDE_u(1,2) = EK
      DSDE_u(1,3) = EK
      DSDE_u(2,3) = EK  

      DO K1 = 1, NTENS
        DO K2 = 1, K1-1
          DSDE_u(K1,K2) = DSDE_u(K2,K1)
        END DO
      END DO
  
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
C Calculate the stiffness and cauchy sttress of elastic fibers   
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
C for deposition stretch
      reference = 0
! Deposition stretch of elastin
      G_e = indentity
      G_e(3,3) = one + t_step*0.1d0                           !Z 
      if (reference.EQ.0) then
        G_e(2,2) = one + t_step*(1.1832064 - one)                     !Theta
      else
        G_e(2,2) = one + t_step*(DFGin(2,2) - one)
      endif
      G_e(1,1) = one + t_step*(one/(G_e(3,3)*G_e(2,2)) - one) !R
      F_e = MATMUL(F,G_e)
!------------------------
!saving deformation gradient of the elastic fiber dominated tissue !!
      STATEV(1)=F_e(1,1)
      STATEV(2)=F_e(1,2)
      STATEV(3)=F_e(1,3)
      STATEV(4)=F_e(2,1)
      STATEV(5)=F_e(2,2)
      STATEV(6)=F_e(2,3)
      STATEV(7)=F_e(3,1)
      STATEV(8)=F_e(3,2)
      STATEV(9)=F_e(3,3)
!------------------------------------
      call determinant(F_e,DET_e)
      Fbar_e = DET_e**(-one/three)*F_e
      Bbar_e = MATMUL(Fbar_e,transpose(Fbar_e))
      BBAR(1) = Bbar_e(1,1)
      BBAR(2) = Bbar_e(2,2)
      BBAR(3) = Bbar_e(3,3)
      BBAR(4) = Bbar_e(1,2)
      BBAR(5) = Bbar_e(1,3)
      BBAR(6) = Bbar_e(2,3)
      Ibar1 = BBAR(1) + BBAR(2) + BBAR(3)
!----------------------
      EG = rho_e*two*C10/DET
      !stress on neo-hook model
      do K1=1,NDI
        STRESS_e(K1)=EG*(BBAR(K1)-Ibar1/three)
      enddo
      do K1=NDI+1,NDI+NSHR
        STRESS_e(K1)=EG*BBAR(K1)
      enddo
      !elasticity tensor on neo-hook model
      EG23=EG*two/three
      DSDE_e(1, 1)= EG23*(BBAR(1) + Ibar1/three)
      DSDE_e(1, 2)=-EG23*(BBAR(1) + BBAR(2) - Ibar1/three)  
      DSDE_e(2, 2)= EG23*(BBAR(2) + Ibar1/three)
      DSDE_e(1, 3)=-EG23*(BBAR(1) + BBAR(3) - Ibar1/three)
      DSDE_e(2, 3)=-EG23*(BBAR(2) + BBAR(3) - Ibar1/three) 
      DSDE_e(3, 3)= EG23*(BBAR(3) + Ibar1/three)
      DSDE_e(1, 4)= EG23*BBAR(4)/two
      DSDE_e(2, 4)= EG23*BBAR(4)/two
      DSDE_e(3, 4)=-EG23*BBAR(4)
      DSDE_e(4, 4)= EG*(BBAR(1) + BBAR(2))/two
      DSDE_e(1, 5)= EG23*BBAR(5)/two
      DSDE_e(2, 5)=-EG23*BBAR(5)
      DSDE_e(3, 5)= EG23*BBAR(5)/two
      DSDE_e(4, 5)= EG*BBAR(6)/two
      DSDE_e(5, 5)= EG*(BBAR(1) + BBAR(3))/two
      DSDE_e(1, 6)=-EG23*BBAR(6)
      DSDE_e(2, 6)= EG23*BBAR(6)/two
      DSDE_e(3, 6)= EG23*BBAR(6)/two
      DSDE_e(4, 6)= EG*BBAR(5)/two
      DSDE_e(5, 6)= EG*BBAR(4)/two
      DSDE_e(6, 6)= EG*(BBAR(2) + BBAR(3))/two
  
      do K1=1, NTENS
        do K2=1, K1-1
          DSDE_e(K1,K2) = DSDE_e(K2,K1)
        enddo
      enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C Calculate the stiffness and cauchy sttress of Collagen fibers                       
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 	  
      alpha(1) = zero                ! axial              ---> z
      alpha(2) = +Pi/two             ! circumferential    ---> teta
      alpha(3) = +Pi/four            ! diagonal
      alpha(4) = -Pi/four            ! diagonal

!------------------------------------------------------	  
      STRESS_c = zero
      DSDE_c = zero
      do i = 1, noff  ! noff is number of callagen fiber family
        c2t = props(4)
        c3t = props(5)
        c2c = props(8) !zero
        c3c = props(9) !zero
        gc = one + t_step*0.1d0     ! gc = lambda_ci

        m0(:) = (/ zero, sin(alpha(i)), cos(alpha(i)) /) ! fiber vector in cylidrical coordinate
        m0_ci = MATMUL(R33,m0)  ! fiber direction in global coordinates

        do k1 = 1,3
          do k2 = 1,3
            m0m0(k1,k2) = m0_ci(k1) * m0_ci(k2)
          enddo
        enddo
        call dd33_product(C,m0m0,lambda_ci)
        lambda_ci = gc*sqrt(lambda_ci)  !multiplicar gc para el CMT
        STATEV(9+i) = lambda_ci
        F_ci = lambda_ci*m0m0+one/sqrt(lambda_ci)*(identity-m0m0)
        Ibar4_c = lambda_ci**two
        mi = MATMUL(F_ci,m0_ci)
C---------------------------------------------------------------------- 
C                          MATERIAL PARAMETERS
C----------------------------------------------------------------------    
        if (Ibar4_c.GE.one) then
          c2 = c2t
          c3 = c3t
        else
          c2 = c2c
          c3 = c3c
        endif  

        coef0=rho_c(i)*c2*(Ibar4_c-one)*exp(c3*(Ibar4_c-one)**two)/DET

        coef1 = coef0*two/three

        coef2 = rho_c(i)*two*c2*exp(c3*(Ibar4_c-one)**two)*(one
     1        + two*c3*(Ibar4_c-one)**two)/(DET*three)

        STRESS_ci(1) = coef0*(mi(1)*mi(1)-Ibar4_c/three)
        STRESS_ci(2) = coef0*(mi(2)*mi(2)-Ibar4_c/three)
        STRESS_ci(3) = coef0*(mi(3)*mi(3)-Ibar4_c/three)
        STRESS_ci(4) = coef0*mi(1)*mi(2)
        STRESS_ci(5) = coef0*mi(1)*mi(3)
        STRESS_ci(6) = coef0*mi(2)*mi(3)

        DSDE_ci(1,1) = coef1*(mi(1)*mi(1)+Ibar4_c/three) 
     1           + coef2*(mi(1)*mi(1)*mi(1)*mi(1)*three
     2           - Ibar4_c*(two*mi(1)*mi(1)-one/three))
        DSDE_ci(1,2) =-coef1*(mi(1)*mi(1)+mi(2)*mi(2)-Ibar4_c/three)
     1           + coef2*(mi(1)*mi(1)*mi(2)*mi(2)*three
     2           - Ibar4_c*(mi(1)*mi(1)+mi(2)*mi(2)-one/three))
        DSDE_ci(1,3) =-coef1*(mi(1)*mi(1)+mi(3)*mi(3)-Ibar4_c/three)
     1           + coef2*(mi(1)*mi(1)*mi(3)*mi(3)*three
     2           - Ibar4_c*(mi(1)*mi(1)+mi(3)*mi(3)-one/three))
        DSDE_ci(1,4) = coef1*mi(1)*mi(2)/two
     1           + coef2*(mi(1)*mi(1)*mi(1)*mi(2)*three
     2           - Ibar4_c*mi(1)*mi(2))
        DSDE_ci(1,5) = coef1*mi(1)*mi(3)/two
     1           + coef2*(mi(1)*mi(1)*mi(1)*mi(3)*three
     2           - Ibar4_c*mi(1)*mi(3))
        DSDE_ci(1,6) =-coef1*mi(2)*mi(3)
     1           + coef2*(mi(1)*mi(1)*mi(2)*mi(3)*three
     2           - Ibar4_c*mi(2)*mi(3))
!----------------
        DSDE_ci(2,2) = coef1*(mi(2)*mi(2)+Ibar4_c/three) 
     1           + coef2*(mi(2)*mi(2)*mi(2)*mi(2)*three
     2           - Ibar4_c*(two*mi(2)*mi(2)-one/three))
        DSDE_ci(2,3) =-coef1*(mi(2)*mi(2)+mi(3)*mi(3)-Ibar4_c/three)
     1           + coef2*(mi(2)*mi(2)*mi(3)*mi(3)*three
     2           - Ibar4_c*(mi(2)*mi(2)+mi(3)*mi(3)-one/three))
        DSDE_ci(2,4) = coef1*mi(1)*mi(2)/two
     1           + coef2*(mi(2)*mi(2)*mi(1)*mi(2)*three
     2           - Ibar4_c*mi(1)*mi(2))
        DSDE_ci(2,5) =-coef1*mi(1)*mi(3)
     1           + coef2*(mi(2)*mi(2)*mi(1)*mi(3)*three
     2           - Ibar4_c*mi(1)*mi(3))
        DSDE_ci(2,6) = coef1*mi(2)*mi(3)/two
     1           + coef2*(mi(2)*mi(2)*mi(2)*mi(3)*three
     2           - Ibar4_c*mi(2)*mi(3))
!----------------
        DSDE_ci(3,3) = coef1*(mi(3)*mi(3)+Ibar4_c/three) 
     1           + coef2*(mi(3)*mi(3)*mi(3)*mi(3)*three
     2           - Ibar4_c*(two*mi(3)*mi(3)-one/three))
        DSDE_ci(3,4) =-coef1*mi(1)*mi(2)
     1           + coef2*(mi(3)*mi(3)*mi(1)*mi(2)*three
     2           - Ibar4_c*mi(1)*mi(2))
        DSDE_ci(3,5) = coef1*mi(1)*mi(3)/two
     1           + coef2*(mi(3)*mi(3)*mi(1)*mi(3)*three
     2           - Ibar4_c*mi(1)*mi(3))
        DSDE_ci(3,6) = coef1*mi(2)*mi(3)/two
     1           + coef2*(mi(3)*mi(3)*mi(2)*mi(3)*three
     2           - Ibar4_c*mi(2)*mi(3))
!----------------
        DSDE_ci(4,4) = coef1*(mi(1)*mi(1)+mi(2)*mi(2))*three/four 
     1           + coef2*mi(1)*mi(2)*mi(1)*mi(2)*three
        DSDE_ci(4,5) = coef1*mi(2)*mi(3)*three/four 
     1           + coef2*mi(1)*mi(2)*mi(1)*mi(3)*three
        DSDE_ci(4,6) = coef1*mi(1)*mi(3)*three/four
     1           + coef2*mi(1)*mi(2)*mi(2)*mi(3)*three
!----------------
        DSDE_ci(5,5) = coef1*(mi(1)*mi(1)+mi(3)*mi(3))*three/four
     1           + coef2*mi(1)*mi(3)*mi(1)*mi(3)*three
        DSDE_ci(5,6) = coef1*mi(1)*mi(2)*three/four
     1           + coef2*mi(1)*mi(3)*mi(2)*mi(3)*three
!----------------
        DSDE_ci(6,6) = coef1*(mi(2)*mi(2)+mi(3)*mi(3))*three/four
     1           + coef2*mi(2)*mi(3)*mi(2)*mi(3)*three
C
        do j=1, 6
          do k=1, j-1
            DSDE_ci(j,k) = DSDE_ci(k,j)
          enddo
        enddo
C fiber effects superposition
        do k1=1,ntens
          STRESS_c(k1) = STRESS_c(k1) + STRESS_ci(k1)
          do k2=1,ntens
            DSDE_c(k1,k2) = DSDE_c(k1,k2) + DSDE_ci(k1,k2)
          enddo
        enddo
      enddo !number of fiber families
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C Calculate the stiffness and cauchy sttress of SMC                      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 	  
      alpha_smc = +Pi/two                                ! circumferential

      c2t_smc = props(6)
      c3t_smc = props(7)
      c2c_smc = props(8)!zero
      c3c_smc = props(9)!zero
      gsmc = one + t_step*0.1d0

      m0(:) = (/ zero, sin(alpha_smc), cos(alpha_smc) /) ! fiber vector
      m0_smc = MATMUL(R33,m0) !fiber vector in global coordinates

      do k1 = 1,3
        do k2 = 1,3
          m0m0(k1,k2) = m0_smc(k1) * m0_smc(k2)
        enddo
      enddo
      call dd33_product(C,m0m0,lambda_smc)
      lambda_smc = gsmc*sqrt(lambda_smc)  !multiplicar gc para el CMT
      STATEV(14) = lambda_smc
      F_smc = lambda_smc*m0m0+one/sqrt(lambda_smc)*(identity-m0m0)
      Ibar4_smc = lambda_smc**two
      mi = MATMUL(F_smc,m0_smc)
C----------------------------------------------------------------------
C                          MATERIAL PARAMETERS
C----------------------------------------------------------------------
      if (Ibar4_smc.GE.one) then
        c2 = c2t_smc
        c3 = c3t_smc
      else
        c2 = c2c_smc
        c3 = c3c_smc
      endif  

      coef0=rho_smc*c2*(Ibar4_smc-one)*exp(c3*(Ibar4_smc-one)**two)/DET

      coef1 = coef0*two/three

      coef2 = rho_smc*two*c2*exp(c3*(Ibar4_smc-one)**two)*(one
     1        + two*c3*(Ibar4_smc-one)**two)/(DET*three)

      STRESS_smc(1) = coef0*(mi(1)*mi(1)-Ibar4_smc/three)
      STRESS_smc(2) = coef0*(mi(2)*mi(2)-Ibar4_smc/three)
      STRESS_smc(3) = coef0*(mi(3)*mi(3)-Ibar4_smc/three)
      STRESS_smc(4) = coef0*mi(1)*mi(2)
      STRESS_smc(5) = coef0*mi(1)*mi(3)
      STRESS_smc(6) = coef0*mi(2)*mi(3)

      DSDE_smc(1,1) = coef1*(mi(1)*mi(1)+Ibar4_smc/three) 
     1           + coef2*(mi(1)*mi(1)*mi(1)*mi(1)*three
     2           - Ibar4_smc*(two*mi(1)*mi(1)-one/three))
      DSDE_smc(1,2) =-coef1*(mi(1)*mi(1)+mi(2)*mi(2)-Ibar4_smc/three)
     1           + coef2*(mi(1)*mi(1)*mi(2)*mi(2)*three
     2           - Ibar4_smc*(mi(1)*mi(1)+mi(2)*mi(2)-one/three))
      DSDE_smc(1,3) =-coef1*(mi(1)*mi(1)+mi(3)*mi(3)-Ibar4_smc/three)
     1           + coef2*(mi(1)*mi(1)*mi(3)*mi(3)*three
     2           - Ibar4_smc*(mi(1)*mi(1)+mi(3)*mi(3)-one/three))
      DSDE_smc(1,4) = coef1*mi(1)*mi(2)/two
     1           + coef2*(mi(1)*mi(1)*mi(1)*mi(2)*three
     2           - Ibar4_smc*mi(1)*mi(2))
      DSDE_smc(1,5) = coef1*mi(1)*mi(3)/two
     1           + coef2*(mi(1)*mi(1)*mi(1)*mi(3)*three
     2           - Ibar4_smc*mi(1)*mi(3))
      DSDE_smc(1,6) =-coef1*mi(2)*mi(3)
     1           + coef2*(mi(1)*mi(1)*mi(2)*mi(3)*three
     2           - Ibar4_smc*mi(2)*mi(3))
!----------------
      DSDE_smc(2,2) = coef1*(mi(2)*mi(2)+Ibar4_smc/three) 
     1           + coef2*(mi(2)*mi(2)*mi(2)*mi(2)*three
     2           - Ibar4_smc*(two*mi(2)*mi(2)-one/three))
      DSDE_smc(2,3) =-coef1*(mi(2)*mi(2)+mi(3)*mi(3)-Ibar4_smc/three)
     1           + coef2*(mi(2)*mi(2)*mi(3)*mi(3)*three
     2           - Ibar4_smc*(mi(2)*mi(2)+mi(3)*mi(3)-one/three))
      DSDE_smc(2,4) = coef1*mi(1)*mi(2)/two
     1           + coef2*(mi(2)*mi(2)*mi(1)*mi(2)*three
     2           - Ibar4_smc*mi(1)*mi(2))
      DSDE_smc(2,5) =-coef1*mi(1)*mi(3)
     1           + coef2*(mi(2)*mi(2)*mi(1)*mi(3)*three
     2           - Ibar4_smc*mi(1)*mi(3))
      DSDE_smc(2,6) = coef1*mi(2)*mi(3)/two
     1           + coef2*(mi(2)*mi(2)*mi(2)*mi(3)*three
     2           - Ibar4_smc*mi(2)*mi(3))
!----------------
      DSDE_smc(3,3) = coef1*(mi(3)*mi(3)+Ibar4_smc/three) 
     1           + coef2*(mi(3)*mi(3)*mi(3)*mi(3)*three
     2           - Ibar4_smc*(two*mi(3)*mi(3)-one/three))
      DSDE_smc(3,4) =-coef1*mi(1)*mi(2)
     1           + coef2*(mi(3)*mi(3)*mi(1)*mi(2)*three
     2           - Ibar4_smc*mi(1)*mi(2))
      DSDE_smc(3,5) = coef1*mi(1)*mi(3)/two
     1           + coef2*(mi(3)*mi(3)*mi(1)*mi(3)*three
     2           - Ibar4_smc*mi(1)*mi(3))
      DSDE_smc(3,6) = coef1*mi(2)*mi(3)/two
     1           + coef2*(mi(3)*mi(3)*mi(2)*mi(3)*three
     2           - Ibar4_smc*mi(2)*mi(3))
!----------------
      DSDE_smc(4,4) = coef1*(mi(1)*mi(1)+mi(2)*mi(2))*three/four 
     1           + coef2*mi(1)*mi(2)*mi(1)*mi(2)*three
      DSDE_smc(4,5) = coef1*mi(2)*mi(3)*three/four 
     1           + coef2*mi(1)*mi(2)*mi(1)*mi(3)*three
      DSDE_smc(4,6) = coef1*mi(1)*mi(3)*three/four
     1           + coef2*mi(1)*mi(2)*mi(2)*mi(3)*three
!----------------
      DSDE_smc(5,5) = coef1*(mi(1)*mi(1)+mi(3)*mi(3))*three/four
     1           + coef2*mi(1)*mi(3)*mi(1)*mi(3)*three
      DSDE_smc(5,6) = coef1*mi(1)*mi(2)*three/four
     1           + coef2*mi(1)*mi(3)*mi(2)*mi(3)*three
!----------------
      DSDE_smc(6,6) = coef1*(mi(2)*mi(2)+mi(3)*mi(3))*three/four
     1           + coef2*mi(2)*mi(3)*mi(2)*mi(3)*three
C
      do j=1, 6
        do k=1, j-1
          DSDE_smc(j,k) = DSDE_smc(k,j)
        enddo
      enddo
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C Calculate TOTAL stiffness and cauchy sttress
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 	  

      do k1=1,ntens
        STRESS(k1) = STRESS_u(k1) + STRESS_e(k1)
     1             + STRESS_c(k1) + STRESS_smc(k1)
        do k2=1,ntens
          DDSDDE(k1,k2) = DSDE_u(k1,k2) + DSDE_e(k1,k2)
     1                  + DSDE_c(k1,k2) + DSDE_smc(k1,k2)
        enddo
      enddo

      RETURN
      END
C================================================================
C Usubroutine to manage Udefined external databases and calculate model-independent history
C================================================================
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

      USE GLOBVAR
      INCLUDE 'ABA_PARAM.INC'

      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR

C  reading decomposion and developmental stretches tensor of elastic fiber dominated tissue (G_e)

      call GETOUTDIR(JOBDIR,LENJOBDIR)
      filename=' '
      filename(1:lenjobdir)=jobdir(1:lenjobdir)
      if (LOP.EQ.0) then
C------------------------------------------------------	Deformation gradient
        filename(lenjobdir+1:lenjobdir+10)='/DFGin.dat'
        open(UNIT=30,file=filename(1:lenjobdir+10), status='unknown')

        do i=1,24*MaxElem
          read(30,*) DFGin(i,1),DFGin(i,2),DFGin(i,3)
        enddo

        close(30)
C------------------------------------------------------	Center of the element 
        filename(lenjobdir+1:lenjobdir+8)='/Xce.dat'
        open(UNIT=70,file=filename(1:lenjobdir+8), status='unknown')

        do i=1,MaxElem
          read(70,*) Xc(i),Yc(i),Zc(i)
        enddo

        close(70)

      endif
  
      RETURN
      END
C================================================================
C User subroutine to manage STATEV=ARRAY
C================================================================
      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)
  
      USE GLOBVAR
      INCLUDE 'ABA_PARAM.INC'
  
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE(ARRAY(1),JRRAY(1,1))
      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR
      REAL*8 U1(MaxNodes), U2(MaxNodes), U3(MaxNodes), avrdisp
      REAL*8 S11c(MaxElem), S22c(MaxElem), S33c(MaxElem), 
     1       av_S11c, av_S22c, av_S33c
      INTEGER m, n, p, i
C 	  	  
      call GETOUTDIR(JOBDIR,LENJOBDIR)
      filename=' '
      filename(1:lenjobdir)=jobdir(1:lenjobdir)

      filename(lenjobdir+1:lenjobdir+10)='/DFGin.dat'
      open(UNIT=30,file=filename(1:lenjobdir+10), status='unknown')
  
      call POSFIL(KSTEP,KINC,ARRAY,JRCD)
  
      p = 0 
      Volume = 0.0
      Landax = 0.0
      Landay = 0.0
      Landaz = 0.0

      S11c = 0.0
      S22c = 0.0
      S33c = 0.0
      kk=0

      do   !cycle do inifinite
c       DBFILE(LOP,ARRAY,JRCD)	  
        call DBFILE(0,ARRAY,JRCD)
        if (JRCD.NE.0) EXIT
        KEY=JRRAY(1,2)
C Record 101 contains U
        if (KEY.EQ.101) then
          p=p+1
          U1(p)=ARRAY(4)
          U2(p)=ARRAY(5)
          U3(p)=ARRAY(6)
        endif
C Record 1 contains element information for subsequent
        if (KEY.EQ.1) then
          m = JRRAY(1,3)   ! element No.
          n = JRRAY(1,4)   !No. of Guass point
        endif

        if (KEY.EQ.5) then
          write(30,100) ARRAY(3),ARRAY(4),ARRAY(5)
          write(30,100) ARRAY(6),ARRAY(7),ARRAY(8)
          write(30,100) ARRAY(9),ARRAY(10),ARRAY(11)
          DFGnew(kk+1,1)=ARRAY(3) - 1.d0  ! we need decimal part
          DFGnew(kk+1,2)=ARRAY(4)
          DFGnew(kk+1,3)=ARRAY(5)

          DFGnew(kk+2,1)=ARRAY(6)
          DFGnew(kk+2,2)=ARRAY(7) - 1.d0
          DFGnew(kk+2,3)=ARRAY(8)

          DFGnew(kk+3,1)=ARRAY(9)
          DFGnew(kk+3,2)=ARRAY(10)
          DFGnew(kk+3,3)=ARRAY(11) - 1.d0

          kk=kk+3

        endif

        if (KEY.EQ.76) then
          Volume(m) = Volume(m) + ARRAY(3)
        endif

        if (KEY.EQ.11) then 
!key=11 stress components -----! key=12 SINV: stress invariants (Starting from ARRAY(3))
          S11c(m) = S11c(m) + ARRAY(3) !direction 1 centroid
          S22c(m) = S22c(m) + ARRAY(4) !direction 2 centroid
          S33c(m) = S33c(m) + ARRAY(5) !direction 3 centroid
        endif

      enddo  !cycle do infinite
C
c110   CONTINUE

      close(30)

100   FORMAT(3(E15.8,1X))

      S11c = S11c/8.0 !direction 1 centroid
      S22c = S22c/8.0 !direction 2 centroid
      S33c = S33c/8.0 !direction 3 centroid

      avrdisp = 0.d0
      do i=1,p
        avrdisp = avrdisp + sqrt(U1(i)**2 + U2(i)**2 + U3(i)**2)
      enddo
      avrdisp = avrdisp/p

      print*,'======================================'
c      print*,'p=', p
      print*,'avrdisp=', avrdisp
      print*,'avrdisp/thickness%=', avrdisp/thickness*100.d0

      totalVolume = 0.0
      av_S11c = 0.0
      av_S22c = 0.0  
      av_S33c = 0.0
      do i=1,m
        totalVolume = totalVolume + Volume(i)
        av_S11c = av_S11c + S11c(i)*Volume(i)
        av_S22c = av_S22c + S22c(i)*Volume(i)
        av_S33c = av_S33c + S33c(i)*Volume(i)
      enddo
      av_S11c = av_S11c/totalVolume
      av_S22c = av_S22c/totalVolume
      av_S33c = av_S33c/totalVolume

      print*,'Total Volume=',totalVolume
      print*,'  '

      print*,'Stresses in three directions'
      print*, 'S11=',av_S11c
      print*, 'S22=',av_S22c
      print*, 'S33=',av_S33c

      print*,'======================================'

120   FORMAT(8(E15.8,1X))

      RETURN
      END
C================================================================
C ORIENT subroutine to manage fiber and material orientaion
C It should be compatible with the mesh and the conectivity
C================================================================  
      SUBROUTINE ORIENT(T,NOEL,NPT,LAYER,KSPT,COORDS,BASIS,ANAME,
     1 NNODES,CNODES,JNNUM)
C
      USE GLOBVAR
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 ANAME
      DIMENSION T(3,3),COORDS(3),BASIS(3,3),CNODES(3,NNODES)
      DIMENSION JNNUM(NNODES)
 
      REAL*8 n1(3),n2(3),n3(3),au(3),r(3),tt(3),z(3)
      REAL*8 norm
 
      n1(:)=CNODES(:,6)
      n2(:)=CNODES(:,2)
      n3(:)=CNODES(:,3)

      if (NNODES.EQ.6) then
        z = n3 - n1
        norm=sqrt(z(1)**2.0 + z(2)**2.0 + z(3)**2.0)
        z = z/norm
        au = n2 - n1
        call crossproduct(z,au,r)
        norm=sqrt(r(1)**2.0 + r(2)**2.0 + r(3)**2.0)
        r = r/norm

      elseif (NNODES.EQ.8) then
        z = n2 - n3
        norm=sqrt(z(1)**2.0 + z(2)**2.0 + z(3)**2.0)
        z = z/norm
        au =  n2 - n1
        call crossproduct(z,au,r)
        norm=sqrt(r(1)**2.0 + r(2)**2.0 + r(3)**2.0)
        r = r/norm
      endif

      call crossproduct(z,r,tt)
      norm=sqrt(tt(1)**2.0 + tt(2)**2.0 + tt(3)**2.0)
      tt = tt/norm

      T(1,1)=r(1); T(1,2)=tt(1); T(1,3)=z(1)
      T(2,1)=r(2); T(2,2)=tt(2); T(2,3)=z(2)
      T(3,1)=r(3); T(3,2)=tt(3); T(3,3)=z(3)

      Rel((NOEL-1)*3+1,:) = T(1,:)
      Rel((NOEL-1)*3+2,:) = T(2,:)
      Rel((NOEL-1)*3+3,:) = T(3,:)

      print*,'hello hello'

      RETURN
      END
C================================================================
C  Elasticity tensor for neo-hook
C================================================================
      SUBROUTINE lr_solver(coef,C,Ce,Fgr,vec,lr)
      !dummy arguments
      real(8),dimension(3,3) :: C,Ce,Fgr,lr
      real(8),dimension(6) :: vec
      real(8) :: coef
      !local arguments
      integer :: i,j,k,l
      real(8),dimension(3,3,3,3) :: special,push,stiff4
      real(8),dimension(3,3) :: A,B,Cgr,ide
      real(8),dimension(3,3) :: inv_C,inv_Ce,inv_Fgr,inv_Cgr,tinv_Fgr
      real(8),dimension(6,6) :: stiff2,invert,aux
      real(8),dimension(6) :: A1,B1,C1,D1,E1,sol
      real(8) :: I1e
!calculate inverses and some tensor of interest
      I1e = Ce(1,1) + Ce(2,2) + Ce(3,3)
      Cgr = MATMUL(transpose(Fgr),Fgr)
      call Inverse33(C,inv_C)
      call Inverse33(Ce,inv_Ce)
      call Inverse33(Cgr,inv_Cgr)
      call Inverse33(Fgr,inv_Fgr)
      tinv_Fgr = transpose(inv_Fgr)
      ide=reshape((/1.d0,0.d0,0.d0,0.d0, 
     1    1.d0,0.d0,0.d0,0.d0,1.d0/), shape(ide)) ! identity matrix
!special tensor product C^-1@C^-1
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          special(i,j,k,l) = 0.5d0*(inv_Ce(i,k)*inv_Ce(j,l)
     1                     + inv_Ce(j,k)*inv_Ce(i,l))
         enddo
        enddo
       enddo
      enddo
! semi push-forward
      call dd393_product(special,inv_Fgr,tinv_Fgr,push)
! stiffness fourth order tensor
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          stiff4(i,j,k,l) = coef*(I1e*inv_Ce(i,j)*inv_C(k,l)/3.d0 
     1 +I1e*push(i,j,k,l)-inv_Ce(i,j)*inv_Cgr(k,l)-ide(i,j)*inv_C(k,l))
         enddo
        enddo
       enddo
      enddo
! turn fourth order tensor to a 6x6 matrix
      call tensor_matrix(stiff4,stiff2)
! calculate of velocity gradient tensor
      call inverse(stiff2,6,invert)
      sol = MATMUL(invert,vec)
      lr(1,1) = sol(1)
      lr(2,2) = sol(2)
      lr(3,3) = sol(3)
      lr(1,2) = sol(4)
      lr(1,3) = sol(5)
      lr(2,3) = sol(6)
      lr(2,1) = lr(1,2)
      lr(3,1) = lr(1,3)
      lr(3,2) = lr(2,3)

      RETURN
      END
C================================================================
C  determinant of a 3*3 matrix
C================================================================
      SUBROUTINE determinant(A,det)

      Real(8), dimension(3,3) :: A
      Real(8) :: det

      det = A(1, 1)*A(2, 2)*A(3, 3)
     1     +A(1, 2)*A(2, 3)*A(3, 1)
     2     +A(1, 3)*A(3, 2)*A(2, 1)
     3     -A(1, 3)*A(3, 1)*A(2, 2)
     4     -A(2, 3)*A(3, 2)*A(1, 1)
     5     -A(1, 2)*A(2, 1)*A(3, 3)

      RETURN
      END
C================================================================
C Cross product (vector product) of two vectrs c = axb
C================================================================
      SUBROUTINE crossproduct(a,b,c)

      real*8 a(3), b(3), c(3)

      c = 0.0
      
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      return
      end
C================================================================
C Tensor Product of two vectors c= a@xb
C================================================================
      SUBROUTINE tensor_product(a,b,c)

      real*8 a(3), b(3), c(3,3)

      do i = 1,3
        do j = 1,3
          c(i,j) = a(i)*b(j)
        enddo
      enddo

      return
      end
C================================================================
C double dot Product of two 3*3 matrix d= A:B
C================================================================
      SUBROUTINE dd33_product(A,B,d)

      Real(8), dimension(3,3) :: A
      Real(8), dimension(3,3) :: B
      real*8 d 

      d = 0.0
      do i = 1,3
        do j = 1,3
          d = d + A(i,j)*B(j,i)
        enddo
      enddo

      return
      end
C================================================================
C double dot Product of two 3*9 matrix C= A:B forth and second order tensors
C================================================================
      SUBROUTINE dd93_product(A,B,C)

      Real(8), dimension(3,3,3,3) :: A
      Real(8), dimension(3,3) :: B, C
      real*8 d 

      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              C(i,j) = C(i,j) + A(i,j,k,l)*B(l,k)
            enddo
          enddo
        enddo
      enddo

      return
      end
C================================================================
C dot Product of three 3*9*3 matrix D= B.A.C forth and second order tensors
C================================================================
      SUBROUTINE dd393_product(A,B,C,E)

      integer :: i,j,k,l,p
      real(8), dimension(3,3,3,3) :: A, D, E
      real(8), dimension(3,3) :: B, C

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          do p = 1,3
            D(i,j,k,l) = D(i,j,k,l) + A(i,j,k,p)*C(p,l)
          enddo
         enddo
        enddo
       enddo
      enddo

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          do p = 1,3
            E(i,j,k,l) = E(i,j,k,l) + B(i,p)*D(p,j,k,l)
          enddo
         enddo
        enddo
       enddo
      enddo

      RETURN
      END
C================================================================
C Inverse of a 3*3 matrix
C================================================================
      SUBROUTINE Inverse33(A,B)

      Real(8) :: D
      Real(8), dimension(3,3) :: A, B

C Determinant of A
      call determinant(A, D)

      if (D.EQ.0.D0) then
        print*,'Error: Determinant is zero for calculation of inverse'
      endif

C inverse of A
      B(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/D
      B(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/D
      B(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))/D
      B(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/D
      B(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/D
      B(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/D
      B(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))/D
      B(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/D
      B(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/D

      RETURN
      END
C================================================================
C  Matrix notation to tensor notation
C================================================================
      SUBROUTINE matrix_tensor(A, B)

      Integer, dimension(2,6) :: indices1,indices2
      Real(8), dimension(6,6) :: A
      Real(8), dimension(3,3,3,3) :: B

      indices1 = reshape((/ 1, 1, 2, 2, 3, 3, 1, 2, 1, 3, 2, 3 /), 
     1           shape(indices1)); ! 2x6
      indices2 = reshape((/ 1, 1, 2, 2, 3, 3, 2, 1, 3, 1, 3, 2 /), 
     1           shape(indices2)); ! 2x6

      B = 0.0
      do m = 1,6
        do n = 1,6
          i1 = indices1(1,m)
          j1 = indices1(2,m)
          k1 = indices1(1,n)
          l1 = indices1(2,n)
          B(i1,j1,k1,l1) = A(m,n) 

          i2 = indices2(1,m)
          j2 = indices2(2,m)
          k2 = indices2(1,n)
          l2 = indices2(2,n)
          B(i2,j2,k2,l2) = A(m,n) 
        enddo
      enddo

      RETURN
      END
C================================================================
C  Fourth order tensor notation to two order tensor
C================================================================
      SUBROUTINE tensor_matrix(A, B)

      Integer, dimension(2,6) :: indices1,indices2
      Real(8), dimension(6,6) :: B
      Real(8), dimension(3,3,3,3) :: A

      indices1 = reshape((/ 1, 1, 2, 2, 3, 3, 1, 2, 1, 3, 2, 3 /), 
     1           shape(indices1)); ! 2x6
      indices2 = reshape((/ 1, 1, 2, 2, 3, 3, 2, 1, 3, 1, 3, 2 /), 
     1           shape(indices2)); ! 2x6

      B = 0.0
      do m = 1,6
        do n = 1,6
          i1 = indices1(1,m)
          j1 = indices1(2,m)
          k1 = indices1(1,n)
          l1 = indices1(2,n)
          B(m,n) = A(i1,j1,k1,l1)

          i2 = indices2(1,m)
          j2 = indices2(2,m)
          k2 = indices2(1,n)
          l2 = indices2(2,n)
          B(m,n) = A(i2,j2,k2,l2)
        enddo
      enddo

      RETURN
      END
C================================================================
C Inverse
C================================================================
      SUBROUTINE inverse(a,n,ainv)
!dummy arguments
      INTEGER :: N
      REAL(8) :: A(N,N)
      REAL(8) :: AINV(N,N)
!local arguments
      INTEGER :: I,J,K,L
      REAL(KIND=8) :: M,MAUM(N,2*N)
!SE CREA LA MATRIZ AUMENTADA
      DO I=1,N
       DO J=1,2*N
        IF (J<=N) THEN
         MAUM(I,J)=A(I,J)
        ELSE IF ((I+N)==J) THEN
         MAUM(I,J)=1.D0
        ELSE
         MAUM(I,J)=0.D0
        END IF
       END DO
      END DO
!SE REDUCE LA MATRIZ A UNA TRIANGULAR superior por eliminacion
      DO K=1,N-1
       DO J=K+1,N
        M=MAUM(J,K)/MAUM(K,K)
        DO I=K,2*N
         MAUM(J,I)=MAUM(J,I)-M*MAUM(K,I)
        END DO
       END DO
      END DO
!elementos de la diagonal iguales a 1
      DO I=1,N
       M=MAUM(I,I)
       DO J=I,2*N
        MAUM(I,J)=MAUM(I,J)/M
       END DO
      END DO
!DIAGONALIZACION
      DO K=N-1,1,-1
       DO I=1,K
        M=MAUM(I,K+1)
        DO J=K,2*N
         MAUM(I,J)=MAUM(I,J)-MAUM(K+1,J)*M
        END DO
       END DO
      END DO
!ALMACENAMIENTO
      DO I=1,N
       DO J=1,N
        AINV(I,J)=MAUM(I,J+N)
       END DO
      END DO

      RETURN
      END
C
