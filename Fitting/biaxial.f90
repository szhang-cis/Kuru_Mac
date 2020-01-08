module constantes
  implicit none
  real(kind=8),parameter :: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0
  real(kind=8),parameter :: pi=3.141592653589d0
end module constantes
!********************************************************
program biaxial
  use constantes
  implicit none
  integer :: inc,counter
  real(kind=8) :: step,tol
  real(kind=8) :: stress_theta1,stress_theta2,stress_z1,stress_z2
  real(kind=8) :: lambda_theta1,lambda_theta2,lambda_z,lambda_r
  real(kind=8) :: stress_theta0,stress_z0,lambda_theta0
!open files for out
  open(unit=10,file='stress.dat',status='replace')
!loading history in the test
  step = 0.01d0
  tol = 1000.d0
  do inc=0,200   !number of steps
    lambda_z = one + step*real(inc)
    lambda_theta1 = lambda_z - step
    lambda_theta2 = lambda_z
    counter = 0
    stress: do !stress cycle
      call stresses(lambda_z,lambda_theta1,stress_z1,stress_theta1)
      call stresses(lambda_z,lambda_theta2,stress_z2,stress_theta2)
      lambda_theta0 = lambda_theta2 - (stress_theta2-stress_z2)*(lambda_theta2-lambda_theta1) &
& /((stress_theta2-stress_z2) - (stress_theta1-stress_z1))
      call stresses(lambda_z,lambda_theta0,stress_z0,stress_theta0)
      if (abs(stress_theta0-stress_z0)<tol) exit stress
      lambda_theta2 = lambda_theta0
      counter = counter + 1
      if (counter>1000) then
        write(*,*) 'error counter greather than ', counter
        exit stress
      end if
    end do stress!stress cycle
!writing
    lambda_r = one/(lambda_theta0*lambda_z)
    write(10,*) lambda_r,lambda_theta0,lambda_z,stress_theta0,stress_z0
  end do !number of steps

  close(10)
  call system ('gnuplot plotting.plt')
end program biaxial
!********************************************************************
subroutine stresses(lambda_z,lambda_theta,stress_z,stress_theta)
  use constantes
  implicit none
  !exteral variables
  real(kind=8),intent(in) :: lambda_z,lambda_theta
  real(kind=8),intent(out) :: stress_z,stress_theta
  !internal variables
  integer :: i
  real(kind=8) :: stress_theta_e,stress_theta_c,stress_theta_ci,stress_theta_m
  real(kind=8) :: stress_z_e,stress_z_c,stress_z_ci,stress_z_m
  real(kind=8) :: lambda_r,lambda_c,lambda_m
  real(kind=8) :: alpha
  real(kind=8) :: g_r,g_theta,g_z,g_c,g_m
  real(kind=8) :: c10,k1c,k2c,k1m,k2m
  real(kind=8) :: fiber_c(4),fiber_m
  real(kind=8) :: den_e,den_c(4),den_m,den_tot,rho_e,rho_c(4),rho_m
!parameter or constant for modeling hyperelastic behavior and deposition tensor
  c10 = 100.0e3
  k1c = 100.0e3
  k2c = 0.01
  k1m = 10.0e3
  k2m = 0.01
  g_z = 1.1d0
  g_theta = 1.18d0
  g_r = one/(g_theta*g_z)
  g_c = 1.1d0
  g_m = 1.1d0
! densities
  den_e = 241.5d0
  den_c(1) = 65.1d0
  den_c(2) = 65.1d0
  den_c(3) = 260.4d0
  den_c(4) = 260.4d0
  den_m = 157.5d0
  !den_e = 146.5d0
  !den_c(1) = 79.1d0
  !den_c(2) = 89.1d0
  !den_c(3) = 333.4d0
  !den_c(4) = 339.4d0
  !den_m = 216.5d0
  den_tot = den_e + den_c(1) + den_c(2) + den_c(3) + den_c(4) + den_m
  rho_e = 0.23
  rho_c(1) = 0.62*0.1
  rho_c(2) = 0.62*0.1
  rho_c(3) = 0.62*0.4
  rho_c(4) = 0.62*0.4
  rho_m = 0.15
!fiber directon for every fiber family
  fiber_c(1) = zero
  fiber_c(2) = pi/two
  fiber_c(3) = pi/four
  fiber_c(4) = -pi/four
  fiber_m = pi/two
!incompresibility
  lambda_r = one/(lambda_theta*lambda_z)
!******************stress in Z direction************************
  !elastin
  stress_z_e = rho_e*two*c10*(g_z**two*lambda_z**two-g_r**two*lambda_r**two)
  !collagen
  stress_z_c = zero
  do i=1,4
   alpha = fiber_c(i)
   lambda_c = g_c*sqrt(lambda_theta**two*sin(alpha)**two+lambda_z**two*cos(alpha)**two)  
   stress_z_ci = rho_c(i)*k1c*lambda_c**two*(lambda_c**two-one) &
& *exp(k2c*(lambda_c**two-one)**two)*cos(alpha)**two
   stress_z_c = stress_z_c + stress_z_ci
  end do
  !smc
  alpha = fiber_m
  lambda_m = g_m*sqrt(lambda_theta**two*sin(alpha)**two+lambda_z**two*cos(alpha)**two)
  stress_z_m = rho_m*k1m*lambda_m**two*(lambda_m**two-one) &
& *exp(k2m*(lambda_m**two-one)**two)*cos(alpha)**two
!stress in longitudinal direction
  stress_z = stress_z_e + stress_z_c + stress_z_m  
!*******************stress in THETA direction***********************
  !elstin
  stress_theta_e = rho_e*two*c10*(g_theta**two*lambda_theta**two-g_r**two*lambda_r**two)
  !collagen
  stress_theta_c = zero
  do i=1,4
   alpha = fiber_c(i)
   lambda_c = g_c*sqrt(lambda_theta**two*sin(alpha)**two+lambda_z**two*cos(alpha)**two)  
   stress_theta_ci = rho_c(i)*k1c*lambda_c**two*(lambda_c**two-one) &
& *exp(k2c*(lambda_c**two-one)**two)*sin(alpha)**two
   stress_theta_c = stress_theta_c + stress_theta_ci
  end do
  !smc
  alpha = fiber_m
  lambda_m = g_m*sqrt(lambda_theta**two*sin(alpha)**two+lambda_z**two*cos(alpha)**two)
  stress_theta_m = rho_m*k1m*lambda_m**two*(lambda_m**two-one) &
& *exp(k2m*(lambda_m**two-one)**two)*sin(alpha)**two
!stress in circumferential direction
  stress_theta = stress_theta_e + stress_theta_c + stress_theta_m
!************************************************************
end subroutine stresses
