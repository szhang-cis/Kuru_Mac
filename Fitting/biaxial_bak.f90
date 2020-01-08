program biaxial
  integer :: inc,i,j,k
  real(kind=8),parameter :: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0
  real(kind=8) :: pi,step
  real(kind=8) :: stress_theta_e,stress_theta_c,stress_theta_ci,stress_theta_m
  real(kind=8) :: stress_z_e,stress_z_c,stress_z_ci,stress_z_m
  real(kind=8) :: lambda_r,lambda_theta,lambda_z,lambda_c,lambda_m
  real(kind=8) :: g_r,g_theta,g_z,g_c,g_m
  real(kind=8) :: c10,k1c,k2c,k1m,k2m
  real(kind=8) :: fiber_c(4),fiber_m,alpha
  real(kind=8) :: den_e,den_c(4),den_m,den_tot,rho_e,rho_c(4),rho_m
!open files for out
  open(unit=10,file='stress.dat',status='replace')
!important constant
  pi = four*atan(one)
  step = 0.01d0
!parameter or constant for modeling hyperelastic behavior and deposition tensor
  c10 = 72.0e3
  k1c = 1136.0e3
  k2c = 11.2
  k1m = 15.2e3
  k2m = 11.4
  g_z = 1.25d0
  g_theta = 1.34d0
  g_c = 1.06d0
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
  rho_e = den_e/den_tot
  rho_c(1) = den_c(1)/den_tot
  rho_c(2) = den_c(2)/den_tot
  rho_c(3) = den_c(3)/den_tot
  rho_c(4) = den_c(4)/den_tot
  rho_m = den_m/den_tot
!loading history in the test
  do inc=0,30   !number of steps
  lambda_theta = one + step*real(inc)
  lambda_z = one + step*real(inc)
!incompresibility
  lambda_r = one/(lambda_theta*lambda_z)
  g_r = one/(g_theta*g_z)
!fiber directon for every fiber family
  fiber_c(1) = zero
  fiber_c(2) = pi/two
  fiber_c(3) = pi/four
  fiber_c(4) = -pi/four
  fiber_m = pi/two
!stress on elastin in both direction
  stress_theta_e = rho_e*two*c10*(g_theta**two*lambda_theta**two-g_r**two*lambda_r**two)
  stress_z_e = rho_e*two*c10*(g_z**two*lambda_z**two-g_r**two*lambda_r**two)
!stress on collagen in both direction
  stress_theta_c = zero
  stress_z_c = zero
  do i=1,4
   alpha = fiber_c(i)
   lambda_c = g_c*sqrt(lambda_theta**two*sin(alpha)**two+lambda_z**two*cos(alpha)**two)  
   stress_theta_ci = rho_c(i)*k1c*lambda_c**two*(lambda_c**two-one)*exp(k2c*(lambda_c**two-one)**two)*sin(alpha)**two
   stress_z_ci = rho_c(i)*k1c*lambda_c**two*(lambda_c**two-one)*exp(k2c*(lambda_c**two-one)**two)*cos(alpha)**two
   stress_theta_c = stress_theta_c + stress_theta_ci
   stress_z_c = stress_z_c + stress_z_ci
  end do
!stress on smooth muscle cells in both direction
  alpha = fiber_m
  lambda_m = g_m*sqrt(lambda_theta**two*sin(alpha)**two+lambda_z**two*cos(alpha)**two)
  stress_theta_m = rho_m*k1m*lambda_m**two*(lambda_m**two-one)*exp(k2m*(lambda_m**two-one)**two)*sin(alpha)**two
  stress_z_m = rho_m*k1m*lambda_m**two*(lambda_m**two-one)*exp(k2m*(lambda_m**two-one)**two)*cos(alpha)**two
!stress in circumferential direction
  stress_theta = stress_theta_e + stress_theta_c + stress_theta_m
!stress in longitudinal direction
  stress_z = stress_z_e + stress_z_c + stress_z_m  

!writing
  write(10,*) lambda_r,lambda_theta,lambda_z,stress_theta,stress_z
  end do !number of steps
  close(10)
  call system ('gnuplot plotting.plt')
end program biaxial
