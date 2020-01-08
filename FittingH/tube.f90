module parameters
  implicit none
  real(kind=8),parameter :: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0
  real(kind=8),parameter :: pi=3.141592653589d0
!geometric parameters
  real(kind=8),parameter :: radius=0.010d0,thickness=0.00141d0
!numeric parameters
  integer,parameter :: tstep=20,maxits=100
  real(kind=8),parameter :: step=0.01d0,tol=0.01d0
!load parameters
  real(kind=8),parameter :: tot_press=13332.2d0
end module parameters
!********************************************************
program tube
  use parameters
  implicit none
  integer :: inc,counter
  real(kind=8) :: pressure,time
  real(kind=8) :: stress_r0,stress_r1,stress_r2,stress_r3
  real(kind=8) :: stress_theta0,stress_theta1,stress_theta2,stress_theta3
  real(kind=8) :: stress_z0,stress_z1,stress_z2,stress_z3
  real(kind=8) :: stress_n0,stress_n1,stress_n2
  real(kind=8) :: stress_circ0,stress_circ1,stress_circ2
  real(kind=8) :: stress_long0,stress_long1,stress_long2
  real(kind=8) :: lambda0,lambda1,lambda2,lambda_r
  real(kind=8) :: multiplier0,multiplier1,multiplier2,det
  real(kind=8),dimension(2,2) :: jacob,jacob_inv
!open files for out
  open(unit=10,file='analstrethin.dat',status='replace')
!GET HOMEOSTATIC STATE
    lambda1 = one
    lambda2 = one + step
    multiplier1 = one
    multiplier2 = one + step
!loading history
  do inc=0,tstep   !number of steps
    time = real(inc)/real(tstep)
    pressure = zero + tot_press*time
    !initial roots
    counter = 0
    stress: do !stress cycle
      !calculate stress in the wall tube
      call walltube(pressure,lambda1,stress_n1,stress_circ1,stress_long1)
      call walltube(pressure,lambda2,stress_n2,stress_circ2,stress_long2)
      !calculate stress from anisotropic model
      call stresses(time,lambda1,multiplier1,stress_r1,stress_theta1,stress_z1)
      call stresses(time,lambda2,multiplier1,stress_r2,stress_theta2,stress_z2)
      call stresses(time,lambda1,multiplier2,stress_r3,stress_theta3,stress_z3)
      !jacobian matrix
      jacob(1,1)=(stress_theta2-stress_circ2-stress_theta1+stress_circ1) &
& /(lambda2-lambda1)
      jacob(1,2)=(stress_theta3-stress_circ1-stress_theta1+stress_circ1) &
& /(multiplier2-multiplier1)
      jacob(2,1)=(stress_r2-stress_n2-stress_r1+stress_n1) &
& /(lambda2-lambda1)
      jacob(2,2)=(stress_r3-stress_n1-stress_r1+stress_n1) &
& /(multiplier2-multiplier1)
      !jacobian inverse
      det = jacob(1,1)*jacob(2,2) - jacob(1,2)*jacob(2,1)
      if (det==zero) then
        write(*,*) 'BAD JACOBIAN'
        exit stress
      end if
      jacob_inv(1,1)=jacob(2,2)/det
      jacob_inv(1,2)=-jacob(1,2)/det
      jacob_inv(2,1)=-jacob(2,1)/det
      jacob_inv(2,2)=jacob(1,1)/det
      !looking for new roots
      lambda0 = lambda1 - jacob_inv(1,1)*(stress_theta1-stress_circ1) &
& - jacob_inv(1,2)*(stress_r1-stress_n1)
      multiplier0 = multiplier1 - jacob_inv(2,1)*(stress_theta1-stress_circ1) &
& - jacob_inv(2,2)*(stress_r1-stress_n1)
      call walltube(pressure,lambda0,stress_n0,stress_circ0,stress_long0)
      call stresses(time,lambda0,multiplier0,stress_r0,stress_theta0,stress_z0)
      !checking convergence
      if (abs((lambda0-lambda1)/lambda0)<tol) exit stress
      !updating roots
      lambda2 = lambda1
      lambda1 = lambda0
      multiplier2 = multiplier1
      multiplier1 = multiplier0
      !limit for de number of iterations
      counter = counter + 1
      if (counter>maxits) then
        write(*,*) 'error: many iterations'
        exit stress
      end if
    end do stress !stress cycle
!writing
    lambda_r = one/(lambda0)
    write(10,100) inc,multiplier0,lambda_r,lambda0,stress_r0,stress_theta0,stress_z0
  end do !number of steps
100 format(i3,f12.3,2f9.5,3f12.3)
  close(10)
!  write(*,*) lambda0*1.3237369272263100
!END OF HOMEOSTATIC PART
!*************************
!  call system ('gnuplot plotting.plt')
end program tube
!================================================================
! Calculate stresses with anisotropic model
!================================================================
subroutine stresses(time,lambda,multiplier,stress_r,stress_theta,stress_z)
  use parameters
  implicit none
  !exteral variables
  real(kind=8),intent(in) :: lambda,multiplier,time
  real(kind=8),intent(out) :: stress_z,stress_theta,stress_r
  !internal variables
  integer :: i,j,k,l
  real(kind=8) :: stress_theta_e,stress_theta_c,stress_theta_ci,stress_theta_m
  real(kind=8) :: stress_z_e,stress_z_c,stress_z_ci,stress_z_m
  real(kind=8) :: lambda_r,lambda_c,lambda_m
  real(kind=8) :: alpha
  real(kind=8) :: g_c,g_m
  real(kind=8) :: c10,k1c,k2c,k1m,k2m
  real(kind=8) :: fiber_c(4),fiber_m
  real(kind=8) :: den_e,den_c(4),den_m,den_tot,rho_e,rho_c(4),rho_m
  real(kind=8),dimension(3) :: m0
  real(kind=8),dimension(3,3) :: identity,F,C,Gh_e,F_e,b_e,m0m0
  real(kind=8),dimension(3,3) :: stress_e,stress_c,stress_ci,stress_m
!identity matrix
  do i=1,3
    do j=1,3
      if (i==j) then
        identity(i,j)=one
      else
        identity(i,j)=zero
      endif
    enddo
  enddo
!parameter or constant for modeling hyperelastic behavior and deposition tensor
  c10 = 72.0*1050.0
  k1c = 1136.0*1050.0
  k2c = 11.2
  k1m = 15.2*1050.0
  k2m = 11.4
! densities
  den_e = 241.5d0
  den_c(1) = 65.1d0
  den_c(2) = 65.1d0
  den_c(3) = 260.4d0
  den_c(4) = 260.4d0
  den_m = 157.5d0
  den_tot = den_e + den_c(1) + den_c(2) + den_c(3) + den_c(4) + den_m
  rho_e = den_e/den_tot
  rho_c(1) = den_c(1)/den_tot
  rho_c(2) = den_c(2)/den_tot
  rho_c(3) = den_c(3)/den_tot
  rho_c(4) = den_c(4)/den_tot
  rho_m = den_m/den_tot
!fiber directon for every fiber family
  fiber_c(1) = zero
  fiber_c(2) = pi/two
  fiber_c(3) = pi/four
  fiber_c(4) = -pi/four
  fiber_m = pi/two
!*******************About deformation gradients***********************
!deformation gradient tensor
  F = identity
  C = identity
  F(1,1) = one/lambda
  F(2,2) = lambda
  F(3,3) = one
  !left cauchy-green tensor
  C = matmul(transpose(F),F)
!Elastin
  Gh_e = identity
  Gh_e(3,3) = one + time*0.25d0
  Gh_e(2,2) = one + time*0.34d0 !3237369272263100
  Gh_e(1,1) = one + time*(one/(Gh_e(2,2)*Gh_e(3,3)) - one)
  !total deposition stretch
  F_e = matmul(F,Gh_e)
  !right cauchy-green tensor for elastin
  b_e = matmul(F_e,transpose(F_e))
  !stresses for elastin
  stress_e = rho_e*c10*b_e
!Collagen
  stress_c = zero
  g_c = one + time*0.06d0
  do i=1,4
    m0(1) = zero
    m0(2) = sin(fiber_c(i))
    m0(3) = cos(fiber_c(i))
    do j=1,3
      do k=1,3
        m0m0(j,k) = m0(j)*m0(k)
      enddo
    enddo
    call dd33_product(C,m0m0,lambda_c)
    lambda_c=g_c*sqrt(lambda_c)
    stress_ci = rho_c(i)*k1c*lambda_c**two*(lambda_c**two-one) &
  & *exp(k2c*(lambda_c**two-one)**two)*m0m0
    stress_c = stress_c + stress_ci
  end do
!SMC
  g_m = one + time*0.1d0
  m0(:) = (/zero,sin(fiber_m),cos(fiber_m)/)
  do j=1,3
    do k=1,3
      m0m0(j,k) = m0(j)*m0(k)
    enddo
  enddo
  call dd33_product(C,m0m0,lambda_m)
  lambda_m=g_m*sqrt(lambda_m)
  stress_m = rho_m*k1m*lambda_m**two*(lambda_m**two-one) &
& *exp(k2m*(lambda_m**two-one)**two)*m0m0
!*******************stress in R direction***********************
  stress_r = stress_e(1,1) + multiplier
!*******************stress in THETA direction***********************
  stress_theta = stress_e(2,2) + stress_c(2,2) + stress_m(2,2) + multiplier
!******************stress in Z direction************************
  stress_z = stress_e(3,3) + stress_c(3,3) + stress_m(3,3) + multiplier
!************************************************************
end subroutine stresses
!================================================================
! Calculate stresses in a thin-walled tube
!================================================================
subroutine walltube(pressure,lambda,stress_r,stress_circ,stress_long)
  use parameters
  implicit none
  real(kind=8),intent(in) :: pressure,lambda
  real(kind=8),intent(out) :: stress_r,stress_circ,stress_long
  real(kind=8) :: ri,ro,t,rm
  real(kind=8) :: lambda_t,lambda_r
!geometry
  lambda_t = lambda
  lambda_r = one/lambda
  ri=lambda_t*radius
  t =lambda_r*thickness
  ro=ri+t
  rm=(ro+ri)/two
!stresses in tube
!  stress_r = pressure*ri**2*(1-ro**2/rm**2)/(ro**2-ri**2)
!  stress_circ = pressure*ri**2*(1+ro**2/rm**2)/(ro**2-ri**2)
!  stress_long = pressure*ri/(two*t)
  stress_r = -(pressure/two)!*(ri/ro)
  stress_circ = pressure*ri/t
  stress_long = pressure*ri/(two*t)

end subroutine walltube
!================================================================
! double dot Product of two 3*3 matrix d= A:B
!================================================================
subroutine dd33_product(A,B,d)
  implicit none
  real(kind=8),dimension(3,3),intent(in) :: A,B
  real(kind=8),intent(out) :: d
  integer :: i,j
  
  d = 0.d0
  do i = 1,3
    do j = 1,3
      d = d + A(i,j)*B(j,i)
    enddo
  enddo

end subroutine
