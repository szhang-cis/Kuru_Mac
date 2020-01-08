# Mathematics libraries
using LinearAlgebra

# Calculate stresses in a thin-walled tube
function ThinWalled(pressure,Variable)
   stretch_t = Variable[1]
   stretch_r = Variable[2]
   R=0.010
   H=0.00141
   #geometry
   r=R #*stretch_t
   h=H #*stretch_r
   #stresses in tube
   stress_t = pressure*r/h
   stress_r = -pressure/2.0
   #stress_z = pressure*r/(2.0*h)
   return stress_t,stress_r
end

# Calculate stresses with anisotropic model
function NeoHookean(Variable)
   # Variable to stretch
   stretch_t = Variable[1]
   stretch_r = Variable[2]
   # parameter of hyperelastic model
   density = 1050.0
   mu = 72.0*density
   kappa = 72.0e3*density
   # avrebation of some parameters
   J = stretch_t*stretch_r
   I1 = stretch_r^2 + stretch_t^2 + 1.0
   # stresses
   stress_t = mu*(stretch_t^2-I1/3.0)*J^(-5.0/3.0) + 2.0*kappa*(J-1.0)
   stress_r = mu*(stretch_r^2-I1/3.0)*J^(-5.0/3.0) + 2.0*kappa*(J-1.0)
   #stress_z = mu*(1.0-I1/3.0)*J^(-5.0/3.0) + 2.0*kappa*(J-1.0)
   return stress_t,stress_r
end

# Compute circumferential residual function
function ResidualFunction(pressure,Variable,fcount)
   #call function to compute the stresses
   ThinStress_t,ThinStress_r = ThinWalled(pressure,Variable)
   HookStress_t,HookStress_r = NeoHookean(Variable)
   #compute the residual
   if fcount==1
      Residual = HookStress_t - ThinStress_t
   else
      Residual = HookStress_r - ThinStress_r
   end
   return Residual
end

# Ridders method for numeric derivatives
function DRidders(pressure,Variable,xcount,fcount)
  """
  Returns the derivative of a function func at a point x by Ridders method of polynomial
  extrapolation. The value h is input as an estimated initial stepsize; it need not be small,
  but rather should be an increment in x over which func changes substantially. An estimate
  of the error in the derivative is returned as err .
  Parameters: Stepsize is decreased by CON at each iteration. Max size of tableau is set by
  NTAB. Return when error is SAFE worse than the best so far.
  """
  h=1.0e-3
  BIG=1.0e30
  NTAB=10
  CON=1.4
  CON2=CON*CON
  SAFE=2.0
  if h==0.0
    println("h must be nonzero in DRidders")
    return
  end
  dfridr=0.0
  a = zeros(Float64,NTAB,NTAB)
  stretch_pos = zeros(Float64,2)
  stretch_neg = zeros(Float64,2)
  hh=h
  for i in 1:2
     stretch_pos[i] = Variable[i]
     stretch_neg[i] = Variable[i]
  end
  stretch_pos[xcount] += hh
  stretch_neg[xcount] -= hh
  func_pos = ResidualFunction(pressure,stretch_pos,fcount)
  func_neg = ResidualFunction(pressure,stretch_neg,fcount)
  a[1,1] = (func_pos-func_neg)/(2.0*hh)
  err=BIG
  # Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation.
  for i in 2:NTAB
    hh=hh/CON
    # Try new, smaller stepsize.
    for i in 1:2
       stretch_pos[i] = Variable[i]
       stretch_neg[i] = Variable[i]
    end
    stretch_pos[xcount] += hh
    stretch_neg[xcount] -= hh
    func_pos = ResidualFunction(pressure,stretch_pos,fcount)
    func_neg = ResidualFunction(pressure,stretch_neg,fcount)
    a[1,i] = (func_pos-func_neg)/(2.0*hh)
    fac=CON2
    # Compute extrapolations of various orders, requiring no new function evaluations.
    for j in 2:i
      a[j,i] = (a[j-1,i]*fac-a[j-1,i-1])/(fac-1.0)
      fac=CON2*fac
      errt=max(abs(a[j,i]-a[j-1,i]),abs(a[j,i]-a[j-1,i-1]))
      #The error strategy is to compare each new extrapolation to one order lower, both at
      #the present stepsize and the previous one.
      #If error is decreased, save the improved answer.
      if errt<=err
        err=errt
        dfridr=a[j,i]
      end
      #If higher order is worse by a significant factor SAFE, then quit early.
      if abs(a[i,i]-a[i-1,i-1])>=SAFE*err
        #print('Early quit in df_ridders function')
        return dfridr
      end
    end
  end
  return dfridr
end

# Evaluation of the function
function FunctionEval(pressure,Variable)
  """
  Function with the aim to get the values of the function and derivative at x.
  """
  func = zeros(Float64,2)
  df = zeros(Float64,2,2)
  for fcount in 1:2
     func[fcount] = ResidualFunction(pressure,Variable,fcount)
     for xcount in 1:2
        df[fcount,xcount] = DRidders(pressure,Variable,xcount,fcount)
     end
  end
  return func,df
end

# Newton-Raphson method
function NewtonRaphson(pressure,Variable)
  """
  Using the Newton-Raphson method, find the root of a function known to lie close to x. 
  The root rtnewt will be refined until its accuracy is known within +/- xacc. funcd
  is a user-supplied subroutine that returns both the function value and the first derivative
  of the function at the point x.
  """
  # Set to maximum number of iterations.
  JMAX=20
  #res = zeros(Float64,2)
  res,_ = FunctionEval(pressure,Variable)
  res_norm = norm(res)
  # Initial guess.
  rtnewt=Variable
  for j in 1:JMAX
    dx = zeros(Float64,2)
    func,df = FunctionEval(pressure,rtnewt)
    df_inv = inv(df)
    mul!(dx,df_inv,func) #func/df
    rtnewt -= dx
    func_norm = norm(func)
    # Convergence.
    if abs(func_norm/res_norm)<1.e-5
      #print('Iterations: '+str(j))
      return rtnewt
    end
  end
  println("newton_raphson exceeded maximum iterations")
  return rtnewt
end

# Test of the Solution method with derivatives from Ridders method
number_of_increments=1
# initial guess
error=1.0
Variable = [1.1,0.9]
# incremental solution
for Inc in 1:number_of_increments
   global Variable
   #global residual_norm
   factor=Inc/number_of_increments
   pressure=13.332e3*factor
   solution = NewtonRaphson(pressure,Variable)
   println("Increment $Inc, Solution= $solution")
end

