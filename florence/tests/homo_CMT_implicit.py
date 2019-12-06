# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
from numpy import einsum
# Build a path for python to Florence
sys.path.append(os.path.join(os.path.expanduser("~"),"florence"))
#import Florence
from Florence import *

#========================== HIPERELASTIC MODEL FUNCTIONS  ======================================
def Directions(mesh,LagrangeGaussCoords):
    """
        Routine dedicated to compute the fibre direction of components in integration point for 
        the Material in Florence and for the auxiliar routines in this script.
    """
    # Array for fibres orientations ['thick','smc','co1','co2','co3','co4']
    ndim = mesh.InferSpatialDimension()
    ngauss = LagrangeGaussCoords.shape[1]
    # Geometric definitions per element
    tangential = np.zeros((mesh.nelem,ngauss,ndim),dtype=np.float64)
    directrix = [0.,1.,0.]
    # Loop throught the element in the mesh
    fibre_direction = np.zeros((6,mesh.nelem,ngauss,ndim),dtype=np.float64)
    for elem in range(mesh.nelem):
        for gcounter in range(ngauss):
            # Geometric definitions per element
            tangential[elem,gcounter,:] = np.cross(directrix,LagrangeGaussCoords[elem,gcounter,:])
            tangential[elem,gcounter,:] = tangential[elem,gcounter,:]/\
                    np.linalg.norm(tangential[elem,gcounter,:])
            fibre_direction[0][elem,gcounter,:] = np.cross(tangential[elem,gcounter,:],directrix)
            fibre_direction[0][elem,gcounter,:] = fibre_direction[0][elem,gcounter,:]/\
                    np.linalg.norm(fibre_direction[0][elem,gcounter,:])
            # Define the anisotropic orientations
            fibre_direction[1][elem,gcounter,:]=tangential[elem,gcounter,:]
            fibre_direction[2][elem,gcounter,:]=tangential[elem,gcounter,:]
            fibre_direction[3][elem,gcounter,:]=np.multiply(directrix,np.cos(np.pi/4)) + \
                    np.multiply(tangential[elem,gcounter,:],np.sin(np.pi/4))
            fibre_direction[4][elem,gcounter,:]=np.multiply(directrix,np.cos(np.pi/4)) - \
                    np.multiply(tangential[elem,gcounter,:],np.sin(np.pi/4))
            fibre_direction[5][elem,gcounter,:]=directrix

    return fibre_direction

def ComputeStress(material,elem,gcounter):
    """
        The purpose of this routine is to compute the Stress by component in the mixture from 
        the Total deformation gradient. This Total deformation gradient come from the equilibrium
        solved in Florence.
    """

    I = np.eye(3,3,dtype=np.float64)
    F = material.TotalDeformation['F'][elem][gcounter]
    J = material.TotalDeformation['J'][elem][gcounter]
    # Directional vector for element
    Axial = material.anisotropic_orientations[5][elem][gcounter][:,None]
    Axial = np.dot(I,Axial)[:,0]
    Tangential = material.anisotropic_orientations[2][elem][gcounter][:,None]
    Tangential = np.dot(I,Tangential)[:,0]
    Normal = material.anisotropic_orientations[0][elem][gcounter][:,None]
    Normal = np.dot(I,Normal)[:,0]
    Rotation = np.eye(3,3,dtype=np.float64)
    for i in range(3):
        Rotation[0,i] = Normal[i]
        Rotation[1,i] = Tangential[i]
        Rotation[2,i] = Axial[i]

    # Total growth gradient deformation
    outerNormal = einsum('i,j',Normal,Normal)
    outerTangential = I - outerNormal
    # Total growth deformation gradient
    F_g = material.Growth[elem][gcounter]*outerNormal + outerTangential
    # Total inelastical deformation gradient
    F_gr_inv = np.linalg.inv(F_g)

    #ELASTIN
    kappa = material.kappa*material.GrowthRemodeling[0][elem][gcounter]
    mu3D = material.mu3D*material.GrowthRemodeling[0][elem][gcounter]
    Gh_ela = material.deposition_stretch['ela'][elem][gcounter]
    Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
    F_ela = np.dot(F,Gh_ela)
    F_ela_e = np.dot(F_ela,F_gr_inv)
    J_ela_e = np.linalg.det(F_ela_e)
    b_ela_e = np.dot(F_ela_e,F_ela_e.T)

    if material.ndim == 3:
        trb_ela_e = trace(b_ela_e)
    elif material.ndim == 2:
        trb_ela_e = trace(b_ela_e) + 1.

    #SMC AND COLLAGEN FIBRES
    Stress_i = np.zeros((5),dtype=np.float64)
    Softness = np.zeros((5),dtype=np.float64)
    for idx in [1,2,3,4,5]:
        fibre_stress = 0.0
        # Fibre direction
        N = material.anisotropic_orientations[idx][elem][gcounter][:,None]
        N = np.dot(I,N)[:,0]
        FN = np.dot(F,N)
        if idx is 1:
            FN = np.dot(material.deposition_stretch['smc'],FN)
            c1 = material.c1m*material.GrowthRemodeling[1][elem][gcounter]
            c2 = material.c2m
        elif idx is not 1:
            FN = np.dot(material.deposition_stretch['col'],FN)
            c1 = material.c1c*material.GrowthRemodeling[idx][elem][gcounter]
            c2 = material.c2c

        # Tensor of reference direction
        outerN = einsum('i,j',N,N)
        # TOTAL deformation
        innerFN = einsum('i,i',FN,FN)
        outerFN = einsum('i,j',FN,FN)
        # Remodeling stretch
        lambda_r = material.GrowthRemodeling[idx+5][elem][gcounter]
        # Elastic deformation
        innerN_e = innerFN/lambda_r**2
        outerN_e = outerFN/lambda_r**2
        # Passive Stress for this fibre
        expo = np.exp(c2*(innerN_e-1.)**2.)
        fibre_stress = c1*(innerN_e-1.)*expo*outerN_e/J
        # Active stress for SMC
        if idx is 1:
            den0_tot = 1050.0
            s_act = 54.0e3*material.GrowthRemodeling[1][elem][gcounter]
            stretch_m = 1.4
            stretch_a = 1.0
            stretch_0 = 0.8
            active_stress = (s_act/(den0_tot*innerFN))*\
                    (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2.)*outerFN/J
            fibre_stress += active_stress

        Stress_i[idx-1] = np.tensordot(fibre_stress,outerN)
        # Fibre softness for remodeling
        Stiffness = innerN_e*(innerN_e-1.)**2.
        Stiffness = (4.*c2*Stiffness+4.*innerN_e-2.)*c1
        Stiffness = Stiffness*np.sqrt(innerN_e)
        Stiffness = Stiffness*expo
        Softness[idx-1] = np.sqrt(innerFN)/(innerN_e*Stiffness)

    return Stress_i,Softness

#===================== ELASTIN DEPOSTION STRETCH FUNCTIONS ===============================
def ThinWalled_Residual(material,R_gauss,elem,gcounter,fcounter):
    """ 
        Thin walled tube formulas for residual to compute deposition stretch
    """
    pressure = 13.3322e3
    R = 0.010
    H = 0.00141
    # some important parameters
    kappa = material.kappa*material.GrowthRemodeling[0][elem][gcounter]
    mu = material.mu3D*material.GrowthRemodeling[0][elem][gcounter]
    c1m = material.c1m*material.GrowthRemodeling[1][elem][gcounter]
    c2m = material.c2m
    c1c1 = material.c1c*material.GrowthRemodeling[2][elem][gcounter]
    c1c2 = material.c1c*material.GrowthRemodeling[3][elem][gcounter]
    c1c3 = material.c1c*material.GrowthRemodeling[4][elem][gcounter]
    c2c = material.c2c
    # deposition stretches
    stretch_r = material.deposition_stretch['ela'][elem][gcounter][0,0]
    stretch_t = material.deposition_stretch['ela'][elem][gcounter][1,1]
    stretch_z = material.deposition_stretch['ela'][elem][gcounter][2,2]
    stretch_m = material.deposition_stretch['smc']
    stretch_c = material.deposition_stretch['col']
    # active smooth muscle cells
    den0_tot = 1050.0
    s_act = 54.0e3*material.GrowthRemodeling[1][elem][gcounter]
    # compute the residual
    if fcounter is 0:
        Residual = pressure*(0.5-(R_gauss-R)/H) + \
        2.*kappa*stretch_r*stretch_t*stretch_z*(stretch_r*stretch_t*stretch_z-1.) + \
        mu*(stretch_r*stretch_t*stretch_z)**(-2./3.)*(stretch_r**2-(stretch_r**2+stretch_t**2+stretch_z**2)/3.)
    else:
        Residual = 2.*kappa*stretch_r*stretch_t*stretch_z*(stretch_r*stretch_t*stretch_z-1.) + \
        mu*(stretch_r*stretch_t*stretch_z)**(-2./3.)*(stretch_t**2-(stretch_r**2+stretch_t**2+stretch_z**2)/3.)+\
        (s_act/den0_tot)*(1.-((1.4-1.0)/(1.4-0.8))**2.) + \
        c1m*(stretch_m**2-1.)*np.exp(c2m*(stretch_m**2-1.)**2)*stretch_m**2 + \
        (c1c1+c1c2/2.+c1c3/2.)*(stretch_c**2-1.)*np.exp(c2c*(stretch_c**2-1.)**2)*stretch_c**2 - \
        pressure*R/H

    return Residual

def ThickWalled_Residual(material,R_gauss,elem,gcounter,fcounter):
    """ 
        Thick walled tube formulas for residual to compute deposition stretch
    """
    pressure = 13.3322e3
    R = 0.010
    H = 0.00141
    Ri = R-0.5*H
    Ro = R+0.5*H
    # some important parameters
    kappa = material.kappa*material.GrowthRemodeling[0][elem][gcounter]
    mu = material.mu3D*material.GrowthRemodeling[0][elem][gcounter]
    c1m = material.c1m*material.GrowthRemodeling[1][elem][gcounter]
    c2m = material.c2m
    c1c1 = material.c1c*material.GrowthRemodeling[2][elem][gcounter]
    c1c2 = material.c1c*material.GrowthRemodeling[3][elem][gcounter]
    c1c3 = material.c1c*material.GrowthRemodeling[4][elem][gcounter]
    c2c = material.c2c
    # deposition stretches
    stretch_r = material.deposition_stretch['ela'][elem][gcounter][0,0]
    stretch_t = material.deposition_stretch['ela'][elem][gcounter][1,1]
    stretch_z = material.deposition_stretch['ela'][elem][gcounter][2,2]
    stretch_m = material.deposition_stretch['smc']
    stretch_c = material.deposition_stretch['col']
    # active smooth muscle cells
    den0_tot = 1050.0
    s_act = 54.0e3*material.GrowthRemodeling[1][elem][gcounter]
    # compute the residual
    Residual = np.zeros(2,dtype=np.float64)
    if fcounter is 0:
        Residual[0] = pressure*Ri**2*(1.-Ro**2/R_gauss**2)/(Ro**2-Ri**2) + \
        2.*kappa*stretch_r*stretch_t*stretch_z*(stretch_r*stretch_t*stretch_z-1.) + \
        mu*(stretch_r*stretch_t*stretch_z)**(-2./3.)*(stretch_r**2-(stretch_r**2+stretch_t**2+stretch_z**2)/3.)
    else:
        Residual[1] = 2.*kappa*stretch_r*stretch_t*stretch_z*(stretch_r*stretch_t*stretch_z-1.) + \
        mu*(stretch_r*stretch_t*stretch_z)**(-2./3.)*(stretch_t**2-(stretch_r**2+stretch_t**2+stretch_z**2)/3.)+\
        (s_act/den0_tot)*(1.-((1.4-1.0)/(1.4-0.8))**2.) + \
        c1m*(stretch_m**2-1.)*np.exp(c2m*(stretch_m**2-1.)**2)*stretch_m**2 + \
        (c1c1+c1c2/2.+c1c3/2.)*(stretch_c**2-1.)*np.exp(c2c*(stretch_c**2-1.)**2)*stretch_c**2 - \
        pressure*Ri**2*(1.+Ro**2/R_gauss**2)/(Ro**2-Ri**2)

    return Residual

def DStretch_Ridders(material,R_gauss,elem,gcounter,xcount,fcount):
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
    if h==0.0:
        print('h must be nonzero in df_ridders')
        return
    a = np.zeros((NTAB,NTAB),dtype=np.float64)
    hh=h
    # store the value of deposition stretch in stretch
    stretch = material.deposition_stretch['ela'][elem][gcounter][xcount,xcount]
    # get the value of the residual function on this deposition stretch
    material.deposition_stretch['ela'][elem][gcounter][xcount,xcount] = stretch + hh
    func_pos = ThinWalled_Residual(material=material,R_gauss=R_gauss,elem=elem,gcounter=gcounter,fcounter=fcount)
    material.deposition_stretch['ela'][elem][gcounter][xcount,xcount] = stretch - hh
    func_neg = ThinWalled_Residual(material=material,R_gauss=R_gauss,elem=elem,gcounter=gcounter,fcounter=fcount)
    a[0,0] = (func_pos-func_neg)/(2.0*hh)
    err=BIG
    # Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation.
    for i in range(1,NTAB):
        hh=hh/CON
        # Try new, smaller stepsize.
        material.deposition_stretch['ela'][elem][gcounter][xcount,xcount] = stretch + hh
        func_pos = ThinWalled_Residual(material=material,R_gauss=R_gauss,
                    elem=elem,gcounter=gcounter,fcounter=fcount)
        material.deposition_stretch['ela'][elem][gcounter][xcount,xcount] = stretch - hh
        func_neg = ThinWalled_Residual(material=material,R_gauss=R_gauss,
                    elem=elem,gcounter=gcounter,fcounter=fcount)
        a[0,i] = (func_pos-func_neg)/(2.0*hh)
        material.deposition_stretch['ela'][elem][gcounter][xcount,xcount] = stretch
        fac=CON2
        # Compute extrapolations of various orders, requiring no new function evaluations.
        for j in range(1,i+1):
            a[j,i] = (a[j-1,i]*fac-a[j-1,i-1])/(fac-1.0)
            fac=CON2*fac
            errt=max(abs(a[j,i]-a[j-1,i]),abs(a[j,i]-a[j-1,i-1]))
            #The error strategy is to compare each new extrapolation to one order lower, both at
            #the present stepsize and the previous one.
            #If error is decreased, save the improved answer.
            if (errt<=err):
                err=errt
                dfridr=a[j,i]
        #If higher order is worse by a significant factor SAFE, then quit early.
        if (abs(a[i,i]-a[i-1,i-1])>=SAFE*err):
            #print('Early quit in df_ridders function '+str(errt))
            return dfridr,err
    #print('Not enough size of the tableau in df_ridders')
    return dfridr,err

def StretchResidual_eval(material,R_gauss,elem,gcounter):
    """
    Return the value of the residual function on the point and its derivative
    """
    func = np.zeros(2,dtype=np.float64)
    df = np.zeros((2,2),dtype=np.float64)
    err = np.zeros((2,2),dtype=np.float64)
    for j in range(2):
        func[j] = ThinWalled_Residual(material=material,R_gauss=R_gauss,elem=elem,gcounter=gcounter,fcounter=j)
        for i in range(2):
            df[j,i],err[j,i] = DStretch_Ridders(material=material,R_gauss=R_gauss,
                                elem=elem,gcounter=gcounter,xcount=i,fcount=j)

    return func,df,err

def Stretch_Newton(material,R_gauss,elem,gcounter):
    """
    Using the Newton-Raphson method, find the root of a function known to lie close to x. 
    The root rtnewt will be refined until its accuracy is known within +/- xacc. funcd
    is a user-supplied subroutine that returns both the function value and the first derivative
    of the function at the point x.
    """
    # Set to maximum number of iterations.
    ITER=0
    Residual = np.zeros(2,dtype=np.float64)
    for j in range(2):
        Residual[j]=ThinWalled_Residual(material=material,R_gauss=R_gauss,elem=elem,gcounter=gcounter,fcounter=j)
    norm_res = np.linalg.norm(Residual)
    if np.isclose(norm_res,0.0):
        norm_res = 1.0e-14
    norm_residual = np.absolute(np.linalg.norm(Residual)/norm_res)
    while norm_residual>1.0e-5 or ITER==0:
        Residual,Derivative,err = StretchResidual_eval(material=material,R_gauss=R_gauss,
                                elem=elem,gcounter=gcounter)
        Derivative_inv = np.linalg.inv(Derivative)
        dx = np.dot(Derivative_inv,Residual)
        material.deposition_stretch['ela'][elem][gcounter][0,0] -= dx[0]
        material.deposition_stretch['ela'][elem][gcounter][1,1] -= dx[1]
        # Convergence.
        ITER += 1
        norm_residual = np.absolute(np.linalg.norm(Residual)/norm_res)
        x = np.zeros(2,dtype=np.float64)
        x[0] = material.deposition_stretch['ela'][elem][gcounter][0,0]
        x[1] = material.deposition_stretch['ela'][elem][gcounter][1,1]
        norm_solution = np.absolute(np.linalg.norm(dx)/np.linalg.norm(x))
        if norm_solution<1.0e-5:
            break
        if np.isnan(norm_residual) or norm_residual>1.0e6 or ITER>20:
            print('newton_raphson diverge, elem= '+str(elem)+' gcounter= '+str(gcounter))
            print('Norm='+str(norm_residual)+' Iter='+str(ITER))
            return True
    return False

def ElastinDepositionStretch0(material,R_gauss,nelem,NoGauss):
    """ 
        This routine aim to get the radial deposition stretch from the boundary conditions
        in radial direction, as stress=-p in inner wall and stress=0 in outer wall. The
        relation for the stresses along the thickness come from thin-walled tube. So radial
        deposition stretch is depending of the position (radius).
    """
    print('START: ElastinDepositionStretch0')
    for elem in range(nelem):
        for gcounter in range(NoGauss):
            Radius = R_gauss[elem][gcounter]
            Newton_Convergence = Stretch_Newton(material=material,R_gauss=Radius,elem=elem,gcounter=gcounter)
            if Newton_Convergence:
                break
        if Newton_Convergence:
            print('ERROR: ElastinDepositionStretch0 did not converge')
            break
    print('END: ElastinDepositionStretch0')

def ElastinDepositionStretch1(material,nelem,NoGauss):
    """ 
        This routine aim to get the radial deposition stretch from the boundary conditions 
        in radial direction, as stress=-p in inner wall and stress=0 in outer wall. The 
        relation for the stresses along the thickness come from thin-walled tube. So radial
        deposition stretch is depending of the position (radius).
    """
    I = np.eye(3,3,dtype=np.float64)        
    print('==> Looking for Elastin Deposition Stretch')
    for elem in range(nelem):
        for gcounter in range(NoGauss):
            # Directional vector for element
            Axial = material.anisotropic_orientations[5][elem][gcounter][:,None]
            Axial = np.dot(I,Axial)[:,0]
            Tangential = material.anisotropic_orientations[2][elem][gcounter][:,None]
            Tangential = np.dot(I,Tangential)[:,0]
            Normal = material.anisotropic_orientations[0][elem][gcounter][:,None]
            Normal = np.dot(I,Normal)[:,0]
            Rotation = np.eye(3,3,dtype=np.float64)
            for i in range(3):
                Rotation[0,i] = Normal[i]
                Rotation[1,i] = Tangential[i]
                Rotation[2,i] = Axial[i]

            # deposition stretches
            F = material.TotalDeformation['F'][elem][gcounter]
            #J = material.TotalDeformation['J'][elem][gcounter]
            Gh_ela = material.deposition_stretch['ela'][elem][gcounter]
            F = np.dot(Rotation,np.dot(F,Rotation.T))
            F_ela = np.dot(F,Gh_ela)
            material.deposition_stretch['ela'][elem][gcounter][0,0] = F_ela[0,0]
            #material.deposition_stretch['ela'][elem][gcounter][1,1] = F_ela[1,1]
            #material.deposition_stretch['ela'][elem][gcounter][2,2] = F_ela[2,2]

    print('==> Elastin Deposition Stretch calculated')

#================== GROWTH AND REMODELING FUNCTIONS   ===================================
def HomeostaticStress_(material,nelem,NoGauss):
    """
    Routine to get the Homeostatic Stress for the reference for Growth and Remodeling
    """
    # Array for homeostatic stress ['smc','co1','co2','co3','co4']
    Stress_H = np.zeros((5,nelem,NoGauss),dtype=np.float64)
    for elem in range(nelem):
        for gcounter in range(NoGauss):
            Stress_i,_ = ComputeStress(material=material, elem=elem, gcounter=gcounter)
            #SMC AND COLLAGEN FIBRES
            for idx in range(5):
                Stress_H[idx][elem][gcounter] = Stress_i[idx]

    return Stress_H

def GrowthRemodeling_Rates(time,material,Stress_H,Y_pos,elem,gcounter):
    """
    This is the Growth and Remodeling consitutive equations.
    """
    #['smc','co1','co2','co3','co4']:
    GR_rate = np.zeros((11),dtype=np.float64)
    Stress_i,Softness = ComputeStress(material=material,elem=elem,gcounter=gcounter)
    # Elastin degradation [days]
    L_dam = 0.010
    t_dam = 40.0
    D_max = 0.5
    T_ela = 101.0*365.25
    den0_e = 241.5
    # elastin density rate
    exponential = np.exp(-0.5*(Y_pos/L_dam)**2 - time/t_dam)
    GR_rate[0]=-material.GrowthRemodeling[0][elem][gcounter]/T_ela-(D_max/t_dam)*np.multiply(exponential,den0_e)
    # Fribres turnover and remodeling
    turnover = 101.0
    gain = 0.05/turnover
    #SMC AND COLLAGEN FIBRES
    for idx in [1,2,3,4,5]:
        # fibres density rates
        DeltaStress = Stress_i[idx-1] - Stress_H[idx-1]
        GR_rate[idx] = gain*np.multiply(material.GrowthRemodeling[idx][elem][gcounter],np.divide(DeltaStress,Stress_H[idx-1]))
        # fibres remodeling rates
        lambda_r_dot = np.divide(GR_rate[idx],material.GrowthRemodeling[idx][elem][gcounter]) + 1.0/turnover
        lambda_r_dot = np.multiply(lambda_r_dot,DeltaStress)
        GR_rate[idx+5] = np.multiply(lambda_r_dot,Softness[idx-1])

    return GR_rate

def GrowthRemodeling_Residual(time,Delta_t,material,Stress_H,Y_pos,elem,gcounter,fcount):
    """ 
    Growth and Remodeling implicit formulas for residual to compute
    """
    # Store some total initial density and elastin
    den0_tot = 1050.0
    # Total Growth magnitud
    den_tot = 0.0
    for idx in range(6):
        den_tot += material.GrowthRemodeling[idx][elem][gcounter]
    material.Growth[elem][gcounter] = den_tot/den0_tot
    GR_rate = GrowthRemodeling_Rates(time=time,material=material,Stress_H=Stress_H,Y_pos=Y_pos,
        elem=elem,gcounter=gcounter)
    # growth and remodeling parameters (densities and remodeling)
    variable0 = material.GrowthRemodeling0[fcount][elem][gcounter]
    variable = material.GrowthRemodeling[fcount][elem][gcounter]
    # compute the residual
    Residual = -variable + variable0 + Delta_t*GR_rate[fcount]
    return Residual

def DGR_Ridders(time,Delta_t,material,Stress_H,Y_pos,elem,gcounter,xcount,fcount):
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
    if h==0.0:
        print('h must be nonzero in df_ridders')
        return
    a = np.zeros((NTAB,NTAB),dtype=np.float64)
    # store the value of deposition stretch in stretch
    parameter = material.GrowthRemodeling[xcount][elem][gcounter]
    hh=h
    # get the value of the residual function on this deposition stretch
    material.GrowthRemodeling[xcount][elem][gcounter] = parameter + hh
    func_pos = GrowthRemodeling_Residual(time=time,Delta_t=Delta_t,material=material,Stress_H=Stress_H,
        Y_pos=Y_pos,elem=elem,gcounter=gcounter,fcount=fcount)
    material.GrowthRemodeling[xcount][elem][gcounter] = parameter - hh
    func_neg = GrowthRemodeling_Residual(time=time,Delta_t=Delta_t,material=material,Stress_H=Stress_H,
        Y_pos=Y_pos,elem=elem,gcounter=gcounter,fcount=fcount)
    a[0,0] = (func_pos-func_neg)/(2.0*hh)
    err=BIG
    # Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation.
    for i in range(1,NTAB):
        hh=hh/CON
        # Try new, smaller stepsize.
        material.GrowthRemodeling[xcount][elem][gcounter] = parameter + hh
        func_pos = GrowthRemodeling_Residual(time=time,Delta_t=Delta_t,material=material,Stress_H=Stress_H,
            Y_pos=Y_pos,elem=elem,gcounter=gcounter,fcount=fcount)
        material.GrowthRemodeling[xcount][elem][gcounter] = parameter - hh
        func_neg = GrowthRemodeling_Residual(time=time,Delta_t=Delta_t,material=material,Stress_H=Stress_H,
            Y_pos=Y_pos,elem=elem,gcounter=gcounter,fcount=fcount)
        a[0,i] = (func_pos-func_neg)/(2.0*hh)
        material.GrowthRemodeling[xcount][elem][gcounter] = parameter
        fac=CON2
        # Compute extrapolations of various orders, requiring no new function evaluations.
        for j in range(1,i+1):
            a[j,i] = (a[j-1,i]*fac-a[j-1,i-1])/(fac-1.0)
            fac=CON2*fac
            errt=max(abs(a[j,i]-a[j-1,i]),abs(a[j,i]-a[j-1,i-1]))
            #The error strategy is to compare each new extrapolation to one order lower, both at
            #the present stepsize and the previous one.
            #If error is decreased, save the improved answer.
            if (errt<=err):
                err=errt
                dfridr=a[j,i]
        #If higher order is worse by a significant factor SAFE, then quit early.
        if (abs(a[i,i]-a[i-1,i-1])>=SAFE*err):
            #print('Early quit in df_ridders function '+str(errt))
            return dfridr,err
    #print('Not enough size of the tableau in df_ridders')
    return dfridr,err

def GrowthRemodelingResidual_eval(time,Delta_t,material,Stress_H,Y,elem,gcounter):
    """
    Return the value of the residual function on the point and its derivative
    """
    func = np.zeros(11,dtype=np.float64)
    df = np.zeros((11,11),dtype=np.float64)
    err = np.zeros((11,11),dtype=np.float64)
    for j in range(11):
        func[j]=GrowthRemodeling_Residual(time=time,Delta_t=Delta_t,material=material,Stress_H=Stress_H,
            Y_pos=Y,elem=elem,gcounter=gcounter,fcount=j)
        for i in range(11):
            df[j,i],err[j,i] = DGR_Ridders(time=time,Delta_t=Delta_t,material=material,Stress_H=Stress_H,
                Y_pos=Y,elem=elem,gcounter=gcounter,xcount=i,fcount=j)
    return func,df,err

def GrowthRemodeling_Newton(time,Delta_t,material,Stress_H,Y_,elem,gcounter):
    """
    Using the Newton-Raphson method, find the root of a function known to lie close to x. 
    The root rtnewt will be refined until its accuracy is known within +/- xacc. funcd
    is a user-supplied subroutine that returns both the function value and the first derivative
    of the function at the point x.
    """
    # Set to maximum number of iterations.
    ITER=0
    Residual = np.zeros(11,dtype=np.float64)
    for j in range(11):
        Residual[j]=GrowthRemodeling_Residual(time=time,Delta_t=Delta_t,material=material,Stress_H=Stress_H,
            Y_pos=Y_,elem=elem,gcounter=gcounter,fcount=j)
    norm_res = np.linalg.norm(Residual)
    if np.isclose(norm_res,0.0):
        norm_res = 1.0e-14
    norm_residual = np.absolute(np.linalg.norm(Residual)/norm_res)
    while norm_residual>1.0e-5 or ITER==0:
        Residual,Derivative,err = GrowthRemodelingResidual_eval(time=time,Delta_t=Delta_t,material=material,
            Stress_H=Stress_H,Y=Y_,elem=elem,gcounter=gcounter)
        determinant = np.linalg.det(Derivative)
        if np.isclose(determinant,0.0):
            print(Derivative)
            print(determinant)
            print('The Jacobian matrix for Newton-Raphson is singular')
            return True
        #print(Derivative)
        Derivative_inv = np.linalg.inv(Derivative)
        dx = np.dot(Derivative_inv,Residual)
        for i in range(11):
            material.GrowthRemodeling[i][elem][gcounter] -= dx[i]
        # Convergence.
        ITER += 1
        norm_residual = np.absolute(np.linalg.norm(Residual)/norm_res)
        x = np.zeros(11,dtype=np.float64)
        for i in range(11):
            x[i] = material.GrowthRemodeling[i][elem][gcounter]
        norm_solution = np.absolute(np.linalg.norm(dx)/np.linalg.norm(x))
        if norm_solution<1.0e-5:
            break
        if np.isnan(norm_residual) or norm_residual>1.0e6 or ITER>20:
            print('newton_raphson diverge, elem= '+str(elem)+' gcounter= '+str(gcounter))
            print('Norm='+str(norm_residual)+' Iter='+str(ITER))
            return True
    return False

def GR_BackwardEuler(Iter,time,Delta_t,material,HomeostaticStress,Y,nelem,NoGauss):
    """
        This is the Growth and Remodeling solver for implicit time integration.
        Backward Euler.
    """
    # Iterative schema to determine Remodeling and Densities
    Newton_failed = False
    print('START: GR_BackwardEuler')
    for elem in range(nelem):
        for gcounter in range(NoGauss):
            Stress_H = np.zeros((5),dtype=np.float64)
            for idx in range(5):
                Stress_H[idx] = HomeostaticStress[idx][elem][gcounter]
            Y_ = Y[elem][gcounter]
            Newton_failed = GrowthRemodeling_Newton(time=time,Delta_t=Delta_t,material=material,
                Stress_H=Stress_H,Y_=Y_,elem=elem,gcounter=gcounter)
            if Newton_failed:
                print('Newton method failed at Gauss point: '+str(gcounter))
                break
        if Newton_failed:
            print('Newton method failed at element: '+str(elem))
            print('ERROR: GR_BackwardEuler failed')
            return True
    print('END: GR_BackwardEuler')
    return False

#=================  MAIN FUNCTION TO SOLVE THE GROWTH AND REMODELING PROBLEM IN ARTERIES===========
def homogenized_CMT():
    """A hyperelastic implicit static example using ArterialMixture model
        of a cylinder under pression with hexahedral elements
    """
    ProblemPath = PWD(__file__)
    mesh_file = ProblemPath + '/Quarter_Cylinder.msh'

#===============  MESH PROCESING  ==========================
    # Build mesh with Florence tools from GMSH mesh
    mesh = Mesh()
    mesh.Read(filename=mesh_file, reader_type="gmsh", element_type="hex",read_surface_info=True)
    ndim = mesh.InferSpatialDimension()
    #Boolean arrays for boundary condition in Dirichlet
    BottomSurface = np.zeros(mesh.nnode,dtype=bool)
    TopSurface = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_Z = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_X = np.zeros(mesh.nnode,dtype=bool)
    #Boolean array for boundary condition in Neumann
    InnerSurface = np.zeros(mesh.faces.shape[0],dtype=bool)
    InnerFaces = []
    for idface in range(mesh.faces.shape[0]):
        if mesh.face_to_surface[idface] == 13:
            for i in range(mesh.faces.shape[1]):
                BottomSurface[mesh.faces[idface,i]] = True
        elif mesh.face_to_surface[idface] == 21:
            for i in range(mesh.faces.shape[1]):
                TopSurface[mesh.faces[idface,i]] = True
        elif mesh.face_to_surface[idface] == 4:
            for i in range(mesh.faces.shape[1]):
                Symmetry_Z[mesh.faces[idface,i]] = True
        elif mesh.face_to_surface[idface] == 26:
            for i in range(mesh.faces.shape[1]):
                Symmetry_X[mesh.faces[idface,i]] = True
        elif mesh.face_to_surface[idface] == 25:
            InnerFaces.append([int(i) for i in mesh.faces[idface,:]])
            InnerSurface[idface] = True
    
    DirichletBoundary = {}
    DirichletBoundary['Bottom'] = BottomSurface
    DirichletBoundary['Top'] = TopSurface
    DirichletBoundary['SymmetryZ'] = Symmetry_Z
    DirichletBoundary['SymmetryX'] = Symmetry_X

    InnerFaces = np.array(InnerFaces,copy=True)
    NeumannBoundary = {}
    NeumannBoundary['InnerLogic'] = InnerSurface
    NeumannBoundary['InnerFaces'] = InnerFaces

#==================  FORMULATION  =========================
    formulation = DisplacementFormulation(mesh)
    NoGauss = formulation.function_spaces[0].AllGauss.shape[0]
    Bases = formulation.function_spaces[0].Bases

    LagrangeGaussCoords = np.zeros((mesh.nelem,NoGauss,3),dtype=np.float64)
    for elem in range(mesh.nelem):
        LagrangeCoords = mesh.points[mesh.elements[elem,:],:]
        ElemLagrangeGaussCoords = np.einsum('ij,ik->jk',Bases,LagrangeCoords)
        for i in range(3):
            LagrangeGaussCoords[elem,:,i] = ElemLagrangeGaussCoords[:,i]

    R_gauss = np.square(LagrangeGaussCoords[:,:,0]) + \
            np.square(LagrangeGaussCoords[:,:,2])
    R_gauss = np.sqrt(R_gauss)
    Y = LagrangeGaussCoords[:,:,1]

#===============  MATERIAL DEFINITION  ====================
    # Kinematic Definiton
    TotalDeformation = {}
    TotalDeformation['F'] = np.zeros((mesh.nelem,NoGauss,3,3),dtype=np.float64)
    TotalDeformation['F'][:,:,0,0] = 1.0
    TotalDeformation['F'][:,:,1,1] = 1.0
    TotalDeformation['F'][:,:,2,2] = 1.0
    TotalDeformation['J'] = np.ones((mesh.nelem,NoGauss),dtype=np.float64)
    # Dictionary for fibres orientations ['thick','smc','co1','co2','co3','co4']
    fibre_direction = np.zeros((6,mesh.nelem,NoGauss,ndim),dtype=np.float64)

    total_density = 1050.0
    # Array for Density mixture and Remodeling
    GrowthRemodeling = np.zeros((11,mesh.nelem,NoGauss),dtype=np.float64)
    GrowthRemodeling0 = np.zeros((11,mesh.nelem,NoGauss),dtype=np.float64)
    # Densities per constituent ['ela','smc','co1','co2','co3','co4']
    GrowthRemodeling[0][:,:] = 0.23*total_density
    GrowthRemodeling[1][:,:] = 0.15*total_density
    GrowthRemodeling[2][:,:] = 0.062*total_density
    GrowthRemodeling[3][:,:] = 0.248*total_density
    GrowthRemodeling[4][:,:] = 0.248*total_density
    GrowthRemodeling[5][:,:] = 0.062*total_density
    # Remodeling directional magnitude ['smc','co1','co2','co3','co4']
    GrowthRemodeling[6:11][:,:] = 1.0
    for i in range(11):
        GrowthRemodeling0[i][:,:] = GrowthRemodeling[i][:,:]
    # Growth magnitude
    Growth = np.ones((mesh.nelem,NoGauss),dtype=np.float64)
    
    # Deposition Stretch
    DepositionStretch = {}
    DepositionStretch['ela'] = np.zeros((mesh.nelem,NoGauss,3,3),dtype=np.float64)
    DepositionStretch['ela'][:,:,2,2] = 1.25
    DepositionStretch['ela'][:,:,1,1] = 1.34
    DepositionStretch['ela'][:,:,0,0] = 1.0/(1.25*1.34)
    DepositionStretch['smc'] = 1.1
    DepositionStretch['col'] = 1.062

    # Directions for fibres and element orientations
    fibre_direction = Directions(mesh,LagrangeGaussCoords)
    
    # Define hyperelastic material for mesh
    material = ArterialWallMixture_im(ndim,
                 mu3D=72.0,
                 c1m=15.2,
                 c2m=11.4,
                 c1c=1136.0,
                 c2c=11.2,
                 kappa=100.0e3,
                 anisotropic_orientations=fibre_direction,
                 deposition_stretch=DepositionStretch,
                 GrowthRemodeling=GrowthRemodeling,
                 GrowthRemodeling0=GrowthRemodeling0,
                 Growth=Growth,
                 TotalDeformation=TotalDeformation
                 )
#    for i in range(11):
#        material.GrowthRemodeling0[i][:,:]=GrowthRemodeling[i][:,:]
#===============  BOUNDARY CONDITIONS  ====================
    # Dirichlet Boundary Conditions
    def Dirichlet_Function(mesh, DirichletBoundary):
        #num_of_points=mesh.point.shape[0], ndim=3 ==> boundary_data(npoin,ndim)
        boundary_data = np.zeros((mesh.nnode, 3))+np.NAN
        # boundary conditions base on BoundarySurface boolean array
        boundary_data[DirichletBoundary['Bottom'],1] = 0.
        boundary_data[DirichletBoundary['Top'],1] = 0.
        boundary_data[DirichletBoundary['SymmetryZ'],2] = 0.
        boundary_data[DirichletBoundary['SymmetryX'],0] = 0.

        return boundary_data

    # Neumann Boundary Conditions
    def Neumann_Function(mesh, NeumannBoundary):
        boundary_flags = np.zeros(mesh.faces.shape[0],dtype=np.uint8)
        boundary_data = np.zeros((mesh.faces.shape[0], 3))
        # Force magnitud
        mag = 13.332e3

        face = 0
        for idf in range(mesh.faces.shape[0]):
            if NeumannBoundary['InnerLogic'][idf]:
                vertex = np.zeros(3)
                distances = np.zeros(mesh.faces.shape[1])
                coord = mesh.points[NeumannBoundary['InnerFaces'][face,:],:]
                # center point in the surface
                midpoint = np.sum(coord,axis=0)/mesh.faces.shape[1]
                midpoint_mag = np.sqrt(np.square(midpoint[0]) + np.square(midpoint[2]))
                # normal to surface
                normal = midpoint/midpoint_mag
                # area of the face
                for i in range(NeumannBoundary['InnerFaces'].shape[1]):
                    vertex = coord[i,:] - midpoint
                    distances[i] = la.norm(vertex)
            
                semi_diagonal = max(distances)
                area = 2.*semi_diagonal**2
                mag0 = normal[0]*mag*area
                mag2 = normal[2]*mag*area
                boundary_data[idf,0] = mag0
                boundary_data[idf,2] = mag2
                face += 1

        boundary_flags[NeumannBoundary['InnerLogic']] = True

        return boundary_flags, boundary_data

#===============  SOLVER DEFINITION  ======================
    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(Dirichlet_Function, mesh, DirichletBoundary)
    boundary_condition.SetNeumannCriteria(Neumann_Function, mesh, NeumannBoundary)

    fem_solver = FEMSolver( number_of_load_increments=1,
                            analysis_nature="nonlinear",
                            analysis_type="static",
                            break_at_stagnation=False,
                            maximum_iteration_for_newton_raphson=50,
                            optimise=False,
                            print_incremental_log=True,
                            memory_store_frequency=1,
                            newton_raphson_tolerance=1.0e-5)

#=================  HOMEOSTATIC SOLUTION  =======================
    # Homeostatic step solution
    print('=====================================')
    print('==  COMPUTE HOMEOSTATIC STATE  ==')
    print('=====================================')
    # Compute the Elastin deposition stretch in radial and circumferential direction
    ElastinDepositionStretch0(material=material,R_gauss=R_gauss,nelem=mesh.nelem,NoGauss=NoGauss)
    deposition_radial = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
    deposition_circum = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
    deposition_longit = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
    for elem in range(mesh.nelem):
        for gcounter in range(NoGauss):
            deposition_radial[elem][gcounter] = material.deposition_stretch['ela'][elem][gcounter][0,0]
            deposition_circum[elem][gcounter] = material.deposition_stretch['ela'][elem][gcounter][1,1]
            deposition_longit[elem][gcounter] = material.deposition_stretch['ela'][elem][gcounter][2,2]

    print(deposition_radial)
    print(deposition_circum)
    print(deposition_longit)
    # Call the FEM solver to compute the homeostatic state of the artery
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh,
        material=material, boundary_condition=boundary_condition)
    # Check the error displacement
    dmesh = Mesh()
    dmesh.points = solution.sol[:,:,0]
    dmesh_bounds = dmesh.Bounds
    error = 100.0*np.sqrt(dmesh_bounds[0,0]**2+dmesh_bounds[0,2]**2)/0.010
    print(error)
    # Write Homeostatic state to paraview
    solution.WriteVTK('GandR_0',quantity=0)
    print('... Homeostatic step finished')
    # Array for Directional Homeostatic Stresses
    Stress_H = HomeostaticStress_(material=material, nelem=mesh.nelem , NoGauss=NoGauss)

#=================  REMODELING SOLUTION  =======================
    print('=====================================')
    print('==  COMPUTE GROWTH AND REMODELING  ==')
    print('=====================================')
    # Growth and Remodeling steps [30 days]
    total_time = 10
    time = 0.0
    # choose a Dt
    Delta_t = 10.0
    step = 0
    fem_solver.newton_raphson_tolerance=1.0e-5
    while time<total_time:
        # Update Euler coordinates
        TotalDisplacements = solution.sol[:,:,0]
        euler_x = mesh.points + TotalDisplacements
        # prepare for next step
        time += Delta_t
        step += 1
        print('==== STEP: '+str(step)+' |8===D| TIME: '+str(time)+' days ====')
        print('*** Compute Solution')
    #**** Iterative loop to find the TotalDeformationGradient
        for Iter in range(1):
            print('Iteration='+str(Iter))
        #**** compute G&R parameters at t_n+1 **** F_{gr}(t_{n+1})
            Euler_failed = GR_BackwardEuler(Iter=Iter,time=time,Delta_t=Delta_t,material=material,
                HomeostaticStress=Stress_H,Y=Y,nelem=mesh.nelem,NoGauss=NoGauss)
            if Euler_failed:
                print('ERROR: Backward Euler method failed')
                break
            print(material.GrowthRemodeling[0][:,:])
            print(material.GrowthRemodeling[1][:,:])
            print(material.GrowthRemodeling[2][:,:])
            print(material.GrowthRemodeling[3][:,:])
            print(material.GrowthRemodeling[4][:,:])
            print(material.GrowthRemodeling[5][:,:])
        #**** compute mechanical equilibrium at t_{n+1} **** F(t_{n+1})
            solution = fem_solver.Solve(formulation=formulation,mesh=mesh,material=material, 
                boundary_condition=boundary_condition,Eulerx=euler_x)
            # Check Residual
            if np.isnan(fem_solver.norm_residual) or fem_solver.norm_residual>1e06:
                fem_solver.newton_raphson_failed_to_converge = True
                print('MODEL DID NOT CONVERGE!')
                break
    #**** OUTPUT ****
        if Euler_failed or fem_solver.newton_raphson_failed_to_converge:
            print('STOP: the method is not working well!')
            break
        solution.WriteVTK('GandR_'+str(step),quantity=0)
        # last step Growth and remodeling parameters
        for i in range(11):
            material.GrowthRemodeling0[i][:,:] = material.GrowthRemodeling[i][:,:]


if __name__ == "__main__":
    homogenized_CMT()

