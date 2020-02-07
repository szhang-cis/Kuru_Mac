# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
from numpy import einsum
# Build a path for python to Kuru
#sys.path.append(os.path.join(os.path.expanduser("~/kuru")))
sys.path.append(os.path.expanduser("~/kuru"))
#import Florence
from Kuru import *

#============================================================
#============= ANISOTROPIC FIBRE DIRECTIONS  ================
#============================================================
def Directions(mesh):
    """
        Routine dedicated to compute the fibre direction of components in integration point for 
        the Material in Florence and for the auxiliar routines in this script. First three directions 
        are taken into the code for Rotation matrix, so always it should be present in this order,
        Normal, Tangential, Axial.
    """
    ndim = mesh.InferSpatialDimension()
    nfibre = 6
    direction = np.zeros((mesh.nelem,nfibre,ndim),dtype=np.float64)
    # Geometric definitions per element
    divider = mesh.elements.shape[1]
    directrix = [0.,1.,0.]
    # Loop throught the element in the mesh
    for elem in range(mesh.nelem):
        # Geometric definitions per element
        center = np.sum(mesh.points[mesh.elements[elem,:],:],axis=0)/divider
        tangential = np.cross(directrix,center)
        tangential = tangential/np.linalg.norm(tangential)
        normal = np.cross(tangential,directrix)
        direction[elem][0][:] = normal/np.linalg.norm(normal)
        # Define the anisotropic orientations
        direction[elem][1][:]=tangential
        direction[elem][2][:]=directrix
        direction[elem][3][:]=np.multiply(directrix,np.cos(np.pi/4.)) + np.multiply(tangential,np.sin(np.pi/4.))
        direction[elem][4][:]=np.multiply(directrix,np.cos(np.pi/4.)) - np.multiply(tangential,np.sin(np.pi/4.))
        direction[elem][5][:]=tangential

    return direction

def GetRadialStretch(mesh,material):
    print("Computing the Radial deposition stretch")
    pressure = 13.3322e3
    R = 0.010
    H = 0.00141
    for node in range(mesh.nnode):
        kappa = material.kappa*material.field_variables[node,11]
        mu3D = material.mu*material.field_variables[node,11]
        stretch_r = material.field_variables[node,0]
        stretch_t = material.field_variables[node,4]
        stretch_z = material.field_variables[node,8]
        radius = np.sqrt(mesh.points[node,0]**2+mesh.points[node,2]**2)
        JMAX=21
        for j in range(JMAX):
            if j==20: exit('Maximum number of iteration reach at Newton method at node: '+str(node))
            func = pressure*(1.-(radius-R)/H)+ \
                    mu3D*(stretch_r*stretch_t*stretch_z)**(-2./3.)*(stretch_r**2-\
                    1./3.*(stretch_r**2+stretch_t**2+stretch_z**2)) +\
                    kappa*stretch_r*stretch_t*stretch_z*(stretch_r*stretch_t*stretch_z-1.)
            df = 4./3.*mu3D*(stretch_r*stretch_t*stretch_z)**(-2./3)*stretch_r-\
                    2./3.*mu3D*stretch_r**(-5./3.)*(stretch_t*stretch_z)**(-2./3.)*(stretch_r**2-\
                    1./3.*(stretch_r**2+stretch_t**2+stretch_z**2)) +\
                    kappa*(2.*stretch_r*stretch_t**2*stretch_z**2-stretch_t*stretch_z)
            dx = func/df
            stretch_r -= dx
            if np.absolute(dx)<1.e-6:
                material.field_variables[node,0] = stretch_r
                break

    print("Radial deposition stretch is ready")

def GetTangentialPenalty(mesh,material):
    print("Computing the Tangential penalty mu2D")
    pressure = 13.3322e3
    R = 0.010
    H = 0.00141
    for node in range(mesh.nnode):
        kappa = material.kappa*material.growth_remodeling[node,0]
        mu3D = material.mu3D*material.growth_remodeling[node,0]
        stretch_r = material.deposition_stretch['Radial'][node]
        stretch_t = material.deposition_stretch['Tangential']
        stretch_z = material.deposition_stretch['Axial']
        k1m = material.k1m*material.growth_remodeling[node,1]
        k2m = material.k2m
        s_act = 54.0e3*growth_remodeling[node,1]
        den0 = 1050.0
        lambda_m = 1.4
        lambda_0 = 0.8
        lambda_a = 1.0
        stretch_m = material.deposition_stretch['Muscle']
        k1c1 = material.k1c*material.growth_remodeling[node,2]
        k1c2 = material.k1c*material.growth_remodeling[node,3]
        k1c3 = material.k1c*material.growth_remodeling[node,4]
        k2c = material.k2c
        stretch_c = material.deposition_stretch['Collagen']
        #radius = np.sqrt(mesh.points[node,0]**2+mesh.points[node,2]**2)
        mu2D = (pressure*R/H - \
                kappa*stretch_r*stretch_t*stretch_z*(stretch_r*stretch_t*stretch_z-1.) - \
                mu3D*(stretch_r*stretch_t*stretch_z)**(-2./3.)*(stretch_t**2 - \
                1./3.*(stretch_r**2+stretch_t**2+stretch_z**2)) - \
                2.*k1m*(stretch_m**2-1.)*np.exp(k2m*(stretch_m**2-1.)**2)*stretch_m**2 - \
                s_act*(1.-((lambda_m-lambda_a)/(lambda_m-lambda_0))**2)/den0 - \
                2.*(k1c1+k1c2/2.+k1c3/2.)*(stretch_c**2-1.)*np.exp(k2c*(stretch_c**2-1.)**2)*stretch_c**2)/\
                (stretch_t**2-1./(stretch_t*stretch_z)**2)
        material.mu2D[node] = mu2D/material.growth_remodeling[node,0]
    print("Tangential penalty mu2D is ready")

#============================================================
#==== GROWTH AND REMODELING RATES AND INTEGRALS  ============
#============================================================
def GrowthRemodelingRates(time,GrowthRemodeling0,mesh,Stress_H,FibreStress,Softness):
    """
        This are the rates of Growth and Remodeling.
    """
    # Elastin degradation [days]
    L_dam = 0.010
    t_dam = 40.0
    D_max = 0.5
    T_ela = 101.0*365.25
    den0_e = 241.5
    # Fribres turnover and remodeling
    turnover = 101.0
    gain = 0.05/turnover
    # Compute GrowthRemodeling
    Rates = np.zeros((mesh.nnode,11),dtype=np.float64)
    for node in range(mesh.nnode):
        AxialCoord = mesh.points[node,1]
        # Update elastin density
        Rates[node,0] = -GrowthRemodeling0[node,0]/T_ela - \
            (D_max/t_dam)*np.exp(-0.5*(AxialCoord/L_dam)**2)*np.exp(-(time)/t_dam)*den0_e
        for key in [1,2,3,4,5]:
            # fibres density rates
            DeltaStress = FibreStress[node,key-1] - Stress_H[node,key-1]
            Rates[node,key] = gain*GrowthRemodeling0[node,key]*DeltaStress/Stress_H[node,key-1]
            # fibres remodeling rates
            Rates[node,key+5] = (gain*DeltaStress/Stress_H[node,key-1] + 1./turnover)*\
                    DeltaStress*Softness[node,key-1]

    return Rates

#============================================================
#======= EXPLICIT FORWARD-EULER METHOD INTEGRAND  ===========
#============================================================
def ForwardEuler_Integrand(time,Delta_t,Stress_H,FibreStress,Softness,
        fem_solver,formulation,mesh,material,boundary_condition):

    den0_tot = 1050.0
    growth_remodeling0 = np.copy(material.field_variables[:,11:23])
    growth_remodeling = np.copy(material.field_variables[:,11:23])
    #**** K1: f(t_n,y_n)
    k1 = GrowthRemodelingRates(time,growth_remodeling0,mesh,Stress_H,FibreStress,Softness)
    growth_remodeling[:,:11] = growth_remodeling0[:,:11] + Delta_t*k1
    den_tot = np.zeros((mesh.nnode),dtype=np.float64)
    for i in range(6):
        den_tot += growth_remodeling[:,i]
    growth_remodeling[:,11] = den_tot/den0_tot
    # solution k1
    material.field_variables[:,11:23] = growth_remodeling
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material,
            boundary_condition=boundary_condition)
    
    # Check Residual
    if np.isnan(fem_solver.norm_residual) or fem_solver.norm_residual>1e06:
        exit('MODEL DID NOT CONVERGE!')

    return solution

#============================================================
#========== EXPLICIT MID-POINT METHOD INTEGRAND  ============
#============================================================
def MidPoint_Integrand(time,Delta_t,Stress_H,FibreStress,Softness,
        fem_solver,formulation,mesh,material,boundary_condition):

    den0_tot = 1050.0
    growth_remodeling0 = np.copy(material.field_variables[:,11:23])
    growth_remodeling = np.copy(material.field_variables[:,11:23])
    #**** K1: f(t_n,y_n)
    k1 = GrowthRemodelingRates(time,growth_remodeling0,mesh,Stress_H,FibreStress,Softness)
    growth_remodeling[:,:11] = growth_remodeling0[:,:11] + Delta_t*k1/2.
    den_tot = np.zeros((mesh.nnode),dtype=np.float64)
    for i in range(6):
        den_tot += growth_remodeling[:,i]
    growth_remodeling[:,11] = den_tot/den0_tot
    # solution k1
    #material.field_variables[:,11:23] = growth_remodeling
    #solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material,
    #        boundary_condition=boundary_condition)
    #solution.StressRecovery()
    #FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
    #Softness = solution.recovered_fields['Softness'][-1,:,:]
    
    #**** K2: f(t_n+h/2,y_n+k1/2)
    k2 = GrowthRemodelingRates(time+Delta_t/2.,growth_remodeling,mesh,Stress_H,FibreStress,Softness)
    growth_remodeling[:,:11] = growth_remodeling0[:,:11] + Delta_t*k2
    den_tot = np.zeros((mesh.nnode),dtype=np.float64)
    for i in range(6):
        den_tot += growth_remodeling[:,i]
    growth_remodeling[:,11] = den_tot/den0_tot
    # solution k2
    material.field_variables[:,11:23] = np.copy(growth_remodeling)
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material,
            boundary_condition=boundary_condition)
    # Check Residual
    #if np.isnan(fem_solver.norm_residual) or fem_solver.norm_residual>1e06:
    #    exit('MODEL DID NOT CONVERGE!')

    return solution

#============================================================
#========== EXPLICIT RUNGE-KUTTA METHOD INTEGRAND  ==========
#============================================================
def RungeKutta_Integrand(time,Delta_t,Stress_H,FibreStress,Softness,
        fem_solver,formulation,mesh,material,boundary_condition):

    den0_tot = 1050.0
    growth_remodeling0 = np.copy(material.field_variables[:,11:23])
    growth_remodeling = np.copy(material.field_variables[:,11:23])
    #**** K1: f(t_n,y_n)
    k1 = GrowthRemodelingRates(time,growth_remodeling0,mesh,Stress_H,FibreStress,Softness)
    growth_remodeling[:,:11] = growth_remodeling0[:,:11] + Delta_t*k1/2.
    den_tot = np.zeros((mesh.nnode),dtype=np.float64)
    for i in range(6):
        den_tot += growth_remodeling[:,i]
    growth_remodeling[:,11] = den_tot/den0_tot
    # solution k1
    """
    material.field_variables[:,11:23] = growth_remodeling
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material,
            boundary_condition=boundary_condition)
    solution.StressRecovery()
    FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
    Softness = solution.recovered_fields['Softness'][-1,:,:]
    """
    #**** K2: f(t_n+h/2,y_n+k1/2)
    k2 = GrowthRemodelingRates(time+Delta_t/2.,growth_remodeling,mesh,Stress_H,FibreStress,Softness)
    growth_remodeling[:,:11] = growth_remodeling0[:,:11] + Delta_t*k2/2.
    den_tot = np.zeros((mesh.nnode),dtype=np.float64)
    for i in range(6):
        den_tot += growth_remodeling[:,i]
    growth_remodeling[:,11] = den_tot/den0_tot
    # solution k2
    """
    material.field_variables[:,11:23] = growth_remodeling
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material,
            boundary_condition=boundary_condition)
    solution.StressRecovery()
    FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
    Softness = solution.recovered_fields['Softness'][-1,:,:]
    """
    #**** K3: f(t_n+h/2,y_n+k2/2)
    k3 = GrowthRemodelingRates(time+Delta_t/2.,growth_remodeling,mesh,Stress_H,FibreStress,Softness)
    growth_remodeling[:,:11] = growth_remodeling0[:,:11] + Delta_t*k3
    den_tot = np.zeros((mesh.nnode),dtype=np.float64)
    for i in range(6):
        den_tot += growth_remodeling[:,i]
    growth_remodeling[:,11] = den_tot/den0_tot
    # solution k3
    """
    material.field_variables[:,11:23] = growth_remodeling
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material,
            boundary_condition=boundary_condition)
    solution.StressRecovery()
    FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
    Softness = solution.recovered_fields['Softness'][-1,:,:]
    """
    #**** K4: f(t_n+h/2,y_n+k2/2)
    k4 = GrowthRemodelingRates(time+Delta_t,growth_remodeling,mesh,Stress_H,FibreStress,Softness)
    growth_remodeling[:,:11] = growth_remodeling0[:,:11] + Delta_t*(k1/6.+k2/3.+k3/3.+k4/6.)
    den_tot = np.zeros((mesh.nnode),dtype=np.float64)
    for i in range(6):
        den_tot += growth_remodeling[:,i]
    growth_remodeling[:,11] = den_tot/den0_tot
    # solution k4
    material.field_variables[:,11:23] = growth_remodeling
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material,
            boundary_condition=boundary_condition)

    # Check Residual
    #if np.isnan(fem_solver.norm_residual) or fem_solver.norm_residual>1e06:
    #    exit('MODEL DID NOT CONVERGE!')

    return solution

#============================================================
#======= IMPLICIT BACKWARD-EULER METHOD INTEGRAND  ==========
#============================================================
def Residual_GR(time,Delta_t,AxialCoord,GrowthRemodeling0,GrowthRemodeling,Stress_H,FibreStress,Softness,ires):
    """
        This are the rates of Growth and Remodeling.
    """
    # Elastin degradation [days]
    L_dam = 0.010
    t_dam = 40.0
    D_max = 0.5
    T_ela = 101.0*365.25
    den0_e = 241.5
    # Fribres turnover and remodeling
    turnover = 101.0
    gain = 0.05/turnover
    # Compute GrowthRemodeling
    Residual = np.zeros((11),dtype=np.float64)
    # Update elastin density
    if ires == 0:
        Residual[0] = -GrowthRemodeling[0] + GrowthRemodeling0[0] - Delta_t*(GrowthRemodeling[0]/T_ela + \
            (D_max/t_dam)*np.exp(-0.5*(AxialCoord/L_dam)**2 - time/t_dam)*den0_e)
    elif ires==1 or ires==2 or ires==3 or ires==4 or ires==5:
        # fibres density rates
        DeltaStress = FibreStress[ires-1] - Stress_H[ires-1]
        Residual[ires] = -GrowthRemodeling[ires] + GrowthRemodeling0[ires] + \
                Delta_t*gain*GrowthRemodeling[ires]*DeltaStress/Stress_H[ires-1]
    elif ires==6 or ires==7 or ires==8 or ires==9 or ires==10:
        # fibres remodeling rates
        DeltaStress = FibreStress[ires-6] - Stress_H[ires-6]
        lambda_r_dot = gain*DeltaStress/Stress_H[ires-6] + 1./turnover
        Residual[ires] = -GrowthRemodeling[ires] + GrowthRemodeling0[ires] + \
                Delta_t*lambda_r_dot*DeltaStress*Softness[ires-6]

    return Residual[ires]

def D_ridders(time,Delta_t,AxialCoord,growth_remodeling0,growth_remodeling,Stress_H,FibreStress,Softness,
        ires,ipar):
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
    x_pos = np.copy(growth_remodeling)
    x_neg = np.copy(growth_remodeling)
    x_pos[ipar] = growth_remodeling[ipar] + hh
    x_neg[ipar] = growth_remodeling[ipar] - hh
    func_pos = Residual_GR(time,Delta_t,AxialCoord,growth_remodeling0,x_pos,Stress_H,FibreStress,Softness,ires)
    func_neg = Residual_GR(time,Delta_t,AxialCoord,growth_remodeling0,x_neg,Stress_H,FibreStress,Softness,ires)
    a[0,0] = (func_pos-func_neg)/(2.0*hh)
    err=BIG
    # Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation.
    for i in range(1,NTAB):
        hh=hh/CON
        # Try new, smaller stepsize.
        x_pos[ipar] = growth_remodeling[ipar] + hh
        x_neg[ipar] = growth_remodeling[ipar] - hh
        func_pos = Residual_GR(time,Delta_t,AxialCoord,growth_remodeling0,x_pos,Stress_H,FibreStress,Softness,ires)
        func_neg = Residual_GR(time,Delta_t,AxialCoord,growth_remodeling0,x_neg,Stress_H,FibreStress,Softness,ires)
        a[0,i] = (func_pos-func_neg)/(2.0*hh)
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
            #print('Early quit in df_ridders function')
            return dfridr
    return dfridr

def ResidualEval_GR(time,Delta_t,AxialCoord,growth_remodeling0,growth_remodeling,Stress_H,FibreStress,Softness):

    nvar = growth_remodeling0.shape[0] - 1
    Residual = np.zeros((nvar),dtype=np.float64)
    DResidual = np.zeros((nvar,nvar),dtype=np.float64)
    for ires in range(nvar):
        Residual[ires] = Residual_GR(time,Delta_t,AxialCoord,growth_remodeling0,growth_remodeling,
                Stress_H,FibreStress,Softness,ires)
        for ipar in range(nvar):
            DResidual[ires,ipar] = D_ridders(time,Delta_t,AxialCoord,growth_remodeling0,growth_remodeling,
                    Stress_H,FibreStress,Softness,ires,ipar)

    return Residual,DResidual

def BackwardEuler_Integrand(time,Delta_t,Stress_H,FibreStress,Softness,
        fem_solver,formulation,mesh,material,boundary_condition):

    JMAX=50
    den0_tot = 1050.0
    growth_remodeling0 = np.copy(material.field_variables[:,11:23])
    growth_remodeling = np.copy(material.field_variables[:,11:23])
    BIGResidual = np.zeros((mesh.nnode,11),dtype=np.float64)
    # Set a first possible solution
    #for node in range(mesh.nnode):
    #    Residual = Residual_GR(time,node,growth_remodeling0[node,:],growth_remodeling[node,:],mesh,
    #            Stress_H[node,:],FibreStress[node,:],Softness[node,:])
    # LOOP ON ITERATIONS
    for Iter in range(JMAX):
        # LOOP ON NODES
        for node in range(mesh.nnode):
            AxialCoord = mesh.points[node,1]
            Residual, DResidual = ResidualEval_GR(time,Delta_t,AxialCoord,growth_remodeling0[node,:],
                    growth_remodeling[node,:],Stress_H[node,:],FibreStress[node,:],Softness[node,:])
            growth_remodeling[node,:11] -= np.dot(np.linalg.inv(DResidual),Residual)
            BIGResidual[node,:] = Residual

        den_tot = np.zeros((mesh.nnode),dtype=np.float64)
        for i in range(6):
            den_tot += growth_remodeling[:,i]
        growth_remodeling[:,11] = den_tot/den0_tot
        # update equilibrium
        material.field_variables[:,11:23] = growth_remodeling
        solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material,
            boundary_condition=boundary_condition)
        solution.StressRecovery()
        FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
        Softness = solution.recovered_fields['Softness'][-1,:,:]
        TotalDisplacements = solution.sol[:,:,-1]
        if np.abs(np.linalg.norm(TotalDisplacements)) < 1.0e-6:
            return solution
        if np.abs(np.linalg.norm(BIGResidual)) < 1.0e-6:
            return solution
        # Check Residual
        if np.isnan(fem_solver.norm_residual) or fem_solver.norm_residual>1e06:
            exit('MODEL DID NOT CONVERGE!')

    exit('Backward Euler exceeded maximum itarations')

#============================================================
#===============  HOMOGENIZED CMT  ==========================
#============================================================
def HomogenizedConstrainedMixtureTheory(p=1,TimeIntegrand="forward_euler"):

    ProblemPath = os.path.dirname(os.getcwd())
    mesh_file = ProblemPath + '/Quarter_Ring.msh'

    #===============  MESH PROCESING  ==========================
    # Build mesh with Florence tools from GMSH mesh
    mesh = Mesh()
    mesh.Read(filename=mesh_file, reader_type="gmsh", element_type="hex",read_surface_info=True)
    ndim = mesh.InferSpatialDimension()
    mesh.GetHighOrderMesh(p=p)
    #Boolean arrays for boundary condition in Dirichlet
    BottomSurface = np.zeros(mesh.nnode,dtype=bool)
    TopSurface = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_Z = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_X = np.zeros(mesh.nnode,dtype=bool)
    #Boolean array for boundary condition in Neumann mesh.faces[id,node]
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
    PressureBoundary = {}
    PressureBoundary['InnerLogic'] = InnerSurface
    PressureBoundary['InnerFaces'] = InnerFaces

    #===============  MATERIAL DEFINITION  ====================
    field_variables = np.zeros((mesh.nnode,23),dtype=np.float64)
    # Deposition Stretches
    field_variables[:,0] = 1.0/(1.34*1.25)
    field_variables[:,4] = 1.34
    field_variables[:,8] = 1.25
    field_variables[:,9] = 1.1
    field_variables[:,10] = 1.062
    # Total initial density
    field_variables[:,11] = 241.5
    field_variables[:,12] = 157.5
    field_variables[:,13] = 65.1
    field_variables[:,14] = 260.4
    field_variables[:,15] = 260.4
    field_variables[:,16] = 65.1
    field_variables[:,17:23] = 1.0

    # fibre directions [thick,smc,co1,co2,co3,co4]
    fibre_direction = Directions(mesh)

    # Total initial density
    #mu2D = np.zeros((mesh.nnode),dtype=np.float64)
    #mu2D = 0.0

    # Define hyperelastic material for mesh
    material = ArterialWallMixture(ndim,
            mu=72.0,
            kappa=72.0*33.0,
            k1m=7.6,
            k2m=11.4,
            k1c=568.0,
            k2c=11.2,
            anisotropic_orientations=fibre_direction,
            field_variables=field_variables)

    # kappa/mu=20  => nu=0.475 (Poisson's ratio)
    # kappa/mu=33  => nu=0.485 (Poisson's ratio)
    # kappa/mu=100 => nu=0.495 (Poisson's ratio)

    # Homogenization of stresses to a thin-walled case
    #GetRadialStretch(mesh,material)
    #print(material.field_variables[:,0])
    #GetTangentialPenalty(mesh,material)
    #print(material.field_variables[:])

    #==================  FORMULATION  =========================
    formulation = DisplacementFormulation(mesh)
    #formulation = DisplacementMixedFormulation(mesh)

    #===============  BOUNDARY CONDITIONS  ====================
    # Dirichlet Boundary Conditions
    def Dirichlet_Function(mesh, DirichletBoundary):
        boundary_data = np.zeros((mesh.nnode, 3))+np.NAN
        # boundary conditions base on BoundarySurface boolean array
        boundary_data[DirichletBoundary['Bottom'],1] = 0.
        boundary_data[DirichletBoundary['Top'],1] = 0.
        boundary_data[DirichletBoundary['SymmetryZ'],2] = 0.
        boundary_data[DirichletBoundary['SymmetryX'],0] = 0.

        return boundary_data

    # Pressure Boundary Conditions
    def Pressure_Function(mesh, PressureBoundary):
        boundary_flags = np.zeros(mesh.faces.shape[0],dtype=np.uint8)
        boundary_data = np.zeros((mesh.faces.shape[0]))
        # Force magnitud
        mag = 13.3322e3

        for idf in range(mesh.faces.shape[0]):
            if PressureBoundary['InnerLogic'][idf]:
                boundary_data[idf] = mag

        boundary_flags[PressureBoundary['InnerLogic']] = True

        return boundary_flags, boundary_data

    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(Dirichlet_Function, mesh, DirichletBoundary)
    boundary_condition.SetPressureCriteria(Pressure_Function, mesh, PressureBoundary)

    #===============  SOLVER DEFINITION  ======================
    fem_solver_h = FEMSolver(analysis_nature="nonlinear",
                       analysis_type="static",
                       break_at_stagnation=False,
                       maximum_iteration_for_newton_raphson=50,
                       optimise=True,
                       print_incremental_log=True,
                       has_moving_boundary=True,
                       number_of_load_increments=1)

    fem_solver_gr = FEMSolver(analysis_nature="nonlinear",
                       analysis_type="static",
                       break_at_stagnation=False,
                       maximum_iteration_for_newton_raphson=50,
                       optimise=True,
                       print_incremental_log=True,
                       has_moving_boundary=True,
                       number_of_load_increments=1)

    #=================  HOMEOSTATIC SOLUTION  =======================
    print('=====================================')
    print('===    COMPUTE HOMEOSTATIC STEP   ===')
    print('=====================================')
    # Call the solver for Homeostatic computation
    solution = fem_solver_h.Solve(formulation=formulation, mesh=mesh,
        material=material, boundary_condition=boundary_condition)
    # Check the error displacement
    dmesh = Mesh()
    dmesh.points = solution.sol[:,:,-1]
    dmesh_bounds = dmesh.Bounds
    distortion = np.sqrt(dmesh_bounds[0,0]**2+dmesh_bounds[0,1]**2+dmesh_bounds[0,2]**2)/0.010
    print('Distortion: '+str(distortion))
    # Write Homeostatic state to paraview
    #solution.WriteVTK('HCMTRing_0',quantity=0)
    print('... Homeostatic step finished')

    # HOMEOSTATIC POSTPROCESS
    solution.StressRecovery()
    DeformationGradient = solution.recovered_fields['F'][-1,:,:,:]
    CauchyStress = solution.recovered_fields['CauchyStress'][-1,:,:]
    FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
    Softness = solution.recovered_fields['Softness'][-1,:,:]
    # Update mesh coordinates
    TotalDisplacements = solution.sol[:,:,-1]
    euler_x = mesh.points + TotalDisplacements
    # Write parameters to files
    file1 = open("growth_remodeling.txt","w+")
    file1.write('%3d %6.3f '%(0,1000.*np.sqrt(euler_x[0,0]**2+euler_x[0,2]**2)))
    for i in range(11,23):
        file1.write('%7.3f '%(material.field_variables[0,i]))
    file1.write('\n')
    file1.close()
    file2 = open("kinematics.txt","w+")
    for i in range(3):
        for j in range(3):
            file2.write('%5.3f '%(DeformationGradient[0,i,j]))
    for i in range(3):
        for j in range(3):
            file2.write('%7.3f '%(CauchyStress[0,i,j]/1000.))
    file2.write('\n')
    file2.close()
    file3 = open("fibre_stress.txt","w+")
    for i in range(5):
        file3.write('%7.3f '%(FibreStress[0,i]/1000.))
    file3.write('\n')
    file3.close()
    # Stock Homeostatic constituent stress
    Stress_H = FibreStress
    
    #=================  REMODELING SOLUTION  =======================
    print('=====================================')
    print('==  COMPUTE GROWTH AND REMODELING  ==')
    print('=====================================')
    # Growth and Remodeling steps [10 days]
    total_time = 5500
    time = 0.0
    Delta_t = 10.0
    step = 0
    while time<total_time:
        # prepare for next step
        time += Delta_t
        step += 1
        print('==== STEP: '+str(step)+' |8===D| TIME: '+str(time)+' days ====')
        print('*** Compute Solution')
        if TimeIntegrand == "mid_point":
            solution = MidPoint_Integrand(time-Delta_t,Delta_t,Stress_H,FibreStress,Softness,
                fem_solver_gr,formulation,mesh,material,boundary_condition)
        elif TimeIntegrand == "runge_kutta":
            solution = RungeKutta_Integrand(time-Delta_t,Delta_t,Stress_H,FibreStress,Softness,
                fem_solver_gr,formulation,mesh,material,boundary_condition)
        elif TimeIntegrand == "backward_euler":
            solution = BackwardEuler_Integrand(time,Delta_t,Stress_H,FibreStress,Softness,
                fem_solver_gr,formulation,mesh,material,boundary_condition)
        elif TimeIntegrand == "forward_euler":
            solution = ForwardEuler_Integrand(time-Delta_t,Delta_t,Stress_H,FibreStress,Softness,
                fem_solver_gr,formulation,mesh,material,boundary_condition)
        else:
            raise ValueError('Integrator not recognized')
        #**** STEPS POSTPROCESS ****
        #solution.WriteVTK('HCMTRing_'+str(step),quantity=0)
        solution.StressRecovery()
        FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
        Softness = solution.recovered_fields['Softness'][-1,:,:]
        DeformationGradient = solution.recovered_fields['F'][-1,:,:,:]
        CauchyStress = solution.recovered_fields['CauchyStress'][-1,:,:]
        # Update mesh coordinates
        TotalDisplacements = solution.sol[:,:,-1]
        euler_x = mesh.points + TotalDisplacements
        # Write parameters to files
        file1 = open("growth_remodeling.txt","a")
        file1.write('%3d %6.3f '%(0,1000.*np.sqrt(euler_x[0,0]**2+euler_x[0,2]**2)))
        for i in range(11,23):
            file1.write('%7.3f '%(material.field_variables[0,i]))
        file1.write('\n')
        file1.close()
        file2 = open("kinematics.txt","a")
        for i in range(3):
            for j in range(3):
                file2.write('%5.3f '%(DeformationGradient[0,i,j]))
        for i in range(3):
            for j in range(3):
                file2.write('%7.3f '%(CauchyStress[0,i,j]/1000.))
        file2.write('\n')
        file2.close()
        file3 = open("fibre_stress.txt","a")
        for i in range(5):
            file3.write('%7.3f '%(FibreStress[0,i]/1000.))
        file3.write('\n')
        file3.close()
    

if __name__ == "__main__":
    HomogenizedConstrainedMixtureTheory(TimeIntegrand="mid_point")
