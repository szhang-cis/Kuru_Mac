# kuru_mac (upadated at 04 Feburary 2021)

This is a Finite-Element code developped by Joan D. Laubrie Soto (joan.laubrie@emse.fr).The code is mainly based on Florence (https://github.com/romeric/florence.git). It can be seen as a reduced version of Florence for the modeling of aortic artery aneurysme's G&R, but with several new features: New follower load (pressure) in Assembly module, New ArterialWallMixture in Material Library and New Explicit Growth and Remodeling time integrator.

First step, download the following codes in your PC: Fastor(https://github.com/romeric/Fastor.git)(master version)(!!!!important: note the path where Fastor is installed), OpenBLAS(https://github.com/xianyi/OpenBLAS.git)(dev version) and Python (https://www.python.org/downloads/)(python 2.7 or later)(!!!!important: note the path where Python is installed).

Second step, download the necessary python modules with pip(tool already installed with Python if its version >2.7.9): Cython      (pip install Cython), NumPy(pip install NumPy)(!!!!important: note the path where it is installed, SciPy(pip install SciPy), pyevtk(pip install pyevtk) and matplotlib(pip install matplotlib).

Third step, compilation of OpenBLAS: open terminal, cd /.../folder/of/your/.../OpenBLAS, tape in the terminal "make" and enter. 
Once it is finished, following the indication of the message and tape in the terminal "make PREFIXE=/.../.../where you want to install OpenBLAS/... intall" (!!!!important: note the path where it is installed (this path may be different to where you have download OpenBLAS)

Fourth step,compilation of Kuru: Go to the each indicated location, tape in the terminal the indicated command and enter ***Important*** Open the makefile in each indicated location where you are going to do the "make" command.Modify the path of Fastor, Numpy and OpenBLAS (if there are in the makefile). These paths should be consistent to where they are installed in your PC.***Important***

	./Tensor --> make
	
	./FunctionSpace/OneDimensional/_OneD --> make
	
	./FiniteElements/Assembly/_Assembly_ --> make cython_assembler_build

	./FiniteElements/LocalAssembly/_KinematicsMeasures_ --> make

	./VariationalPrinciple/_GeometricStiffness_ --> make

	./VariationalPrinciple/_VolumetricStiffness_ --> make

	./VariationalPrinciple/_ConstitutiveStiffness_ --> make

	./MaterialLibrary/LLDispatch --> make MATERIAL=ArteryWallMixture

	./FiniteElements/Assembly/_Assembly_ --> make ASSEMBLY_NAME=RobinForces
	./FiniteElements/Assembly/_Assembly_ --> make ASSEMBLY_NAME=SparseAssemblyNative
	./FiniteElements/Assembly/_Assembly_ --> make ASSEMBLY_NAME=RHSAssemblyNative
	./FiniteElements/Assembly/_Assembly_ --> make ASSEMBLY_NAME=ComputeSparsityPattern

Fifth step: modify the path of Kuru in the python launch scripts and run 

Sixth step: modify the path of Kuru in the python postprocess scripts to convert the generated results to vtk files.


