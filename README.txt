Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D

-- Simulation Package --

The Simulation Package of Metos3D is a software for computation and simulation of steady periodic cycles of passive biogeochemical tracers. It complies to a general programming interface, embedding the biogeochemical modelling into an optimization [Piwonski and Slawig, 2012]. Metos3D is founded on the PETSc library [Balay et al., 2011] as a basis for efficient, parallelized numerical methods. Moreover it's portable and relies only on software, which is freely and public available.
The application of marine transport is based on the Transport Matrix Method, whereas the idea and its use are presented in [Khatiwala et al., 2005], [Khatiwala, 2007] and [Khatiwala, 2008].

-- Quick start --

Assuming PETSc 3.3 is installed, your Metos3D working directory is stored in a variable named $METOS3D, do the following:

1. Prepare the data, i.e. clone the data archive, cd into it and follow the instructions in README.txt:

$>
cd $METOS3D
git clone https://github.com/metos3d/data.git
cd data
cat README.txt

2. Prepare the models, i.e. clone the model archive:

$>
cd $METOS3D
git clone https://github.com/metos3d/model.git

3. Prepare the simulation package, i.e. clone the simpack archive:

$>
cd $METOS3D
git clone https://github.com/metos3d/simpack.git

Assuming now your Metos3D working directory looks like as follows:

$METOS3D/
	data/
	model/
	simpack/

To compile and run the software, do the following:

1. Change into the simulation package directory:

$>
cd simpack

2. Link the data and model directories. For example:

$>
ln -s ../data
ln -s ../model

3. Set the PETSc environment variables. For example:

$>
. petsc/de.uni-kiel.rz.nesh-fe.petsc-3.3-p5.opt.txt

4. Compile the software with a chosen model. For example:

$>
make BGC=model/MITgcm-PO4-DOP/

5. Run executable with chosen options. For example:

$>
./metos3d-simpack-MITgcm-PO4-DOP.exe model/MITgcm-PO4-DOP/option/de.uni-kiel.rz.nesh-fe.option.MITgcm-PO4-DOP.simpack.txt

-- Documentation --

Inspect the 'doc' directory for more detailed information.

-- References --

[Piwonski and Slawig, 2012]
Piwonski, J., Slawig, T., 2012.
Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D -- Simulation Package --.
in preparation.

[Balay et al., 2011]
Balay, S., Buschelman, K., Gropp, W.D., Kaushik, D., Knepley, M.G., McInnes, L.C., Smith, B.F., Zhang, H., 2009.
PETSc web page. http://www.mcs.anl.gov/petsc.

[Khatiwala et al., 2005]
Khatiwala, S., Visbeck, M., Cane, M., 2005.
Accelerated simulation of passive tracers in ocean circulation models.
Ocean Modelling 9, 51– 69.

[Khatiwala, 2007]
Khatiwala, S., 2007.
A computational framework for simulation of biogeochemical tracers in the ocean.
Global Biogeochemical Cycles 21.

[Khatiwala, 2008]
Khatiwala, S., 2008.
Fast spin up of ocean biogeochemical models using matrix-free newton-krylov.
Ocean Modelling 23, 121–129.
