Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D

-- Simulation Package --

The Simulation Package of Metos3D is a software for computation and simulation of steady periodic cycles of passive biogeochemical tracers. It complies to a general programming interface, embedding the biogeochemical modelling into an optimization or optimal control context [Piwonski and Slawig, 2012]. Metos3D is founded on the PETSc library [Balay et al., 2011] as a basis for efficient, parallelized numerical methods. Moreover it's portable and relies only on software, which is freely and public available.
The application of marine transport is based on the Transport Matrix Method, whereas the idea and its use are presented in [Khatiwala et al., 2005], [Khatiwala, 2007] and [Khatiwala, 2008].

-- Quick start --

Assuming PETSc, version 3.1, is installed do the following:

1. Prepare the data, i.e. download an archive from https://github.com/metos3d/data, extract it and follow the instructions.

2. Prepare the model, i.e. download an archive from https://github.com/metos3d/model and extract it.

3. Prepare the simulation package, i.e. download an archive from https://github.com/metos3d/simpack, extract it and change into the package directory.

4. Link the data and model directories. For example:

$>
ln -s ../metos3d-data-v0.2 data
ln -s ../metos3d-model-v0.2 model

5. Set the PETSc environment variables. For example:

$>
. petsc/local.jserver.petsc.opt.txt

6. Compile the software with a chosen model. For example:

$>
make BGC=model/I-Cs/

7. Run executable with chosen options. For example:

$>
./metos3d-simpack-I-Cs.exe model/I-Cs/option/local.jserver.option.I-Cs.simpack.txt

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
