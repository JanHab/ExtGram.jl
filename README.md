# TrixiHyQMOM

*A discontinuous Galerkin implementation of the extended Gramian closure for moment equations*

## Physical Background

- To-Do

## Main Features

- Solving the 1D(space) - 1D(velocity) moment equations on a Riemann problem
- Numerical grid convergence
- Solving the 1D(space) - 1D(velocity) moment equations for the Vlasov-Poisson system with initial conditions:
    - Landau Damping
    - Two-stream instability
- Extending the equations to 1D(space) - 3D(velocity) for *M=4*

## Installation

First, clone the git repository

``` bash
git clone git@git.rwth-aachen.de:JanHab/trixihyqmom.git
cd trixihyqmom
```

We need the `Revise` julia package to execute most of the scripts.
In principle, this is not necessary for the numerical backend, but was used to accelerate implementations.

``` julia
using Pkg;
Pkg.add(Revise);
Pkg.instantiate()
```

Now we can install the `HyQMOM` package locally with:

``` bash

julia --project
using Pkg;
Pkg.instantiate()
```

- ToDo: Test this out!

### Testing

Clone the repository and execute the testing pipeline with:

``` bash
julia test/tests.jl
```

## Usage

Find the package source code in [src](https://git.rwth-aachen.de/JanHab/trixihyqmom/-/tree/main/src?ref_type=heads), which implements the momentum closure and the system of PDEs.
The subfolder [src/TrixiTree2Triangulation](https://git.rwth-aachen.de/JanHab/trixihyqmom/-/tree/main/src/TrixiTree2Triangulation?ref_type=heads) is a subpackage developed by Matthias Geratz for storing the solution in a *.tsv file.

The folder [examples](https://git.rwth-aachen.de/JanHab/trixihyqmom/-/tree/main/examples?ref_type=heads) contains some main files for physical relevant benchmark problems
[slurm](https://git.rwth-aachen.de/JanHab/trixihyqmom/-/tree/main/slurm?ref_type=heads) gives shell scripts for executing the benchmark problems either on your local laptop or on a HPC server with a slurm script. You need to either comment the `os.system(command)` or `createsbatch` file, depending on your system.

[miscellaneous](https://git.rwth-aachen.de/JanHab/trixihyqmom/-/tree/main/miscellaneous?ref_type=heads) contains miscellaneous files used for implementation and checking.

[notebooks](https://git.rwth-aachen.de/JanHab/trixihyqmom/-/tree/main/notebooks?ref_type=heads) has all the visualizations.

Finally, [test](https://git.rwth-aachen.de/JanHab/trixihyqmom/-/tree/main/test?ref_type=heads) has the testing pipeline.

## Contact

- **Jan Habscheid** (Thesis)
  - [Jan.Habscheid@rwth-aachen.de](mailto:Jan.Habscheid@rwth-aachen.de)
- **Eda Yilmaz** (Supervisor)
  - ACoM - Applied and Computational Mathematics
  - RWTH Aachen University
  - [yilmaz@acom.rwth-aachen.de](mailto:yilmaz@acom.rwth-aachen.de)
- **Matthias Geratz** (Initial Development)
  - ACoM - Applied and Computational Mathematics
  - RWTH Aachen University
  - [geratz@acom.rwth-aachen.de](mailto:geratz@acom.rwth-aachen.de)
- **Prof. Dr. Manuel Torrilhon** (Supervising Professor)
  - ACoM - Applied and Computational Mathematics
  - RWTH Aachen University
  - [mt@acom.rwth-aachen.de](mailto:mt@acom.rwth-aachen.de)