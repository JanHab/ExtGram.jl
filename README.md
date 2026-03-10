# ExtGram.jl

*A discontinuous Galerkin implementation of the extended Gramian closure for moment equations with Trixi.jl*

## Abstract (from the thesis)

The classical continuum mechanical equations of Navier-Stokes and Fourier (NSF) accurately describe gases with a high collision frequency, measured by a small Knudsen number (Kn).
However, the NSF equations fail for gases with a high Knudsen number and, consequently, a small collision frequency.
Directly solving Boltzmann's equation is an alternative approach of high accuracy, but is computationally very demanding due to its high dimensionality.
Another method is to approximate the Boltzmann equation by directly solving for macroscopic properties, such as density, velocity, and temperature, which can be described by moments.
Therefore, the so-called moment system is solved, which is an infinitely large system of hierarchically coupled partial differential equations (PDEs), with the moments as unknowns.
The i-th PDE depends on the i-th moment and is coupled to the (i+1)-st PDE, again depending on the (i+1)-st moment.
Directly solving the moment system is impossible, as it is infinitely large and always has one unknown more than equations.
Therefore, the moment system is truncated after M equations, requiring a moment closure for the (M+1)-st moment.

This thesis discusses the extended Gramian closure, a recently invented moment closure based on orthogonal polynomials.
The extended Gramian closure has attractive structure-preserving properties and overcomes shortcomings of classical closures, such as the lack of global hyperbolicity in Grad's closure and the computational effort required by the maximum entropy closure.
This work introduces the theoretical framework of the extended Gramian closure and suggests an updated closure under the condition that the equations are truncated with M being odd.
Rigorous mathematical proofs for the structure-preserving properties of this updated closure are discussed. 
While the closure is constructed in one dimension, a projection method for applying it to a three-dimensional velocity phase, based on matrix rotations, is introduced.
Further, the extended Gramian closure is numerically solved with the discontinuous Galerkin method and its accuracy is analyzed and compared with that of the Gramian closure and Grad's closure.
Therefore, various benchmark problems, including the shock tube, shock structure, Landau damping, and two-stream instability, are analyzed.

## Main Features

- Solving the 1D(space) - 1D(velocity) moment equations 
  - Shock tube problem
  - Shock structure problem
  - Two beam problem was tested, but not yet included
- Discrete velocity method (DVM) as reference solution for 1D(space) - 1D(velocity)
  - Shock tube
  - Shock structure
- Extending the equations to 1D(space) - 3D(velocity) for *M=4*
  - Tested for the shock tube problem
- Solving the 1D(space) - 1D(velocity) moment equations, coupled to the Vlasov-Poisson equation with initial conditions for:
  - Landau Damping
  - Two-stream instability
- Numerical grid convergence
  - On periodic boundary conditions with initial conditions from Landau damping (no shock)

## Installation

First, clone the git repository

``` bash
git clone git@git.rwth-aachen.de:JanHab/ExtGram.jl.git
cd ExtGram.jl
```

<!-- We need the `Revise` julia package to execute most of the scripts.
In principle, this is not necessary for the numerical backend, but was used to accelerate implementations.

``` julia
using Pkg;
Pkg.add(Revise);
Pkg.instantiate()
``` 
-->

Now we can install the `ExtGram` package locally with:

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

### Documentation

Generate the docs with:

``` bash
julia --project docs/make.jl
```

## Usage

Find the package source code in [src](https://git.rwth-aachen.de/JanHab/ExtGram.jl/-/tree/main/src?ref_type=heads), which implements the momentum closure and the system of PDEs.
The subfolder [src/TrixiTree2Triangulation](https://git.rwth-aachen.de/JanHab/ExtGram.jl/-/tree/main/src/TrixiTree2Triangulation?ref_type=heads) is a subpackage developed by Matthias Geratz for storing the solution in a *.tsv file.

The folder [examples](https://git.rwth-aachen.de/JanHab/ExtGram.jl/-/tree/main/examples?ref_type=heads) contains some main files for physical relevant benchmark problems.
[slurm](https://git.rwth-aachen.de/JanHab/ExtGram.jl/-/tree/main/slurm?ref_type=heads) gives shell scripts for executing the benchmark problems either on your local laptop or on a HPC server with a slurm script. You need to either comment the `os.system(command)` or `createsbatch` file, depending on your system.

[miscellaneous](https://git.rwth-aachen.de/JanHab/ExtGram.jl/-/tree/main/miscellaneous?ref_type=heads) contains miscellaneous files used for implementation and checking.

[notebooks](https://git.rwth-aachen.de/JanHab/ExtGram.jl/-/tree/main/notebooks?ref_type=heads) has all the visualizations.

Finally, [test](https://git.rwth-aachen.de/JanHab/ExtGram.jl/-/tree/main/test?ref_type=heads) contains a small testing pipeline.

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