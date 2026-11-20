# ExtGram.jl

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://JanHab.github.io/ExtGram.jl/)
[![Build](https://img.shields.io/github/actions/workflow/status/JanHab/ExtGram.jl/documentation.yml?label=build)](https://github.com/JanHab/ExtGram.jl/actions)
[![GitLab Version](https://img.shields.io/badge/version-1.0.1-blue.svg)](https://github.com/JanHab/ExtGram.jl/releases/tag/v1.0.1)
[![License](https://img.shields.io/badge/license-CCBY4.0-blue)](https://github.com/JanHab/ExtGram.jl/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18998361.svg)](https://doi.org/10.5281/zenodo.18998361)

*A discontinuous Galerkin implementation of the extended Gramian closure for moment equations with Trixi.jl.*

## Abstract
The classical continuum mechanical equations of Navier-Stokes and Fourier (NSF) accurately describe gases with a high collision frequency, measured by a small Knudsen number (Kn).
However, the NSF equations fail for gases with a high Knudsen number and, consequently, a small collision frequency.
Directly solving Boltzmann's equation is an alternative approach of high accuracy, but is computationally very demanding due to its high dimensionality.
Another method is to approximate the Boltzmann equation by directly solving for macroscopic properties, such as density, velocity, and temperature, which can be described by moments.
Therefore, the so-called moment system is solved, which is an infinitely large system of hierarchically coupled partial differential equations (PDEs), with the moments as unknowns.
The i-th PDE depends on the i-th moment and is coupled to the (i+1)-st PDE, again depending on the (i+1)-st moment.
Directly solving the moment system is impossible, as it is infinitely large and always has one unknown more than equations.
Therefore, the moment system is truncated after M equations, requiring a moment closure for the (M+1)-st moment.

This is an implementation of the extended Gramian closure, a recently invented moment closure based on orthogonal polynomials.
For comparison reasons, it also has implementations of the Morin-McDonald closure, and Grad's closure.
The extended Gramian closure has attractive structure-preserving properties and overcomes shortcomings of classical closures, such as the lack of global hyperbolicity in Grad's closure and the computational effort required by the maximum entropy closure.
While the closure is constructed in one dimension, a projection method for applying it to a three-dimensional velocity phase, based on matrix rotations, is implemented.
This allows to go from a 1D-1V setting to 1D-3V.

## Main Features

- Solving the 1D(space) - 1V(velocity) moment equations 
  - Shock tube
  - Shock structure
  - Two beam
- Discrete velocity method (DVM) as reference solution for 1D(space) - 1V(velocity)
  - Shock tube
  - Shock structure
- Extending the equations to 1D(space) - 3V(velocity) in slab-geometry and full three velocity dimensions
- Solving the 1D(space) - 1V(velocity) moment equations, coupled to the Vlasov-Poisson equation with initial conditions for:
  - Landau Damping
  - Two-stream instability
- Numerical grid convergence
  - On periodic boundary conditions with initial conditions from Landau damping (no shock)

## Installation

Clone the git repository and install the `ExtGram.jl` package locally with:

``` bash
git clone git@github.com:JanHab/ExtGram.jl.git
cd ExtGram.jl
julia --project -e 'import Pkg; Pkg.instantiate()'
```

The code was developed with Julia version 1.12.6.
The solutions are stored in the `out` folder (in .gitignore).
Before executing the example files, make sure that you create all necessary subfolders, where the solutions are stored.
The dependency KitBase.jl (used only for the kinetic reference solutions) builds a small Python environment during installation and requires `pip3` version 23 or newer on the `PATH`; if its build of `FiniteMesh` fails, the package still works and the error can be ignored.

```bash
mkdir -p out/Shock/{Moments,DVM,GaugeInvestigation} out/Reference out/CompareClosure \
  out/GridConvergencePeriodic out/VlasovPoisson/{LandauDamping,TwoStreamInstability} \
  out/1D3V/{1D3V_ShockTube,1D3V_TwoBeam,1D3Vfull_ShockTube} \
  out/Figures/{1D3V/slab,1D3V/full,CompareClosure,EigenvaluesEquilibrium,GaugeInvestigation,LandauDamping,ShockStructure,ShockTube,TwoStreamInstability}
```
Further, Python (with Jupyter notebooks) is used for visualizations.
Install the requirements from the `requirements.txt` file with (tested for Python version 3.12.2) `pip install -r requirements.txt`.
The visualization notebooks use the `science` plot package, which requires a local LaTeX installation.
You can always comment the Python import.

## Testing

Clone the repository and execute the testing pipeline with:

``` bash
julia --project test/tests.jl
```

Further, you can test the installation by running a shock tube example with:

``` bash
julia --project examples/ShockTubeExample.jl
```

## Usage

Find the package source code in [src](https://github.com/JanHab/ExtGram.jl/tree/main/src), which implements the momentum closure and the system of PDEs.
The DVM reference curves (for 1D-3V) are computed with [KitBase.jl](https://github.com/vavrines/KitBase.jl).

The folder [examples](https://github.com/JanHab/ExtGram.jl/tree/main/examples) contains the main files for physically relevant benchmark problems.
[slurm](https://github.com/JanHab/ExtGram.jl/tree/main/slurm) gives shell scripts for executing the benchmark problems either on your local laptop or on a HPC cluster with a slurm script. You need to either comment out the `os.system(command)` or the `createsbatch` file, depending on your system.

[miscellaneous](https://github.com/JanHab/ExtGram.jl/tree/main/miscellaneous) contains miscellaneous files used for implementation and checking.
Here, you can also find a symbolic proof, that the [Morin-McDonald](https://doi.org/10.1016/j.jcp.2024.113659) closure is hyperbolic for M=3.

[notebooks](https://github.com/JanHab/ExtGram.jl/tree/main/notebooks) has all the visualizations.
The results are stored in `out`.
Note that this folder is in `.gitignore`. 
The visualizations will only work if you have the necessary data.
<!-- You can find minimal working data archived with Zenodo here:  [https://doi.org/10.5281/zenodo.18998361](https://doi.org/10.5281/zenodo.18998361)   -->

Finally, [test](https://github.com/JanHab/ExtGram.jl/tree/main/test) contains a small testing pipeline.

## Referencing 

If you find these results useful or use the implementations provided here, please cite the corresponding article.

```bibtex
@thesis{habscheid_modeling_2026,
  title = {Modeling and analysis of nonlinear closures in hierarchical moment equations using discontinuous {Galerkin} methods},
  type = {mathesis},
  institution = {RWTH Aachen University},
  url = {https://publications.rwth-aachen.de/record/1035026},
  doi = {10.18154/RWTH-2026-04820},
  author = {Habscheid, Jan},
  month = mar,
  year = {2026},
}
```

If you use the implementations provided here, please also cite this repository as

```bibtex
@software{habscheid_2026_18998361,
  author       = {Habscheid, Jan and
                  Yilmaz, Eda and
                  Geratz, Matthias and
                  Torrilhon, Manuel},
  title        = {ExtGram.jl},
  month        = mar,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {1.0.0},
  doi          = {10.5281/zenodo.18998361},
  url          = {https://doi.org/10.5281/zenodo.18998361},
}
```

## Contact

- **Jan Habscheid**
  - TIMA - Applied Mathematics
  - Linköping University
  - [Jan.Habscheid@liu.se](mailto:Jan.Habscheid@liu.se)
- **Eda Yilmaz**
  - ACoM - Applied and Computational Mathematics
  - RWTH Aachen University
  - [yilmaz@acom.rwth-aachen.de](mailto:yilmaz@acom.rwth-aachen.de)
- **Matthias Geratz**
  - ACoM - Applied and Computational Mathematics
  - RWTH Aachen University
  - [geratz@acom.rwth-aachen.de](mailto:geratz@acom.rwth-aachen.de)
- **Prof. Dr. Manuel Torrilhon**
  - ACoM - Applied and Computational Mathematics
  - RWTH Aachen University
  - [mt@acom.rwth-aachen.de](mailto:mt@acom.rwth-aachen.de)