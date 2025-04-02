# UVLJscatter

## About

UVLJscatter is a Monte-Carlo simulation program that models electron transport in liquid microjets.

*[Detailed description of code]*


## Installation

### Cloning the repository
This repository may be cloned to your local machine using git:

```
git clone https://github.com/edosim/UVLJscatter.git
```

It can also be downloaded from <https://github.com/edosim/UVLJscatter.git>


### Julia setup
The Monte-Carlo simulation in [basisFunctionCode](basisFunctionCode) uses Julia and has only been tested with version 1.7.0.
Julia 1.7.0 may be downloaded and installed from [here](https://julialang.org/downloads/oldreleases/).
The installation may be verified by running `julia --version` in the terminal or command prompt.

The required dependencies are included in `Project.toml` and `Manifest.toml` and may be installed by running 

```
julia
]
activate .
instantiate
```
in [basisFunctionCode](basisFunctionCode) from the terminal or command prompt.


### Python setup

The retrieval code in [Analysis](Analysis) uses Python and Jupyter notebooks. One way to install and manage the packages in Python is to use [Anaconda Distribution](https://www.anaconda.com/download/success).

[`requirements.txt`](Analysis/requirements.txt) may be used to create a conda environment with all the required dependencies: 

```
conda create --name <env> --file <this file>
conda activate <env>
```


## Usage

### Monte-Carlo simulation

The code in [basisFunctionCode](basisFunctionCode) is used to run the Monte-Carlo simulation and apply concentration profiles to the results.

The file [input_scattering.csv](basisFunctionCode\input_scattering.csv) includes the input required to run the code and metadata input by the user.

The file [ScatteringCode.jl](basisFunctionCode\ScatteringCode.jl) runs the simulation and returns an `h5` file with histograms of final kinetic energies as a function of initial kinetic energy and initial depth. The output is saved in the folder [bases](basisFunctionCode\bases) as `allData.h5`. Histograms obtained by integrating over equally weighed depths, depths weighed with a Gaussian with mean 1 nm and standard deviation 1 nm and depths weighed with an exponential distribution with mean 1 are also saved as `basis_uni.txt`, `basis_gau.txt`, `basis_exp.txt`.

The file [applyConcProfile.jl](basisFunctionCode\applyConcProfile.jl) scales the  output of [ScatteringCode.jl](basisFunctionCode\ScatteringCode.jl) by a depth profile and integrates it over depth. This creates histograms of final kinetic energies as a function of initial kinetic energy. The depth profiles used are described by the sum of a vertical offset and a gaussian centred at a depth of 0.1 nm and a FWHM of 0.4 nm. The ratio between the height of the Guassian and the offset is varied to obtain a series of histograms. The mean, FWHM and ratios employed in the concentration profiles can be changed in the file. The output is saved in the folder [bases](basisFunctionCode\bases) as `basisConcProfiles.h5`.

The file [ScatteringCode_angles.jl](basisFunctionCode\ScatteringCode_angles.jl) is similar to [ScatteringCode.jl](basisFunctionCode\ScatteringCode.jl) but returns histograms of incidence angles of escaped electrons as a function of initial kinetic energy and initial depth. The output is saved in the folder [bases](basisFunctionCode\bases) as `allAngles.h5`. Histogram obtained by integrating over equally weighed depths is also saved as `basis_angles.txt`.

### Retrieval of true photoelectron spectra

The code in [Analysis](Analysis) is used to retrive true photoelectron spectra from experimental ones. The notebook [SpectralRetrieval.ipynb](Analysis\SpectralRetrieval.ipynb) contains the code and instructions for the retrieval. [ShowH5structure.ipynb](Analysis\ShowH5structure.ipynb) shows the structure and metadata of the simulation results.
