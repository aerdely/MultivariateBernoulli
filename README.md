# MultivariateBernoulli
A subcopula characterization of dependence for the Multivariate Bernoulli Distribution.

> Author: Arturo Erdely

> Preprint: [arXiv](https://drive.google.com/file/d/1bIAhVZOcmEsMejxec8j1luH5kCUCCcXy/view?usp=sharing)

### Instructions for reproducibility

1. Download and install the [Julia](https://julialang.org/downloads/) programming language.
2. Download the code files clicking in the green button `<> Code` of this GitHub repository and `DownloadZIP`. Unzip the downloaded file and move the following files into a directory of your choice: `Packages.jl`, `MultivariateBernoulli.jl`, `Examples.jl`, and `covid2020.csv`
3. Open the `Julia` terminal and change to the working directory where you unzipped the files. You may do this by defining a string variable `path` with the path to the files directory and then execute in the terminal `cd(path)`. For example, in the operating system *Windows* it may look something like:
   ```julia
   path = "D:/MyFiles/mychoice"
   cd(path)
   readdir()
   ```
5. Install the required packages by executing the following command in the `Julia` terminal:
   ```julia
   include("Packages.jl")
   ```
6. Generate figures and example calculations by executing:
   ```julia
   include("Examples.jl")
   ```
