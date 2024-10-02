### Packages for Multivariate Bernoulli Distribution
### Author: Dr. Arturo Erdely
### Version: 2024-09-29

# Install required Julia packages

using Pkg

begin
    paquete = ["Distributions", "Plots", "LaTeXStrings", "CSV", "DataFrames"]
    for p in paquete
        println("*** Installing package: ", p)
        Pkg.add(p)
    end
    println("*** End of package list. \n")
end
