# SCCUCSupplement
Repository containing supplementary data and code for the paper "Unit Commitment with N-1 Security and Wind Uncertainty". [arxiv](http://arxiv.org/abs/1602.00079)

*Installation instructions:*

The optimization model was implemented by using the [JuMPChance](https://github.com/mlubin/JuMPChance.jl) extension to [JuMP](https://github.com/JuliaOpt/JuMP.jl) in the [Julia](http://julialang.org/downloads/) programming language.
Additionally, we used CPLEX 12.2 in our numerical experiments. [CPLEX](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) is a commercial solver which must be installed and licensed separately (one may easily use a different solver if CPLEX is not available, see the JuMP documentation).

The experiments require Julia 0.4.2 or later, and the following Julia packages:
- [JuMP](https://github.com/JuliaOpt/JuMP.jl) 0.11.1
- [JuMPChance](https://github.com/mlubin/JuMPChance.jl) 0.2.2
- [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl) 0.1.0
- [MAT](https://github.com/simonster/MAT.jl)
- [JLD](https://github.com/JuliaIO/JLD.jl)
- [Graphs](https://github.com/JuliaLang/Graphs.jl)
- [MatpowerCases.jl](https://github.com/kersulis/MatpowerCases.jl) 0.3.2

You should force the use of particular versions of these Julia packages with
```
julia> Pkg.pin("JuMP", v"0.11.1")
julia> Pkg.pin("JuMPChance", v"0.2.2")
julia> Pkg.pin("CPLEX", v"0.1.0")
julia> Pkg.pin("MatpowerCases", v"0.3.2")
```

*Running the code:*

The code for the experiments is contained in the ``codejl`` directory. The file ``input.jl`` contains routines and data structures for processing the input, and the file ``sccuc_simulation.jl`` contains the main simulation logic and optimization model. The ``scenarios.jl`` and ``matrix.jl`` contain certain additional functions to generate the scenarios and bus admittance matrices, respectively.

You can run the model by entering the ``run`` directory and executing:
```
julia ../codejl/sccuc_simulation.jl case96.dat
```

The output of the simulation is a ``.jld`` file (``case96.jld``) which can be opened directly loaded in Julia via the JLD package. Contained in the ``.jld`` file are the optimal objective values, solution times, solution status, and optimal values for the decision variables. 

The file ``case96.dat`` specifies all of the input paths and parameters for the simulation. In particular, one can modify the parameters for different runs.


Further documentation is available on request.