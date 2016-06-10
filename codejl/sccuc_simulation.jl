# This model refers to the paper: {place holder for the reference}.
# Interested readers are referred to the paper for further description of the model.

# Define packages to be used:
using JuMP # Optimization package
using JuMPChance # Chance constraints package
using CPLEX # Solver, could be altered if needed
using MAT # Package to interface with .mat files

include("input.jl")
include("model.jl")
include("matrix.jl")
include("scenarios.jl")

# Read the input data:
casefilename, refbus, loadscale, windscale, thermalLimitscale, extras, ϵij, ϵg, linebufferamt, mvaBase = readconfig(ARGS[1])

# Read specific data for the test case:
generators, generatorlist, sucost, buses, lines, farms, hydro_idx, nuclear_idx = readcase(casefilename, extras, loadscale, windscale, thermalLimitscale, mvaBase)

sucost = sucost'

# Define the set of indices for hydro generators. Will be used to enforce the uniform participation factors.
fixed_alpha_subsets = collect(hydro_idx)

# Extract the name of the  out file from the outputs of readconfig. (line 142):
output_file = extras["output_file"]

numbuses = length(buses)
numlines = length(lines)
numfarms = length(farms)
numgenerators = length(generators)

println(">>>> Buses: $(numbuses), Lines: $(numlines), Farms: $(numfarms)")
println(">>>> Generators: $(numgenerators)")
println(">>>> # Hydro generators: $(length(hydro_idx)), genIDs: $(hydro_idx)")
println(">>>> # Nuclear generators: $(length(nuclear_idx)), genIDs: $(nuclear_idx)")

scenariofilename = "none"

if scenariofilename == "none"
    CreateScenarios(generators, lines, buses, farms, hydro_idx, nuclear_idx, refbus)
end

scenarios = Scenario[]
scenariodata = matread("scenariodata.mat")
kind = scenariodata["kind"]
IDs = scenariodata["IDs"]
heads = scenariodata["heads"]
tails = scenariodata["tails"]
Bhats = scenariodata["Bhats"]
Btildes = scenariodata["Btildes"]
numscenarios = length(kind)
@assert size(heads) == size(kind)
@assert size(IDs) == size(kind)
@assert size(Bhats) == size(kind)
@assert size(Btildes) == size(kind)
@assert size(tails) == size(kind)

for i in 1:numscenarios
   s = Scenario(kind[i], IDs[i], tails[i], heads[i], Bhats[i], Btildes[i])
   push!(scenarios, s)
end

println(">>>> Number of contigencies: $(numscenarios)")

# Factorization, this only needs to be done once
B = genB(buses, lines)

#=
sumvar = 0.0
for f in farms
    sumvar += f.stddev^2
end
println("-> sumvar ", sumvar, " stddev ", sqrt(sumvar))
=#

violations = Set{Int}()
solvemodel(refbus, ϵij, ϵg, linebufferamt, mvaBase, loadscale,
    thermalLimitscale, generators, buses, lines, farms,
    scenarios, generatorlist, sucost, B, hydro_idx, nuclear_idx, extras, violations)

#createandsolvemodel(refbus, ϵij, ϵg, linebufferamt, mvaBase, loadscale,
#                    thermalLimitscale, generators, buses, lines, farms,
#                    scenarios, generatorlist, sucost, B, hydro_idx, nuclear_idx, extras)
