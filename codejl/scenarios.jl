using Graphs
using MAT
include("matrix.jl")

function CreateScenarios(generators, lines, buses, farms, hydro_idx, nuclear_idx, refbus)
    numbuses = length(buses)
    numlines = length(lines)
    numfarms = length(farms)
    numgenerators = length(generators)
    kind = Int[]
    IDs = Int[]
    tails = Int[]
    heads = Int[]
    Bhats = SparseMatrixCSC{Float64, Int64}[]
    Btildes = Matrix{Float64}[]
    
    # Creating line scenario B and Btilde
    count = 0
    for line in lines
        g = simple_graph(numbuses, is_directed=false)
        for i in lines
            if (i==line)
                continue
            end
            add_edge!(g, i.head, i.tail)
        end
        if (length(connected_components(g)) != 1)
            continue
        end

        push!(kind, 1)
        push!(IDs, line.arcID)
        push!(tails, line.tail)
        push!(heads, line.head)

        B = genBscenario(buses, lines, line)
        Bhat = remove_col_and_row(B,refbus)
        Bhatsp = sparse(Bhat)
        Blu = lufact(Bhatsp)
        Btilde = zeros(numbuses, numbuses)
        rhs = zeros(numbuses-1)
        for b in buses
            if (b.nodeID == refbus) 
                colbarinv = Array(Float64, numbuses)
                colbarinv = zeros(numbuses)
                Btilde[:,refbus] = colbarinv
                continue
            end
            rhsindex = b.nodeID
            if (b.nodeID > refbus)
                rhsindex -= 1
            end
            rhs[rhsindex] = 1
            colbarinv = Array(Float64, numbuses)
            colbarinv[1:numbuses .!= refbus] = Blu\rhs
            colbarinv[refbus] = 0
            rhs[rhsindex] = 0
            Btilde[:,b.nodeID] = colbarinv
        end
        push!(Bhats, Bhatsp)
        push!(Btildes, Btilde)
        count = count + 1
    end
    
    # creating generator scenario B and Btilde
    B = genB(buses, lines)
    Bhat = remove_col_and_row(B, refbus)
    Bhatsp = sparse(Bhat)
    Blu = lufact(Bhatsp)
    Btilde = zeros(numbuses, numbuses)
    rhs = zeros(numbuses-1)
    for b in buses
        if (b.nodeID == refbus) 
            colbarinv = Array(Float64, numbuses)
            colbarinv = zeros(numbuses)
            Btilde[:,refbus] = colbarinv
            continue
        end
        rhsindex = b.nodeID
        if (b.nodeID > refbus)
            rhsindex -= 1
        end
        rhs[rhsindex] = 1
        colbarinv = Array(Float64, numbuses)
        colbarinv[1:numbuses .!= refbus] = Blu\rhs
        colbarinv[refbus] = 0
        rhs[rhsindex] = 0
        Btilde[:,b.nodeID] = colbarinv
    end

    for g in generators
        if (g.genID in hydro_idx || g.genID in nuclear_idx)
            continue
        end
        push!(kind, 2)
        push!(IDs, g.genID)
        push!(tails, 0)
        push!(heads, 0)
        push!(Bhats, Bhatsp)
        push!(Btildes, Btilde)
        count = count + 1
    end
   
    println(">>>> Number of contigencies: $(count)")
    mfile = matopen("scenariodata.mat", "w")
    write(mfile, "kind", kind)
    write(mfile, "IDs", IDs)
    write(mfile, "tails", tails)
    write(mfile, "heads", heads)
    write(mfile, "Bhats", Bhats)
    write(mfile, "Btildes", Btildes)
    close(mfile)

end

