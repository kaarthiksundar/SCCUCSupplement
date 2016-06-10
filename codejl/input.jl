using MAT
using MatpowerCases

type Bus
    nodeID::Int
    kind::Int
    Pd::Vector{Float64}               # For unit commitment
    Qd::Float64
    Pg::Float64
    Qg::Float64
    Pgmax::Float64
    Qgmax::Float64
    Pgmin::Float64
    Qgmin::Float64
    pi2::Float64                 # Objective coefficient           
    pi1::Float64                 # Objective coefficient
    Farmids::Vector{Int}
    meanwindminusload::Vector{Float64}
    qobjcoeff::Float64
    genids::Vector{Int}          # Bus can contain more than one generator
    outlist::Vector{Int}         # outgoing line indices
    inlist::Vector{Int}          # incoming line indices
    function Bus(nodeID, kind, Pd, Qd)
        b = new(nodeID, kind, Pd, Qd)
        b.Pg = 0
        b.Qg = 0
        b.Pgmax = 0
        b.Qgmax = 0
        b.pi1 = 0
        b.pi2 = 0
        b.Farmids = Int[]
        b.meanwindminusload = Float64[]
        b.genids = Int[]
        b.outlist = Int[]
        b.inlist = Int[]
        return b
    end
end

# Functions for bus
function setg(b::Bus, genidx, Pg, Qg, Pgmax, Pgmin)
    b.Pg += Pg
    b.Qg += Qg
    b.Pgmax += Pgmax
    b.Pgmin += Pgmin
    if b.kind == 1
        warn("Generator $genidx was assigned to bus $(b.nodeID), but this bus has type 1")
    end
    b.kind = 2
    push!(b.genids,genidx)
end

function getBus(buses, nodeID) 
    for b in buses
        if (b.nodeID == nodeID)
            return b
        end
    end
end

getLoad(b::Bus) = b.Pd
getPg(b::Bus) = b.Pg
isgen(b::Bus) = b.kind == 2
getb(b::Bus, t) = b.Pg - b.Pd[t]
getLoadMinusMeanWind(b::Bus) = -b.meanwindminusload

function setfarm(b::Bus, farmid)
    println(">> node ", b.nodeID, " is farm ", farmid)
    push!(b.Farmids, farmid)
end

function setMeanWindMinusLoad(b::Bus, farms)
    b.meanwindminusload = -b.Pd
    for farmid in b.Farmids
        b.meanwindminusload += farms[farmid].mean
    end
end

setqobjcoeff(b::Bus, coeff) = (b.qobjcoeff = coeff)

type Generator
    genID::Int
    busidx::Int
    Pg::Float64
    Qg::Float64
    Pgmax::Float64
    Pgmin::Float64
    prodslope::Vector{Float64}
    b0::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    pi1::Float64
    pi2::Float64
    countinitoff::Int
    countiniton::Int
    noloadcost::Float64
    rampdown::Float64
    rampup::Float64
    downtime::Int
    uptime::Int
    outputatinit::Float64
    function Generator(genID, busidx, Pg, Qg, Pgmax, Pgmin) 
        g = new(genID, busidx, Pg, Qg, Pgmax, Pgmin)
        g.prodslope = Float64[]
        g.b0 = 0.0
        g.b1 = Pgmax/3
        g.b2 = 2*Pgmax/3
        g.b3 = Pgmax
        g.pi1 = 1.0
        g.pi2 = 1.0
        g.countinitoff = 0
        g.countiniton = 0
        g.noloadcost = 0.0
        g.rampdown = 0.0
        g.rampup = 0.0
        g.downtime = 0
        g.uptime = 0
        g.outputatinit = 0.0
        return g
    end
end

setPgmax(g::Generator, Pgmax) = (g.Pgmax = Pgmax)
setcountinitoff(g::Generator, countinitoff) = (g.countinitoff = countinitoff)
setcountiniton(g::Generator, countiniton) = (g.countiniton = countiniton)
setnoloadcost(g::Generator, noloadcost) = (g.noloadcost = noloadcost)
setrampdown(g::Generator, rampdown) = (g.rampdown = rampdown)
setrampup(g::Generator, rampup) = (g.rampup = rampup)
setdowntime(g::Generator, downtime) = (g.downtime = downtime)
setuptime(g::Generator, uptime) = (g.uptime = uptime)
setPgmin(g::Generator, Pgmin) = (g.Pgmin = Pgmin)
setoutputatinit(g::Generator, outputatinit) = (g.outputatinit = outputatinit)

function setprodcost(generators, prodcost, mvaBase)
    for i in 1:length(generators)
        push!(generators[i].prodslope, prodcost[1,i]*mvaBase)
        push!(generators[i].prodslope, prodcost[2,i]*mvaBase)
        push!(generators[i].prodslope, prodcost[3,i]*mvaBase)
    end
    return
end

type Line
    arcID::Int
    tail::Int # the "to" node
    head::Int # the "from" node
    y::Float64 # the susceptance value
    x::Float64 # the reactance value
    u::Float64 # the capacity of the line
    distance_scale::Float64 # this will be used to scale u
    Line(arcID, tail, head, y, x, u, d) = new(arcID, tail, head, y, x, u, d)
end

Line(arcID, tail, head, y, distance_scale) = Line(arcId, tail, head, y, 1/y, 0.0, distance_scale)

getThermalCapacity(l::Line, mvaBase) = l.u#/mvaBase  # line limits
getSyncCapacity(l::Line, mvaBase) = l.y 

type Farm
    node::Int
    mean::Vector{Float64}
    stddev::Vector{Float64}
    var::Vector{Float64}
    Farm(idx, mean, stddev) = new(idx, mean, stddev, stddev.^2)
end

type Scenario
    kind::Int # 1 for line scenario and 2 for generator scenario
    ID::Int # if kind is 1 then lineID else genID
    tail::Int # if kind is 1
    head::Int # if kind is 1
    Bhat::SparseMatrixCSC{Float64,Int64}
    Btilde::Matrix{Float64}
    Scenario(kind, ID, tail, head, Bhat, Btilde) = new(kind, ID, tail, head, Bhat, Btilde)
end


function readcase(casefilename, extras, loadscale, windscale, thermalLimitscale, mvaBase)

    # case = loadcase("case24_ieee_rts")
    case = loadcase("case96")

    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    busmat = case["bus"]
    busIDmap = Dict()
    
    buses = Bus[]
    for i in 1:size(busmat,1)
        nodeID = i
        busIDmap[busmat[i,1]] = i
        bustype = round(Int, busmat[i,2])
        Pd = busmat[i,3]
        Qd = busmat[i,4]
        Gs = busmat[i,5]
        Bs = busmat[i,6]
        area = busmat[i,7]
        Vm = busmat[i,8]
        Va = busmat[i,9]
        baseKV = busmat[i,10]
        zone = busmat[i,11]
        Vmax = busmat[i,12]
        Vmin = busmat[i,13]
        
        b = Bus(nodeID, bustype, Float64[], Qd./mvaBase)
        push!(buses, b)
    end

    ## generator data
    # bus Pg Qg	Qmax Qmin Vg mBase status Pmax Pmin Pc1 Pc2	Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ramp_q apf
    generatorlist = Int[]
    generators = Generator[]
    genmat = case["gen"]
    genucdata = readcsv(extras["gen_file"])
    removegenindices = readcsv(extras["remove_indices"]) # only when using DC power flow model
    gencounter = 1
    for i in 1:size(genmat,1)
        #if genmat[i,1] == removegenindices[1]
        if (genmat[i,1] == removegenindices[1] || genmat[i,1] == removegenindices[2] || genmat[i,1] == removegenindices[3])
            continue
        end
        
        busidx = busIDmap[genmat[i,1]]
        Pg = genmat[i,2]
        Qg = genmat[i,3]
        Pgmax = genmat[i,9]
        Pgmin = genmat[i,10]
        g = Generator(gencounter, busidx, Pg/mvaBase, Qg/mvaBase, Pgmax/mvaBase, Pgmin/mvaBase)
        # g = Generator(gencounter, busidx, Pg, Qg, Pgmax, Pgmin)
        setPgmax(g, genucdata[gencounter,2]./mvaBase)
        # setPgmax(g, genucdata[gencounter,2]./mvaBase)
        setcountinitoff(g, genucdata[gencounter,3])
        setcountiniton(g, genucdata[gencounter,4])
        setnoloadcost(g, genucdata[gencounter,5])
        setrampdown(g, genucdata[gencounter,6]./mvaBase)
        setrampup(g, genucdata[gencounter,7]./mvaBase)
        setdowntime(g, genucdata[gencounter,8])
        setuptime(g, genucdata[gencounter,9])
        setPgmin(g, genucdata[gencounter,10]./mvaBase)
        setoutputatinit(g, genucdata[gencounter,11]./mvaBase)
        push!(generators, g)
        push!(generatorlist, busidx)
        setg(buses[busidx], gencounter, Pg/mvaBase, Qg/mvaBase, g.Pgmax, g.Pgmin)
        # setg(buses[busidx], gencounter, Pg, Qg, g.Pgmax, g.Pgmin)
        gencounter = gencounter + 1
    end
    
    for i in 1:length(buses)
        if buses[i].kind == 2 && length(buses[i].genids) == 0
            warn("Bus $i is listed as a generator, but no corresponding generator information found, changing to kind 1")
            buses[i].kind = 1
        end
    end

    ## branch data
    # fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax
    branchmat = case["branch"]
    lines = Line[]
    for i in 1:size(branchmat,1)
        fbus = busIDmap[branchmat[i,1]]
        tbus = busIDmap[branchmat[i,2]]
        x = branchmat[i,4]
        y = 1/x
        u = branchmat[i,6]/mvaBase
        push!(buses[fbus].outlist, i)
        push!(buses[tbus].inlist, i)
        l = Line(i, tbus, fbus, y, x, u*thermalLimitscale, thermalLimitscale)
        push!(lines,l)
    end

    generatorlist = unique(generatorlist)
    
    costsfilename = "none"
    if costsfilename != "none"
        costsmat = readcsv(costsfilename, Float64)
        for i in 1:size(costsmat,1)
            busid = round(Int, costsmat[i,1])
            @assert length(buses[busid].genids) == 1
            buses[busid].pi2 = costsmat[i,2] # quadratic coefficient
            buses[busid].pi1 = costsmat[i,3] # linear coefficient
        end
    else
        gencost = case["gencost"]
        for g in 1:length(generators)
            generators[g].pi1 = gencost[g,5]
            generators[g].pi2 = gencost[g,6]*mvaBase
        end
    end
    
    numbuses = length(buses)

    #loads
    loads = readcsv(extras["load_file"])
    for i in 1:length(buses)
        buses[i].Pd = (loadscale.*loads[:,i+1]./mvaBase)
    end

    # Setting up windfarms
    farms = Farm[]
    winddata = readcsv(extras["wind_cap"])
    windmean = readcsv(extras["wind_mean"])
    windstddev = readcsv(extras["wind_stddev"])
    for i in 1:size(winddata,1)
        busidx = busIDmap[winddata[i,1]]
        mean = windscale.*windmean[:,i]./mvaBase
        stddev = windstddev[:,i]./mvaBase
        f = Farm(busidx, mean, stddev)
        push!(farms, f)
        setfarm(buses[busidx],i)
    end

    for b in buses
        setMeanWindMinusLoad(b, farms)
    end
    
    # startup cost for generators
    sucost = readcsv(extras["gen_startup"])
    prodcost = readcsv(extras["gen_prod"])
    setprodcost(generators, prodcost, mvaBase)
    
    # hydro_idx and nuclear_idx
    hydro_idx = Int[]
    nuclear_idx = Int[]
    hydro = readcsv(extras["hydro_indices"])
    nuclear = readcsv(extras["nuclear_indices"])

    for i in hydro
        push!(hydro_idx, busIDmap[i])
    end

    for i in nuclear
        push!(nuclear_idx, busIDmap[i])
    end
    
    @assert issubset(hydro_idx, generatorlist)
    @assert issubset(nuclear_idx, generatorlist)
    
    return generators, generatorlist, sucost, buses, lines, farms, hydro_idx, nuclear_idx
end


function readconfig(configfilename)
    println("\nreading config $configfilename")
    refbus = 0
    uniformalphas = false
    
    lines = readlines(open(configfilename,"r"))

    numlines = length(lines)

    windfilename = "none"
    costsfilename = "none"
    logfilename = "none"
    casefilename = "none"

    lookingforend = true
    line_probability_threshold = 0.0
    mvaBase = 100
    gen_probability_threshold = 0.0
    linebufferamt = 0
    loadscale = 1.0
    thermalLimitscale = 1.0
    windscale = 1.0

    extras = Dict()
    
    for l in lines
        startswith(l,'#') && continue
        
        thisline = split(l)
        length(thisline) > 0 || continue
        if thisline[1] == "END"
            break
        elseif thisline[1] == "case"
            casefilename = thisline[2]
        elseif thisline[1] == "wind"
            windfilename = thisline[2]
        elseif thisline[1] == "costs"
            costsfilename = thisline[2]
        elseif thisline[1] == "refbus"
            refbus = parse(Int,thisline[2])
        elseif thisline[1] == "line_probability_threshold"
            line_probability_threshold = float(thisline[2])
            println(">>>> line_probability_threshold = $line_probability_threshold")
        elseif thisline[1] == "gen_probability_threshold"
            gen_probability_threshold = float(thisline[2])
            println(">>>> gen_probability_threshold = $gen_probability_threshold")
        elseif thisline[1] == "linebufferamt"
            linebufferamt = float(thisline[2])
            println(">>>> linebufferamt = $linebufferamt")
        elseif thisline[1] == "uniformalphas"
            uniformalphas = true
            println(">>>> uniform alphas")
        elseif thisline[1] == "mvaBase"
            mvaBase = float(thisline[2])
            println(">>>> mvaBase = $mvaBase")
        elseif thisline[1] == "loadscale"
            loadscale = float(thisline[2])
            println(">>>> loadscale = $loadscale")
        elseif thisline[1] == "thermalLimitscale"
            thermalLimitscale = float(thisline[2])
            println(">>>> thermalLimitscale = $thermalLimitscale")
        elseif thisline[1] == "windscale"
            windscale = float(thisline[2])
            println(">>>> windscale = $windscale")
        elseif thisline[1] == "logfile"
            logfilename = thisline[2]
            println(">>>> logfilename = $logfilename")
        elseif thisline[1] == "factor" || thisline[1] == "rmatrix" || thisline[1] == "line_sync_nu"
            println(">>>> ignoring $(thisline[1])")
        else
            extras[thisline[1]] = chomp(thisline[2])
            println(">>>> $(thisline[1]) = $(thisline[2])")
        end
    end


    return casefilename, refbus, loadscale, windscale, thermalLimitscale, extras, line_probability_threshold, gen_probability_threshold, linebufferamt, mvaBase
end


