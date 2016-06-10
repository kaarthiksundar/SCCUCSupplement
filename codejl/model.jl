using MAT
using JLD
using Distributions

function solvemodel(refbus, ϵij, ϵg, linebufferamt, mvaBase, loadscale, thermalLimitscale, generators, buses, lines, farms, scenarios, generatorlist, 
    sucost, B, hydro_idx, nuclear_idx, extras, violations)

    lineviolations = solveinnerproblem(refbus, ϵij, ϵg, linebufferamt, mvaBase, loadscale, thermalLimitscale, generators, buses, lines, farms, scenarios, 
    generatorlist, sucost, B, hydro_idx, nuclear_idx, extras, violations)

    while length(lineviolations) != 0
        violations = union(violations, lineviolations)
		println(" Before next solve: $lineviolations")
	    lineviolations = solveinnerproblem(refbus, ϵij, ϵg, linebufferamt, mvaBase, loadscale, thermalLimitscale, generators, buses, lines, farms, scenarios, 
	    generatorlist, sucost, B, hydro_idx, nuclear_idx, extras, violations)
		println("After solve: $lineviolations")
        # break
    end


end

function solveinnerproblem(refbus, ϵij, ϵg, linebufferamt, mvaBase, loadscale, thermalLimitscale,
    generators, buses, lines, farms, scenarios, generatorlist, sucost, B, hydro_idx, nuclear_idx, extras, violations)

    numbuses = length(buses)
    numlines = length(lines)
    numfarms = length(farms)
    numscenarios = length(scenarios)
    numgenerators = length(generators)
    Bhat = remove_col_and_row(B,refbus)
    Bhatsp = sparse(Bhat)

    numlinescenarios = 0
    numgenscenarios = 0
    for k in 1:numscenarios
        if (scenarios[k].kind == 1)
            numlinescenarios += 1
        else
            numgenscenarios += 1
        end
    end

    println(">>>> Line contigencies : $(numlinescenarios)")
    println(">>>> Generator contigencies : $(numgenscenarios)")

    @assert scenarios[numscenarios].kind == 2

    Btilde = scenarios[numscenarios].Btilde 

    # UC paramter calculations
    T = 24
    numblocks = 8
    onoffinit = zeros(numgenerators) # Initially assume every generator is off at t=0 and update
    Lupmin = zeros(numgenerators)
    Ldownmin = zeros(numgenerators)
    for i in 1:numgenerators
        onoffinit[i] = ((generators[i].countiniton > 0) ? 1 : 0)
        Lupmin[i] = min(T, (generators[i].uptime)*onoffinit[i]) # set the intial up time based on onoffinit
        Ldownmin[i] = min(T, (generators[i].downtime)*(1 - onoffinit[i])) # set the intial down time based on onoffinit
    end

    #
    m = ChanceModel(solver=CplexSolver(CPX_PARAM_EPGAP=0.01, CPX_PARAM_TILIM=3600, CPX_PARAM_VARSEL=4, CPX_PARAM_FLOWCOVERS=1,
    CPX_PARAM_LPMETHOD=4, CPX_PARAM_MIPDISPLAY=4))

    tic()
    # Binary Variables
    @defVar(m, x[1:numgenerators, 1:T], Bin)
    @defVar(m, y[1:numgenerators, 1:T], Bin)
    @defVar(m, z[1:numgenerators, 1:T], Bin)
    @defVar(m, suc[1:numgenerators, 1:numblocks, 1:T], Bin)

    # Continuous variables
    @defVar(m, xp[1:numgenerators, 1:T] >= 0)
    @defVar(m, xr[1:numgenerators, 1:T] >= 0)
    @defVar(m, pseg[1:numgenerators, 1:T, 1:3])
    @defVar(m, p[1:numgenerators, 1:T])
    @defVar(m, r⁺[1:numgenerators, 1:T] >= 0)
    @defVar(m, r⁻[1:numgenerators, 1:T] >= 0)
    @defVar(m, surrogate[1:T] >= 0)
    @defVar(m, α[1:numgenerators, 1:T] >= 0)
    @defVar(m, SC[1:numgenerators, 1:T])
    @defVar(m, pᵇ[i=1:numbuses, 1:T])
    @defVar(m, αᵇ[i=1:numbuses, 1:T] >= 0)
    @defVar(m, pcurtail[1:T] >= 0)
    @defVar(m, delta[1:numbuses, 1:T])
    @defVar(m, theta[1:numbuses, 1:T])
    #

    # Random variables
    @defIndepNormal(m, ω[i=1:numfarms, t=1:T], mean=0, var=farms[i].var[t])

    # Mccormick for xp, xr⁺ variables
    @addConstraint(m, xp1[i=1:numgenerators, t=1:T], xp[i,t] <= p[i,t])
    @addConstraint(m, xp2[i=1:numgenerators, t=1:T], xp[i,t] <= generators[i].Pgmax*x[i,t])
    @addConstraint(m, xp3[i=1:numgenerators, t=1:T], xp[i,t] >= p[i,t] - generators[i].Pgmax*(1-x[i,t]))
    @addConstraint(m, xr1[i=1:numgenerators, t=1:T], xr[i,t] <= r⁺[i,t])
    @addConstraint(m, xr2[i=1:numgenerators, t=1:T], xr[i,t] <= generators[i].Pgmax*x[i,t]/2)
    @addConstraint(m, xr3[i=1:numgenerators, t=1:T], xr[i,t] >= r⁺[i,t] - generators[i].Pgmax*(1-x[i,t])/2)

    # Reserve bounds
    @addConstraint(m, resbounds1[i=1:numgenerators, t=1:T], sum{x[j,t]*generators[j].Pgmax - xp[j,t] - xr[j,t] , j=1:numgenerators; j!=i} >= p[i,t])
    @addConstraint(m, resbounds2[t=1:T], sum{x[i,t], i=1:numgenerators} >= 2)
    # Moved to sub-problem
    # @addConstraint(m, resbounds2[i=1:numgenerators, t=1:T], sum{rup[j,t], j=1:numgenerators} >= p[i,t])

    # Other parameters calculation
    Ω = Array(Any, T)
    for t=1:T
        Ω[t] = sum([ω[i,t] for i in 1:numfarms])
    end

    println(">>>> Variable creation time: $(toq()) seconds")
    @assert refbus == numbuses
    @assert !(refbus in generatorlist)

    tic()

    # δref and θref values to zero 
    @addConstraint(m, δref[t=1:T], delta[refbus,t] == 0)
    @addConstraint(m, θref[t=1:T], theta[refbus,t] == 0)

    # piecewise linear production cost
    @addConstraint(m, pseg1lb[i=1:numgenerators, t=1:T], pseg[i,t,1] >= 0)
    @addConstraint(m, pseg1ub[i=1:numgenerators, t=1:T], pseg[i,t,1] <= generators[i].b1*x[i,t])
    @addConstraint(m, pseg2lb[i=1:numgenerators, t=1:T], pseg[i,t,2] >= 0)
    @addConstraint(m, pseg2ub[i=1:numgenerators, t=1:T], pseg[i,t,2] <= generators[i].b1*x[i,t])
    @addConstraint(m, pseg3lb[i=1:numgenerators, t=1:T], pseg[i,t,3] >= 0)
    @addConstraint(m, pseg3ub[i=1:numgenerators, t=1:T], pseg[i,t,3] <= generators[i].b1*x[i,t])
    @addConstraint(m, totalpower[i=1:numgenerators, t=1:T], p[i,t] == sum(pseg[i,t,:]))

    # Binary variable logic 
    @addConstraint(m, binlogic1[i=1:numgenerators, t=1:T], y[i,t] - z[i,t] == x[i,t] - ((t == 1) ? onoffinit[i] : x[i,t-1]))
    @addConstraint(m, binlogic2[i=1:numgenerators, t=1:T], y[i,t] + z[i,t] <= 1)

    # Generator limits and chance constraints on generator limits - base case
    @addConstraint(m, genlimits1a[i=1:numgenerators, t=1:T], p[i,t] + r⁺[i,t] <= generators[i].Pgmax*x[i,t])
    # @addConstraint(m, genlimits1b[i=1:numgenerators, t=1:T], p[i,t] <= generators[i].Pgmax*x[i,t])
    @addConstraint(m, genlimits1c[i=1:numgenerators, t=1:T], p[i,t] - r⁻[i,t] >= generators[i].Pgmin*x[i,t])
    @addConstraint(m, reslimits1[i=1:numgenerators, t=1:T], r⁺[i,t] <= (generators[i].Pgmax)*x[i,t]/2)
    @addConstraint(m, reslimits2[i=1:numgenerators, t=1:T], r⁻[i,t] <= (generators[i].Pgmax)*x[i,t]/2)
    @addConstraint(m, ccgenlimits1[i=1:numgenerators, t=1:T], -r⁻[i,t] <= -Ω[t]*α[i,t] <= r⁺[i,t], with_probability=(1-ϵg), approx="2.0")

    println(">>> Added binary logic, gen limits, gen cc")

    # Participation factor settings 
    @addConstraint(m, αcommited[i=1:numgenerators, t=1:T], α[i,t] <= x[i,t])
    @addConstraint(m, sumα[t=1:T], sum(α[:,t]) == 1)

    # Aggregare the power produced and alpha values for the buses
    @addConstraint(m, pb[b=1:numbuses, t=1:T], pᵇ[b,t] == ((length(buses[b].genids) > 0) ? sum([p[k,t] for k in buses[b].genids]) : 0))
    @addConstraint(m, αb[b=1:numbuses, t=1:T], αᵇ[b,t] == ((length(buses[b].genids) > 0) ? sum([α[k,t] for k in buses[b].genids]) : 0))

    # Minimum up time and downtime constraints
    @addConstraint(m, minupdown1[i=1:numgenerators, t=1:(Lupmin[i] + Ldownmin[i])], x[i,t] == onoffinit[i])
    for i in 1:numgenerators
        for t in max(1, Lupmin[i]):T
            LB = max(t-generators[i].uptime+1, 1)
            @addConstraint(m, sum{y[i,r], r=LB:t} <= x[i,t])
        end
        for t in max(1, Ldownmin[i]):T
            LB = max(t-generators[i].downtime+1, 1)
            @addConstraint(m, sum{z[i,r], r=LB:t} <= (1-x[i,t]))
        end
    end

    println(">>> Added α settings, mismatch value setup, up time and down time constraints")
    # Ramping constraints
    @addConstraint(m, rampdown[i=1:numgenerators, t=1:T], p[i,t] - ((t==1) ? generators[i].outputatinit : p[i,t-1]) >= -generators[i].rampdown)
    @addConstraint(m, rampup[i=1:numgenerators, t=1:T], p[i,t] - ((t==1) ? generators[i].outputatinit : p[i,t-1]) <= generators[i].rampup)

    # Startup cost setup \cite{simouglou2010}
    @addConstraint(m, startupcost2[i=1:numgenerators, t=1:T], sum{suc[i,j,t], j=1:numblocks} == y[i,t])
    @addConstraint(m, startupcost3[i=1:numgenerators, t=1:T], SC[i,t] == sum{sucost[i,j]*suc[i,j,t], j=1:8})
    @addConstraint(m, startupcost1[i=1:numgenerators, j=1:numblocks, t=1:T], suc[i,j,t] <= sum{z[i,t-j], k=j:min(t-1,j)} +
    ((j<8 && j<=(generators[i].countinitoff+t-1) && (generators[i].countinitoff+t-1)<(j+1)) ? 1 : 0) + 
    ((j==8 && j<=(generators[i].countinitoff+t-1)) ? 1 : 0))


    # Power balance for base case and scenarios
    LoadMinusMeanWind = zeros(numbuses, T)
    for b in 1:numbuses
        LoadMinusMeanWind[b,:] = getLoadMinusMeanWind(buses[b])
    end
    @addConstraint(m, powerbalance[t=1:T], sum(p[:,t]) - pcurtail[t] == sum(LoadMinusMeanWind[:,t]))

    # B̂δ = α and B̂θ = p + μ - d constraints for base case and line scenarios
    Brow = Bhatsp'
    for b in 1:(numbuses-1)
        for t in 1:T
            @addConstraint(m, sum{Brow.nzval[idx]*delta[Brow.rowval[idx],t], idx in Brow.colptr[b]:(Brow.colptr[b+1]-1)} - (isgen(buses[b]) ? αᵇ[b,t] : 0 ) == 0)
            @addConstraint(m, sum{Brow.nzval[idx]*theta[Brow.rowval[idx],t], idx in Brow.colptr[b]:(Brow.colptr[b+1]-1)} - pᵇ[b,t] + LoadMinusMeanWind[b,t] == 0)
        end
    end

    # Ref bus alpha values
    @addConstraint(m, αrefbus[t=1:T], αᵇ[refbus,t] == 0)

    println(">>> Added ramping, startup cost, power balance and ref settings")
    # Line limit constraints for base case
    Blu = lufact(Bhatsp)
    A = zeros(numbuses, numfarms)
    rhs = zeros(numbuses-1)
    for f in 1:numfarms
        b = buses[farms[f].node]
        if (b.nodeID == refbus) 
            colbarinv = Array(Float64, numbuses)
            colbarinv = zeros(numbuses)
            A[:,f] = colbarinv
            continue
        end
        rhsindex = b.nodeID
        rhs[rhsindex] = 1
        colbarinv = Array(Float64, numbuses)
        colbarinv[1:numbuses .!= refbus] = Blu\rhs
        colbarinv[refbus] = 0
        rhs[rhsindex] = 0
        A[:,f] = colbarinv
    end
    fexpr = Array(Any, (numbuses,T))
    fexpr[numbuses,:] = zeros(1,T)
    for b in 1:numbuses-1
        for t in 1:T
            deterministic_part = theta[b,t]
            chance_part = @defExpr(-delta[b,t])*Ω[t] + @defExpr(sum{A[b,f]*ω[f,t], f=1:numfarms})
            fexpr[b,t] = deterministic_part + chance_part
        end
    end

    @addConstraint(m, cclinelimits1[i=1:numlines, t=1:T], -getThermalCapacity(lines[i], mvaBase) <= 
    lines[i].y*(fexpr[lines[i].head,t] - fexpr[lines[i].tail,t]) <= getThermalCapacity(lines[i], mvaBase), with_probability=(1-ϵij), approx="2.0")

    println(">>> Added base case line limit cc")

    #
    busTofarmmap = Dict() # Make sure each bus has only one farm - if not, aggregate them
    for b in buses
        busTofarmmap[b.nodeID] = ((length(b.Farmids) == 1) ? b.Farmids[1] : 0)
    end

    if length(violations) != 0
        for s in 1:numlinescenarios
            fexprv = Array(Any, (numbuses, T))
            Btildeˢ = scenarios[s].Btilde
            Bhatsp = scenarios[s].Bhat
            Blu = lufact(Bhatsp)
            A = zeros(numbuses, numfarms)
            rhs = zeros(numbuses-1)
            for f in 1:numfarms
                b = buses[farms[f].node]
                if (b.nodeID == refbus) 
                    colbarinv = Array(Float64, numbuses)
                    colbarinv = zeros(numbuses)
                    A[:,f] = colbarinv
                    continue
                end
                rhsindex = b.nodeID
                rhs[rhsindex] = 1
                colbarinv = Array(Float64, numbuses)
                colbarinv[1:numbuses .!= refbus] = Blu\rhs
                colbarinv[refbus] = 0
                rhs[rhsindex] = 0
                A[:,f] = colbarinv
            end

            for t in 1:T
                for b in 1:numbuses-1
                    deterministic_part = @defExpr(sum{Btildeˢ[b,k]*(pᵇ[k,t] - LoadMinusMeanWind[k,t]), k=1:numbuses})
                    chance_part = @defExpr(-sum{Btildeˢ[b,k]*αᵇ[k,t], k=1:numbuses})*Ω[t] + @defExpr(sum{A[b,f]*ω[f,t], f=1:numfarms})
                    fexprv[b,t] = deterministic_part + chance_part
                end
                fexprv[numbuses,t] = 0
            end

            for l in violations
                for t in 1:T
                    if (scenarios[s].kind == 1 && lines[l].arcID == scenarios[s].ID)
                        continue
                    end
                    @addConstraint(m, -getThermalCapacity(lines[l], mvaBase) <= lines[l].y*(fexprv[lines[l].head,t] - fexprv[lines[l].tail,t]) <= 
                    getThermalCapacity(lines[l], mvaBase), with_probability=(1-2*ϵij), approx="2.0")
                end
                println(">>> Creating line $l cc for scenario : $s")
            end    
        end
    end	
    
    for s=1:numlinescenarios
        fexpr = Array(Any, (numbuses,T))
        fexpr[numbuses,:] = zeros(1,T)
        Btildeˢ = scenarios[s].Btilde
        for b in 1:numbuses-1
            for t in 1:T
                deterministic_part = @defExpr(sum{Btildeˢ[b,k]*(pᵇ[k,t] - LoadMinusMeanWind[k,t]), k=1:numbuses})
                fexpr[b,t] = deterministic_part
            end
        end

        for i in 1:numlines
            if (scenarios[s].kind == 1 && lines[i].arcID == scenarios[s].ID)
                continue
            end
            for t=1:T
                @addConstraint(m, -getThermalCapacity(lines[i], mvaBase) <= lines[i].y*(fexpr[lines[i].head,t] - fexpr[lines[i].tail,t]) <= getThermalCapacity(lines[i], mvaBase))
            end
        end
        println(">>> Creating hard limit for line : $s")
    end
    #
    println(">>>> Constraint add time: $(toq()) seconds")
    # Objective 
    sumvar = Array(Float64, T)
    for t=1:T
        sumvar[t] = sum([(farms[i].var)[t] for i in 1:numfarms])
    end
    #
    @setObjective(m, Min, sum{sum{generators[i].noloadcost*x[i,t] + sum{generators[i].prodslope[k]*pseg[i,t,k], k=1:3} + 
    generators[i].pi2*2*(r⁺[i,t]+r⁻[i,t]), i=1:numgenerators}, t=1:T} + 100000*sum{pcurtail[t], t=1:T} + sum{surrogate[t], t=1:T})

    function lazyconstraintgenerator(cb)
        pvals = getValue(p)
        xpvals = getValue(xp)
        xrvals⁺ = getValue(xr)
        pvalsᵇ = getValue(pᵇ)
        αvals = getValue(α)
        αvalsᵇ = getValue(αᵇ)
        xvals = getValue(x)
        rvals⁺ = getValue(r⁺)
        rvals⁻ = getValue(r⁻)

        infcount = 0
        feasibilityCuts = Any[]
        optimalityCuts = Any[]
        UBObjective = zeros(T)
        for t in 1:T
            msub = Model(solver = CplexSolver(CPX_PARAM_PREIND=0))

            # Subproblem Variables
            @defVar(msub, rup[1:numgenerators] >= 0)
            @defVar(msub, δ[1:numgenerators, 1:numgenscenarios])
            @defVar(msub, θ[1:numbuses, 1:numgenscenarios])
            @defVar(msub, δᵇ[1:numbuses, 1:numgenscenarios])

            @addConstraint(msub, ub1sub[i=1:numgenerators], rup[i] <= generators[i].Pgmax*xvals[i,t] - xpvals[i,t] - xrvals⁺[i,t])
            @defExpr(ub1subexpr[i=1:numgenerators], generators[i].Pgmax*x[i,t]-xp[i,t]-xr[i,t])

            @addConstraint(msub, ub2sub[i=1:numgenerators], rup[i] <= (generators[i].Pgmax)*xvals[i,t]/2)
            @defExpr(ub2subexpr[i=1:numgenerators], generators[i].Pgmax*x[i,t]/2)

            @addConstraint(msub, lbsub[i=1:numgenerators, k=1:numgenscenarios], -rup[i] <= -δ[i,k])

            @addConstraint(msub, resboundsub[i=1:numgenerators], sum{-rup[j], j=1:numgenerators; j!=i} <= -pvals[i,t])
            @defExpr(resboundsubexpr[i=1:numgenerators], -p[i,t])

            @addConstraint(msub, δlbsub[i=1:numgenerators, k=1:numgenscenarios], -δ[i,k] <= generators[i].Pgmax*xvals[scenarios[k+numlinescenarios].ID,t])
            @defExpr(δlbsubexpr[i=1:numgenerators, k=1:numgenscenarios], generators[i].Pgmax*x[scenarios[k+numlinescenarios].ID,t])

            @addConstraint(msub, δubsub[i=1:numgenerators, k=1:numgenscenarios], δ[i,k] <= generators[i].Pgmax*xvals[scenarios[k+numlinescenarios].ID,t])
            @defExpr(δubsubexpr[i=1:numgenerators, k=1:numgenscenarios], generators[i].Pgmax*x[scenarios[k+numlinescenarios].ID,t])

            @addConstraint(msub, δbussub[b=1:numbuses, k=1:numgenscenarios], δᵇ[b,k] == ((length(buses[b].genids) > 0) ? sum([δ[i,k] for i in buses[b].genids]) : 0))

            @addConstraint(msub, δpfsub[k=1:numgenscenarios], sum(δ[:,k]) == 0)

            @addConstraint(msub, δcsub[k=1:numgenscenarios], δ[scenarios[k+numlinescenarios].ID,k] == -pvals[scenarios[k+numlinescenarios].ID,t])
            @defExpr(δcsubexpr[k=1:numgenscenarios], -p[scenarios[k+numlinescenarios].ID,t])

            @addConstraint(msub, θrefbussub[k=1:numgenscenarios], θ[refbus,k] == 0)      

            committedgenerators = 0

            for k in 1:numgenscenarios
                if (xvals[scenarios[k+numlinescenarios].ID,t] >= 0.9)
                    committedgenerators += 1
                end
            end
            @defConstrRef pfll[1:committedgenerators,1:numbuses-1+numlines*2]

            generatorcount = 0
            for k in 1:numgenscenarios
                if (xvals[scenarios[k+numlinescenarios].ID,t] >= 0.9)
                    generatorcount += 1
                    for b in 1:(numbuses-1)
                        pfll[generatorcount, b] = @addConstraint(msub, sum{Brow.nzval[idx]*θ[Brow.rowval[idx],k], idx in Brow.colptr[b]:(Brow.colptr[b+1]-1)} 
                        - δᵇ[b,k] == pvalsᵇ[b,t] - LoadMinusMeanWind[b,t])
                    end

                    for i in 1:numlines
                        pfll[generatorcount, i+numbuses-1] = @addConstraint(msub, lines[i].y*(θ[lines[i].head,k] - θ[lines[i].tail,k]) <= getThermalCapacity(lines[i], mvaBase))
                        pfll[generatorcount, i+numbuses-1+numlines] = @addConstraint(msub, -lines[i].y*(θ[lines[i].head,k] - θ[lines[i].tail,k]) <= getThermalCapacity(lines[i], mvaBase))
                    end
                else
                    @addConstraint(msub, θ[:,k] .== 0)
                end
            end

            @setObjective(msub, Min, sum{generators[i].pi2*2*rup[i], i=1:numgenerators})
            substatus = solve(msub)
            UBObjective[t] = getObjectiveValue(msub)

            # Bender's cut
            cut =  dot(getDual(ub1sub), ub1subexpr)
            append!(cut,dot(getDual(ub2sub), ub2subexpr))
            append!(cut,dot(getDual(resboundsub), resboundsubexpr))
            append!(cut,sum(getDual(δlbsub).*δlbsubexpr))
            append!(cut,sum(getDual(δubsub).*δubsubexpr))
            append!(cut,dot(getDual(δcsub), δcsubexpr))
            generatorcount = 0
            for k in 1:numgenscenarios
                if (xvals[scenarios[k+numlinescenarios].ID,t] >= 0.9)
                    generatorcount += 1
                    for b in 1:(numbuses-1)
                        append!(cut, getDual(pfll[generatorcount, b])*(pᵇ[b,t]-LoadMinusMeanWind[b,t]))
                    end

                    for i in 1:numlines
                        append!(cut, getDual(pfll[generatorcount, i+numbuses-1])*getThermalCapacity(lines[i], mvaBase))
                        append!(cut, getDual(pfll[generatorcount, i+numbuses-1+numlines])*getThermalCapacity(lines[i], mvaBase))
                    end
                end
            end

            if substatus != :Optimal
                push!(feasibilityCuts, @LinearConstraint(cut <= 0))
                infcount += 1
            end

            if substatus == :Optimal
                push!(optimalityCuts, @LinearConstraint(surrogate[t] >= cut))
            end

        end
        #=       
        if infcount == 0
            upperbound = MathProgBase.cbgetbestbound(cb) - sum(getValue(surrogate)) + sum(UBObjective)
            lowerbound = MathProgBase.cbgetbestbound(cb)

            println(">>> Upper bound = $upperbound")
            println(">>> Lower bound = $lowerbound")
            if (abs(upperbound-lowerbound) < 1e-6 && lowerbound > 0 && upperbound < 1e20 )
                println(" ------------------------- Returning ---------------------- ")
                return
            end

            if (lowerbound < upperbound)
                for cut in optimalityCuts
                    addLazyConstraint(cb, cut)
                end
            end
            =#
            #else
            for cut in optimalityCuts
                addLazyConstraint(cb, cut)
            end

            for cut in feasibilityCuts
                addLazyConstraint(cb, cut)
            end
            #end


            println(">>> Infeasibilities: $(infcount), Bender's Feasibility Cuts added: $infcount")
            println(">>> Bender's Optimality Cuts added: $(24-infcount)")

        end

        addLazyCallback(m, lazyconstraintgenerator)

        status = solve(m, method=:Cuts, debug=false, lazy_constraints=true)

        if status != :Optimal
            println(status)
            quit()
        end
		#=
        extras["output_file"] = chomp(extras["output_file"])
        jldopen("case96.jld", "w") do file
            write(file, "x", getValue(x))
            write(file, "y", getValue(y))
            write(file, "z", getValue(z))
            write(file, "pseg", getValue(pseg))
            write(file, "α", getValue(α))
            write(file, "p", getValue(p))
            write(file, "r⁺", getValue(r⁺))
            write(file, "r⁻", getValue(r⁻))
            write(file, "surrogate", getValue(surrogate))
            write(file, "pcurtail", getValue(pcurtail))
            write(file, "obj", getObjectiveValue(m))
        end
		=#

        Violations = Set{Int}()
        αvals = getValue(αᵇ)

        for s in 1:numlinescenarios
            fbus = zeros(numbuses, T)
            Btildeˢ = scenarios[s].Btilde
            pvals = getValue(pᵇ)
            fbus = Btildeˢ*(pvals - LoadMinusMeanWind)
            δvals = Btildeˢ*(αvals)

            Blu = lufact(scenarios[s].Bhat)
            A = zeros(numbuses, numfarms)
            rhs = zeros(numbuses-1)
            for f in 1:numfarms
                b = buses[farms[f].node]
                if (b.nodeID == refbus) 
                    colbarinv = Array(Float64, numbuses)
                    colbarinv = zeros(numbuses)
                    A[:,f] = colbarinv
                    continue
                end
                rhsindex = b.nodeID
                rhs[rhsindex] = 1
                colbarinv = Array(Float64, numbuses)
                colbarinv[1:numbuses .!= refbus] = Blu\rhs
                colbarinv[refbus] = 0
                rhs[rhsindex] = 0
                A[:,f] = colbarinv
            end

            with_prob = (1-ϵij)
            nu = quantile(Normal(0,1), 1-2*ϵij)
            for i in 1:numlines
                if (scenarios[s].kind == 1 && lines[i].arcID == scenarios[s].ID)
                    continue
                end
                for t=1:T
                    mean = lines[i].y*(fbus[lines[i].head, t] - fbus[lines[i].tail, t])
                    var = 0
                    for f in 1:numfarms
                        var += farms[f].var[t]*(A[lines[i].head, f] - A[lines[i].tail, f] - δvals[lines[i].head, t] + δvals[lines[i].tail, t])^2
                    end
                    stddevline = lines[i].y*sqrt(var)
                    #println("$(mean + nu*stddevline) > $(getThermalCapacity(lines[i], mvaBase) + 1e-5) || $(mean - nu*stddevline) < $(-getThermalCapacity(lines[i], mvaBase) - 1e-5)")
                    if (mean + nu*stddevline > getThermalCapacity(lines[i], mvaBase) + 1e-3 || mean - nu*stddevline < -getThermalCapacity(lines[i], mvaBase) - 1e-3)
                        println("($i, $t, $s)")
                        println("$(mean + nu*stddevline) > $(getThermalCapacity(lines[i], mvaBase) + 1e-5) || $(mean - nu*stddevline) < $(-getThermalCapacity(lines[i], mvaBase) - 1e-5)")
                        #push!(violations, (i,t,s))
                        push!(Violations, i)
                    end
                end
            end

        end
        println(">>> No of violations = $(length(Violations)) ")

        return Violations

    end

    #=
    status = solve(m, method=:Cuts, debug=false, lazy_constraints=true)
    if status != :Optimal
        println(status)
        quit()
    end

    extras["output_file"] = chomp(extras["output_file"])
    jldopen("case96.jld", "w") do file
        write(file, "x", getValue(x))
        write(file, "y", getValue(y))
        write(file, "z", getValue(z))
        write(file, "pseg", getValue(pseg))
        write(file, "α", getValue(α))
        write(file, "p", getValue(p))
        write(file, "r⁺", getValue(r⁺))
        write(file, "r⁻", getValue(r⁻))
        write(file, "surrogate", getValue(surrogate))
        write(file, "pcurtail", getValue(pcurtail))
        write(file, "obj", getObjectiveValue(m))
    end

    quit()
    #
    #    violations = Tuple{Int64, Int64, Int64}[] # Line, time, scenario
    violations = Set{Int}()
    αvals = getValue(αᵇ)

    for s in 1:numlinescenarios
        fbus = zeros(numbuses, T)
        Btildeˢ = scenarios[s].Btilde
        pvals = getValue(pᵇ)
        fbus = Btildeˢ*(pvals - LoadMinusMeanWind)
        δvals = Btildeˢ*(αvals)

        Blu = lufact(scenarios[s].Bhat)
        A = zeros(numbuses, numfarms)
        rhs = zeros(numbuses-1)
        for f in 1:numfarms
            b = buses[farms[f].node]
            if (b.nodeID == refbus) 
                colbarinv = Array(Float64, numbuses)
                colbarinv = zeros(numbuses)
                A[:,f] = colbarinv
                continue
            end
            rhsindex = b.nodeID
            rhs[rhsindex] = 1
            colbarinv = Array(Float64, numbuses)
            colbarinv[1:numbuses .!= refbus] = Blu\rhs
            colbarinv[refbus] = 0
            rhs[rhsindex] = 0
            A[:,f] = colbarinv
        end

        with_prob = (1-ϵij)
        nu = quantile(Normal(0,1), 1-2*ϵij)
        for i in 1:numlines
            if (scenarios[s].kind == 1 && lines[i].arcID == scenarios[s].ID)
                continue
            end
            for t=1:T
                mean = lines[i].y*(fbus[lines[i].head, t] - fbus[lines[i].tail, t])
                var = 0
                for f in 1:numfarms
                    var += farms[f].var[t]*(A[lines[i].head, f] - A[lines[i].tail, f] - δvals[lines[i].head, t] + δvals[lines[i].tail, t])^2
                end
                stddevline = lines[i].y*sqrt(var)
                #println("$(mean + nu*stddevline) > $(getThermalCapacity(lines[i], mvaBase) + 1e-5) || $(mean - nu*stddevline) < $(-getThermalCapacity(lines[i], mvaBase) - 1e-5)")
                if (mean + nu*stddevline > getThermalCapacity(lines[i], mvaBase) + 1e-5 || mean - nu*stddevline < -getThermalCapacity(lines[i], mvaBase) - 1e-5)
                    println("($i, $t, $s)")
                    println("$(mean + nu*stddevline) > $(getThermalCapacity(lines[i], mvaBase) + 1e-5) || $(mean - nu*stddevline) < $(-getThermalCapacity(lines[i], mvaBase) - 1e-5)")
                    #push!(violations, (i,t,s))
                    push!(violations, i)
                end
            end
        end

    end
    println(">>> No of violations = $(length(violations)) ")

    while (length(violations) != 0)

        while (length(violations) != 0)

            for s in 1:numlinescenarios
                # l, t, s = pop!(violations)
                fexprv = Array(Any, (numbuses, T))
                Btildeˢ = scenarios[s].Btilde
                Bhatsp = scenarios[s].Bhat
                Blu = lufact(Bhatsp)
                A = zeros(numbuses, numfarms)
                rhs = zeros(numbuses-1)
                for f in 1:numfarms
                    b = buses[farms[f].node]
                    if (b.nodeID == refbus) 
                        colbarinv = Array(Float64, numbuses)
                        colbarinv = zeros(numbuses)
                        A[:,f] = colbarinv
                        continue
                    end
                    rhsindex = b.nodeID
                    rhs[rhsindex] = 1
                    colbarinv = Array(Float64, numbuses)
                    colbarinv[1:numbuses .!= refbus] = Blu\rhs
                    colbarinv[refbus] = 0
                    rhs[rhsindex] = 0
                    A[:,f] = colbarinv
                end

                for t in 1:T
                    for b in 1:numbuses-1
                        deterministic_part = @defExpr(sum{Btildeˢ[b,k]*(pᵇ[k,t] - LoadMinusMeanWind[k,t]), k=1:numbuses})
                        chance_part = @defExpr(-sum{Btildeˢ[b,k]*αᵇ[k,t], k=1:numbuses})*Ω[t] + @defExpr(sum{A[b,f]*ω[f,t], f=1:numfarms})
                        fexprv[b,t] = deterministic_part + chance_part
                    end
                    fexprv[numbuses,t] = 0
                end

                for l in violations
                    for t in 1:T
                        if (scenarios[s].kind == 1 && lines[l].arcID == scenarios[s].ID)
                            continue
                        end
                        @addConstraint(m, -getThermalCapacity(lines[l], mvaBase) <= lines[l].y*(fexprv[lines[l].head,t] - fexprv[lines[l].tail,t]) <= 
                        getThermalCapacity(lines[l], mvaBase), with_probability=(1-2*ϵij), approx="2.0")
                    end
                end    
            end # End of for loop
                empty!(violations)
        end # end of inner while loop


            println("------------------------------- Resolving --------------------------")
            status = solve(m, method=:Cuts, debug=false, lazy_constraints=true)
            if status != :Optimal
                println(status)
                quit()
            end

            αvals = getValue(αᵇ)

            for s in 1:numlinescenarios
                fbus = zeros(numbuses, T)
                Btildeˢ = scenarios[s].Btilde
                pvals = getValue(pᵇ)
                fbus = Btildeˢ*(pvals - LoadMinusMeanWind)
                δvals = Btildeˢ*(αvals)

                Blu = lufact(scenarios[s].Bhat)
                A = zeros(numbuses, numfarms)
                rhs = zeros(numbuses-1)
                for f in 1:numfarms
                    b = buses[farms[f].node]
                    if (b.nodeID == refbus) 
                        colbarinv = Array(Float64, numbuses)
                        colbarinv = zeros(numbuses)
                        A[:,f] = colbarinv
                        continue
                    end
                    rhsindex = b.nodeID
                    rhs[rhsindex] = 1
                    colbarinv = Array(Float64, numbuses)
                    colbarinv[1:numbuses .!= refbus] = Blu\rhs
                    colbarinv[refbus] = 0
                    rhs[rhsindex] = 0
                    A[:,f] = colbarinv
                end

                with_prob = (1-ϵij)
                nu = quantile(Normal(0,1), 1-2*ϵij)
                for i in 1:numlines
                    if (scenarios[s].kind == 1 && lines[i].arcID == scenarios[s].ID)
                        continue
                    end
                    for t=1:T
                        mean = lines[i].y*(fbus[lines[i].head, t] - fbus[lines[i].tail, t])
                        var = 0
                        for f in 1:numfarms
                            var += farms[f].var[t]*(A[lines[i].head, f] - A[lines[i].tail, f] - δvals[lines[i].head, t] + δvals[lines[i].tail, t])^2
                        end
                        stddevline = lines[i].y*sqrt(var)
                        #println("$(mean + nu*stddevline) > $(getThermalCapacity(lines[i], mvaBase) + 1e-5) || $(mean - nu*stddevline) < $(-getThermalCapacity(lines[i], mvaBase) - 1e-5)")
                        if (mean + nu*stddevline > getThermalCapacity(lines[i], mvaBase) + 1e-3 || mean - nu*stddevline < -getThermalCapacity(lines[i], mvaBase) - 1e-3)
                            println("$(mean + nu*stddevline) > $(getThermalCapacity(lines[i], mvaBase) + 1e-5) || $(mean - nu*stddevline) < $(-getThermalCapacity(lines[i], mvaBase) - 1e-5)")
                            println("($i, $t, $s)")
                            push!(violations, i)
                        end
                    end
                end
            end

    end # end of outer while loop

    end
    =#
