using StatsBase # Weights() and sample()
using Random    # randexp()
using Plots
using Interpolations
using Statistics # mean()
Random.seed!(2022)

#=
Stochastic chemical reaction: Gillespie Algorithm (direct method)
Adapted from: Chemical and Biomedical Enginnering Calculations Using Python Ch.4-3
=#
function ssa_direct(model, u0::AbstractArray, tend, p, stoich; tstart=zero(tend))
    t = tstart   # Current time
    ts = [t]     # Time points
    u = copy(u0) # Current state
    us = copy(u) # States over time
    step = 0
    while step<20000 # t < tend
        step += 1
        a = model(u, p, t)               # propensities
        dt = randexp() / sum(a)          # Time step for the direct method
        du = sample(stoich, Weights(a))  # Choose the stoichiometry for the next reaction
        u .+= du  # Update state
        t += dt   # Update time
        
        us = [us u]  # Append state
        push!(ts, t) # Append time point
    end
    # Trasnpose to make columns as variables, rows as observations
    us = collect(us')
    return (t = ts, u = us)
end

#=
Stochastic chemical reaction: Gillespie Algorithm (first reaction method)
Adapted from: Chemical and Biomedical Enginnering Calculations Using Python Ch.4-3
=#
function ssa_first(model, u0, tend, p, stoich; tstart=zero(tend))
    t = tstart   # Current time
    ts = [t]     # Time points
    u = copy(u0) # Current state
    us = copy(u) # States over time
    step = 0
    while step<20000 # t < tend
        step += 1
        a = model(u, p, t)              # propensities
        dts = randexp(length(a)) ./ a   # dts from all reactions
        # Choose the reaction 
        i = argmin(dts)
        dt = dts[i]
        du = stoich[i]
        # Update state and time
        u .+= du
        t += dt
        us = [us u]  # Append state variable to record
        push!(ts, t) # Append time point to record
    end
    # Make column as variables, rows as observations
    us = collect(us')
    return (t = ts, u = us)
end

"""
Propensity model for this reaction.
Reaction of A <-> B with rate constants k1 & k2
"""
model(u, p, t) = [p.alpha/(1 + u[2]^p.beta),  p.alpha/(1+ u[1]^p.beta), p.delta*u[1], p.delta*u[2]]
parameters = (alpha=5.0, beta=4.0, delta=1.0, stoich=[[1, 0], [0, 1], [-1, 0], [0, -1]])
u0 = [75,75]
tend = 100.0



soldirect = ssa_direct(model, u0, tend, parameters, parameters.stoich)
solfirst = ssa_first(model, u0, tend, parameters, parameters.stoich)

plot(soldirect.t, soldirect.u, 
    xlabel="time", ylabel="# of molecules", 
    title = "SSA (direct method)", label=["A" "B"])

plot(solfirst.t, solfirst.u, 
xlabel="time", ylabel="# of molecules", 
title = "SSA (1st reaction method)", label=["A" "B"])
















