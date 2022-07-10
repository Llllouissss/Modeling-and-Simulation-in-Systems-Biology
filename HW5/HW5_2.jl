using Agents, Random

mutable struct Agent <: AbstractAgent
    id::Int                 # Mandatory Agent identifier
    pos::NTuple{2,Float64}  # Position, required for agents in the ContinuousSpace
    vel::NTuple{2,Float64}  # Moving speeds
    mass::Float64           # Can move or not
end

function ball_model(; speed = 0.002)
    space2d = ContinuousSpace((1, 1), 0.02)
    model = ABM(Agent, space2d, properties = Dict(:dt => 1.0), rng = MersenneTwister(42))

    # Add agents to the model
    for ind in 1:500
        pos = Tuple(rand(model.rng, 2))
        vel = sincos(2π * rand(model.rng)) .* speed
        mass = 1.0
        add_agent!(pos, model, vel, mass)
    end
    return model
end

model = ball_model()

# Agents.move_agent!()
agent_step!(agent, model) = move_agent!(agent, model, model.dt)

function model_step!(model)
    for (a1, a2) in interacting_pairs(model, 0.012, :nearest)
        elastic_collision!(a1, a2, :mass)
    end
end

mutable struct PoorSoul <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    mass::Float64
    days_infected::Int  # number of days since is infected
    status::Symbol  # :S, :I or :R
    β::Float64
end

const steps_per_day = 24 # One tick per hour

function sir_initiation(;
    infection_period = 30 * steps_per_day,
    detection_time = 14 * steps_per_day,
    reinfection_probability = 0.02,
    isolated = 0.0, # in percentage
    interaction_radius = 0.010,
    dt = 1.0,
    speed = 0.001,
    death_rate = 0.044,
    N = 1000,
    initial_infected = 2,
    seed = 42,
    βmin = 0.01,
    βmax = 0.04,
)

    properties = (;
        infection_period,
        reinfection_probability,
        detection_time,
        death_rate,
        interaction_radius,
        dt,
    )
    space = ContinuousSpace((1,1), 0.02)
    model = ABM(PoorSoul, space, properties = Dict(pairs(properties)), rng = MersenneTwister(seed))

    # Add initial individual agents
    for ind in 1:N
        pos = Tuple(rand(model.rng, 2))
        status = ind ≤ N - initial_infected ? :S : :I
        isisolated = ind ≤ isolated * N
        mass = isisolated ? Inf : 1.0
        vel = isisolated ? (0.0, 0.0) : sincos(2π * rand(model.rng)) .* speed

        β = (βmax - βmin) * rand(model.rng) + βmin
        add_agent!(pos, model, vel, mass, 0, status, β)
    end

    return model
end

sir_model = sir_initiation()

sir_colors(a) = a.status == :S ? "#2b2b33" : a.status == :I ? "#bf2642" : "#338c54"


function transmit!(a1, a2, reinfection_probability)

    infected, other = a1.status == :I ? (a1, a2) : (a2, a1)

    # No one infected
    if infected.status != :I
        return
    end

    # if the other agent is already infected, do nothing
    if other.status == :I
        return
    end

    # Lucky and not infected
    if rand(model.rng) <= infected.β
        return
    end

    # Risk of reinfection
    if other.status == :R && rand(model.rng) > reinfection_probability
        return
    end

    # You got virus
    other.status = :I
end

function sir_model_step!(model)
    r = model.interaction_radius
    for (a1, a2) in interacting_pairs(model, r, :nearest)
        transmit!(a1, a2, model.reinfection_probability)
        elastic_collision!(a1, a2, :mass)
    end
end

# Agent-specific functions
function update!(agent) 
    if agent.status == :I
        agent.days_infected += 1
    end
end

function recover_or_die!(agent, model)
    if agent.days_infected ≥ model.infection_period
        if rand(model.rng) ≤ model.death_rate
            kill_agent!(agent, model)
        else
            agent.status = :R
            agent.days_infected = 0
        end
    end
end

function sir_agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    update!(agent)
    recover_or_die!(agent, model)
end

infected(x) = count(i == :I for i in x)
recovered(x) = count(i == :R for i in x)
# Aggregated data for number of infected and recovered indivisuals
adata = [(:status, infected), (:status, recovered)]

# Try different parameters
β1, β2, β3,β4,β5 = 0.01, 0.02,0.03,0.04,0.05
sir_model1 = sir_initiation(βmax = β1,isolated = 0)
sir_model2 = sir_initiation(isolated = 0, βmax = β2)
sir_model3 = sir_initiation(isolated = 0, βmax = β3)
sir_model4 = sir_initiation(isolated = 0, βmax = β4)
sir_model5 = sir_initiation(isolated = 0, βmax = β5)

data1, _ = run!(sir_model1, sir_agent_step!, sir_model_step!, 5000; adata)
data2, _ = run!(sir_model2, sir_agent_step!, sir_model_step!, 5000; adata)
data3, _ = run!(sir_model3, sir_agent_step!, sir_model_step!, 5000; adata)
data4, _ = run!(sir_model4, sir_agent_step!, sir_model_step!, 5000; adata)
data5, _ = run!(sir_model5, sir_agent_step!, sir_model_step!, 5000; adata)

data1[(end-10):end, :]

using CairoMakie

figure = Figure()
ax = figure[1, 1] = Axis(figure; ylabel = "Infected", xlabel="Steps")
l1 = lines!(ax, data1[:, dataname((:status, infected))], color = :orange)
l2 = lines!(ax, data2[:, dataname((:status, infected))], color = :blue)
l3 = lines!(ax, data3[:, dataname((:status, infected))], color = :green)
l4 = lines!(ax, data4[:, dataname((:status, infected))], color = :black)
l5 = lines!(ax, data5[:, dataname((:status, infected))], color = :red)
figure[1, 2] = Legend(figure, [l1, l2, l3,l4,l5], ["isolated = 0, beta=$β1", "isolated = 0, beta=$β2", "isolated = 0, beta=$β3", "isolated = 0, beta=$β4","isolated = 0, beta=$β5"])

figure






