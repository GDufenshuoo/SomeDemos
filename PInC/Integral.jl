
include("Depend.jl")

include("EI.jl")
include("Unit.jl")
include("Model.jl")

import ForwardDiff

Problem = K_Modelly(1,1000000.0,Atomicᵁ)

dim = 3*Problem.N
T = as(Array, dim);
P = TransformedLogDensity(T, Problem);  
∇P = ADgradient(:ForwardDiff, P);

using Pathfinder

result_pf = pathfinder(∇P)

init_params = result_pf.draws[:, 1]
inv_metric = result_pf.fit_distribution.Σ

result_dhmc2 = mcmc_with_warmup(
    Random.GLOBAL_RNG,
    ∇P,
    10000;
    initialization=(; q=init_params, κ=GaussianKineticEnergy(inv_metric)),
    warmup_stages=default_warmup_stages(),
    reporter=NoProgressReport(),
)

chains = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, 10)

chains = mcmc_with_warmup(
    Random.GLOBAL_RNG, ∇P, 100; 
    initialization = (ϵ = 1e-2,q=rand(3) ),  
    warmup_stages = default_warmup_stages(; stepsize_search = nothing))

function Rdistribution(Problem,P)
    @unpack N,B = Problem
    t = size(P,2)
    box = Vector{Vector{Float64}}(undef,N)
    for i in eachindex(box)
        box[i] = []
    end
    for i in 1:t
        Rdistribution!(box,reshape(P[:,i],3,B,N),N,B)
    end
    return box
end

function Rdistribution!(box,P,N,B)
    for i in 1:N, b in 1:B
        push!(box[i],√(sum(abs2,P[:,b,i])))
    end
end

Rd = Rdistribution(Problem,chains.posterior_matrix)

using Plots
histogram(Rd)