using MuseInference
using LinearAlgebra
using Zygote
using ForwardDiff

# 512-dimensional noisy funnel
prob = SimpleMuseProblem(
    rand(10),
    function sample_x_z(rng, θ)
        z = rand(rng, MvNormal(zeros(10), exp(θ) * I))
        x = rand(rng, MvNormal(z, I))
        (;x, z)
    end,
    function logLike(x, z, θ)
        -(1//2) * (sum((x .- z).^2) + sum(z.^2) / exp(θ) + 512*θ)
    end, 
    function logPrior(θ)
        -θ^2/(2*3^2)
    end;
)

# get solution
muse(prob, (θ=1,))

