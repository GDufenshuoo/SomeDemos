#The simple test
include("Fermions.jl")

using Distributions
using Random
using LinearAlgebra
using Plots
using ProgressBars
using Printf
using StatsFuns

begin #"SI2019"
	ħ = 1.0545718176461565e-34
	E_hatree = 4.359744722207211e-18
	E_epr = 8.238723498265451e-8
	k_e = 4.359744722207211e-18
	k_B = 1.380649e-23
	𝘦 = 1.602176634e-19
	μ_B = 9.274010078302855e-24
	m_e = 9.109383701558256e-31
	m_p = 1.6726219236940502e-27
	m_u = 1.6605390666030467e-27
	a_0 = 5.291772109022829e-11

    m_he⁴ = 4.00260325413*m_u
end

# make it reproducible
Random.seed!(1)

# define distribution
N=2 #电子数量
B=256   #beads数量
ndim = N*B*3    #总的维度，不需要调整
T=10    #温度
beta=1/(k_B*T)  #β
k_05mw2=0.5*m_e*B/sqrt(beta*ħ)^2    #路径积分谐振子势能
# P=rand(N,3,B)

# W_f = W_bf(P)

function V_r(P) #势能 只有一个正电中心
    V = 0.0
    for b in 1:B
        for i in 1:N    #+
            V += -3/sqrt(sum(abs2,P[i,:,b])) *E_epr
        end
        for i in 2:N    #-
            for j in 1:i-1
                V += 1/sqrt(sum(abs2,P[i,:,b]-P[j,:,b])) *E_epr
            end
        end
    end
    V
end

function lp(p)  #总的势能
    nc = size(p,2)
    lp = zeros(nc)
    for i in 1:nc
        P = reshape(p[:,i],N,3,B)
        W_f = W_bf(P)   #交换虚拟势能
        V = V_r(P)/B
        lp[i] = -1 *((W_f <0 ? -10e2*log(-W_f) : -log(W_f)) + V)
    end
    return lp
end

function lp_debug(p)
    nc = size(p,2)
    lp = zeros(nc,2)
    for i in 1:nc
        P = reshape(p[:,i],N,3,B)
        W_f = W_bf(P)
        V = V_r(P)/B
        lp[i,1] = -1 *((W_f <0 ? -10e2*log(-W_f) : -log(W_f)))
        println(W_f," ",lp[i,1]," ",V)
    end
end

“”“
一个实验采样器
”“”
function RunDIME(init::Array, niter::Int; sigma::Float64=1e-5, gamma=nothing, aimh_prob::Float64=0.1, df_proposal_dist::Int=10, progress::Bool=true)

    ndim,nchain = size(init)

    # get some default values
    dft = df_proposal_dist

    if gamma == nothing 
        g0 = 2.38 / sqrt(2 * ndim)
    else
        g0 = gamma
    end

    # fix that MvTDist does not accept positive demi-definite covariance matrices
    fixPSD = Matrix(1e-7I, ndim, ndim)

    # initialize
    ccov = Matrix(1.0I, ndim, ndim)
    cmean = zeros(ndim)
    dist = MvTDist(dft, cmean, ccov + fixPSD)
    accepted = ones(nchain)
    cumlweight = -Inf

    # calculate intial values
    x = copy(init)
    lprob = lp(x)

    # preallocate
    lprobs = Array{Float64,2}(undef, niter, nchain)
    lprobs = fill!(lprobs, 0.0)

    chains = Array{Float64,3}(undef, niter, nchain, ndim)
    chains = fill!(chains, 0.0)

    # optional progress bar
    if progress
        iter = ProgressBar(1:niter)
    else
        iter = 1:niter
    end

    @inbounds for i in iter

        # get differential evolution proposal
        # draw the indices of the complementary chains
        i1 = collect(0:nchain-1) .+ rand(1:nchain-1, nchain)
        i2 = collect(0:nchain-1) .+ rand(1:nchain-2, nchain)
        i2[i2 .>= i1] .+= 1
        # add small noise and calculate proposal
        f = sigma * rand(Normal(0,1), (1,nchain))
        q = x + g0 * (x[:,(i1 .% nchain) .+ 1] - x[:,(i2 .% nchain) .+ 1]) .+ f
        factors = zeros(nchain)

        # log weight of current ensemble
        lweight = logsumexp(lprobs) + log(sum(accepted)) - log(nchain)

        # calculate stats for current ensemble
        ncov = cov(transpose(x))
        nmean = mean(x, dims=2)

        # update AIMH proposal distribution
        newcumlweight = logaddexp(cumlweight, lweight)
        statelweight = cumlweight - newcumlweight
        ccov = exp(statelweight) * ccov + exp(lweight - newcumlweight) * ncov
        cmean = exp(statelweight) * cmean + exp(lweight - newcumlweight) * nmean
        cumlweight = newcumlweight
        # get AIMH proposals if any chain is drawn
        xchnge = rand(Uniform(0,1), nchain) .<= aimh_prob

        # lp_debug(q)
        if sum(xchnge) > 0
            # draw alternative candidates and calculate their proposal density
            # println(ccov*(dft - 2)/dft + fixPSD)
            dist = MvTDist(dft, cmean[:], ccov*(dft - 2)/dft + fixPSD)
            xcand = rand(dist, sum(xchnge))
            lprop_old = logpdf(dist, x[:, xchnge])
            lprop_new = logpdf(dist, xcand)

            # update proposals and factors
            q[:,xchnge] = xcand
            factors[xchnge] = lprop_old - lprop_new
        end

        # Metropolis-Hasings 
        newlprob = lp(q)
        lnpdiff = factors + newlprob - lprob
        accepted = lnpdiff .> log.(rand(Uniform(0,1), nchain))
        naccepted = sum(accepted)
        # update chains
        x[:,accepted] = q[:,accepted]
        lprob[accepted] = newlprob[accepted]

        # store
        chains[i,:,:] = transpose(x)
        lprobs[i,:] = lprob

        if progress
            set_description(iter, string(@sprintf("[ll/MAF: %7.3f(%1.0e)/%2.0d%% | %1.0e]", maximum(lprob), std(lprob), 100*naccepted/nchain, statelweight)))
        end
    end

    return chains, lprobs, dist
end

initvar = 1
nchain = N*5#ndim*5 # a sane default
initcov = I(ndim)*initvar
initmean = zeros(ndim)
initchain = rand(MvNormal(initmean, initcov), nchain)

niter = 1000
chains, lprobs, propdist = RunDIME(initchain, niter, progress=true, aimh_prob=0.1)