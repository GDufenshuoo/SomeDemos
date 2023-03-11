using LinearAlgebra,Colors
difussion(x) = exp(-100*x^2/2)

function o(i,l)
    if i < 1
        return l+i
    elseif i > l
        return i-l
    else
        return i
    end
end

function zop(i)
    return abs(i)+1
end

function op!(ψ,operater,N)
    lx = length(ψ)
    m = length(operater)-1
    ψ⁺ = zeros(lx)

    nom = zeros(2)
    for i in 1:lx
        for dx in -m:m
            ri = o(i+dx,lx)
            ro = ψ[ri] * operater[zop(dx)]
            if N[i] == N[ri] || N[i] == 0 || N[ri] == 0
                ψ⁺[i] += ro
                N[i] = N[ri]
            elseif ψ⁺[i] > ro
                ψ⁺[i] = ψ⁺[i] - ro
            elseif  ψ⁺[i] < ro
                ψ⁺[i] = ro - ψ⁺[i]
                N[i] = N[ri]
    end end end 
    for i in eachindex(ψ⁺)
        if N[i] != 0
            nom[N[i]]+=ψ⁺[i] 
        end
    end
    for i in eachindex(ψ⁺)
        if N[i] != 0
            ψ⁺[i] /= nom[N[i]]
        end
    end
    for i in eachindex(ψ⁺)
        if N[i] == 0
            ψ⁺[i] = 0
        end
    end
    ψ = ψ⁺
end

function gaussian_mask_1d(px::Int)
    [difussion((i - px)/px) for i in px:2px]
end
mask = gaussian_mask_1d(5)
mask ./ maximum(mask) .|> ColorTypes.Gray{Float32}

girds = zeros(100)
N = zeros(Int,100)

using Plots
heatmap(girds)
heatmap((girds⊗mask)[1])

plot(girds)


ans = girds
ans = op(ans,mask)
a = ans
plot(a)

girds = zeros(100)
N = zeros(Int,100)
girds[20:30].=10
N[20:30] .= 1
girds[77:87].=20
N[77:87] .= 2
@gif for i in 1:100
    global girds
    girds = op!(girds,mask,N)
    plot(girds)
end

girds = op!(girds,mask,N)

A = copy(girds)
B = copy(girds)
for i in eachindex(N)
    if N[i] != 1
        A[i] = 0
    end
end
for i in eachindex(N)
    if N[i] != 2
        B[i] = 0
    end
end
plot(A)
plot(B)
plot(N)

