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

function ⊗(ψ, operater)
    lx,ly = size(ψ)
    m = length(operater)-1
    gx = similar(ψ)
    for i in 1:lx, j in 1:ly
        for dx in -m:m
            ri = o(i+dx,lx)
            ro = ψ[ri, j] * operater[zop(dx)]
            if N[i,j] == N[ri,j]
                gx[i, j] += ro
            elseif gx[i, j] < ro
                gx[i, j] = ro
                N[i,j] = N[ri,j]
    end end end 
    return gx
end

function gaussian_mask_1d(px::Int)
    ans = [difussion((i - px)/px) for i in px:2px]
end
mask = gaussian_mask_1d(5)
mask ./ maximum(mask) .|> ColorTypes.Gray{Float32}

girds = zeros(100)
N = zeros(Int,100)

girds[50:60].=10

using Plots
heatmap(girds)
heatmap((girds⊗mask)[1])

plot(girds)

plot(girds⊗mask)
