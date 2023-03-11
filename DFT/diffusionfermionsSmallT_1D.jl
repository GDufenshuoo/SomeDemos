using LinearAlgebra,Colors

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

function op(ψ, operater)
    lx = size(ψ)
    m = length(operater)-1
    𝜓 = zeros(lx)
    for k in eachindex(ψ[1,:])
        for i in eachindex(ψ[:,1])
            for dx in -m:m
                ri = o(i+dx,lx[1])
                𝜓[i,k] += ψ[ri,k] * operater[zop(dx)]
    end end end
    return 𝜓
end

freemove(x) = exp(-100*x^2/2)
function freemove(px::Int)
    return [freemove((i - px)/px) for i in px:2px]
end

function fermions_a!(ψ)
    for i in eachindex(ψ[:,1])
        if !iszero(ψ[i,:])
            max = findmax(ψ[i,:])
            for j in eachindex(ψ[i,:])
                ψ[i,j] = (j == max[2]) ? max[1] : 0
            end
    end end
    for i in eachindex(ψ[1,:])
        ψ[:,i] ./= sum(ψ[:,i])
    end
end

function fermions_b!(ψ)
    for i in eachindex(ψ[:,1])
        if !in(0,ψ[i,:])
            U = sum(ψ[i,:])
            for j in eachindex(ψ[i,:])
                a = 2*ψ[i,j] .- U
                ψ[i,j] = a < 0 ?  0 : a
    end end end
    for i in eachindex(ψ[1,:])
        ψ[:,i] ./= sum(ψ[:,i])
    end
end

function fermions_c!(ψ)
    ψ ./= sum(x->(x>0) ? x : 0,ψ)
    ϕ = copy(ψ)
    ϕ = op(ϕ,mask)
    ϕ = (ϕ.-ψ)./sum(x->(x>0) ? x : 0,ψ)
    for i in eachindex(ϕ)
        ψ[i] += (ψ[i] == 0) ? ϕ[i] : -ϕ[i]
    end
    ψ ./= sum(x->(x>0) ? x : 0,ψ)
end

function fermions!(ψ)
    for i in eachindex(ψ[:,1])
        if !iszero(ψ[i,:])
            U = sum(ψ[i,:])
            for j in eachindex(ψ[i,:])
                a = 2*ψ[i,j] .- U
                ψ[i,j] = a < 0 ?  0 : a
    end end end
    for i in eachindex(ψ[1,:])
        ψ[:,i] ./= sum(ψ[:,i])
    end
end

function bosons_c!(ψ)
    ϕ = copy(ψ)
    ϕ = op(ϕ,mask)
    ϕ_o = sum(ψ[:])
    ϕ .-= ψ
    ϕ ./= ϕ_o
    ψ .+= ϕ
    ψ ./= sum(ψ)
end

function bosons!(ψ)
    for i in eachindex(ψ[:,1])
        if !iszero(ψ[i,:])
            a = sum(ψ[i,:])
            ψ[i,:] .= a
    end end
    for i in eachindex(ψ[1,:])
        ψ[:,i] ./= sum(ψ[:,i])
    end
end

N = 2
mask = freemove(100)
girds = zeros(100,N)

using Plots
plot(girds)

N = 2
ϕ = zeros(100,N)
ϕ[45,1] = 1
ϕ[55:57,2] .= 1
ϕ = op(ϕ,mask)
fermions!(ϕ)
plot([ϕ[:,i] for i in 1:N])

@gif for i in 1:500
    global ϕ,mask
    ϕ = op(ϕ,mask)
    fermions!(ϕ)
    plot(ϕ)
end

N = 2
ϕ = zeros(100,N)
ϕ[45,1] = 1
ϕ[55,2] = 1
@gif for i in 1:200
    global ϕ,mask
    ϕ = op(ϕ,mask)
    fermions_b!(ϕ)
    plot(ϕ)
end

N = 2
ϕ = zeros(100,2)
ϕ[45,1] = 1
ϕ[55,2] = 1
ϕ1 = zeros(100,2)
ϕ1[45,1] = 1
ϕ1[55,2] = 1
@gif for i in 1:200
    global ϕ,ϕ1,mask
    ϕ = op(ϕ,mask)
    fermions!(ϕ)
    ϕ1 = op(ϕ1,mask)
    fermions_a!(ϕ1)
    plot([ϕ,ϕ1])
end

N = 2
ϕ = zeros(100,2)
ϕ[45,1] = 1
ϕ[55,2] = 1
ϕ1 = zeros(100,2)
ϕ1[45,1] = 1
ϕ1[55,2] = 1
@gif for i in 1:200
    global ϕ,ϕ1,mask
    ϕ = op(ϕ,mask)
    fermions!(ϕ)
    ϕ1 = op(ϕ1,mask)
    fermions_b!(ϕ1)
    plot([ϕ,ϕ1])
end

ϕ = zeros(500)
ϕ[[230,270]] .= 1
ϕ1 = copy(ϕ)
@gif for i in 1:2000
    global ϕ,mask
    fermions_c!(ϕ)
    bosons_c!(ϕ1)
    plot([ϕ./sum(abs,ϕ),abs.(ϕ),ϕ1])
end

plot([ϕ,ϕ1])

N = 2
ϕ = zeros(100,N)
ϕ[45,1] = 1
ϕ[55,2] = 1
@gif for i in 1:200
    global ϕ,mask
    ϕ = op(ϕ,mask)
    bosons!(ϕ)
    plot(ϕ)
end


N = 2
ϕ = zeros(100,2)
ϕ[45,1] = 1
ϕ[55,2] = 1
ϕ1 = zeros(100,1)
ϕ1[45,1] = 1
ϕ1[55,1] = 1
@gif for i in 1:200
    global ϕ,ϕ1,mask
    ϕ = op(ϕ,mask)
    bosons!(ϕ)
    ϕ1 = op(ϕ1,mask)
    bosons!(ϕ1)
    plot(ϕ.-ϕ1)
end
