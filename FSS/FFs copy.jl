using Random
using LinearAlgebra

function markov_chain_sampling(x::AbstractVector, m::AbstractVector, δx::Real, Mcond::Integer)
    N = length(x)
    k = findfirst(isequal(N), m)
    ux = view(U(x, m[1:k]), :, k) # get k-th row of U(x, m)
    h = nullspace(transpose(ux))[:, 1] # get normal vector to ux1, ..., uxk-1
    Uxk = ux * (ux' * h) # project h onto ux to get U(x1, ..., xk+δx, mk)
    Uxkδ = U(x, m)[:, 1:k] # get U(x1, ..., xk, mk)
    Uxkδ[:, k] += δx # make the trial update xk → xk + δx
    Uxkδ = U(Uxkδ, m[1:k]) # get U(x1, ..., xk+δx, mk)
    α = det(Uxk) / det(ux) # compute acceptance probability
    α *= α # compute squared modulus of α
    for i in 1:Mcond # iterate to let the Markov chain equilibrate
        if rand() < α
            x[k] += δx # accept trial update
            ux = view(U(x, m[1:k]), :, k)
            h = nullspace(transpose(ux))[:, 1]
            Uxk = ux * (ux' * h)
            Uxkδ = U(x, m)[:, 1:k]
            Uxkδ[:, k] += δx
            Uxkδ = U(Uxkδ, m[1:k])
            α = det(Uxk) / det(ux)
            α *= α
        end
    end
end

function U(x::AbstractVector, m::AbstractVector)
    N = length(x)
    U = Matrix{ComplexF64}(I, N, N)
    for k in m
        ux = view(U, :, k)
        for j in k+1:N
            vj = randn(ComplexF64, N-j+1)
            vj -= ux' * vj / (ux' * ux) * ux
            vj /= norm(vj)
            ux[j:end] = vj
        end
    end
    for k in m[end:-1:1]
        ux = view(U, :, k)
        for j in 1:k-1
            vj = randn(ComplexF64, j)
            vj -= vj' * ux[1:j] / (ux[1:j]' * ux[1:j]) * ux[1:j]
            vj /= norm(vj)
            ux[1:j-1] = vj[1:j-1]
        end
    end
    for k in m
        ux = view(U, :, k)
        h = nullspace(transpose(ux))[:, 1]
        U *= diagm(exp.(im .* x[k] * h))
        U = ux * (ux' * U)
    end
    U
end

function FFS_continuous(N::Integer, Mcond::Integer, δx::Real)
    m = randperm(N)
    x = zeros(N)
    for k in 1:N
        markov_chain_sampling(x, m, δx, Mcond)
    end
    x
end



FFS_continuous(10, 10, 0.1)


function get_normal_vector(U::Matrix{Float64}, k::Int)
    # Initialize normal vector as zeros
    h = zeros(k)

    # Perform iterative Gaussian elimination
    for j = 1:k-1
        for i = j+1:k
            if U[i,j] != 0
                multiplier = U[i,j]/U[j,j]
                U[i,j+1:end] -= multiplier * U[j,j+1:end]
                h[j] += multiplier^2
            end
        end
    end

    # Normalize the normal vector
    h[k] = 1.0
    for j = 1:k-1
        h[j] = sqrt(h[j])
    end

    return h
end


