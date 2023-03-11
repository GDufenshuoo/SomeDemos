using Random
function FFs(N)
    for i in 1:Mcond
        for k in 1:N
            m = shuffle!(collect(1:N))
            
        end
    end
end




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


