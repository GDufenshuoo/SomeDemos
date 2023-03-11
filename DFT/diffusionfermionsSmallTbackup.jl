using LinearAlgebra,Colors
gaussian(x) = 1 / sqrt(2π) * exp(-x^2 / 2)
gaussian(x, y) = 1 / sqrt(2π) * exp(-(x^2 + y^2) / 2)

function gaussian_mask(px::Int)
    result = [gaussian((i - px - 1) * 3 / px, (j - px - 1) * 3 / px) for i in 1:(2px + 1), j in 1:(2px + 1)]
    result /= sum(result)
end
mask = gaussian_mask(30)
mask ./ maximum(mask) .|> ColorTypes.Gray{Float32}

function ⊗(matrix::AbstractMatrix{T}, mask::AbstractMatrix) where {T}
    matrix_size = size(matrix)
    mask_size = size(mask)
    px = size(mask) .÷ 2
    result = zeros(T, matrix_size .- 2 .* px)
    for i in 1:size(result)[1], j in 1:size(result)[2]
        for idx in 1:mask_size[1], jdx in 1:mask_size[2]
            result[i, j] += matrix[i + idx - 1, j + jdx - 1] * mask[idx, jdx]
        end
    end
    return result
end

function gaussian_mask_1d(px::Int)
    result = [gaussian((i - px - 1) * 3 / px, (j - px - 1) * 3 / px) for i in 1:1, j in 1:(2px + 1)]
    # result /= sum(result)
end
mask = gaussian_mask_1d(20)
mask ./ maximum(mask) .|> ColorTypes.Gray{Float32}

girds = zeros(100,100)

girds[50:51,50:51].=0.25

ans = girds⊗mask

using Plots
heatmap(girds⊗mask)
