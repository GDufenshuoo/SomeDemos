using Random
using LinearAlgebra
using Distributions

# 设置随机数生成器的种子
Random.seed!()

# 定义分子的能量
function energy(r)
    return -1.0 / (norm(r))
end

# 计算移动后分子的能量
function energy_new(r, dr)
    r_new = r .+ dr
    return energy(r_new)
end

# 定义 DMC 的主循环
function run_dmc(walkers, steps, Δt, E_guess)
    # 存储所有步骤的平均位置
    box = []
    sum
    E_o = E_guess
    E_n = 0.0
    # 迭代所有步骤
    for step in 1:steps
        # 迭代所有步行者
        for i in 1:length(walkers)
            # 随机生成步行者的位移
            dr = randn(3).* sqrt(Δt/3)

            # 计算步行者的新能量
            E_new = energy_new(walkers[i], dr)

            # 计算接受概率
            p_accept = exp(-Δt * (E_new - E_o))

            # 使用随机数决定是否接受此位移
            if rand() < p_accept
                E_n += E_new
                walkers[i] .+= dr
                push!(box,walkers[i] .+ dr)
            else
                walkers[i] = sample(box[])
            end
        end
        E_o = E_n/length(box)
        println(E_o)
    end     
    return box
end

# 定义初始条件
N_walkers = 1000
steps = 1000
Δt = 0.01

# 生成初始的步行者
walkers = [rand(3).-0.5 for i in 1:1000]

# 运行 DMC
box = run_dmc(walkers, steps, Δt,-1)

using Plots

l = zeros(length(box))
for i in 1:length(box)
    l[i] = norm(box[i])
end

lx = zeros(length(box))
ly = similar(lx)
for i in 1:length(box)
    lx[i] = box[i][1]
    ly[i] = box[i][2]
end

histogram(l)
histogram2d(lx,ly)

