using Random
using LinearAlgebra
using Distributions

# 设置随机数生成器的种子
Random.seed!()

# 定义分子的能量
function energy(r)
    V = 0.0
    for i in 1:3
        V -= 3/norm(r[i])
    end
    for i in 2:3, j in 1:i-1
        V += 1/norm(r[i].-r[j])
    end
    return V
end

# 定义 DMC 的主循环
function run_dmc(walkers, steps, Δt, E_guess)
    # 存储所有步骤的平均位置
    # box = [[] for i in 1:3]
    box = []
    sum
    E_o = E_guess
    E_n = 0.0
    N_walkers = length(walkers)
    N_sample = 0
    # 迭代所有步骤
    for step in 1:steps
        # 迭代所有步行者
        for w in 1:N_walkers
            for i in 1:3
                # 随机生成步行者的位移
                dr = randn(2).* sqrt(Δt/2)
                walkers[w][i] .+= dr

                # 计算步行者的新能量
                E_new = energy(walkers[w])

                # 计算接受概率
                p_accept = exp(-Δt * (E_new - E_o))

                # 使用随机数决定是否接受此位移
                if rand() < p_accept
                    E_n += E_new
                    push!(box,walkers[w][i])
                    N_sample += 1
                else
                    if length(box[i]) != 0
                        walkers[w][i] = sample(box) 
                    else
                        walkers[w][i] .-= dr
                    end
                end
            end
        end
        E_o = E_n/N_sample
        println("E: ",E_o," Acc: ", N_sample/(N_walkers*3*step))
    end     
    return box
end

# 定义初始条件
N_walkers = 1000
steps = 10000
Δt = 0.01

# 生成初始的步行者
N = 3
walkers = [[3*(randn(2).-0.5) for i in 1:N] for i in 1:1000]

# 运行 DMC
box = run_dmc(walkers, steps, Δt,-5)
ln = length(box)
ep = copy(box[div(ln,2):ln])
walkers = [[sample(ep) for i in 1:N] for i in 1:1000]
box = run_dmc(walkers, steps, Δt,-5)

using Plots

l=[]#[[] for i in 1:N]
for i in 1:3,j in 1:length(box)
    push!(l, norm(box[j]))
end 

lx = zeros(length(box))
ly = similar(lx)
for i in 1:length(box)
    lx[i] = box[i][1]
    ly[i] = box[i][2]
end

histogram(l)
histogram2d(lx,ly)

