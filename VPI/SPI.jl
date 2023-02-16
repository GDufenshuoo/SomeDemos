"""

Atomic Unit 
"""
𝑐 = 1/0.0072973525628

U(r) = -1/r
F(r) = -1/r^2
Lₘ(r) = 1
Pₘ(r) = 3/r

""" 
fr = mv2
    p = L/r = nħ/r
fr = p^2/m
√kₑm = p = mv
∴ v = 1/sqrt(m)

∵ √kₑm = p = L/r = nħ/r
∴ r√kₑm = nħ

check: ħ = 1.0545718176461565e-34
        r√m = 5.291772109022829e-11*√(9.109383701558256e-31)
            = 5.050630892238351e-26
        kₑ = 4.359744722207211e-18
        E = (1.0545718176461565e-34/5.291772109022829e-11)^2
            /(2*9.109383701558256e-31)
            + 4.359744722207211e-18/5.291772109022829e-11


fr³ = L² 
f = 1/r
L = n/r
1/r2 = n/r^2
"""

using Plots
Ex(r) =  max(Pₘ(r)^2,1)/2 + U(r)

plot(0:0.01:1,[Pᵣ,Pₘ])

"""
a = fr/m
"""
function Tᵣ(r)
    a = F(r)*r/m
    a² = a^2
    v² = (sqrt(a²*(4+a²))-a²)/2
    return 0.5*m*v²/√(1-v²)
end

"""
a = fr/m
"""
function Pᵣ(r)
    a = F(r)*r/m
    a² = a^2
    v² = (sqrt(a²*(4+a²))-a²)/2
    return m/√(1/v²-1)
end

f(r) = Pᵣ(r)/Pₘ(r)

plot(3:0.01:10,[Pᵣ,Pₘ])

sample_size =3000
dataset = randn(sample_size)
using Distributions
new_sample = sample(dataset,10)

using Plots
plot(dataset)
plot(new_sample)
