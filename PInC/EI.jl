

function H(P,N,mₑ)
    return ΔE(P) - 1/sqrt(sum(abs2,P[1:3,1]))
end

"""
f R = m v² 

f R = m v² / √(1-v²)
a = f R / m
-a + a √(4 + a²)
"""
function ΔE(P)
    m = 1
    f = 1/sum(abs2,P[1:3,1])
    R = √(sum(abs2,P[1:3,1]))
    a = f*R / m
    v² = 0.5*(√(a*(4+a))-a)
    return v²/2
end

using Plots

f(x) = ΔE([x,0,0])-1 / x
g(x) = H([x,0,0],1,1)

plot(0.0:0.01:2,g)
66