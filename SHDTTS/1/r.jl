function Tr(r)
    n = 1000
    a = zeros(n)

    a[1] = 2-r
    for i in 2:n
        a[i] = (-1)^i * (r-1)^i + a[i-1]
    end
    return a[n]
end

using Plots
f(x) = 1/x
plot(0.1:0.01:2,[Tr,f])
g(x) = f(x) - Tr(x)
plot(0.0001:0.00001:0.001,g)
g(big(0.01))
g(x) = (f(x) - Tr(x))*x
plot(1:0.1:2,g)
g(big(0.01))

using BenchmarkTools
@benchmark Tr(rand())
@benchmark 1/rand()
