# using("Physics.jl")
include("CFinddistence.jl")

function W_bf(P)
    _E_N = E_N(P)
    W_bf=zeros(N+1)
    W_bf[1]=-1.0    #文献上 =1 我觉得不对
    for n=1:N
        for k=1:n
            W_bf[n+1] += (-exp(-beta*_E_N[k])* W_bf[n+1-k])/n
            #debug println(W_bf[n+1]," ",-exp(-beta*_E_N[k])," ",W_bf[n+1-k],"\n",W_bf)
        end
    end
    return W_bf[N+1]
end

function E_N(P)
    _E_N = zeros(N)
    _Rd2 = Rd2(P)
    for b in 1:B
        for k in 1:N
            for l=(N-k+1):N
                _E_N[k]+=_Rd2[l,b]
            end
        end
    end
    _E_N[1]+=_Rd2[1,B+1]-_Rd2[N,B]
    return _E_N.*k_05mw2.*a_0^2
end


