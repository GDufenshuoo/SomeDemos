function o(i, N)
	n = 0
	if i == 1
		n = N
	else
		n = i-1
	end
	return n
end

function fo(i, N)
	n = 0
	if i == N
		n = 1
	else
		n = i+1
	end
	return n
end

function bosons()
end

"""
Forward fermions sampling
"""
function fermions(P)
    d2,d2_s = d2(P)
    db2,db2_s = db2(P)
end

function d2(P)
    N = size(P)[1]
    B = size(P)[3]

    d2 = zeros(N,N,B)
    for i=1:N
        d2[i,i,:].=1e33
    end
    for b=1:B
        for i=2:N
            for j=1:i-1
                d = sum(abs2,P[i,:,b].-P[j,:,b])
                d2[i,j,b] = d
                d2[j,i,b] = d
            end
        end
    end
    return d2
end


function d2_Reshape(P)
    N = size(P)[1]
    B = size(P)[3]

    d2 = zeros(N,N,B)
    for i=1:N
        d2[i,i,:].=1e33
    end
    for b=1:B
        for i=2:N
            for j=1:i-1
                d = sum(abs2,P[i,:,b].-P[j,:,b])
                d2[i,j,b] = d
                d2[j,i,b] = d
            end
        end
        d2[:,:,b] .= sort(d2[:,:,b], dims=2)
    end
    return d2
end

"""
According to 
    Hirshberg, Invernizzi, and Parrinello, 
    “Path Integral Molecular Dynamics for Fermions...”

    It's not necessary to consider all permutation.
"""
function Rd2(P)
    N = size(P)[1]
    B = size(P)[3]

    Rd2 = zeros(N,B+1)
    Rd2.= 1e33
    for b in 1:B
        for i=1:N
            Rd2[i,b] = sum(abs2,P[i,:,b].-P[fo(i,N),:,fo(b,B)])
        end
    end
    for i=1:N
        Rd2[i,B+1] = sum(abs2,P[N-i+1,:,B].-P[N,:,1])
    end
    return Rd2
end

function db2(P)
    N = size(P)[1]
    B = size(P)[3]

    db2_sp = zeros(N,N,B)
    db2_sf = zeros(N,N,B)
    db2 = zeros(N,N,B)

    for b=1:B
        for j=1:N
            for i=1:N
                d = sum(abs2,P[i,:,b].-P[j,:,o(b,B)])
                db2[i,j,b] = d
            end
            db2_sp[:,j,b] = sort(db2[:,j,b])
        end
        for i in 1:N
            db2_sf[i,:,b] = sort(db2[i,:,b])
        end
    end
    return db2,db2_sp,db2_sf
end
