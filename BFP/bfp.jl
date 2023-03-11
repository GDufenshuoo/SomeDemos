
struct Fermi{T<:AbstractFloat} 
    Value::T
end

struct Boson{T<:AbstractFloat} 
    Value::T
end

function +(a::Fermi,b::Fermi)
    if a.Value >= b.Value
        return a
    else
        return b
    end
end

function -(a::Fermi,b::Fermi)
    return Fermi(a.Value - b.Value)
end

function +(a::Boson,b::Boson)
    return Boson(a.Value + b.Value)
end

function -(a::Boson,b::Boson)
    return a
end 

Fermi(rand()) + Fermi(rand())
Fermi(rand()) - Fermi(rand())

Boson(rand()) + Boson(rand())
Boson(rand()) - Boson(rand())


