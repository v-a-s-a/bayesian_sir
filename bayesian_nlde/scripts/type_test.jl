
using Turing

@model tmodel1(x, y) = begin
    p,n = size(x)
    
    params = [Vector{Float64}(undef, p) for _ in 1:n]
    params[1] ~ MvNormal(zeros(p), ones(p))

    y[1] ~ MvNormal(params[i], 1.0)
    for i = 2:n
        y[i] ~ MvNormal(params[i], 1.0)
    end

end


model = tmodel1(ones(10), ones(10))
varinfo = Turing.VarInfo(tmodel1)
spl = Turing.SampleFromPrior()

@code_warntype model.f(varinfo, spl, Turing.DefaultContext(), tmodel1)

@model tmodel2(x, y, ::Type{T}=Array{Float64, (72, 3)}) where {T} = begin
    p,n = size(x)

    params = T(undef)
    params[1, :] ~ MvNormal(zeros(p), ones(p))
    y[1] ~ MvNormal(params[i], 1.0)

    for i = 2:n
        y[i] ~ MvNormal(x[i] * params[i-1], 1.0)
    end

end

model = tmodel2(ones(10), ones(10))
varinfo = Turing.VarInfo(tmodel1)
spl = Turing.SampleFromPrior()

@code_warntype model.f(varinfo, spl, Turing.DefaultContext(), tmodel2)