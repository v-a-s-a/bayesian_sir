
using MacroTools    

# write a macro that generates numbered variables

ex = quote
    θ_t[1][1]
end

@capture(ex, var__[idx1_][idx2_])

## rewrite: θ_t = [T(undef, 3) for _ in 1:timepoints]
## to:      θ_1 = T(undef, 3); θ_2 = T(undef, 3); ...
function init_vars(basevname, vlength)
    expr = ""
    for i in 1:vlength
        exprs[i] =
        quote
            $basevname_$i = T(undef, 3)
        end
    end
    return exprs
end

println(init_vars(:θ, 10))

## rewrite: θ_t[1][1] ~ Beta(1.0, 1.0 / S[1])
## to:      θ_1[1] ~ Beta(1.0, 1.0 / S[1])

## rewrite: θ_[t][3] = 1.0 - θ_t[1][1] - θ_t[1][2] 
## to:      θ_1[3] = 1.0 - θ_1[1] - θ_1[2] 

## rewrite: sir_prob = ODEProblem(sir_ode, θ_t[t-1], tspan)
## to:      sir_prob = ODEProblem(sir_ode, θ_$(t-1), tspan)

## rewrite: θ_t[t] ~ Dirichlet(κ * sol[2])

## rewrite: θ_t[t][1] = θ_t[t][1] - eps()

## rewrite: I[t] ~ Beta(λ_I * θ_t[t][2], λ_I * (1.0 - θ_t[t][2]))

createvar = macro (varname)
    @capture(ex, varname__[idx1_][idx2_])
end

for t in 1:5
    @createvar (θ_t) 
end