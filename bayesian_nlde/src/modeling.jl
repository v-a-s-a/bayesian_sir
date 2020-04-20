
# Turing.setadbackend(:zygote)
Turing.setadbackend(:reversediff)

sir_ode = @ode_def SIRModel begin
    dS = -β * S * I
    dI = β * S * I - γ * I
    dR = γ * I
    end β γ


@model hierarchical_sir_model(grouped_sir, m_countries, country_pops, tspan = [0.0, 365.0], ::Type{T}=Vector{Float64}) where {T} = begin

    let T = T, grouped_sir = grouped_sir, m_countries = m_countries, country_pops = country_pops, tspan = tspan

    # hierarchical priors across countries
    #   locations
    β_loc ~ truncated(Normal(0, 0.5), 0, Inf) 
    γ_loc ~ truncated(Normal(0, 0.5), 0, Inf) 
    #   scales
    β_scale ~ truncated(Normal(0, 0.5), 0, Inf) 
    γ_scale ~ truncated(Normal(0, 0.5), 0, Inf)

    β = T(undef, m_countries)
    γ = T(undef, m_countries)
    for j in 1:m_countries
        # prior on infection rate
        β[j] ~ LogNormal(β_loc, β_scale)
        # prior on recovery rate
        γ[j] ~ LogNormal(γ_loc, γ_scale)
    end

    # global prior over measurement noise
    σ ~ LogNormal(0, 2)

    init = [[(country_pops[j] - 100.0) / country_pops[j], 100.0 / country_pops[j], 0.0] for j in 1:m_countries]
    init = [convert(T, init[j]) for j in 1:m_countries]

    sir_probs = [ODEProblem(sir_ode, init[j], tspan) for j in 1:m_countries]
    timepoints = [size(grouped_sir[j], 2) for j in 1:m_countries]

    # likelihood
    for j = 1:m_countries

        sol = concrete_solve(sir_probs[j], Tsit5(), init[j], [β[j], γ[j]]; saveat = 1:timepoints[j])

        if size(sol, 2) < timepoints[j]
            @logpdf() = zero(1) * sum(grouped_sir[j]) + -Inf
            return
        end

        for i in 1:timepoints[j]
            grouped_sir[j][:, i] ~ MvNormal(sol[:, i], σ)
        end

    end
    end
end

