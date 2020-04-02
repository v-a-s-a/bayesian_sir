
# Turing.setadbackend(:zygote)
Turing.setadbackend(:forward_diff)

sir_ode = @ode_def SIRModel begin
    dS = -β * S * I
    dI = β * S * I - γ * I
    dR = γ * I
    end β γ

@model hierarchical_sir_model(grouped_sir, m_countries, country_pops) = begin

    # hierarchical priors across countries
    #   locations
    β_loc ~ truncated(Normal(0, 0.5), 0, Inf) 
    γ_loc ~ truncated(Normal(0, 0.5), 0, Inf) 
    #   scales
    β_scale ~ truncated(Normal(0, 0.5), 0, Inf) 
    γ_scale ~ truncated(Normal(0, 0.5), 0, Inf)

    β = Vector{Real}(undef, m_countries)
    γ = Vector{Real}(undef, m_countries)
    for j in 1:m_countries
        # prior on infection rate
        β[j] ~ LogNormal(β_loc, β_scale)
        # prior on recovery rate
        γ[j] ~ LogNormal(γ_loc, γ_scale)
    end

    # global prior over measurement noise
    σ ~ LogNormal(0, 2)

    # likelihood
    for j = 1:m_countries
        timepoints = size(grouped_sir[j], 2)
        init = [(country_pops[j] - 100.0) / country_pops[j], 100.0 / country_pops[j], 0.0]
        sol = concrete_solve(sir_prob, Tsit5(), init, [β[j], γ[j]]; saveat = 1:timepoints)

        if size(sol, 2) < timepoints
            @logpdf() = zero(1) * sum(grouped_sir[j]) + -Inf
            return
        end

        for i in 1:timepoints
            grouped_sir[j][:, i] ~ MvNormal(sol[:, i], σ)
        end

    end
end

