

sir_ode = @ode_def SIRModel begin
    dS = -β * S * I
    dI = β * S * I - γ * I
    dR = γ * I
    end β γ


@model hierarchical_sir_model(grouped_sir, m_countries, country_pops, tspan = [0.0, 365.0], ::Type{T}=Float64) where {T<:Real} = begin
    let T = T, grouped_sir = grouped_sir, m_countries = m_countries, country_pops = country_pops, tspan = tspan

    timepoints = [size(grouped_sir[j], 2) for j in 1:m_countries]
    longest = maximum(timepoints)

    # hierarchical priors across countries
    #   locations
    β_loc ~ truncated(Normal(0, 0.5), 0, Inf) 
    γ_loc ~ truncated(Normal(0, 0.5), 0, Inf) 
    #   scales
    β_scale ~ truncated(Normal(0, 0.5), 0, Inf) 
    γ_scale ~ truncated(Normal(0, 0.5), 0, Inf)

    β = Vector{T}(undef, m_countries)
    γ = Vector{T}(undef, m_countries)
    for j in 1:m_countries
        # prior on infection rate
        β[j] ~ LogNormal(β_loc, β_scale)
        # prior on recovery rate
        γ[j] ~ LogNormal(γ_loc, γ_scale)
    end

    # global prior over measurement noise
    σ ~ LogNormal(0, 2)

    init = Vector{T}(undef, 3)
    sol = Array{T, 2}(undef, 3, longest)

    # likelihood
    for j = 1:m_countries
        init = [(country_pops[j] - 100.0) / country_pops[j], 100.0 / country_pops[j], 0.0]
        sir_prob = ODEProblem(sir_ode, init, tspan)

        sol = concrete_solve(sir_prob, Tsit5(), init, [β[j], γ[j]]; saveat = 1:longest)[:]

        if size(sol, 2) < timepoints[j]
            @logpdf() = zero(1) * sum(grouped_sir[j]) + -Inf
            return
        end

        grouped_sir[j] ~ arraydist([MvNormal(sol[:, i], σ) for i in 1:timepoints[j]])

    end
    end
end

