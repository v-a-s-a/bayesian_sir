
using bayesian_nlde
using DifferentialEquations
using DataFrames
using Plots, StatsPlots
using Turing
using Random
# using DiffEqSensitivities
# using DistributionsAD
using SpecialFunctions

Random.seed!(12035)

# Turing.setadbackend(:forwarddiff)
Turing.setadbackend(:reversediff)
# Turing.setadbackend(:zygote)


@model eSIR_model(sir, population, tspan = [0.0, 365.0], ::Type{T}=Vector{Float64}) where {T} = begin

    timepoints = size(sir, 1)

    # reproduction number
    R_0 ~ LogNormal(1.099, 0.096)
    
    # recovery rate
    γ ~ LogNormal(-2.955, 0.910)

    # infection rate
    β = R_0 * γ

    # time dependent infection rate modifier
    # π_t = Vector{Real}(undef, timepoints)

    # priors over observed rate variances 
    λ_I ~ Gamma(2.0, 1.0 / 0.001)
    λ_R ~ Gamma(2.0, 1.0 / 0.001)

    # prior over dirichlet variance
    κ ~ Gamma(2.0, 1.0 / 0.001)

    # priors over initial proportions
    θ_t = [T(undef, 3) for _ in 1:timepoints]
    θ_t[1][1] ~ Beta(1.0, 1.0 / sir[1][1])   # S
    θ_t[1][2] ~ Beta(1.0, 1.0 / sir[1][2])   # I
    θ_t[1][3] = 1.0 - θ_t[1][1] - θ_t[1][2]  # R

    # likelihood
    for t in 2:timepoints

        # disease dynamics at time t
        sir_prob = ODEProblem(sir_ode, θ_t[t-1], tspan)
        # println(θ_t[t-1])
        sol = concrete_solve(sir_prob, Tsit5(), θ_t[t-1], [β, γ]; saveat = [0, t-1])
    
        if size(sol, 2) < 2
            @logpdf() = zero(1) * sum(sir[t]) + -Inf
            return
        end

        # probably not a good idea
        sol[2][sol[2] .< 0] .= eps()

        # latent infection and recovery proportions
        θ_t[t] ~ Dirichlet(κ * sol[2])

        # also probably not a good idea        
        θ_t[t][1] = θ_t[t][1] - eps()
        θ_t[t][2] = θ_t[t][2] + eps()
        θ_t[t][3] = θ_t[t][3] + eps()

        # observed infection proportion
        sir[t][2] ~ Beta(λ_I * θ_t[t][2], λ_I * (1.0 - θ_t[t][2]))
        # observed recovery proportion
        sir[t][3] ~ Beta(λ_R * θ_t[t][3], λ_R * (1.0 - θ_t[t][3]))
    end

end


@model eSIR_model2(S, I, R, population, tspan = [0.0, 365.0], ::Type{T}=Vector{Float64}) where {T} = begin

    timepoints = size(I, 1)

    # reproduction number
    R_0 ~ LogNormal(1.099, 0.096)
    
    # recovery rate
    γ ~ LogNormal(-2.955, 0.910)

    # infection rate
    β = R_0 * γ

    # time dependent infection rate modifier
    # π_t = Vector{Real}(undef, timepoints)

    # priors over observed rate variances 
    λ_I ~ Gamma(2.0, 1.0 / 0.0001)
    λ_R ~ Gamma(2.0, 1.0 / 0.0001)

    # prior over dirichlet variance
    κ ~ Gamma(2.0, 1.0 / 0.0001)

    # priors over initial proportions
    # θ_t = [Vector{Float64}(undef, 3) for _ in 1:timepoints]

    # θ_t = [ones(Float64, 3) for _ in 1:timepoints]
    # θ_t = T(undef, timepoints)
    let T=T

        θ_t = [T(undef, 3) for _ in 1:timepoints]
        # println(typeof(θ_t))

        # @show θ_t[1:2]
        θ_t[1][1] ~ Beta(1.0, 1.0 / S[1]) # S
        # println("hit")

        θ_t[1][2] ~ Beta(1.0, 1.0 / I[1])   # I
        # println("hit")

        θ_t[1][3] = 1.0 - θ_t[1][1] - θ_t[1][2]  # R
        # println("hit")
        # @show θ_t[1:2]


        # likelihood
        for t in 2:timepoints

            # disease dynamics at time t
            sir_prob = ODEProblem(sir_ode, θ_t[t-1], tspan)
            # println(θ_t[t-1])
            sol = concrete_solve(sir_prob, Tsit5(), θ_t[t-1], [β, γ]; saveat = [0, t-1])

            if size(sol, 2) < 2
                @logpdf() = zero(1) * sum(sir[t]) + -Inf
                return
            end

            # probably not a good idea
            sol[2][sol[2] .< 0] .= eps()

            # println("sol at $t: ", sol[2])

            # latent infection and recovery proportions
            θ_t[t] ~ Dirichlet(κ * sol[2])

            # also probably not a good idea        
            θ_t[t][1] = θ_t[t][1] - eps()
            θ_t[t][2] = θ_t[t][2] + eps()
            θ_t[t][3] = θ_t[t][3] + eps()

            # observed infection proportion
            I[t] ~ Beta(λ_I * θ_t[t][2], λ_I * (1.0 - θ_t[t][2]))
            # observed recovery proportion
            R[t] ~ Beta(λ_R * θ_t[t][3], λ_R * (1.0 - θ_t[t][3]))
        end
    end
    # @show θ_t[1:2]


end



@model eSIR_model3(S, I, R, tspan = [0.0, 72.0], ::Type{T}=Matrix{Float64}) where {T} = begin

    timepoints = size(I, 1)
    

    # reproduction number
    R_0 ~ LogNormal(1.099, 0.096)
    
    # recovery rate
    γ ~ LogNormal(-2.955, 0.910)

    # infection rate
    β = R_0 * γ

    # time dependent infection rate modifier
    # π_t = Vector{Real}(undef, timepoints)

    # priors over observed rate variances 
    λ_I ~ Gamma(2.0, 1.0 / 0.0001)
    λ_R ~ Gamma(2.0, 1.0 / 0.0001)

    # prior over dirichlet variance
    κ ~ Gamma(2.0, 1.0 / 0.0001)

    # priors over initial proportions
    # θ_t = [Vector{Float64}(undef, 3) for _ in 1:timepoints]

    # θ_t = [ones(Float64, 3) for _ in 1:timepoints]
    # θ_t = T(undef, timepoints)

    θ_t = T(undef, timepoints, 3)
    # println(typeof(θ_t))

    # @show θ_t[1:2]
    θ_t[1, 1] ~ Beta(1.0, 1.0 / S[1]) # S
    # println("hit")

    θ_t[1, 2] ~ Beta(1.0, 1.0 / I[1])   # I
    # println("hit")

    θ_t[1, 3] = 1.0 - (θ_t[1, 1] - θ_t[1, 2])  # R
    # println("hit")
    # @show θ_t[1:2]


    # likelihood
    for t in 2:timepoints

        # disease dynamics at time t
        sir_prob = ODEProblem(sir_ode, θ_t[t-1, :], tspan)

        # @show θ_t[1, :]

        @show θ_t[t-1, :]
        @show [β, γ]

        sol = concrete_solve(sir_prob, Tsit5(), θ_t[t-1, :], [β, γ]; saveat = [0, 1])

        @show sol[2]
        @show κ
        # @show θ_t[t-1, :]
        # @show [β, γ]
        @show κ * sol[2]
        # sleep(0.1)

        if size(sol, 2) < 2
            @logpdf() = zero(1) * (S[t] + I[t] + R[t]) + -Inf
            return
        end

        # probably not a good idea
        # sol[2][sol[2] .<= 0.0] .= eps()

        # println("sol at $t: ", sol[2])
        # @show sol[2]
        # latent infection and recovery proportions
        # θ_t[t, :] ~ Dirichlet(κ * (sol[2] .+ eps()))
        θ_t[t, :] ~ Dirichlet(κ * sol[2])

        # @show κ * (sol[2] .+ eps())

        # observed infection proportion
        I[t] ~ Beta((λ_I * θ_t[t, 2]) + eps(), (λ_I * (1.0 - θ_t[t, 2])) + eps())
        # observed recovery proportion
        R[t] ~ Beta((λ_R * θ_t[t, 3]) + eps(), (λ_R * (1.0 - θ_t[t, 3])) + eps())

    end

# @show θ_t[1:2]


end


# # load data and format for inference
dat = import_data(download_data = false)
grouped_dat = [DataFrame(df) for df in groupby(dat, [:Province_State, :Country_Region])];
grouped_sir_data = [Matrix((df[:, [:susc_prop, :case_prop, :recov_prop]]))' for df in groupby(dat, [:Province_State, :Country_Region])];
country_pops = [df.pop[1] for df in grouped_dat]
country_names = [df.Country_Region[1] for df in grouped_dat]

# look at countries with more than a months worth of reporting
sub_idx = [i for (i,x) in enumerate(grouped_sir_data) if size(x, 2) > 28]
sub_grouped_sir_data = grouped_sir_data[sub_idx]
sub_country_pops = country_pops[sub_idx]

hubei_ir = Matrix(sub_grouped_sir_data[1][[2,3], :]')


# china_sir = [sub_grouped_sir_data[1][:, t] for t in 1:size(sub_grouped_sir_data[1], 2)]

model = eSIR_model3(S = sub_grouped_sir_data[1][1, :],
                    I = sub_grouped_sir_data[1][2, :],
                    R = sub_grouped_sir_data[1][3, :])

varinfo = Turing.VarInfo(eSIR_model3)
spl = Turing.SampleFromPrior()
# @code_warntype model.f(varinfo, spl, Turing.DefaultContext(), model)




# infer paramters
# @time chains = sample(eSIR_model2(sir = china_sir, population = sub_country_pops[1]),
    # NUTS(50, 0.80), 500; progress = true, verbose = true)
@time chains = sample(eSIR_model3(S = sub_grouped_sir_data[1][1, :],
    I = sub_grouped_sir_data[1][2, :],
    R = sub_grouped_sir_data[1][3, :]),
    # NUTS(5, 0.80), 10; progress = true, verbose = true)
    SMC(), 10; progress = true, verbose = true)
write(string(dirname(@__DIR__)) * "/data/chains/eSIR_china_posterior_samples.jls", chains) 

chains_dat = DataFrame(chains)

θ_S = chains_dat[:, [Symbol("θ_t[$i,Colon()][1]") for i in 2:size(sub_grouped_sir_data[1], 2)]]
θ_I = chains_dat[:, [Symbol("θ_t[$i,Colon()][2]") for i in 2:size(sub_grouped_sir_data[1], 2)]]
θ_R = chains_dat[:, [Symbol("θ_t[$i,Colon()][3]") for i in 2:size(sub_grouped_sir_data[1], 2)]]

θ_S_mean = Array(aggregate(θ_S, mean))
θ_I_mean = Array(aggregate(θ_I, mean))
θ_R_mean = Array(aggregate(θ_R, mean))

plotly()
scatter(θ_S_mean', label="S")
scatter!(θ_I_mean', label="I")
scatter!(θ_R_mean', label="R")


mean(chains_dat[:, :γ])
mean(chains_dat[:, :R_0])

# get mean latent proportions


# # plot distribution over trajectories
# for i in sub_idx
#     plt = plot_country_posterior(chains, grouped_dat, i; subsample_size = 200)
#     plt_path = string(dirname(@__DIR__)) * "/data/plots/" * string(country_names[i]) * ".svg"
#     savefig(plt, plt_path)
# end


