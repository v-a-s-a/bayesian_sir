using bayesian_nlde
using DataFrames
using Plots
using Turing
using ReverseDiff
using DiffEqSensitivity

Turing.setadbackend(:reversediff)

# load data and format for inference
dat = import_data(download_data = false)
grouped_dat = [DataFrame(df) for df in groupby(dat, :Country_Region)];
grouped_sir_data = [Matrix((df[:, [:susc_prop, :case_prop, :recov_prop]]))' for df in groupby(dat, :Country_Region)];
country_pops = [df.pop[1] for df in grouped_dat]
country_names = [df.Country_Region[1] for df in grouped_dat]

# look at countries with more than a months worth of reporting
sub_idx = [i for (i,x) in enumerate(grouped_sir_data) if size(x, 2) > 28]
sub_grouped_sir_data = grouped_sir_data[sub_idx]
sub_country_pops = country_pops[sub_idx]

model = hierarchical_sir_model(grouped_sir = sub_grouped_sir_data,
    m_countries = length(sub_grouped_sir_data),
    country_pops = sub_country_pops)


varinfo = Turing.VarInfo(hierarchical_sir_model)
spl = Turing.SampleFromPrior()
@code_warntype model.f(varinfo, spl, Turing.DefaultContext(), model)


# infer paramters
@time chains = sample(hierarchical_sir_model(grouped_sir = sub_grouped_sir_data,
    m_countries = length(sub_grouped_sir_data), country_pops = sub_country_pops),
    NUTS(100, 0.8), 500; progress = true, verbose = true)
write(string(dirname(@__DIR__)) * "/data/chains/posterior_samples.jls", chains) 

#plotlyjs()

# plot distribution over trajectories
for i in sub_idx
    plt = plot_country_posterior(chains, grouped_dat, i; subsample_size = 200)
    plt_path = string(dirname(@__DIR__)) * "/data/plots/" * string(country_names[i]) * ".svg"
    savefig(plt, plt_path)
end

