import bayesian_nlde

# load data and format for inference
dat = import_data()
grouped_sir_data = [Matrix((df[:, [:susc_prop, :case_prop, :recov_prop]]))' for df in groupby(dat, :Country_Region)];
country_pops = [df.pop[1] for df in groupby(dat, :Country_Region)]
country_names = [df.Country_Region[1] for df in groupby(dat, :Country_Region)]

# look at countries with about a months worth of reporting
sub_idx = [i for (i,x) in enumerate(grouped_sir_data) if size(x, 2) > 28]
sub_grouped_sir_data = grouped_sir_data[sub_idx]
sub_country_pops = country_pops[sub_idx]

# infer paramters
@time chains = sample(hierarchical_sir_model(grouped_sir = biggest_grouped_sir_data,
    m_countries = length(biggest_grouped_sir_data), country_pops = biggest_country_pops),
    NUTS(1000, 0.65), 5000; progress = true, verbose = true)

# plot distribution over trajectories
for i in biggest_idx
    plt = plot_country_posterior(chains, 1, biggest_country_pops[1], biggest_country_names[1]; subsample_size = 100)
    plt_path = string(@__DIR__) * "/data/plots/" * string(country_names[i]) * ".svg"
    savefig(plt, plt_path)
end

