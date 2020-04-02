

function plot_country_posterior(chains, grouped_dat, country; tspan = [0.0, 365.0], subsample_size=100)

    country_pop = country_dat.pop[1]
    country_name = country_dat.pop[1]

    post_dat = DataFrame(chains);

    init = [(country_pop - 100.0) / country_pop, 100.0 / country_pop, 0.0]
    country_param_cols = [Symbol("β[" * string(country) * "]"), Symbol("γ[" * string(country) * "]")]

    posterior_mean_params = mean.(eachcol(post_dat[:, country_param_cols]))
    posterior_mean_prob = ODEProblem(sir_ode, init, tspan, post_dat[i, country_param_cols])
    posterior_mean_sol = solve(posterior_mean_prob);

    plt = plot(mean_sol;
        title = country_name * "\n",
        title_location = :left,
        top_margin = 5mm,
        vars = [1],
        xguidefontsize = 10, yguidefontsize = 10,
        xlabel = "Days since 100 confirmed cases",
        ylabel = "Proportion of population",
        label = "Susceptible",
        color = "blue",
        linewidth = 2)
    plot!(mean_sol; label="Infected", vars = [2], color = "red", linewidth = 2)
    plot!(mean_sol; label="Recovered", vars = [3], color = "green", linewidth = 2)

    scatter!()

    # subsample to avoid a crowded plot
    sub_idx = sample(axes(post_dat, 1), subsample_size; replace = false);
    for i in sub_idx
        prob = ODEProblem(sir_ode, init, tspan, post_dat[i, country_param_cols])
        sol = solve(prob);
        plot!(sol; vars = [1], label="", color = "blue", alpha = 0.1, xlabel = "Days since 100 confirmed cases")
        plot!(sol; vars = [2], label="", color = "red", alpha = 0.1, xlabel = "Days since 100 confirmed cases")
        plot!(sol; vars = [3], label="", color = "green", alpha = 0.1, xlabel = "Days since 100 confirmed cases")
    end

    return plt
end