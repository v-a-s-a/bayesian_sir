module bayesian_nlde

    using DifferentialEquations
    using Turing
    using CSV
    using DataFrames
    using Dates
    using Plots
    using StatsPlots
    using Statistics
    using Plots.PlotMeasures
    using Pipe
    import NamedTupleTools.namedtuple

    include("preprocessing.jl");
    include("modeling.jl");
    include("criticism.jl");

    export sir_ode, hierarchical_sir_model
    export import_data, import_hopkins_data, import_edcd_data
    export plot_country_posterior

end # module
