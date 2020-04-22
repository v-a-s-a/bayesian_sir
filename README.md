# Bayesian parameter estimates for non-linear systems of differential equations

**Please note, the following work is for educational purposes only. If you'd like to draw valid inferences, please consult an infectious disease epidemiologist. Modeling the outbreak is [hard](https://fivethirtyeight.com/features/why-its-so-freaking-hard-to-make-a-good-covid-19-model/), and this data must be used [responsibly](https://medium.com/nightingale/ten-considerations-before-you-create-another-chart-about-covid-19-27d3bd691be8) (thanks to [@simonw](https://github.com/simonw/covid-19-datasette/blob/master/README.md) for raising these points).**

This was inspired by Thomas Wiecki's [work](https://github.com/twiecki/covid19) on fitting PyMC3 models to the COVID data, as well as the elegant [`DiffEqBayes.jl`](https://github.com/JuliaDiffEq/DiffEqBayes.jl) package.

## Data

Case data comes from [Johns Hopkins Center for Systems Science and Engineering](https://github.com/CSSEGISandData/COVID-19). Country population data comes from the [European Centre for Disease Prevention and Control ](https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide).

## Model

Currently, the model is hierarchical SIR model, with infection and recovery rates grouped by country. The system is specified and solved using [`DifferentialEquations.jl`](https://pkg.julialang.org/docs/DifferentialEquations/UQdwS/6.6.0/) and parameters are infered using [`Turing.jl`](https://turing.ml/dev/).

## Results

The following show distributions of disease trajectories over time. Thick lines represent the expected posterior distribution of disease trajectory.

Data/fit from 20.04.2020:

![US](bayesian_nlde/data/plots/US.svg)
![Germany](bayesian_nlde/data/plots/Germany.svg)
![China](bayesian_nlde/data/plots/China.svg)
![Italy](bayesian_nlde/data/plots/Italy.svg)
![Switzerland](bayesian_nlde/data/plots/Switzerland.svg)
![Spain](bayesian_nlde/data/plots/Spain.svg)
![Iran](bayesian_nlde/data/plots/Iran.svg)
![Japan](bayesian_nlde/data/plots/Japan.svg)
![Singapore](bayesian_nlde/data/plots/Singapore.svg)

