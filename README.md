# VaccinationDistribution

# Research questions

The research of this repository raises the question whether the European Union could allocate the purchased vaccines more efficiently by more flexible vaccination strategies than constant rates. Attempting to answer this question, a deterministic Susceptible, Infectious, Recovered, Deceased (SIRD) model with two countries, two vaccines, and two virus types is developed and calibrated using parameters from the literature. Real-world vaccine purchase numbers of the EU are used as exogenous vaccine inflow and scaled down to the population size of the two-country model. Within the model, vaccines have varying efficacies with respect to the variants. Variants are distributed heterogeneously across countries. Piecewise constant functions and logistically transformed cubic Hermite splines are used as vaccination channels that determine the fractions of the vaccine doses each country receives. Both channels allow the fraction to be non-constant over the time course of the pandemic yielding potentially more complex vaccination strategies. As objective, the number of deaths caused by the pandemic is minimized over the parameters of the channel functions, once with additional Pareto constraints and once without the additional constraints. Subsequently, the optimized strategies are benchmarked against the current EU vaccination strategy. The deterministically derived optimal strategies are validated using a stochastic SIRD model.

# Tools

We use Python and mainly its libraries [libSBML](http://sbml.org/Main_Page) (Bornstein, 2008) and [AMICI](https://amici.readthedocs.io/en/latest/) (Fröhlich, 2021) to implement our models. [pyPESTO](https://pypesto.readthedocs.io/en/latest/) is our main tool for optimization.

## To-Do:
- Finish implementation of neural nets
- sensitivity analysis (hyperbolic tangent vs logistic)
- Implement inter-spatial infections by individuals that are in another country (new model)
- Incorporate reinfections
- Incorporate tests (Can we substitute vaccination by testing? Optimal distribution of vaccines and tests)
