# VaccinationDistribution

# Research question

How can a supranational institution allocate scarce vaccines across countries in tims of a pandemic?


### Description
The research of this repository raises the question whether the European Union could allocate the purchased vaccines more efficiently by more flexible vaccination strategies than constant rates. Attempting to answer this question, a deterministic Susceptible, Infectious, Recovered, Deceased (SIRD) model with two countries, two vaccines, and two virus types is developed and calibrated using parameters from the literature. Real-world vaccine purchase numbers of the EU are used as exogenous vaccine inflow and scaled down to the population size of the two-country model. Within the model, vaccines have varying efficacies with respect to the variants. Variants are distributed heterogeneously across countries. Piecewise constant functions and logistically transformed cubic Hermite splines are used as vaccination channels that determine the fractions of the vaccine doses each country receives. Both channels allow the fraction to be non-constant over the time course of the pandemic yielding potentially more complex vaccination strategies. As objective, the number of deaths caused by the pandemic is minimized over the parameters of the channel functions, once with additional Pareto constraints and once without the additional constraints. Subsequently, the optimized strategies are benchmarked against the current EU vaccination strategy. The deterministically derived optimal strategies are validated using a stochastic SIRD model.

### Findings 
Our results show that the optimal derived vaccination strategies differ from the current EU strategy, which could indicate that more complex vaccination policies can lower the number of death cases caused by the pandemic. Enhancing the robustness of the results, we qualitatively find highly similar results using piecewise constant and spline vaccination channels. Leaving aside country-specific interests, it appears that assigning most vaccine doses to one country can be as much as 17\% beneficial, leading to a reduction of approximately 80\% fewer deaths within this country. The other country then experiences a substantial increase in deaths of about 23\%. But since the decrease in the number of cases in the first country is substantially lower, this leads to an overall decrease in the number of deaths. However, policymakers of the second country would not be willing to agree to this policy due to the higher case numbers. We find that imposing additional Pareto constraints yield an overall improvement by 5\% in comparison to the current strategy but an overall deterioration in comparison to the unrestricted optimal strategy. The Pareto optimal strategy derived by us assigns just enough vaccine doses to the second country such that it has the same case numbers as for the baseline EU strategy. However, the first country is better off, leading to an overall improvement of the Pareto strategy in comparison to the current strategy. Since both countries are not worse off, they have no incentive to vote against the Pareto optimal strategy, making it more likely to be implemented in practice. Stochastic simulations yielded higher death numbers but underlined the general findings with respect to the efficiency of the strategies.


# Tools

We use Python and mainly its libraries [libSBML](http://sbml.org/Main_Page) (Bornstein, 2008) and [AMICI](https://amici.readthedocs.io/en/latest/) (Fröhlich, 2021) to implement our models. [pyPESTO](https://pypesto.readthedocs.io/en/latest/) is our main tool for optimization.

## To-Do:
- Finish implementation of neural nets
- sensitivity analysis (hyperbolic tangent vs logistic)
- Implement inter-spatial infections by individuals that are in another country (new model)
- Incorporate reinfections
- Incorporate tests (Can we substitute vaccination by testing? Optimal distribution of vaccines and tests)
