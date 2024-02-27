# EpiEKF

Code for the EpiEKF method (SIRS model coupled with EKF) used for the [Respicast](https://respicast.ecdc.europa.eu/) hub. It is adapted from a method used for midterm projections within the Covid-19 task force in Luxembourg ([Aalto et al., 2022](#references)).

## Instructions

 - Download files [_latest-ILI_incidence.csv_](https://github.com/european-modelling-hubs/flu-forecast-hub/blob/main/target-data/latest-ILI_incidence.csv) and [_forecasting_weeks.csv_](https://github.com/european-modelling-hubs/flu-forecast-hub/blob/main/supporting-files/forecasting_weeks.csv) from the [flu forecast hub](https://github.com/european-modelling-hubs/flu-forecast-hub) github page.
 - Run the file `SIRS_EKF_main.m`
 - Output csv file in the required format will be generated and the most recent projections plotted

### Other notes

 - To inspect past projections, simply truncate the data vector `Y` while the full data is stored in the vector `Yfull`.
 - Note that if a new country is added, an error will be generated.

## Model details

The model is a stochastic SIRS model (discrete-time with time step of 1 day). The stochastication is done by accompanying each state transition by a stochastic process. For example, in the deterministic model, $\beta(t)S(t)I(t)/N$ infections occur on day $t$. In the stochastic version, the number of new infections is normally distributed with mean $\beta(t)S(t)I(t)/N$ and variance $\beta(t)S(t)I(t)/N$. This stochastic formalisation arises from the assumption that on day $t$, each susceptible person has the probability $p=\beta(t)I(t)/N$ to become infected. Then, the number of new infections is binomially distributed. Since $p$ is rather small, the term $(1-p)$ is ignored from the variance. Moreover, with high enough number of new infections, the binomial distribution can be well approximated by the normal distribution. More details on the EKF implementation on such model can be found in ([Proverbio et al., 2022](#references), Appendices A-B).

### Adaptive hyperparameter estimation

A key parameter affecting the quality of projections is the ratio of detected and total cases. In particular, the link between the exponential growth rate at the onset of a new epidemic wave and the amplitude of the wave is strongly dependent on this parameter. With too large parameter, the projections will overshoot the wave, and conversely, with too small parameter, the projections are too optimistic. Unfortunately, one parameter rarely works for every epidemic wave in a region. Therefore we have implemented an adaptive scheme to retroactively adjust this hyperparameter. This estimation works so that if $10\beta(t)-S(t)/0.45 > 1$, that is, $\beta$-parameter becomes too large and/or the susceptible-pool becomes too small, the state estimate is deemed unrealistic. The simulation then jumps backwards in time for 12 weeks, and re-runs the last 12 weeks with an increased ratio.

### Method parameters

To be added...

## References

 * A. Aalto, S. Martina, D. Proverbio, F. Kemp, P. Wilmes, J. Goncalves, A. Skupin. “Covid-19 report: Update on the current epidemic status in Luxembourg” (30 June 2022) [Link](https://www.researchluxembourg.org/en/covid-19-task-force/publications/).

 * D. Proverbio, F. Kemp, S. Magni, L. Ogorzaly, H.-M. Cauchie, J. Goncalves, A. Skupin, A. Aalto. “Model-based assessment of COVID-19 epidemic dynamics by wastewater analysis”, _Science of the Total Environment_ __827__, 154235, (2022), [doi.org/10.1016/j.scitotenv.2022.154235](doi.org/10.1016/j.scitotenv.2022.154235).
