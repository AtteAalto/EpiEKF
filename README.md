# EpiEKF

__Atte Aalto__<sup>1</sup>, __Daniele Proverbio__<sup>2</sup>, __Giulia Giordano__<sup>2</sup>, and __Jorge Goncalves__<sup>1</sup>

  1. Luxembourg Centre for Systems Biomedicine, University of Luxembourg  
  2. University of Trento

#

This is the code for the EpiEKF method (SIRS model coupled with EKF) used for ILI incidence projections for the [Respicast](https://respicast.ecdc.europa.eu/) hub. It is adapted from a method used for midterm projections within the Covid-19 task force in Luxembourg ([Aalto et al., 2022](#references)).


## Instructions

 - Download files [_latest-ILI_incidence.csv_](https://github.com/european-modelling-hubs/flu-forecast-hub/blob/main/target-data/latest-ILI_incidence.csv) and [_forecasting_weeks.csv_](https://github.com/european-modelling-hubs/flu-forecast-hub/blob/main/supporting-files/forecasting_weeks.csv) from the [flu forecast hub](https://github.com/european-modelling-hubs/flu-forecast-hub) github page.
 - Run the file `SIRS_EKF_main.m`
 - Output csv file in the required format will be generated and the most recent projections plotted

### Other notes

 - To inspect past projections, simply truncate the data vector `Y` while the full data is stored in the vector `Yfull`.
 - Note that if a new country is added, an error will be generated.

## Model details

The model is a stochastic SIRS model (discrete-time with time step of 1 day) 

$S(t+1) = S(t) - \frac{\beta(t)I(t)}{N}S(t) + \varphi R(t) - \sqrt{\frac{\beta(t)I(t)}{N}S(t)}w_1(t) + \sqrt{\varphi R(t)} w_3(t)$

$I(t+1) = I(t) + \frac{\beta(t)I(t)}{N}S(t) - \mu I(t) + \sqrt{\frac{\beta(t)I(t)}{N}S(t)}w_1(t) - \sqrt{\mu I(t)} w_2(t)$

$R(t+1) = R(t) + \mu I(t) - \varphi R(t) + \sqrt{\mu I(t)} w_2(t) - \sqrt{\varphi R(t)} w_3(t)$

where $w_1$, $w_2$ and $w_3$ are mutually independent discrete-time white noise processes with variance 1. The stochastic model is derived by accompanying each state transition by a stochastic process. For example, in the deterministic model, $\beta(t)S(t)I(t)/N$ infections occur on day $t$. In the stochastic version, the number of new infections is normally distributed with mean $\beta(t)S(t)I(t)/N$ and variance $\beta(t)S(t)I(t)/N$. This stochastic formalisation arises from the assumption that on day $t$, each susceptible person has the probability $p(t)=\beta(t)I(t)/N$ to become infected. Then, the number of new infections is binomially distributed with mean $S(t)p(t)$ and variance $S(t)p(t)(1-p(t))$. Since $p(t)$ is rather small, the term $(1-p(t))$ is omitted from the variance. Moreover, with high enough number of new infections, the binomial distribution can be well approximated by the normal distribution. 

The weekly number of detected cases according to the deterministic part of the model (which is used for the Kalman filtering) is
$y(t) = C(t) \sum_{\tau = t-6,...,t} \frac{\beta(t)I(t)}{N}S(t)$
where $C(t)$ is the time-varying (as explained below) ratio of detected and total infections.

### State estimation and future projections

Details on the EKF implementation on such a model can be found in ([Proverbio et al., 2022](#references), Appendices A-B). Every week, the deterministic part of the model is simulated forward in time for one week. The model-predicted number of cases is then compared to the real data, and the difference is used to adjust the model's state after multiplication by the so-called Kalman gain. The stochastic part of the model is used to determine the gain. To generate the future projections, the stochastic model is simulated forward in time to generate 1000 stochastic trajectories. The medians and required quantiles are calculated from the replicates. The EKF returns also an error covariance matrix for the model state. The initial state for the stochastic replicates is drawn from the normal distribution whose mean is the current state estimate and the covariance is the error covariance from the EKF.

### Adaptive hyperparameter estimation

A key parameter affecting the quality of projections is the ratio of detected and total cases $C(t)$. In particular, the link between the exponential growth rate at the onset of a new epidemic wave and the amplitude of the wave is strongly dependent on this parameter. With too large parameter, the projections will overshoot the wave, and conversely, with too small parameter, the projections are too optimistic. Unfortunately, one parameter rarely works for every epidemic wave in a region. Therefore we have implemented an adaptive scheme to retroactively adjust this hyperparameter. This estimation works so that if $10\beta(t)-S(t)/(0.45N) > 1$, that is, if $\beta$-parameter becomes too large and/or the susceptible-pool becomes too small, then the state estimate is deemed unrealistic. The simulation then jumps backwards in time for 12 weeks, and re-runs the last 12 weeks with $C(t)$ increased by 20% (in a piecewise linear way). 26 weeks after the adjustment was made, $C(t)$ has again returned to the baseline value. 

### Method parameters

Here we list the parameters of the model, and region-dependent hyperparameters. In this section, we denote by $Y_j$ the vector of length $m$ containing the weekly case numbers four country $j$. In the method, there are three global parameters $a$, $b$, and $c$ that modulate local parameters that depend on the historical data. These three parameters are fitted by optimising over all 24 regions available at the end of Jan 2024. The cost function for the optimisation is

$J(a,b,c) = \sum_{j=1,...,24} \sum_{t=1,...,m-4} \sum_{\tau = 1,2,3,4} \frac{|\sqrt{Y_j(t+\tau)} - \sqrt{\hat Y_j(t+\tau|t)}|}{\sqrt{Y_j(t+\tau)+1}}$

where $\tau$ is the projection horizon, that is, $\hat Y_j(t+\tau|t)$ stands for the modelled case numbers for day $t+\tau$ using data until day $t$. The square root is used as a variance-stabilising transformation. Without this transformation, too much emphasis is put on the overshooting tendency of the projections, which results in rather optimistic projections for the majority of time.

 - Rate for transition $I \to R$, $\mu = 0.06$. It is unrealistically small, but to guarantee stable state estimation, the time scale of the model dynamics should not be too far from the time scale of data collection (one week).

 - Rate for transition $R \to S$,  $\varphi = log(2)/60$, corresponding to immunity half life of 60 days. Influenza-like-illnesses may be caused by many different viruses, with varying cross-immunity properties. This parameter should be considered as an effective parameter rather than a reflection of the human immune system.

- The rate for transition $S \to I$, denoted by $\beta(t)$ is time-varying, and it is also part of the system’s state estimated by EKF. The stochastic model for it is $\beta(t+1) = \beta(t) + w_{\beta}(t)$ where $w_{\beta}(t)$ is normally distributed with mean 0 and variance $0.012^2$. Here the number 0.012 is chosen as one fifth of the rate I to R.

 - The stochastic component in the model is rather small. To account for uncertainty due to modelling errors, the covariance matrix for the state noise in the EKF is multiplied by the tuning parameter $a$.

 - The baseline for the ratio of detected and total cases is $b\frac{52 \sum_{t=1,...,m} Y(t)}{mN}$, that is, number of detected cases per capita per year multiplied by the tuning parameter $b$.

- The measurement noise variance for region $j$ on day $t$ is $c K_j Y_j(t)$ where $c$ is the global tuning parameter, and $K_j$ is obtained by $K_j = 1/m \sum_{t=1,...,m} (Y(t) - Y_s(t))^2 / Y_s(t)$ where $Y_s$ is a moving window average of $Y$ (over $t-2,…,t+2$).



## References

 * A. Aalto, S. Martina, D. Proverbio, F. Kemp, P. Wilmes, J. Goncalves, A. Skupin. “Covid-19 report: Update on the current epidemic status in Luxembourg” (30 June 2022) [Link](https://www.researchluxembourg.org/en/covid-19-task-force/publications/).

 * D. Proverbio, F. Kemp, S. Magni, L. Ogorzaly, H.-M. Cauchie, J. Goncalves, A. Skupin, A. Aalto. “Model-based assessment of COVID-19 epidemic dynamics by wastewater analysis”, _Science of the Total Environment_ __827__, 154235, (2022), [doi.org/10.1016/j.scitotenv.2022.154235](https://doi.org/10.1016/j.scitotenv.2022.154235).
