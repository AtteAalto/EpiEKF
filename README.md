# EpiEKF

Code for the EpiEKF method (SIRS model coupled with EKF) used for the [Respicast](https://respicast.ecdc.europa.eu/) hub.

## Instructions

 - Download files [_latest-ILI_incidence.csv_](https://github.com/european-modelling-hubs/flu-forecast-hub/blob/main/target-data/latest-ILI_incidence.csv) and [_forecasting_weeks.csv_](https://github.com/european-modelling-hubs/flu-forecast-hub/blob/main/supporting-files/forecasting_weeks.csv) from the [flu forecast hub](https://github.com/european-modelling-hubs/flu-forecast-hub) github page.
 - Run the file `SIRS_EKF_main.m`
 - Output csv file in the required format will be generated and the most recent projections plotted

## Other notes

 - To inspect past projections, simply truncate the data vector `Y` while the full data is stored in the vector `Yfull`.
 - Note that if a new country is added, an error will be generated.
