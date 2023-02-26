# COVID-19-data-analysis
The aim of our work is to (1) explore the lagged dependence between the time series of case counts and the time series of death counts; and (2) utilize such a relationship for prediction. Time-varying coefficient models are considered: local polynomial regression (death_predict_Github_localpoly.R) and piece-wise linear regression (death_predict_Github_pwl.R).

## An example for using the two R functions
<sub>
library(tidyverse)

# Data in Canada
# https://github.com/ccodwg/Covid19Canada
# https://github.com/ccodwg/CovidTimelineCanada (new data since May 2022)

urldeath <- "https://raw.githubusercontent.com/ccodwg/Covid19Canada/master/timeseries_prov/mortality_timeseries_prov.csv"
death = read.csv(urldeath)
death$Date <- as.Date(death$date_death_report, format="%d-%m-%Y")


urlcase <- 'https://raw.githubusercontent.com/ccodwg/Covid19Canada/master/timeseries_prov/cases_timeseries_prov.csv'
cases = read.csv(urlcase)
cases$Date <- as.Date (cases$date_report, format="%d-%m-%Y")

prov_name <-          c("Alberta", "BC", "Manitoba", "Ontario", "New Brunswick", 
                        "NL", "Nova Scotia", "PEI", "Quebec", "Saskatchewan")
omicron_wave_start <- c("2021-10-05", "2021-12-05", "2021-12-05", "2021-10-31", "2021-11-28", "2021-12-05", "2021-12-05", "2021-11-28", "2021-10-31", "2021-12-05")
data_end_date <- c("2023-01-01", "2022-03-20", "2022-03-26", "2023-01-01", "2022-03-13", "2023-01-01", "2022-03-06", "2022-04-06", "2023-01-01", "2022-02-08")

omicron_waves <- data.frame(prov_name, omicron_wave_start, data_end_date)
omicron_waves$omicron_wave_start <- as.Date(omicron_waves$omicron_wave_start)
omicron_waves$data_end_date <- as.Date(omicron_waves$data_end_date)

region <- "Ontario"

cutoff_start <- omicron_waves$omicron_wave_start[prov_name == region]
cutoff_end <- omicron_waves$data_end_date[prov_name == region]

death_long <- death %>% filter (Date >= cutoff_start & province %in% region & Date <= cutoff_end)
cases_long <- cases %>% filter (Date >= cutoff_start & province %in% region & Date <= cutoff_end)

df <-  merge (death_long[,c( "deaths", "cumulative_deaths","Date" )],
              cases_long[,c( "cases", "cumulative_cases", "Date" ) ],
              by = "Date")

#### combining all the output from (1) pwl vs localpoly (2) forecast vs predict (3) ll vs lc

## Ontario
pwl_ON <- pwl_pred(df, bp_index = c(70,85,105))
ll_forecast <- localpoly(df)
lc_forecast <- localpoly(df, method = "lc")
</sub>
