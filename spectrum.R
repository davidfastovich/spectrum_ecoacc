library(astrochron)
library(tidyverse)
library(patchwork)
library(scales)

# ##############################################################################
# ACCLIMATION AND FUNCTION MODEL
# ##############################################################################

# functions - all from Michael Stemkovski
acclimation.func <- function(time, climate, lambda, E_0, kappa=0){
  E <- rep(NA, length(time))
  E[1] <- E_0
  for(i in seq_along(time)[-1]){
    t <- time[i]
    E[i] <- E[i-1] - (1/lambda) * (E[i-1] - climate[i-1])
  }
  return(E)
}

response.func <- function(climate, ecosystem, gamma = 0.02, alpha=1000, beta=1, disequil_abs = FALSE, disequil_sqr = FALSE, form=1){
  disequil <- ecosystem - climate
  if(disequil_abs) disequil <- abs(disequil)
  if(disequil_sqr) disequil <- disequil^2
  
  if(form == 0) response <- alpha + beta*climate # null linear model - no effect of disequilibrium
  # response <- alpha + (beta - gamma*disequilibrium)*climate
  if(form == 1) response <- alpha + beta*climate - gamma*disequil*climate # disequilibrium affects slope
  if(form == 2) response <- alpha + beta*climate - gamma*disequil # disequilibrium affects intercept
  if(form == 3) response <- alpha - gamma*disequil # there is no climate sensitivity, baseline response is a result of disequilibrium - to be used with multiple acclimation processes
  if(form == 4) response <- alpha - gamma*disequil*climate # there is no direct climate sensitivity, climate sensitivity is a result of disequilibrium - to be used with multiple acclimation processes
  
  if(form == 5) response <- alpha*ecosystem - gamma*disequil # baseline response is defined by ecosystem structure
  if(form == 6) response <- alpha*ecosystem + beta*climate - gamma*disequil*climate # baseline response is result of structure, and there is baseline climate sensitivity, and disequilibrium affects climate sensitivity
  
  if(form == 7) response <- alpha + (beta*climate)/(gamma*(1+disequil)) # climate sensitivity is reduced by disequilibrium. +1 so that sensitivity doesn't approach infinity if disequil is small
  if(form == 8) response <- alpha/(gamma*(1+disequil)) + beta*climate # baselines function is reduced by disequilibrium, sensitivity is unaffected
  
  if(form == 9) response <- (alpha + beta*climate) * exp(-(disequil)/gamma) # photosynthesis efficiency acclimation from Friend (2010)
  if(form == 10) response <- alpha + (beta*climate) * exp(-(disequil)/gamma) # photosynthesis efficiency acclimation from Friend (2010)
  
  if(form == 11) response <- alpha + (beta*climate) * -exp(disequil/gamma)^2 # squares the whole exp term - to be used with disequil_sqr=F - more similar to Friend 2010
  
  return(response)
}

# ##############################################################################
# SET OUR PARAMETERS FOR THE ACCLIMATION FUNCTION, RESPONSE FUNCTION, AND INPUT
# CLIMATE DATA 
# ##############################################################################

# For acclimation function
lambda <- 2

# For ecosystem response function 
form <- 1
response_func_beta <- 10
alpha <- 100
gamma <- 0.02

# Incput climate data
input_clim_beta <- 2

# ##############################################################################
# SIMULATE CLIMATE WITH KNOWN SPECTRAL PROPERTIES AND CALCULATE ITS POWER 
# SPECTRUM
# ##############################################################################

# Generate the input climate data with that returns a power spectrum with a
# beta of 2
input_climate <- pwrLaw(
  npts = 1000,
  mean = 20,
  sdev = 1,
  beta = 2,
  genplot = FALSE,
  verbose = FALSE # Set this to TRUE if you would like to see statistics about the input climate data
)

clim_spectra <- mtmPL(
  input_climate,
  verbose = FALSE,
  genplot = FALSE,
  output = 1
)

# Need to normalize to account for a fourier transform extending from positive
# to negative infinity which would produce infinite variance.
clim_spectra$Power <- clim_spectra$Power * (diff(input_climate[,1])[1] * (dim(input_climate)[1])^3)

# Some data wrangling to plot data nicely using ggplot2
clim_spectra$id <- "Climate"
input_climate$id <- "Climate"
input_climate <- rename(input_climate, "value" = "V2")

# ##############################################################################
# CALCULATE ECOSYSTEM ACCLIMATION AND ITS POWER SPECTRUM
# ##############################################################################

acclimation <- acclimation.func(
  time = input_climate[, 1],
  climate = input_climate[, 2],
  lambda = lambda,
  E_0 = mean(input_climate[1:100, 2]),
  kappa = 0
)

acclimation_spectra <- mtmPL(
  dat = data.frame(
    time = input_climate[,1],
    equil = acclimation
  ),
  verbose = FALSE,
  genplot = FALSE,
  output = 1
)

# Need to normalize to account for a fourier transform extending from positive
# to negative infinity which would produce infinite variance.
acclimation_spectra$Power <- acclimation_spectra$Power * (diff(input_climate[,1])[1] * (dim(input_climate)[1])^3)

# Some data wrangling to plot data nicely using ggplot2
acclimation_spectra$id <- "Ecosystem Acclimation"

# Data frame that holds the acclimation time series results
acclimation_timeseries <- data.frame(
  time = input_climate[, 1],
  value = acclimation,
  id = "Ecosystem Acclimation"
)

# ##############################################################################
# CALCULATE ECOSYSTEM FUNCTION AND ITS POWER SPECTRUM
# ##############################################################################

ecosystem_function <- response.func(
  input_climate[,2],
  acclimation,
  gamma = gamma,
  alpha = alpha,
  beta = response_func_beta,
  disequil_sqr = TRUE,
  form = form
)

ecosystem_function_spectra <- mtmPL(
  dat = data.frame(
    time = input_climate[,1],
    equil = ecosystem_function
  ),
  verbose = FALSE,
  genplot = FALSE,
  output = 1
)

# Need to normalize to account for a fourier transform extending from positive
# to negative infinity which would produce infinite variance.
ecosystem_function_spectra$Power <- ecosystem_function_spectra$Power * (diff(input_climate[,1])[1] * (dim(input_climate)[1])^3)

# Some data wrangling to plot data nicely using ggplot2
ecosystem_function_spectra$id <- "Ecosystem Function"

# Data frame that holds the acclimation time series results
ecosystem_function_timeseries <- data.frame(
  time = input_climate[, 1],
  value = ecosystem_function,
  id = "Ecosystem Function"
)

# ##############################################################################
# DATA WRANGLING TO GET EVERYITHING TO PLOT TOGETHER
# ##############################################################################

result_time_series <- rbind(
  input_climate,
  acclimation_timeseries,
  ecosystem_function_timeseries
)

result_power_spectra <- rbind(
  clim_spectra,
  acclimation_spectra,
  ecosystem_function_spectra
)

# ##############################################################################
# PLOT
# ##############################################################################

# Results of the models

timeseries_plot <- result_time_series %>% 
  ggplot(aes(x = time, y = value, color = id)) + 
  geom_line(linewidth = 0.25) + 
  xlab("Time") + 
  ylab("Temperature and Community Niche") + 
  facet_wrap(~ id, scales = "free_y", nrow = 3) + 
  theme_bw()

spectra_plot <- result_power_spectra %>%
  ggplot(aes(x = Frequency, y = Power, color = id)) +
  geom_line(linewidth = 0.25) + 
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), breaks = 10^seq(-2.5, 0.25, 0.25)) + 
  geom_smooth(method = 'lm', formula = y ~ x) + 
  xlab(expression(paste("Frequency (cycles ", yr^{-1}, ")"))) +
  ylab("Power Spectral Density") + 
  theme_bw()

combined_plot <- timeseries_plot + spectra_plot +
  plot_layout(guides = 'collect')

combined_plot
