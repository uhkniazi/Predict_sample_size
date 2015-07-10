# Predict_sample_size
Possible sample sizes for PREDICT study calculated using Simulation

# Simulate_posterior_predictive.R
Reads data from the population, and finds possible outliers using a maximum likelihood ratio cutoff of 0.2. Using the data to calculate the hyperparameters for mean and variance, it simulates the variance from an inverse chisquare distribution, and a prior mean from the normal distribution given a sample from the prior distribution for variance. Finally it calculates a likelihood matrix for the data. Using the maximum likelihood vector for the variance and mean, it samples (similar to rejection sampling) from the 
prior distributions of variances and means (with replacement) to create a posterior distribution. Then it uses the posterior 
to simulate the data for the posterior predictive distribution. NOTE: we do not use n0 to divide the variance by n0 while 
calculating a prior mean (using the prior variance sample) as we would like to keep the variance large. Similarly in the likelihood
calculations we are also not dividing the variance by n.

We perform some power calculations, by calculating the parameters for the function using the simulated data.






