This project utilizes maximum likelihood (ML) estimators to predict the
value of the alpha and lambda parameters in the rayleigh and exponential
distributions, respectively.  In part 1, data is generated from randomly
drawing samples from ground truth distributions. The formula to find the
ML estimator for each distribution was derived by hand before applying it
to the data.  Once the estimator was calculated, it's corresponding mean
squared error (MSE), bias, and variance was plotted for differing numbers of
observations.  Note that for each number of observations 100000 trials
was completed. In part 2, data was provided from either a rayleigh or
exponential distribution.  Using the calculations in part 1, the
parameter was estimated for each potential distribution. The log
likelihood of these corresponding distributions was then found.  Since
the rayleigh distribution had a greater log likelihood, the data
most likely came from sampling a rayleigh pdf.
