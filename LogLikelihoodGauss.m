function log_sum = LogLikelihoodGauss(trait_values, mean, sigma)
log_sum = sum(log((1 / (sigma * sqrt(2 * pi))) * exp(-((trait_values - mean) .* (trait_values - mean)) / (2 * sigma * sigma))));
end