function [x_values, y_values] = KernelDensityFunction(bins, fill_levels, normalise)
% Function to smooth over a histogram by making every point in the histogram
% a Gaussian curve and finding the sum of these curves.
sigma_0 = bins(end) - bins(1);
sigma_min = sigma_0 / (length(bins) * 2);
x_values = linspace(bins(1), bins(end), 1000);
y_values = zeros(1, 1000);
Gaussian = @(mean, var) ((1 / sqrt(2 * pi * var)) * exp(-(x_values - mean).^2 / (2 * var)));
for i = 1:length(bins)
    y_values = y_values + Gaussian(bins(i), sigma_min) * fill_levels(i);
end
    
if normalise == true
    y_values = y_values / NumericalIntegrator(x_values, y_values);
end
end
