function [integral_value] = NumericalIntegrator(x_values, y_values)
integral_value = 0;
for i = 1:(length(x_values) - 1)
    integral_value = integral_value + (abs((y_values(i+1) + y_values(i)) / 2) * (x_values(i+1) - x_values(i)));
end
end