function [solution_value] = SSConAb(x, a, D, L, x0, t)
% Function to return the value of the solution to the partial
% differential equation dy/dt = -a dy/dx + D/2 d^2y/dx^2 (y = y(x,t))
% given a delta function initial condition: y(x,0) = delta(x - x0)
% and Dirichlet boundaries at x=0,L: y(0,t) = y(L,t) = 0.
% The code works for single or multiple values of time (t).
% Since the analytical solution is an infinite series solution, this
% calculation is truncated at k_max terms.
% The solution is normalised so that the area under the solution
% curve is one.
if isscalar(t)
    k_max = 100;
    
    k = (1:k_max).';
    % k_sum = sum(exp(a * x / D) .* sin(k * x0 * pi / L) .* sin(k .* x * pi / L) .* exp(-(k .* k * pi * pi) * t * D / (2 * L * L)));
    k_cumsum = cumsum(exp(a * x / D) .* sin(k * x0 * pi / L) .* sin(k .* x * pi / L) .* exp(-(k .* k * pi * pi) * t * D / (2 * L * L)));
    k_sum = zeros(1, length(x));
    for i = 1:length(x)
        try
            k_sum(i) = k_cumsum(find(k_cumsum(:, i), 1, 'last'), i);
        catch
        end
    end
    k_sum(logical((x <= 0) + (x >= L))) = 0;

    integral_value = sum(sum(sin(k * x0 * pi / L) .* exp(-(k .* k * pi * pi) * t * D / (2 * L * L)) .* (pi * D * D * k * L .* (1 - (((-1).^k) * exp(a * L / D))) ./ ((D * D * pi * pi * k .* k) + (a * a * L * L)))));

    solution_value = k_sum / integral_value;
else
    solution_value = zeros(length(t), length(x));
    for i = 1:length(t)
        solution_value(i, :) = SSConAb(x, a, D, L, x0, t(i));
    end
end
end
