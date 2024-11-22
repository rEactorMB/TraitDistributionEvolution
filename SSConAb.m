% Old method
% function [solution_value] = SSConAb(x, a, D, L, x0, t)
% if isscalar(t)
%     k_max = 100;
% 
%     outside_sum = ((2 / L) * exp((a * (x - x0) / D) - ((a * a) * t / (2 * D))));
%     k = (1:k_max).';
%     k_sum = sum(sin(k * x0 * pi / L) .* sin(k .* x * pi / L) .* exp(-(k .* k * pi * pi) * t * D / (2 * L * L)));
% 
%     temp_x = linspace(0, L, 100);
%     temp_outside_sum = ((2 / L) * exp((a * (temp_x - x0) / D) - ((a * a) * t / (2 * D))));
%     temp_k_sum = sum(sin(k * x0 * pi / L) .* sin(k .* temp_x * pi / L) .* exp(-(k .* k * pi * pi) * t * D / (2 * L * L)));
%     temp_y = temp_outside_sum .* temp_k_sum;
%     temp_y = temp_y / min(temp_y(temp_y ~= 0));
% 
%     integral_value = sum((2 / L) * exp(-(a * x0 / D) - ((a * a) * t / (2 * D))) * sum(sin(k * x0 * pi / L) .* exp(-(k .* k * pi * pi) * t * D / (2 * L * L)) .* (pi * D * D * k * L .* (1 - (((-1).^k) * exp(a * L / D))) ./ ((D * D * pi * pi * k .* k) + (a * a * L * L)))));
%     disp(integral_value);
%     solution_value = (outside_sum .* k_sum) / NumericalIntegrator(temp_x, temp_y);
% else
%     solution_value = zeros(length(t), length(x));
%     for i = 1:length(t)
%         solution_value(i, :) = SSConAb(x, a, D, L, x0, t(i));
%     end
% end
% end

function [solution_value] = SSConAb(x, a, D, L, x0, t)
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