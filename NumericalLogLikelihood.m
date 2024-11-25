function [LogLikelihood] = NumericalLogLikelihood(data_t, data_x, a, D, x0, L, method)
% Function to calculate the log likelihood when the probability is given by the
% numerical solution to the partial differential equation 
% dy/dt = -a dy/dx + D/2 d^2y/dx^2, (y=y(x,t))
% subject to a delta function initial condition: y(x,0) = delta(x - x0),
% and a Dirichlet boundaries at x=0: y(0,t) = 0.
% The type of boundary at the upper boundary, x=L, and the method of solution
% is set by the method argument which corresponds to the methods in 
% NumericalPDESolve.m
data_count = length(data_t);
t_counts = ceil(100 * max(data_t));
x_counts = 200;
t_points = linspace(0, max(data_t), t_counts);
[~, num_x, num_y] = NumericalPDESolve(a, D, x0, 0, L, t_points, x_counts, false, method, false);

% Option to visualise the numerical solution to ensure it looks correct.
% Good for test cases before running optimisations with this function.
% ResultsBar(t_points(1:10000:t_counts), num_x, num_y(1:10000:t_counts ,:), false, 0.000001);

t_indices = zeros(1, data_count);
x_indices = zeros(size(data_x));
LogLikelihood = 0;
for i = 1:data_count
    t_indices(i) = find(abs(data_t(i) - t_points) == min(abs(data_t(i) - t_points)), 1);
    for j = 1:length(data_x(i, ~isnan(data_x(i, :))))
        x_indices(i, j) = find(abs(data_x(i, j) - num_x(2:x_counts-1)) == min(abs(data_x(i, j) - num_x(2:x_counts-1))), 1);
        LogLikelihood = LogLikelihood + log(num_y(t_indices(i), x_indices(i, j)) / NumericalIntegrator(num_x, num_y(t_indices(i), :)));
    end
end
end
