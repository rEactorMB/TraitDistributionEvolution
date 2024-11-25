function [extract_times, x_values, y_values] = NumericalPDESolve(a, D, x0, lower_boundary, upper_boundary, extract_times, traitspace_splits, wtf, method, large_time)
% Various methods to numerically solve the differential equation
% dy/dt = -a dy/dx + D/2 d^2y/dx^2, (y = y(x,t)).
% The initial state is assumed to be a single unit in a single bin
% (delta function like initial condtion).  The lower boundary is
% always a Dirichlet boundary: y(0,t) = 0.
% The upper boundary can either be a Dirichlet boundary: y(L,t) = 0,
% a Neumann boundary: dy/dx|x=L = 0,
% or a no flux (Robin type) boundary: dy/dx|x=L = 2a/D y(L,t)

% Forward Euler method with Dirichlet boundary at x=L
if method == "forwardEulerUpperA"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = linspace(lower_boundary, upper_boundary, traitspace_splits);
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');

    y_values = zeros(length(extract_times), length(x_values));
    y_values(1, ix0) = 1 / dx;
    for k = 1:(length(extract_times) - 1)
        for idx = 1:length(x_values)
            if idx == 1
                y_values(k+1, idx) = 0;
            elseif idx == length(x_values)
                y_values(k+1, idx) = 0;
            else
                y_values(k+1, idx) = (((-(a * dt) / dx) * (y_values(k, idx + 1) - y_values(k, idx))) + (((D * dt) / (2 * dx * dx)) * (y_values(k, idx+1) - 2 * y_values(k, idx) + y_values(k, idx-1))) + y_values(k, idx));
            end
        end
    end

% Crank Nicolson method with a Dirichlet boundary at x=L
elseif method == "CrankNicolsonUpperA"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = linspace(lower_boundary, upper_boundary, traitspace_splits);
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');
    if ix0 == 1 || ix0 == length(x_values)
        y_values = zeros(length(extract_times), length(x_values));
    else
        y_values = zeros(length(extract_times), length(x_values));
        y_values(1, ix0) = 1 / dx;
    end

    for k = 1:(length(extract_times) - 1)
        A = zeros(length(x_values));
        B = zeros(length(x_values));
        for idx = 1:length(x_values)
            if idx == 1
                A(idx, idx) = 1;
                B(idx, idx) = 1;
            elseif idx == length(x_values)
                A(idx, idx) = 1;
                B(idx, idx) = 1;
            else
                f1 = (a * dt) / (2 * dx);
                f2 = (D * dt) / (4 * dx * dx);
                A(idx, idx) = 1 - f1 + (2 * f2);
                A(idx, idx + 1) = f1 - f2;
                A(idx, idx - 1) = -f2;
                B(idx, idx) = 1 + f1 - (2 * f2);
                B(idx, idx + 1) = f2 - f1;
                B(idx, idx - 1) = f2;
            end
        end
        b = B * y_values(k, :).';
        y_values(k+1, :) = A \ b;
    end

% Crank Nicolson method with a Neumann boundary at x=L
elseif method == "CrankNicolsonUpperRneu"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = linspace(lower_boundary, upper_boundary, traitspace_splits);
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');
    y_values = zeros(length(extract_times), length(x_values));
    y_values(1, ix0) = 1 / dx;

    for k = 1:(length(extract_times) - 1)
        A = zeros(length(x_values));
        B = zeros(length(x_values));
        for idx = 1:length(x_values)
            if idx == 1
                A(idx, idx) = 1;
                B(idx, idx) = 1;
            elseif idx == length(x_values)
                f2 = (D * dt) / (4 * dx * dx);
                A(idx, idx) = 1 + f2;
                A(idx, idx - 1) = -f2;
                B(idx, idx) = 1 - f2;
                B(idx, idx - 1) = f2;
            else
                f1 = (a * dt) / (2 * dx);
                f2 = (D * dt) / (4 * dx * dx);
                A(idx, idx) = 1 - f1 + (2 * f2);
                A(idx, idx + 1) = f1 - f2;
                A(idx, idx - 1) = -f2;
                B(idx, idx) = 1 + f1 - (2 * f2);
                B(idx, idx + 1) = f2 - f1;
                B(idx, idx - 1) = f2;
            end
        end
        b = B * y_values(k, :).';
        y_values(k+1, :) = A \ b;
    end

% Crank Nicolson method with a no flux boundary at x=L
elseif method == "CrankNicolsonUpperRrob"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = linspace(lower_boundary, upper_boundary, traitspace_splits);
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');
    y_values = zeros(length(extract_times), length(x_values));
    y_values(1, ix0) = 1 / dx;

    for k = 1:(length(extract_times) - 1)
        grad_vec = zeros(1, length(x_values));
        upper_G = (2 * a / D) * y_values(k, end);
        A = zeros(length(x_values));
        B = zeros(length(x_values));
        for idx = 1:length(x_values)
            if idx == 1
                A(idx, idx) = 1;
                B(idx, idx) = 1;
            elseif idx == length(x_values)
                f1 = (a * dt) / (2 * dx);
                f2 = (D * dt) / (4 * dx * dx);
                A(idx, idx) = 1 + f2;
                A(idx, idx - 1) = -f2;
                B(idx, idx) = 1 - f2;
                B(idx, idx - 1) = f2;
                grad_vec(idx) = (f2 - f1) * 2 * dx * upper_G;
            else
                f1 = (a * dt) / (2 * dx);
                f2 = (D * dt) / (4 * dx * dx);
                A(idx, idx) = 1 - f1 + (2 * f2);
                A(idx, idx + 1) = f1 - f2;
                A(idx, idx - 1) = -f2;
                B(idx, idx) = 1 + f1 - (2 * f2);
                B(idx, idx + 1) = f2 - f1;
                B(idx, idx - 1) = f2;
            end
        end
        b = B * y_values(k, :).' + grad_vec.';
        y_values(k+1, :) = (A \ b).';
    end

% Throw an error if the method is not one of those above
else
    error('Error : method "' + method + '" not supported.');
end

end
