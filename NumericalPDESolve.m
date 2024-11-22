function [extract_times, x_values, y_values] = NumericalPDESolve(a, D, x0, lower_boundary, upper_boundary, extract_times, traitspace_splits, wtf, method, large_time)
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
elseif method == "forwardEulerUpperA_var_a"
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
                y_values(k+1, idx) = (((-(a(idx) * dt) / dx) * (y_values(k, idx + 1) - y_values(k, idx))) + (((D * dt) / (2 * dx * dx)) * (y_values(k, idx+1) - 2 * y_values(k, idx) + y_values(k, idx-1))) + y_values(k, idx));
            end
        end
    end
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
elseif method == "DiscreteODEUpperA"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = zeros(1, traitspace_splits);
    for i = 1:traitspace_splits
        if i == 1
            x_values(i) = lower_boundary + (0.5 * dx);
        else
            x_values(i) = x_values(i-1) + dx;
        end
    end
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');
    y_values = zeros(length(extract_times), length(x_values));
    y_values(1, ix0) = 1;

    D_mat = zeros(length(x_values));
    a_mat = zeros(length(x_values));
    D = D / (dx * dx);
    a = a / dx;
    for i = 1:length(x_values)
        if i == 1
            D_mat(i, i+1) = 0.5 * D;
            D_mat(i, i) = -D;
            a_mat(i, i+1) = abs(a);
            a_mat(i, i) = -abs(a);
        elseif i == length(x_values)
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -D;
            a_mat(i, i) = -abs(a);
        else
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -D;
            D_mat(i, i+1) = 0.5 * D;
            a_mat(i, i+1) = abs(a);
            a_mat(i, i) = -abs(a);
        end
    end
    if a < 0
        a_mat = a_mat.';
    end
    for i = 1:(length(extract_times) - 1)
        y_values(i+1, :) = y_values(i, :) + (dt * (y_values(i, :) * (D_mat + a_mat)));
    end
    if large_time == true
        [~, eig_val, eig_vecL] = eig(a_mat + D_mat);
        eig_vals = zeros(1, length(x_values));
        for i = 1:length(x_values)
            eig_vals(i) = eig_val(i, i);
        end
        eig_maxi = eig_vals == max(eig_vals);
        plot(x_values, eig_vecL(:, eig_maxi));
        title("(greatest) eigenvalue = " + eig_vals(eig_maxi));
        pause(2);

        [~, eig_val, eig_vecL] = eig(eye(length(x_values)) + dt * (a_mat + D_mat));
        eig_vals = zeros(1, length(x_values));
        for i = 1:length(x_values)
            eig_vals(i) = eig_val(i, i);
        end
        eig_onei = abs(1 - eig_vals) == min(abs(1 - eig_vals));
        plot(x_values, eig_vecL(:, eig_onei));
        title("eigenvector = s = " + eig_vals(eig_onei));
    end
elseif method == "DiscreteODEUpperA_var_a"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = zeros(1, traitspace_splits);
    for i = 1:traitspace_splits
        if i == 1
            x_values(i) = lower_boundary + (0.5 * dx);
        else
            x_values(i) = x_values(i-1) + dx;
        end
    end
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');
    y_values = zeros(length(extract_times), length(x_values));
    y_values(1, ix0) = 1;

    D_mat = zeros(length(x_values));
    a_mat = zeros(length(x_values));
    D = D / (dx * dx);
    a = a / dx;
    for i = 1:length(x_values)
        if i == 1
            D_mat(i, i+1) = 0.5 * D;
            D_mat(i, i) = -D;
            a_mat(i, i+1) = a(i);
            a_mat(i, i) = -a(i);
        elseif i == length(x_values)
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -D;
            a_mat(i, i) = -abs(a(i-1));
        else
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -D;
            D_mat(i, i+1) = 0.5 * D;
            a_mat(i, i+1) = a(i);
            a_mat(i, i) = -a(i);
        end
    end
    a_T = a_mat.';
    a_mat(101:200, 101:200) = a_T(101:200, 101:200);
    % a_mat(101, 100) = a(1);
    a_mat(100, 101) = 0;
    for i = 1:(length(extract_times) - 1)
        y_values(i+1, :) = y_values(i, :) + (dt * (y_values(i, :) * (D_mat + a_mat)));
    end
elseif method == "DiscreteODEUpperA_var_aD"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = zeros(1, traitspace_splits);
    for i = 1:traitspace_splits
        if i == 1
            x_values(i) = lower_boundary + (0.5 * dx);
        else
            x_values(i) = x_values(i-1) + dx;
        end
    end
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');
    y_values = zeros(length(extract_times), length(x_values));
    y_values(1, ix0) = 1;

    D_mat = zeros(length(x_values));
    a_mat = zeros(length(x_values));
    D = D / (dx * dx);
    a = a / dx;
    for i = 1:length(x_values)
        if i == 1
            D_mat(i, i+1) = 0.5 * D(i);
            D_mat(i, i) = -D(i);
            a_mat(i, i+1) = a(i);
            a_mat(i, i) = -a(i);
        elseif i == length(x_values)
            D_mat(i, i-1) = 0.5 * D(i);
            D_mat(i, i) = -D(i);
            a_mat(i, i) = -abs(a(i-1));
        else
            D_mat(i, i-1) = 0.5 * D(i);
            D_mat(i, i) = -D(i);
            D_mat(i, i+1) = 0.5 * D(i);
            a_mat(i, i+1) = a(i);
            a_mat(i, i) = -a(i);
        end
    end
    a_T = a_mat.';
    % a_mat(101:200, 101:200) = a_T(101:200, 101:200);
    % % a_mat(101, 100) = a(1);
    % a_mat(100, 101) = 0;
    for i = 1:(length(extract_times) - 1)
        y_values(i+1, :) = y_values(i, :) + (dt * (y_values(i, :) * (D_mat + a_mat)));
    end
elseif method == "DiscreteODEUpperRrob"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = zeros(1, traitspace_splits);
    for i = 1:traitspace_splits
        if i == 1
            x_values(i) = lower_boundary + (0.5 * dx);
        else
            x_values(i) = x_values(i-1) + dx;
        end
    end
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');
    y_values = zeros(length(extract_times), length(x_values));
    y_values(1, ix0) = 1;

    D_mat = zeros(length(x_values));
    a_mat = zeros(length(x_values));
    D = D / (dx * dx);
    a = a / dx;
    for i = 1:length(x_values)
        if i == 1
            D_mat(i, i+1) = 0.5 * D;
            D_mat(i, i) = -D;
            a_mat(i, i+1) = abs(a);
            a_mat(i, i) = -abs(a);
        elseif i == length(x_values)
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = - 0.5 * D;
        else
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -D;
            D_mat(i, i+1) = 0.5 * D;
            a_mat(i, i+1) = abs(a);
            a_mat(i, i) = -abs(a);
        end
    end
    if a < 0
        a_mat = a_mat.';
        a_mat(length(x_values), length(x_values)) = -abs(a);
    end
    for i = 1:(length(extract_times) - 1)
        y_values(i+1, :) = y_values(i, :) + (dt * (y_values(i, :) * (D_mat + a_mat)));
    end
    if large_time == true
        [~, eig_val, eig_vecL] = eig(a_mat + D_mat);
        eig_vals = zeros(1, length(x_values));
        for i = 1:length(x_values)
            eig_vals(i) = eig_val(i, i);
        end
        eig_maxi = eig_vals == max(eig_vals);
        plot(x_values, eig_vecL(:, eig_maxi));
        title("(greatest) eigenvalue = " + eig_vals(eig_maxi));
        pause(2);

        [~, eig_val, eig_vecL] = eig(eye(length(x_values)) + dt * (a_mat + D_mat));
        eig_vals = zeros(1, length(x_values));
        for i = 1:length(x_values)
            eig_vals(i) = eig_val(i, i);
        end
        eig_onei = abs(1 - eig_vals) == min(abs(1 - eig_vals));
        plot(x_values, eig_vecL(:, eig_onei));
        title("eigenvector = s = " + eig_vals(eig_onei));
    end
elseif method == "DiscreteODEUpperRneu"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = zeros(1, traitspace_splits);
    for i = 1:traitspace_splits
        if i == 1
            x_values(i) = lower_boundary + (0.5 * dx);
        else
            x_values(i) = x_values(i-1) + dx;
        end
    end
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');
    y_values = zeros(length(extract_times), length(x_values));
    y_values(1, ix0) = 1 / dx;

    D_mat = zeros(length(x_values));
    a_mat = zeros(length(x_values));
    D = D / (dx * dx);
    a = a / dx;
    for i = 1:length(x_values)
        if i == 1
            D_mat(i, i+1) = 0.5 * D;
            D_mat(i, i) = -D;
            a_mat(i, i+1) = abs(a);
            a_mat(i, i) = -abs(a);
        elseif i == length(x_values)
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -0.5 * D;
            a_mat(i, i) = -abs(a);
        else
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -D;
            D_mat(i, i+1) = 0.5 * D;
            a_mat(i, i+1) = abs(a);
            a_mat(i, i) = -abs(a);
        end
    end
    if a < 0
        a_mat = a_mat.';
        a_mat(length(x_values), length(x_values)) = 0;
    end
    for i = 1:(length(extract_times) - 1)
        y_values(i+1, :) = y_values(i, :) + (dt * (y_values(i, :) * (D_mat + a_mat)));
    end
    if large_time == true
        [~, eig_val, eig_vecL] = eig(a_mat + D_mat);
        eig_vals = zeros(1, length(x_values));
        for i = 1:length(x_values)
            eig_vals(i) = eig_val(i, i);
        end
        eig_maxi = eig_vals == max(eig_vals);
        plot(x_values, eig_vecL(:, eig_maxi));
        title("(greatest) eigenvalue = " + eig_vals(eig_maxi));
        pause(2);

        [~, eig_val, eig_vecL] = eig(eye(length(x_values)) + dt * (a_mat + D_mat));
        eig_vals = zeros(1, length(x_values));
        for i = 1:length(x_values)
            eig_vals(i) = eig_val(i, i);
        end
        eig_onei = abs(1 - eig_vals) == min(abs(1 - eig_vals));
        plot(x_values, eig_vecL(:, eig_onei));
        title("eigenvector = s = " + eig_vals(eig_onei));
    end
elseif method == "DiscreteODEUpperRneu2"
    dx = (upper_boundary - lower_boundary) / traitspace_splits;
    x_values = zeros(1, traitspace_splits);
    for i = 1:traitspace_splits
        if i == 1
            x_values(i) = lower_boundary + (0.5 * dx);
        else
            x_values(i) = x_values(i-1) + dx;
        end
    end
    dt = extract_times(2) - extract_times(1);
    truth_array = min(abs(x_values - x0)) == abs(x_values - x0);
    ix0 = find(truth_array, 1, 'last');
    y_values = zeros(length(extract_times), length(x_values));
    y_values(1, ix0) = 1 / dx;

    D_mat = zeros(length(x_values));
    a_mat = zeros(length(x_values));
    D = D / (dx * dx);
    a = a / dx;
    for i = 1:length(x_values)
        if i == 1
            D_mat(i, i+1) = 0.5 * D;
            D_mat(i, i) = -D;
            a_mat(i, i+1) = abs(a);
            a_mat(i, i) = -abs(a);
        elseif i == length(x_values)
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -0.5 * D;
        elseif i > (length(x_values) - 20) && i < length(x_values)
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -D;
            D_mat(i, i+1) = 0.5 * D;
        elseif i == (length(x_values) - 20)
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -D;
            D_mat(i, i+1) = 0.5 * D;
            a_mat(i, i) = -abs(a);
        else
            D_mat(i, i-1) = 0.5 * D;
            D_mat(i, i) = -D;
            D_mat(i, i+1) = 0.5 * D;
            a_mat(i, i+1) = abs(a);
            a_mat(i, i) = -abs(a);
        end
    end
    if a < 0
        a_mat = a_mat.';
        a_mat(length(x_values), length(x_values) - 1) = abs(a);
    end
    for i = 1:(length(extract_times) - 1)
        y_values(i+1, :) = y_values(i, :) + (dt * (y_values(i, :) * (D_mat + a_mat)));
    end
    if large_time == true
        [~, eig_val, eig_vecL] = eig(a_mat + D_mat);
        eig_vals = zeros(1, length(x_values));
        for i = 1:length(x_values)
            eig_vals(i) = eig_val(i, i);
        end
        eig_maxi = eig_vals == max(eig_vals);
        plot(x_values, eig_vecL(:, eig_maxi));
        title("(greatest) eigenvalue = " + eig_vals(eig_maxi));
        pause(2);

        [~, eig_val, eig_vecL] = eig(eye(length(x_values)) + dt * (a_mat + D_mat));
        eig_vals = zeros(1, length(x_values));
        for i = 1:length(x_values)
            eig_vals(i) = eig_val(i, i);
        end
        eig_onei = abs(1 - eig_vals) == min(abs(1 - eig_vals));
        plot(x_values, eig_vecL(:, eig_onei));
        title("eigenvector = s = " + eig_vals(eig_onei));
    end
elseif method == "CrankNicolsonUpperA_var_a"
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
                f1 = (a(idx) * dt) / (2 * dx);
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
else
    error('Error : method "' + method + '" not supported.');
end

if wtf == true
    info = "a = " + a + "; D = " + D + "; x0 = " + x0 + "; lower boundary : " + lower_boundary + "; upper boundary : " + upper_boundary + "; method : " + method;
    file = fopen('NumericalIntegrationResults.txt', 'w');
    fprintf(file, info + '\n');
    fprintf(file, strjoin(compose('%g', extract_times), ','));
    fprintf(file, '\n');
    fprintf(file', strjoin(compose('%g', x_values), ','));
    for i = 1:length(extract_times)
        fprintf(file, '\n');
        fprintf(file, strjoin(compose('%g', y_values(i, :)), ','));
    end
    fclose(file);
    disp("'NumericalIntegrationResults.txt' written.")
end

end

