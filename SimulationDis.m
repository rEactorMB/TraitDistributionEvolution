function [extract_times, bin_values, population_history] = SimulationDis(a, D, b, d, x0, a_params, D_params, b_params, d_params, lower_boundary, lower_boundary_type, upper_boundary, upper_boundary_type, extract_times, traitspace_splits, runs)
% Gillespie style algorithm so simulate the evolution of a population
% where each population member gives birth with rate b, dies with rate
% d, mutates symmetrically with rate D and mutates in a deterministic
% direction with rate |a|.  The direction of the deterministic mutation
% is given by the sign of a.  The space in which the population lives
% is a discrete, one dimensional space.

% a, D, b and d are function inputs such whose arguments are x (the
% space) and a set of parameters: a_params, D_params, b_params and d_params.
% lower_boundary and upper_boundary set the position of the boundaries on
% the space.  x0 sets the position the initial population member is placed
% in the space.  The boundary types (lower_boundary_type and
% upper_boundary_type) are 'A' (Dirichlet), 'N' (Neumann), or 'R' (zero-
% flux (Robin-type)).
% extract_times are the times the algorithm will output the state of the
% population for.
% traitspace_splits is the number of bins the space is split into.
% runs is an integer which can be used to perform multiple simulations
% and average the results.

if upper_boundary_type ~= 'A' && upper_boundary_type ~= 'R' && upper_boundary_type ~= 'N'
    error('Error : Not a valid upper boundary type.  Use "A" (absorbing (Dirichlet)), "R" (zero-flux (Robin-type)) or "N" (reflecting (Neumann))')
end
if lower_boundary_type ~= 'A' && upper_boundary_type ~= 'R' && upper_boundary_type ~= 'N'
    error('Error : Not a valid upper boundary type.  Use "A" (absorbing (Dirichlet)), "R" (zero-flux (Robin-type)) or "N" (reflecting (Neumann))')
end

if extract_times < 0
    error('Error : entry in extract_times is less than 0')
end

if traitspace_splits <= 0
    error('Error : traitspace_splits must be a positive integer')
elseif traitspace_splits ~= round(traitspace_splits)
    error('Error : traitspace_splits must be a positive integer')
end

edges = linspace(lower_boundary, upper_boundary, (traitspace_splits + 1));
bin_separation = edges(2) - edges(1);
bin_values = edges(2:end) - (bin_separation * 0.5);
truth_array = edges <= x0;
start_bin = find(truth_array, 1, 'last');
population_array = zeros(1, traitspace_splits);
population_array(start_bin) = 1;
event_time = 0;
extract_times = sort(extract_times);
t_stop = extract_times(end);

if extract_times(1) == 0
    remaining_times = extract_times(2:end);
    population_history = population_array;
else
    remaining_times = extract_times;
    population_history = [];
end
if remaining_times(1) == 0
    warning('Warning : extract_times has 0 in it twice which may cause problems later on')
end

extinction = false;

birth_rate = zeros(1, traitspace_splits);
death_rate = zeros(1, traitspace_splits);
drift_rate = zeros(1, traitspace_splits);
diffusion_rate = zeros(1, traitspace_splits);
for i = 1:traitspace_splits
    x = bin_values(i);
    birth_rate(i) = b(x, b_params);
    death_rate(i) = d(x, d_params);
    drift_rate(i) = abs(a(x, a_params)) / bin_separation;
    diffusion_rate(i) = D(x, D_params) / (bin_separation * bin_separation);
end
summed_rate = birth_rate + death_rate + drift_rate + diffusion_rate;

while event_time < t_stop && extinction == false
    total_rate = 0;
    for i = 1:traitspace_splits
        total_rate = total_rate + (population_array(i) * summed_rate(i));
    end
    event_time = event_time + ExpRand(total_rate);
    
    random_choice = rand * total_rate;
    chosen_index = 1;
    cummulative_rate = summed_rate(chosen_index) * population_array(chosen_index);
    while random_choice > cummulative_rate
        chosen_index = chosen_index + 1;
        cummulative_rate = cummulative_rate + (summed_rate(chosen_index) * population_array(chosen_index));
    end

    if population_array(chosen_index) <= 0
        error('Error : bin chosen has no members - something has gone wrong')
    end

    bin_total_rate = summed_rate(chosen_index);
    bin_random_choice = rand * bin_total_rate;
    birth_value = birth_rate(chosen_index);
    death_value = birth_value + death_rate(chosen_index);
    drift_value = death_value + drift_rate(chosen_index);
    diffusion_value = drift_value + diffusion_rate(chosen_index);
    if bin_random_choice <= birth_value
        population_array(chosen_index) = population_array(chosen_index) + 1;
    elseif bin_random_choice <= death_value
        population_array(chosen_index) = population_array(chosen_index) - 1;
    elseif bin_random_choice <= drift_value
        if a(bin_values(chosen_index), a_params) > 0
            sign_a = 1;
        elseif a(bin_values(chosen_index), a_params) < 0
            sign_a = -1;
        end
        if chosen_index == traitspace_splits && sign_a == 1
            if upper_boundary_type == 'A' || upper_boundary_type == 'N'
                population_array(chosen_index) = population_array(chosen_index) - 1;
            end
        elseif chosen_index == traitspace_splits && sign_a == -1 && upper_boundary_type == 'N'
            population_array(chosen_index + sign_a) = population_array(chosen_index + sign_a) + 1;
        elseif chosen_index == 1 && sign_a == -1
            if lower_boundary_type == 'A'
                population_array(chosen_index) = population_array(chosen_index) - 1;
            end
        else
            population_array(chosen_index + sign_a) = population_array(chosen_index + sign_a) + 1;
            population_array(chosen_index) = population_array(chosen_index) - 1;
        end
    elseif bin_random_choice <= diffusion_value
        if rand < 0.5
            direction = 1;
        else
            direction = -1;
        end
        if chosen_index == traitspace_splits && direction == 1
            if upper_boundary_type == 'A'
                population_array(chosen_index) = population_array(chosen_index) - 1;
            end
        elseif chosen_index == 1 && direction == -1
            if lower_boundary_type == 'A'
                population_array(chosen_index) = population_array(chosen_index) - 1;
            end
        else
            population_array(chosen_index + direction) = population_array(chosen_index + direction) + 1;
            population_array(chosen_index) = population_array(chosen_index) - 1;
        end
    end

    try
        if remaining_times(1) <= event_time
            remaining_times(1) = [];
            population_history = cat(1, population_history, population_array);
            % disp("Current event time : " + event_time);
        end
    catch
    end
    
    if sum(population_array) == 0
        extinction = true;
        if ~isempty(remaining_times)
            zero_array = zeros(1, traitspace_splits);
            for i = 1:length(remaining_times)
                population_history = cat(1, population_history, zero_array);
            end
        end
    end
end

if runs > 1
    runs_completed = 1;
    while runs_completed < runs
        [~, ~, p_h] = SimulationDis(a, D, b, d, x0, a_params, D_params, b_params, d_params, lower_boundary, lower_boundary_type, upper_boundary, upper_boundary_type, extract_times, traitspace_splits, 1, false);
        if isscalar(extract_times)
            population_history = population_history + p_h;
        else
            for i = 1:length(extract_times)
                population_history(i, :) = population_history(i, :) + p_h(i, :);
            end
        end
        runs_completed = runs_completed + 1;
    end
    population_history = population_history / runs;
end
end
