function ResultsBar(extract_times, bins, population_history, normalise, t_pause, x_lim)
for i = 1:length(extract_times)
    population = population_history(i, :);
    if normalise == true
        bin_width = bins(2) - bins(1);
        population_size = sum(population * bin_width);
        population = population / population_size;
    end
    bar(bins, population);
    % title("Pop size = " + sum(population_history(i, :)));
    title("Time = " + extract_times(i));
    if nargin == 6
        xlim(x_lim);
    end
    pause(t_pause);
end
end
