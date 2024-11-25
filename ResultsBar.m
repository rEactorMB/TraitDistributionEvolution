function ResultsBar(extract_times, bins, population_history, normalise, t_pause, x_lim)
% Function to display data which are sequential histograms.  If extract_times
% is of length n and bins is of length m, population_history is an nxm matrix
% where the ith row is a the bin fill levels for bin values given by bins for
% the histogram at time = extract_times(i).
% normalise allows you to normalise the area under each histogram to one.
% t_pause sets the pause between displaying each subsequent histogram in seconds.
% x_lim sets the x-axis domain for the entire plotting duration.
for i = 1:length(extract_times)
    population = population_history(i, :);
    if normalise == true
        bin_width = bins(2) - bins(1);
        population_size = sum(population * bin_width);
        population = population / population_size;
    end
    bar(bins, population);
    title("Time = " + extract_times(i));
    if nargin == 6
        xlim(x_lim);
    end
    pause(t_pause);
end
end
