% Import data from BivalveSummary.csv and BrachiopodSummary.csv
Bivalvia_summary = readtable("BivalveSummary.csv");
Brachiopoda_summary = readtable("BrachiopodSummary.csv");

% Get out size, first appearance date (fad) and last appearance date (lad)
% fad is the maximum value (mya), i.e. furthest back in time
% lad is the minimum (mya), i.e. cloest in time to current day
BiSum_size = table2array(Bivalvia_summary(:, "size"));
BiSum_fad = table2array(Bivalvia_summary(:, "firstapp_max_ma"));
BiSum_lad = table2array(Bivalvia_summary(:, "lastapp_min_ma"));

BrSum_size = table2array(Brachiopoda_summary(:, "size"));
BrSum_fad = table2array(Brachiopoda_summary(:, "firstapp_max_ma"));
BrSum_lad = table2array(Brachiopoda_summary(:, "lastapp_min_ma"));

% Strip out NaNs from datasets
BiTruth = logical(~isnan(BiSum_size) .* ~isnan(BiSum_fad) .* ~isnan(BiSum_lad));
BrTruth = logical(~isnan(BrSum_size) .* ~isnan(BrSum_fad) .* ~isnan(BrSum_lad));

BiSum_size = BiSum_size(BiTruth);
BiSum_fad = BiSum_fad(BiTruth);
BiSum_lad = BiSum_lad(BiTruth);

BrSum_size = BrSum_size(BrTruth);
BrSum_fad = BrSum_fad(BrTruth);
BrSum_lad = BrSum_lad(BrTruth);

% Use this code to remove bivalve outliers: two outliers (i = 1:2) are
% Tridacna and Pinna.
warning("data removed from dataset bivalve (line 29).");
for i = 1:2
    BiSum_lad = BiSum_lad(BiSum_size ~= max(BiSum_size));
    BiSum_fad = BiSum_fad(BiSum_size ~= max(BiSum_size));
    BiSum_size = BiSum_size(BiSum_size ~= max(BiSum_size));
end

% Make sure the minimum size is not less than 1 so no negative values
% appear when you look at log10(size)
Bi_10_scale = 0;
while min(BiSum_size) <= 1
    BiSum_size = BiSum_size * 10;
    Bi_10_scale = Bi_10_scale + 1;
end
Br_10_scale = 0;
while min(BrSum_size) <= 1
    BrSum_size = BrSum_size * 10;
    Br_10_scale = Br_10_scale + 1;
end
disp("Bivalve size data scaled by 10^" + Bi_10_scale);
disp("Brachiopod size data scaled by 10^" + Br_10_scale);

% If you only want to look at data before or after the Permian-Triassic
% extinction event then use these boolean objects to set which time period
% you wish to look at (true -> include data).
prePTE = true;
postPTE = true;

%%
% Bivalve analysis
BiSum_fad = -BiSum_fad;
BiSum_lad = -BiSum_lad;
Bi_min_fad = min(BiSum_fad);
Bi_max_lad = max(BiSum_lad);
BiSum_log_size = log10(BiSum_size);
Bi_max_size = max(BiSum_log_size);
Bi_min_size = min(BiSum_log_size);
Bi_dates = linspace(Bi_min_fad, Bi_max_lad, 101);
Bi_dates = Bi_dates(1:100) + (0.5 * (Bi_dates(2) - Bi_dates(1)));

% Add dates you want for the histograms (for plotting specific
% dates later)
Bi_plot_dates = [-485.4, -359.2, -145.5, -2.588];
Bi_dates(end+1:end+length(Bi_plot_dates)) = Bi_plot_dates;
Bi_dates = sort(Bi_dates);

%%
% Bivalve log likelihood analysis V1 (analytical solution)
% --> Dirichlet boundaries only.
% Find all the time intervals which have a unique histogram (i.e.
% a unique set of species) and the midpoints of these intervals.
Bi_unique_dates = zeros(1, 2 * length(BiSum_lad));
next_fill = 1;
for i = 1:length(BiSum_fad)
    if ~any(Bi_unique_dates == BiSum_fad(i))
        Bi_unique_dates(next_fill) = BiSum_fad(i);
        next_fill = next_fill + 1;
    end
    if ~any(Bi_unique_dates == BiSum_lad(i))
        Bi_unique_dates(next_fill) = BiSum_lad(i);
        next_fill = next_fill + 1;
    end
end
Bi_unique_dates = Bi_unique_dates(1:next_fill-1);
Bi_unique_dates = sort(Bi_unique_dates);

Bi_unique_im = 0.5 * (Bi_unique_dates(2:end) + Bi_unique_dates(1:end-1));

if prePTE && ~postPTE
    Bi_unique_im = Bi_unique_im(Bi_unique_im < -252);
    next_fill = length(Bi_unique_im) + 2;
elseif ~prePTE && postPTE
    Bi_unique_im = Bi_unique_im(Bi_unique_im > -252);
    next_fill = length(Bi_unique_im) + 2;
end

% Create the datapoints which are the time (unique interval midpoint) and
% logsize value for each Bivalve present in each unique interval.

Bi_sample_size = zeros(1, next_fill-2);
for i = 1:(next_fill-2)
    Bi_sample_size(i) = sum((BiSum_fad <= Bi_unique_im(i)) .* (BiSum_lad > Bi_unique_im(i)));
end

% New faster method
Bi_times = Bi_unique_im;
Bi_log_sizes = zeros(next_fill-2, max(Bi_sample_size));
for i = 1:(next_fill-2)
    Bi_log_sizes(i, 1:Bi_sample_size(i)) = BiSum_log_size(logical((BiSum_fad <= Bi_unique_im(i)) .* (BiSum_lad > Bi_unique_im(i))));
end

% Set the furthest point back in time you will allow the search algorithm
% to go.  e.g. -514 -> 514 mya
Bi_t_offset = -514;
Bi_log_sizes(Bi_log_sizes == 0) = NaN;

func2min = @(theta) -LogLikelihoodSSConAb(Bi_log_sizes - theta(5), theta(1), theta(2), theta(4), theta(3), Bi_times - theta(6));

options = optimoptions('fmincon', 'PlotFcn','optimplotfval', 'Algorithm', 'interior-point', 'EnableFeasibilityMode', true, 'SubproblemAlgorithm', 'cg', 'MaxFunctionEvaluations', 5000);
Bi_OptVals = fmincon(func2min, [0.0001, 0.0005, 1.4, 3, 0.2, -510], [0, 0, 1, -1, 0, 0; 0, -1, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0; 0, 0, -1, 0, 0, 0; 0, 0, 0, -1, -1, 0; 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, -1], [0; 0; Bi_min_size; 0; -Bi_max_size; min(Bi_times); -Bi_t_offset], [], [], [], [], [], options);

% Warn user if search algorithm returns values close to the constraints.
% Proximity can be set to determine how close the value needs to be
% before a warning is issued.
proximity = 0.001;
if abs(Bi_OptVals(5) - Bi_min_size) < abs(proximity * Bi_min_size)
    disp("WARNING: Bi_OptVals(5) (size offset) is sitting near the constraing boundary.");
end
if abs(Bi_OptVals(4) - Bi_max_size + Bi_OptVals(5)) < abs(proximity * (Bi_max_size - Bi_OptVals(5)))
    disp("WARNING: Bi_OptVals(4) (L) is sitting near the constraint boundary.");
end
if abs(Bi_OptVals(6) - min(Bi_times)) < abs(proximity * min(Bi_times)) || abs(Bi_OptVals(6) - Bi_t_offset) < abs(proximity * Bi_t_offset)
    disp("WARNING: Bi_OptVals(6) (time offset) is sitting near the constraint boundary.");
end

%%
% Bivalve log likelihood analysis V2 (numerical solution)
% --> Dirichlet, Neumann or no flux (Robin type) boundaries.
% Find all the time intervals which have a unique histogram and the
% midpoints of these intervals.
Bi_unique_dates = zeros(1, 2 * length(BiSum_lad));
next_fill = 1;
for i = 1:length(BiSum_fad)
    if ~any(Bi_unique_dates == BiSum_fad(i))
        Bi_unique_dates(next_fill) = BiSum_fad(i);
        next_fill = next_fill + 1;
    end
    if ~any(Bi_unique_dates == BiSum_lad(i))
        Bi_unique_dates(next_fill) = BiSum_lad(i);
        next_fill = next_fill + 1;
    end
end
Bi_unique_dates = Bi_unique_dates(1:next_fill-1);
Bi_unique_dates = sort(Bi_unique_dates);

Bi_unique_im = 0.5 * (Bi_unique_dates(2:end) + Bi_unique_dates(1:end-1));

if prePTE && ~postPTE
    Bi_unique_im = Bi_unique_im(Bi_unique_im < -252);
    next_fill = length(Bi_unique_im) + 2;
elseif ~prePTE && postPTE
    Bi_unique_im = Bi_unique_im(Bi_unique_im > -252);
    next_fill = length(Bi_unique_im) + 2;
end

% Create the datapoints which are the time (unique interval midpoint) and
% logsize value for each Bivalve present in each unique interval.

Bi_sample_size = zeros(1, next_fill-2);
for i = 1:(next_fill-2)
    Bi_sample_size(i) = sum((BiSum_fad <= Bi_unique_im(i)) .* (BiSum_lad > Bi_unique_im(i)));
end

% New faster method
Bi_times = Bi_unique_im;
Bi_log_sizes = zeros(next_fill-2, max(Bi_sample_size));
for i = 1:(next_fill-2)
    Bi_log_sizes(i, 1:Bi_sample_size(i)) = BiSum_log_size(logical((BiSum_fad <= Bi_unique_im(i)) .* (BiSum_lad > Bi_unique_im(i))));
end

% Set the furthest point back in time you will allow the search algorithm
% to go.  e.g. -514 -> 514 mya
Bi_t_offset = -514;
Bi_log_sizes(Bi_log_sizes == 0) = NaN;

% Set the upper boundary type by specifying the method
% (Corresponds to methods in NumericalPDESolve.m)
method = "CrankNicolsonUpperRrob";
func2min = @(theta) -NumericalLogLikelihood(Bi_times - theta(6), Bi_log_sizes - theta(5), theta(1), theta(2), theta(3), theta(4), method);

options = optimoptions('fmincon', 'PlotFcn','optimplotfval', 'Algorithm', 'interior-point', 'EnableFeasibilityMode', true, 'SubproblemAlgorithm', 'cg', 'MaxFunctionEvaluations', 5000);
Bi_OptVals = fmincon(func2min, [0.0005628, 0.0004676, 1.254, 3.032, 0.205, -514], [0, 0, 1, -1, 0, 0; 0, -1, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0; 0, 0, -1, 0, 0, 0; 0, 0, 0, -1, -1, 0; 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, -1], [0; 0; Bi_min_size; 0; -Bi_max_size; min(Bi_times); -Bi_t_offset], [], [], [], [], [], options);

% Warn user if search algorithm returns values close to the constraints.
% Proximity can be set to determine how close the value needs to be
% before a warning is issued.
proximity = 0.001;
if abs(Bi_OptVals(5) - Bi_min_size) < abs(proximity * Bi_min_size)
    disp("WARNING: Bi_OptVals(5) (size offset) is sitting near the constraing boundary.");
end
if abs(Bi_OptVals(4) - Bi_max_size + Bi_OptVals(5)) < abs(proximity * (Bi_max_size - Bi_OptVals(5)))
    disp("WARNING: Bi_OptVals(4) (L) is sitting near the constraint boundary.");
end
if abs(Bi_OptVals(6) - min(Bi_times)) < abs(proximity * min(Bi_times)) || abs(Bi_OptVals(6) - Bi_t_offset) < abs(proximity * Bi_t_offset)
    disp("WARNING: Bi_OptVals(6) (time offset) is sitting near the constraint boundary.");
end

%%

% Bivalve analysis plotting
% Calculate the numerical solution to the PDE to overlay on
% histograms of the data.
extract_times = linspace(0, -Bi_OptVals(6), 100 * -round(Bi_OptVals(6)));
method = "CrankNicolsonUpperRrob";
[~, NS_x_values, NS_y_values] = NumericalPDESolve(Bi_OptVals(1), Bi_OptVals(2), Bi_OptVals(3), 0, Bi_OptVals(4), extract_times, 200, false, method, false);
ti_arr = zeros(1, length(Bi_dates));
% pause_count allows you to observe the plots corresponding to the dates in
% Bi_plot_dates for longer.
pause_count = 1;

% Plot sequentially the hisotgrams of size data with the solution to the
% PDE with maximally likely parameters
figure;
set(gcf, 'Position', [600, 400, 600, 300]);
for i = 1:length(Bi_dates)
    date_in_interval = logical((BiSum_fad <= Bi_dates(i)) .* (BiSum_lad > Bi_dates(i)));
    histogram(BiSum_log_size(date_in_interval) - Bi_OptVals(5), 20, Normalization="pdf");
    hold on;
    
    ti = find(abs(extract_times - Bi_dates(i) + Bi_OptVals(6)) == min(abs(extract_times - Bi_dates(i) + Bi_OptVals(6))));
    ti_arr(i) = ti;
    plot(NS_x_values, NS_y_values(ti, :) / NumericalIntegrator(NS_x_values, NS_y_values(ti, :)), "LineWidth", 2, "Color", "g");
    
    hold off;
    xlim([0, Bi_max_size]);
    ylim([0, 2]);
    xlabel("$\log_{10}(\mathrm{size})$", Interpreter='latex');
    ylabel("Density", Interpreter='latex');
    title("Time = " + -Bi_dates(i) + " million years ago.  (" + sum(date_in_interval) + " units)");
    fontsize(gcf, 25, 'pixels');
    pause(0.1);

    try
        if Bi_plot_dates(pause_count) <= Bi_dates(i)
            title("");
            pause_count = pause_count + 1;
            pause(2);
        end
    catch
    end
end

%%
% Bivalve peak and bulk motion analysis
% Find the poition of the peak (smoothed using KernelDensityFunction)
% and the motion of the mean (bulk).
bin_number = 50;
edges = linspace(0, Bi_max_size, bin_number + 1);
bins = edges(1:end-1) + (0.5 * (edges(2) - edges(1)));
distribution_history = zeros(length(Bi_dates), 1000);
for i = 1:length(Bi_dates)
    date_in_interval = logical((BiSum_fad <= Bi_dates(i)) .* (BiSum_lad > Bi_dates(i)));
    [number, ~] = histcounts(BiSum_log_size(date_in_interval) - Bi_OptVals(5), edges);
    [x_values, distribution_history(i, :)] = KernelDensityFunction(bins, number, false);
end
% Compare data and numerical solution with maximally likely parameters.
[t_data, p_data, b_data, t_sim, p_sim, b_sim] = PeakandBulkMotion(Bi_dates, x_values, distribution_history, 'Bivalves', Bi_dates, NS_x_values, NS_y_values(ti_arr, :), 'Numerical solution');

% Calculating the diversity at each time in t_data
diversity = zeros(1, length(t_data));
for i = 1:length(t_data)
    diversity(i) = sum((t_data(i) >= BiSum_fad) .* (t_data(i) <= BiSum_lad));
end

% Plot peak and bulk motion with diversity and colours for geological periods.
% Ensure that the dates of the periods are only the dates you are plotting
% otherwise the scale will be wrong.  E.g. if you want to have a date range
% from 0 to 500 mya then you need to set:
% Cambrian = int64(scale*(500-488.3)).
% scale should be used to make sure there are no decimals so the widths of
% each colour are not rounded up or down.
% Dates from https://ucmp.berkeley.edu/help/timeform.php
% Last accessed: 25/11/2024
% Colours from https://timescalefoundation.org/charts/RGB.pdf
% Last accessed: 25/11/2024
scale = 1000;
% Proterozoic = int64(scale*(2500-542.0));
Proterozoic = 0;
% Cambrian = int64(scale*(542.0-488.3));
Cambrian = int64(scale*(500-488.3));
Ordovician = int64(scale*(488.3-443.7));
Silurian = int64(scale*(443.7-416.0));
Devonian = int64(scale*(416.0-359.2));
Carboniferous = int64(scale*(359.2-299.0));
Permian = int64(scale*(299.0-251.0));
Triassic = int64(scale*(251.0-199.6));
Jurassic = int64(scale*(199.6-145.5));
Cretaceous = int64(scale*(145.5-65.5));
Paleogene = int64(scale*(65.5-23.03));
Neogene = int64(scale*(23.03-2.588));
Quaternary = int64(scale*(2.588-0));
cmap = [repmat([247,53,99], [Proterozoic, 1]); repmat([127,160,86], [Cambrian, 1]); repmat([0,146,112], [Ordovician, 1]); repmat([179,225,182], [Silurian,1]); repmat([203,140,55], [Devonian,1]); repmat([103,165,153], [Carboniferous,1]); repmat([240,64,40], [Permian,1]); repmat([129,43,146], [Triassic,1]); repmat([52,178,201], [Jurassic,1]); repmat([127,198,78], [Cretaceous,1]); repmat([253,154,82], [Paleogene,1]); repmat([255,230,25], [Neogene,1]); repmat([249,249,127], [Quaternary,1])];
% RGB colours on scale from 0 to 255
cmap = cmap / 255;
cmap = flip(cmap);

figure;
tiledlayout(4, 1, TileSpacing='none');
set(gcf, 'Position', [600, 400, 500, 400]);
nexttile;
plot(-t_data, diversity, color='black');
ylim([0,510]);
xlim([0, 500]);
set(gca, 'xdir', 'reverse');
title("Peak", Interpreter='latex');
ylabel("Diversity", Interpreter='latex');
xticks([]);
yticks([0,250,500]);
nexttile(2, [3,1]);
plot(-t_data, p_data, -t_sim, p_sim);
legend('Data', 'Maximally likely model', Location='northwest', Interpreter='latex');
ylabel("$\log_{10}(\mathrm{size})$", Interpreter='latex');
xlabel("Time (mya)", Interpreter='latex');
ylim([1.15,1.7]);
xlim([0, 500]);
yticks(1.2:0.1:1.6);
fontsize(gcf, 20, 'pixel');

ax = gca;
ax.Colormap = colormap(cmap);
cb = colorbar;
cb.Location = 'south';
cb.Position(1) = ax.Position(1);
cb.Position(2) = ax.Position(2);
cb.Position(3) = ax.Position(3);
cb.Ticks = [];
set(gca, 'xdir', 'reverse');
cb.Direction = 'reverse';

figure;
tiledlayout(4, 1, TileSpacing='none');
set(gcf, 'Position', [600, 400, 500, 400]);
nexttile;
plot(-t_data, diversity, color='black');
ylim([0, 510]);
xlim([0, 500]);
set(gca, 'xdir', 'reverse');
title("Mean", Interpreter='latex');
ylabel("Diversity", Interpreter='latex');
xticks([]);
yticks([0,250,500]);
nexttile(2, [3,1]);
plot(-t_data, b_data, -t_sim, b_sim);
legend('Data', 'Maximally likely model', Location='northwest', Interpreter='latex');
ylabel("$\log_{10}(\mathrm{size})$", Interpreter='latex');
xlabel("Time (mya)", Interpreter='latex');
ylim([1.1,1.7]);
xlim([0, 500]);
fontsize(gcf, 20, 'pixel');

ax = gca;
ax.Colormap = colormap(cmap);
cb = colorbar;
cb.Location = 'south';
cb.Position(1) = ax.Position(1);
cb.Position(2) = ax.Position(2);
cb.Position(3) = ax.Position(3);
cb.Ticks = [];
set(gca, 'xdir', 'reverse');
cb.Direction = 'reverse';

%%
% Brachiopod analysis
BrSum_fad = -BrSum_fad;
BrSum_lad = -BrSum_lad;
Br_min_fad = min(BrSum_fad);
Br_max_lad = max(BrSum_lad);
BrSum_log_size = log10(BrSum_size);
Br_max_size = max(BrSum_log_size);
Br_min_size = min(BrSum_log_size);
Br_dates = linspace(Br_min_fad, Br_max_lad, 101);
Br_dates = Br_dates(1:100) + (0.5 * (Br_dates(2) - Br_dates(1)));

% Add dates you want for the histograms
Br_plot_dates = [-541, -488.3, -433.7, -359.2, -145.5, -2.588];
Br_dates(end+1:end+length(Br_plot_dates)) = Br_plot_dates;
Br_dates = sort(Br_dates);

%%
% Brachiopod log likelihood analysis V1 (analytical solution)
% --> Dirichlet upper boundary
% Find all the time intervals which have a unique histogram and the
% midpoints of these intervals.
Br_unique_dates = zeros(1, 2 * length(BrSum_lad));
next_fill = 1;
for i = 1:length(BrSum_fad)
    if ~any(Br_unique_dates == BrSum_fad(i))
        Br_unique_dates(next_fill) = BrSum_fad(i);
        next_fill = next_fill + 1;
    end
    if ~any(Br_unique_dates == BrSum_lad(i))
        Br_unique_dates(next_fill) = BrSum_lad(i);
        next_fill = next_fill + 1;
    end
end
Br_unique_dates = Br_unique_dates(1:next_fill-1);
Br_unique_dates = sort(Br_unique_dates);

Br_unique_im = 0.5 * (Br_unique_dates(2:end) + Br_unique_dates(1:end-1));

if prePTE && ~postPTE
    Br_unique_im = Br_unique_im(Br_unique_im < -252);
    next_fill = length(Br_unique_im) + 2;
elseif ~prePTE && postPTE
    Br_unique_im = Br_unique_im(Br_unique_im > -252);
    next_fill = length(Br_unique_im) + 2;
end

% Create the datapoints which are the time (unique interval midpoint) and
% logsize value for each Bivalve present in each unique interval.
Br_sample_size = zeros(1, next_fill-2);
for i = 1:(next_fill-2)
    Br_sample_size(i) = sum((BrSum_fad <= Br_unique_im(i)) .* (BrSum_lad > Br_unique_im(i)));
end

% New faster method
Br_times = Br_unique_im;
Br_log_sizes = zeros(next_fill-2, max(Br_sample_size));
for i = 1:(next_fill-2)
    Br_log_sizes(i, 1:Br_sample_size(i)) = BrSum_log_size(logical((BrSum_fad <= Br_unique_im(i)) .* (BrSum_lad > Br_unique_im(i))));
end

% Set the furthest point back in time you will allow the search algorithm
% to go.  e.g. -607 -> 607 mya
Br_t_offset = -607;
Br_log_sizes(Br_log_sizes == 0) = NaN;

func2min = @(theta) -LogLikelihoodSSConAb(Br_log_sizes - theta(5), theta(1), theta(2), theta(4), theta(3), Br_times - theta(6));

options = optimoptions('fmincon', 'PlotFcn','optimplotfval', 'Algorithm', 'interior-point', 'EnableFeasibilityMode', true, 'SubproblemAlgorithm', 'cg', 'MaxFunctionEvaluations', 5000);
Br_OptVals = fmincon(func2min, [0.002, 0.0007, 1.4, 2.8, 0.5, -600], [0, 0, 1, -1, 0, 0; 0, -1, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0; 0, 0, -1, 0, 0, 0; 0, 0, 0, -1, -1, 0; 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, -1], [0; 0; Br_min_size; 0; -Br_max_size; min(Br_times); -Br_t_offset], [], [], [], [], [], options);

% Warn user if search algorithm returns values close to the constraints.
% Proximity can be set to determine how close the value needs to be
% before a warning is issued.
proximity = 0.001;
if abs(Br_OptVals(5) - Br_min_size) < abs(proximity * Br_min_size)
    disp("WARNING: Br_OptVals(5) (size offset) is sitting near the constraing boundary.");
end
if abs(Br_OptVals(4) - Br_max_size + Br_OptVals(5)) < abs(proximity * (Br_max_size - Br_OptVals(5)))
    disp("WARNING: Br_OptVals(4) (L) is sitting near the constraint boundary.");
end
if abs(Br_OptVals(6) - min(Br_times)) < abs(proximity * min(Br_times)) || abs(Br_OptVals(6) - Br_t_offset) < abs(proximity * Br_t_offset)
    disp("WARNING: Br_OptVals(6) (time offset) is sitting near the constraint boundary.");
end

%%
% Brachiopod log likelihood analysis V2 (numerical solution)
% --> Dirichlet, Neumann and no flux (Robin type) upper boundary.
% Find all the time intervals which have a unique histogram and the
% midpoints of these intervals.
Br_unique_dates = zeros(1, 2 * length(BrSum_lad));
next_fill = 1;
for i = 1:length(BrSum_fad)
    if ~any(Br_unique_dates == BrSum_fad(i))
        Br_unique_dates(next_fill) = BrSum_fad(i);
        next_fill = next_fill + 1;
    end
    if ~any(Br_unique_dates == BrSum_lad(i))
        Br_unique_dates(next_fill) = BrSum_lad(i);
        next_fill = next_fill + 1;
    end
end
Br_unique_dates = Br_unique_dates(1:next_fill-1);
Br_unique_dates = sort(Br_unique_dates);

Br_unique_im = 0.5 * (Br_unique_dates(2:end) + Br_unique_dates(1:end-1));

if prePTE && ~postPTE
    Br_unique_im = Br_unique_im(Br_unique_im < -252);
    next_fill = length(Br_unique_im) + 2;
elseif ~prePTE && postPTE
    Br_unique_im = Br_unique_im(Br_unique_im > -252);
    next_fill = length(Br_unique_im) + 2;
end

% Create the datapoints which are the time (unique interval midpoint) and
% logsize value for each Bivalve present in each unique interval.
Br_sample_size = zeros(1, next_fill-2);
for i = 1:(next_fill-2)
    Br_sample_size(i) = sum((BrSum_fad <= Br_unique_im(i)) .* (BrSum_lad > Br_unique_im(i)));
end

% New faster method
Br_times = Br_unique_im;
Br_log_sizes = zeros(next_fill-2, max(Br_sample_size));
for i = 1:(next_fill-2)
    Br_log_sizes(i, 1:Br_sample_size(i)) = BrSum_log_size(logical((BrSum_fad <= Br_unique_im(i)) .* (BrSum_lad > Br_unique_im(i))));
end

% Set the furthest point back in time you will allow the search algorithm
% to go.  e.g. -607 -> 607 mya
Br_t_offset = -607;
Br_log_sizes(Br_log_sizes == 0) = NaN;

% Set the upper boundary type by specifying the method
% (Corresponds to methods in NumericalPDESolve.m)
method = "CrankNicolsonUpperRrob";
func2min = @(theta) -NumericalLogLikelihood(Br_times - theta(6), Br_log_sizes - theta(5), theta(1), theta(2), theta(3), theta(4), method);

options = optimoptions('fmincon', 'PlotFcn','optimplotfval', 'Algorithm', 'interior-point', 'EnableFeasibilityMode', true, 'SubproblemAlgorithm', 'cg', 'MaxFunctionEvaluations', 5000);
Br_OptVals = fmincon(func2min, [0.001493, 0.000741, 1.411, 2.797, 0.4759, -607], [0, 0, 1, -1, 0, 0; 0, -1, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0; 0, 0, -1, 0, 0, 0; 0, 0, 0, -1, -1, 0; 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, -1], [0; 0; Br_min_size; 0; -Br_max_size; min(Br_times); -Br_t_offset], [], [], [], [], [], options);

% Warn user if search algorithm returns values close to the constraints.
% Proximity can be set to determine how close the value needs to be
% before a warning is issued.
proximity = 0.001;
if abs(Br_OptVals(5) - Br_min_size) < abs(proximity * Br_min_size)
    disp("WARNING: Br_OptVals(5) (size offset) is sitting near the constraing boundary.");
end
if abs(Br_OptVals(4) - Br_max_size + Br_OptVals(5)) < abs(proximity * (Br_max_size - Br_OptVals(5)))
    disp("WARNING: Br_OptVals(4) (L) is sitting near the constraint boundary.");
end
if abs(Br_OptVals(6) - min(Br_times)) < abs(proximity * min(Br_times)) || abs(Br_OptVals(6) - Br_t_offset) < abs(proximity * Br_t_offset)
    disp("WARNING: Br_OptVals(6) (time offset) is sitting near the constraint boundary.");
end

%%
% Brachiopod analysis plotting
extract_times = linspace(0, -Br_OptVals(6), 1000 * -round(Br_OptVals(6)));
method = "CrankNicolsonUpperRrob";
[~, NS_x_values, NS_y_values] = NumericalPDESolve(Br_OptVals(1), Br_OptVals(2), Br_OptVals(3), 0, Br_OptVals(4), extract_times, 1000, false, method, false);
ti_arr = zeros(1, length(Br_dates));
pause_count = 1;

figure;
set(gcf, 'Position', [600, 400, 600, 300]);
for i = 1:length(Br_dates)
    date_in_interval = logical((BrSum_fad <= Br_dates(i)) .* (BrSum_lad > Br_dates(i)));
    histogram(BrSum_log_size(date_in_interval) - Br_OptVals(5), Normalization="pdf");
    hold on;
    
    ti = find(abs(extract_times - Br_dates(i) + Br_OptVals(6)) == min(abs(extract_times - Br_dates(i) + Br_OptVals(6))));
    ti_arr(i) = ti;
    plot(NS_x_values, NS_y_values(ti, :) / NumericalIntegrator(NS_x_values, NS_y_values(ti, :)), "LineWidth", 2, "Color", "r");

    hold off;
    xlim([0, Br_OptVals(4)]);
    ylim([0, 1.5]);
    xlabel("$\log_{10}(\mathrm{size})$", Interpreter='latex');
    ylabel("Density", Interpreter='latex');
    title("Time = " + -Br_dates(i) + " million years ago.  (" + sum(date_in_interval) + " units)")
    fontsize(gcf, 25, 'pixels');
    pause(0.1);

    try
        if Br_plot_dates(pause_count) <= Br_dates(i)
            title("");
            pause_count = pause_count + 1;
            pause(2);
        end
    catch
    end
end
%%

% Brachiopod peak and bulk motion analysis
% Calculate the numerical solution to the PDE to overlay on
% histograms of the data.
bin_number = 20;
edges = linspace(0, Br_OptVals(4), bin_number + 1);
bins = edges(1:end-1) + (0.5 * (edges(2) - edges(1)));
distribution_history = zeros(length(Br_dates), 1000);
for i = 1:length(Br_dates)
    date_in_interval = logical((BrSum_fad <= Br_dates(i)) .* (BrSum_lad > Br_dates(i)));
    [number, ~] = histcounts(BrSum_log_size(date_in_interval) - Br_OptVals(5), edges);
    [x_values, distribution_history(i, :)] = KernelDensityFunction(bins, number, false);
end
% Compare data and numerical solution with maximally likely parameters.
[t_data, p_data, b_data, t_sim, p_sim, b_sim] = PeakandBulkMotion(Br_dates, x_values, distribution_history, 'Brachiopods', Br_dates, NS_x_values, NS_y_values(ti_arr, :), 'Numerical solution');

% Calculating the diversity at each time in t_data
diversity = zeros(1, length(t_data));
for i = 1:length(t_data)
    diversity(i) = sum((t_data(i) >= BrSum_fad) .* (t_data(i) <= BrSum_lad));
end

% Plot peak and bulk motion with diversity and colours for geological periods.
% There are may commented lines in the following section of code.  These are
% for plotting the results for three different scenarios: all data; data for
% pre-PTME species; data for post-PTME species.  The variables to be adjusted
% for each scenario are the color-date scalings, the y-limits and the x-limits
% for the plotting region.
% Ensure that the dates of the periods are only the dates you are plotting
% otherwise the scale will be wrong.  E.g. if you want to have a date range
% from 0 to 500 mya then you need to set:
% Cambrian = int64(scale*(500-488.3)).
% scale should be used to make sure there are no decimals so the widths of
% each colour are not rounded up or down.
% Dates from https://ucmp.berkeley.edu/help/timeform.php
% Last accessed 25/11/2024
% Colours from https://timescalefoundation.org/charts/RGB.pdf
% Last accessed 25/11/2024
scale = 1000;
% Proterozoic = int64(scale*(2500-542.0));
Proterozoic = int64(scale*(550-542.0));
% Proterozoic = 0;
Cambrian = int64(scale*(542.0-488.3));
% Cambrian = 0;
Ordovician = int64(scale*(488.3-443.7));
% Ordovician = 0;
Silurian = int64(scale*(443.7-416.0));
% Silurian = 0;
Devonian = int64(scale*(416.0-359.2));
% Devonian = 0;
Carboniferous = int64(scale*(359.2-299.0));
% Carboniferous = 0;
Permian = int64(scale*(299.0-251.0));
% Permian = int64(scale*(260-251.0));
Triassic = int64(scale*(251.0-199.6));
% Triassic = int64(scale*(251.0-240));
Jurassic = int64(scale*(199.6-145.5));
% Jurassic = 0;
Cretaceous = int64(scale*(145.5-65.5));
% Cretaceous = 0;
Paleogene = int64(scale*(65.5-23.03));
% Paleogene = 0;
Neogene = int64(scale*(23.03-2.588));
% Neogene = 0;
Quaternary = int64(scale*(2.588-0));
% Quaternary = 0;
cmap = [repmat([247,53,99], [Proterozoic, 1]); repmat([127,160,86], [Cambrian, 1]); repmat([0,146,112], [Ordovician, 1]); repmat([179,225,182], [Silurian,1]); repmat([203,140,55], [Devonian,1]); repmat([103,165,153], [Carboniferous,1]); repmat([240,64,40], [Permian,1]); repmat([129,43,146], [Triassic,1]); repmat([52,178,201], [Jurassic,1]); repmat([127,198,78], [Cretaceous,1]); repmat([253,154,82], [Paleogene,1]); repmat([255,230,25], [Neogene,1]); repmat([249,249,127], [Quaternary,1])];
% RGB colours on scale from 0 to 255
cmap = cmap / 255;
cmap = flip(cmap);

figure;
tiledlayout(4, 1, TileSpacing='none');
set(gcf, 'Position', [600, 400, 500, 400]);
nexttile;
plot(-t_data, diversity, color='black');
ylim([0,510]);
xlim([0,550]);
set(gca, 'xdir', 'reverse');
title("Peak");
ylabel("Diversity", Interpreter='latex');
xticks([]);
yticks([0,250,500]);
nexttile(2, [3,1]);
plot(-t_data, p_data, -t_sim, p_sim);
ylabel("$\log_{10}(\mathrm{size})$", Interpreter='latex');
xlabel("Time (mya)", Interpreter='latex');
ylim([1.1,2.3]);
% ylim([1.15,2.4]);
% ylim([0.8658,2.1158]);
xlim([0,550]);
% xlim([240, 550]);
% xlim([0, 260]);
legend('Data', 'Maximally likely model', Location='northwest', Interpreter='latex');
% xline(252, Color='black', LineStyle='--');
% legend('Data', 'Maximally likely model', 'PTME', Location='northwest', Interpreter='latex');
fontsize(gcf, 20, 'pixel');

ax = gca;
ax.Colormap = colormap(cmap);
cb = colorbar;
cb.Location = 'south';
cb.Position(1) = ax.Position(1);
cb.Position(2) = ax.Position(2);
cb.Position(3) = ax.Position(3);
cb.Ticks = [];
set(gca, 'xdir', 'reverse');
cb.Direction = 'reverse';

figure;
tiledlayout(4, 1, TileSpacing='none');
set(gcf, 'Position', [600, 400, 500, 400]);
nexttile;
plot(-t_data, diversity, color='black');
ylim([0,510]);
xlim([0,550]);
set(gca, 'xdir', 'reverse');
title("Mean");
ylabel("Diversity", Interpreter='latex');
xticks([]);
yticks([0,250,500]);
nexttile(2, [3,1]);
plot(-t_data, b_data, -t_sim, b_sim);
ylabel("$\log_{10}(\mathrm{size})$", Interpreter='latex');
xlabel("Time (mya)", Interpreter='latex');
ylim([1.1,2.3]);
% ylim([1.15,2.4]);
% ylim([0.8658,2.1158]);
xlim([0,550]);
% xlim([240, 550]);
% xlim([0, 260]);
legend('Data', 'Maximally likely model', Location='northwest', Interpreter='latex');
% xline(252, Color='black', LineStyle='--');
% legend('Data', 'Maximally likely model', 'PTME', Location='northwest', Interpreter='latex');
fontsize(gcf, 20, 'pixel');

ax = gca;
ax.Colormap = colormap(cmap);
cb = colorbar;
cb.Location = 'south';
cb.Position(1) = ax.Position(1);
cb.Position(2) = ax.Position(2);
cb.Position(3) = ax.Position(3);
cb.Ticks = [];
set(gca, 'xdir', 'reverse');
cb.Direction = 'reverse';
%%

% Plot solutions to the PDE for Dirichlet, Neumann and zero-flux
% (Robin type) upper boundaries on the same axis overlaid on
% the corresponding data (Brachiopod size data), for the date
% given below:
date = -145.5;
date_in_interval = logical((BrSum_fad <= date) .* (BrSum_lad > date));
histogram(BrSum_log_size(date_in_interval), Normalization="pdf", HandleVisibility='off');
hold on;

% Specify the maximally likely parameters for each solution
Br_OptVals_Dir = [0.1370, 0.1352, 1.281, 2.642, 0.6283, -606.9];
Br_OptVals_Neu = [9.273e-4, 6.380e-4, 1.404, 2.767, 0.5609, -606.9];
Br_OptVals_Rob = [1.089e-3, 6.444e-4, 1.430, 2.874, 0.4909, -606.9];

extract_times_Dir = linspace(0, -Br_OptVals_Dir(6), 200*-round(Br_OptVals_Dir(6)));
[~, xv_Dir, yv_Dir] = NumericalPDESolve(Br_OptVals_Dir(1), Br_OptVals_Dir(2), Br_OptVals_Dir(3), 0, Br_OptVals_Dir(4), extract_times_Dir, 200, false, "CrankNicolsonUpperA", false);
xv_Dir = xv_Dir + Br_OptVals_Dir(5);

extract_times_Neu = linspace(0, -Br_OptVals_Neu(6), 200*-round(Br_OptVals_Neu(6)));
[~, xv_Neu, yv_Neu] = NumericalPDESolve(Br_OptVals_Neu(1), Br_OptVals_Neu(2), Br_OptVals_Neu(3), 0, Br_OptVals_Neu(4), extract_times_Neu, 200, false, "CrankNicolsonUpperRneu", false);
xv_Neu = xv_Neu + Br_OptVals_Neu(5);

extract_times_Rob = linspace(0, -Br_OptVals_Rob(6), 200*-round(Br_OptVals_Rob(6)));
[~, xv_Rob, yv_Rob] = NumericalPDESolve(Br_OptVals_Rob(1), Br_OptVals_Rob(2), Br_OptVals_Rob(3), 0, Br_OptVals_Rob(4), extract_times_Rob, 200, false, "CrankNicolsonUpperRrob", false);
xv_Rob = xv_Rob + Br_OptVals_Rob(5);

ti_Dir = find(abs(extract_times_Dir - date + Br_OptVals_Dir(6)) == min(abs(extract_times_Dir - date + Br_OptVals_Dir(6))));
plot(xv_Dir, yv_Dir(ti_Dir, :) / NumericalIntegrator(xv_Dir, yv_Dir(ti_Dir, :)), "LineWidth", 2, "Color", "k", "LineStyle", "-", "DisplayName", "Dirichlet");
xline(max(xv_Dir), "LineWidth", 2, "Color", "k", "LineStyle", "-", HandleVisibility='off');

ti_Neu = find(abs(extract_times_Neu - date + Br_OptVals_Neu(6)) == min(abs(extract_times_Neu - date + Br_OptVals_Neu(6))));
plot(xv_Neu, yv_Neu(ti_Neu, :) / NumericalIntegrator(xv_Neu, yv_Neu(ti_Neu, :)), "LineWidth", 2, "Color", "r", "LineStyle", "--", "DisplayName", "Neumann");
xline(max(xv_Neu), "LineWidth", 2, "Color", "r", "LineStyle", "--", HandleVisibility='off');

ti_Rob = find(abs(extract_times_Rob - date + Br_OptVals_Rob(6)) == min(abs(extract_times_Rob - date + Br_OptVals_Rob(6))));
plot(xv_Rob, yv_Rob(ti_Rob, :) / NumericalIntegrator(xv_Rob, yv_Rob(ti_Rob, :)), "LineWidth", 2, "Color", "b", "LineStyle", "-.", "DisplayName", "Zero-flux");
xline(max(xv_Rob), "LineWidth", 2, "Color", "b", "LineStyle", "-.", HandleVisibility='off');

xlabel("$\log_{10}(\mathrm{size})$", Interpreter='latex');
ylabel("Density", Interpreter='latex');
ylim([0,1.25]);
xlim([0.4,3.4]);

set(gcf, 'Position', [600, 400, 600, 300]);
fontsize(gcf, 20, 'pixel');
legend(Interpreter='latex', Location='best');
