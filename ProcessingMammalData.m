% Import data from MOMv3.3.txt and extract the relevant columns
Data = readtable("MOMv3.3.txt");
continent = Data(:, "Var1");
status = Data(:, "Var2");
% mass = Data(:, "Var7");
mass = Data(:, "Var8");
status = table2cell(status);
continent = table2cell(continent);
mass = table2cell(mass);
mass = cell2mat(mass);

% Separate out the extant and extinct mammal data
extinct_mass = mass;
extant_mass = mass;
for i = 1:length(status)
    if strcmp(status(i), 'extinct')
        extant_mass(i) = NaN;
    elseif strcmp(status(i), 'extant')
        extinct_mass(i) = NaN;
    else
        extinct_mass(i) = NaN;
        extant_mass(i) = NaN;
    end
    if extinct_mass(i) == -999
        extinct_mass(i) = NaN;
    end
    if extant_mass(i) == -999
        extant_mass(i) = NaN;
    end
end
extinct_mass = extinct_mass(~isnan(extinct_mass));
extant_mass = extant_mass(~isnan(extant_mass));

% Take log of all mass values with minimum mass mapping to 0
extant_mass = log(extant_mass / min(extant_mass));

% Bin the data variously
[number, edges] = histcounts(extant_mass, 20, "Normalization", "pdf");
[number2, edges2] = histcounts(extant_mass, 50, "Normalization", "pdf");
[number3, edges3] = histcounts(extant_mass, 100, "Normalization", "pdf");
bins = edges(1:end-1) + (0.5 * (edges(2) - edges(1)));
bins2 = edges2(1:end-1) + (0.5 * (edges2(2) - edges2(1)));
bins3 = edges3(1:end-1) + (0.5 * (edges3(2) - edges3(1)));

x_range = linspace(0, max([max(edges), max(edges2), max(edges3)]), 1000);
%%

% Minimise the log likelihood function for the series solution to
% the partial differential equation:
% dy/dt = -a dy/dx + D/2 d^2y/dx^2 (y = y(x,t))
% subject to a delta function initial condition:
% y(x,0) = delta(x - x0)
% and Dirichlet boundaries at x = 0,L:
% y(0,t) = y(L,t) = 0
% in full and for the special case of a = 0 (purely diffusive case).
% Also minimise the log likelihood for a Gaussian function (solution
% in the case of no boundaries on the domain).

% Set the first appearance date of the single common ancestor
% species (million years ago).
timeRecent = 225;

x_ml2_i = [-0.1, 0.3, 20, 1, 0.01];
func2min2 = @(theta) -LogLikelihoodSSConAb(extant_mass + theta(5), theta(1), theta(2), theta(3), theta(4), timeRecent);
x_ml2 = fmincon(func2min2, x_ml2_i, [0,-1,0,0,0; 0,0,0,-1,0; 0,0,-1,1,0; 0,0,-1,0,1; 0,0,0,0,-1], [0, 0, 0, -max(extant_mass), 0]);

func2mina0 = @(theta) -LogLikelihoodSSConAb(extant_mass + theta(4), 0, theta(1), theta(2), theta(3), timeRecent);
x_mla0 = fmincon(func2mina0, [0.3, 20, 1, 0.01], [-1,0,0,0; 0,0,-1,0; 0,-1,1,0; 0,-1,0,1; 0,0,0,-1; 0,1,0,-1], [0, 0, 0, -max(extant_mass), 0, 1.5*max(extant_mass)]);

func2minG = @(theta) -LogLikelihoodGauss(extant_mass, theta(1), theta(2));
x_mlG = fminsearch(func2minG, [2.8, 4]);

% Calculate the log likelihood values
LL2 = -func2min2(x_ml2);
LLa0 = -func2mina0(x_mla0);
LLG = -func2minG(x_mlG);
%%

% Plot the data with the different maximally likely functions alongside
% comparrison of the cummulative distribution for each function.
figure;
set(gcf, 'Position', [600, 400, 600, 500]);
tiledlayout(4, 1, "TileSpacing", "none");
nexttile([3, 1]);
lw = 2;
bar(bins2, number2, 'DisplayName', 'Mammal mass');
xlim([0, max(extant_mass)]);
legend(Interpreter='latex', Location='northeast');
hold on;
plot(x_range, SSConAb(x_range + x_ml2(5), x_ml2(1), x_ml2(2), x_ml2(3), x_ml2(4), timeRecent), 'DisplayName', 'Dirichlet boundaries', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'red');
plot(x_range, SSConAb(x_range + x_mla0(4), 0, x_mla0(1), x_mla0(2), x_mla0(3), timeRecent), 'DisplayName', 'Dirichlet boundaries, $a=0$', 'LineWidth', lw, 'LineStyle', '-.', 'Color', 'blue');
Gaussian = @(x, xdata) (1 / (x(2) * sqrt(2 * pi))) * exp(-((xdata - x(1)) .* (xdata - x(1))) / (2 * x(2) * x(2)));
plot(x_range, Gaussian(x_mlG, x_range), 'DisplayName', "No boundaries", 'LineWidth', lw, 'Color', 'black', 'LineStyle', "--");
xticks([]);
yticks([0.05, 0.1, 0.15, 0.2]);
ylabel('Density', Interpreter='latex');

cummulative_data = CummulativeFunction(extant_mass, extant_mass);
cummulative_2 = CummulativeFunction(@(x) SSConAb(x + x_ml2(5), x_ml2(1), x_ml2(2), x_ml2(3), x_ml2(4), timeRecent), extant_mass, -x_ml2(5));
cummulative_a0 = CummulativeFunction(@(x) SSConAb(x + x_mla0(4), 0, x_mla0(1), x_mla0(2), x_mla0(3), timeRecent), extant_mass, -x_mla0(4));
cummulative_G = CummulativeFunction(@(x) Gaussian(x_mlG, x), extant_mass, -20);
L2_2 = sqrt(sum((cummulative_data - cummulative_2).^2));
L2_a0 = sqrt(sum((cummulative_data - cummulative_a0).^2));
L2_G = sqrt(sum((cummulative_data - cummulative_G).^2));

nexttile;
masses = linspace(0, 20, 1000);
cummulative_data = CummulativeFunction(extant_mass, masses);
bar(masses, cummulative_data, 'HandleVisibility', 'off');
hold on;
cummulative_2 = CummulativeFunction(@(x) SSConAb(x + x_ml2(5), x_ml2(1), x_ml2(2), x_ml2(3), x_ml2(4), timeRecent), masses, -x_ml2(5));
cummulative_a0 = CummulativeFunction(@(x) SSConAb(x + x_mla0(4), 0, x_mla0(1), x_mla0(2), x_mla0(3), timeRecent), masses, -x_mla0(4));
cummulative_G = CummulativeFunction(@(x) Gaussian(x_mlG, x), masses, -20);
plot(masses, cummulative_2, 'DisplayName', "$l^2$-error="+L2_2, 'LineWidth', lw, 'LineStyle', '-', 'Color', 'red');
plot(masses, cummulative_a0, 'DisplayName', "$l^2$-error="+L2_a0, 'LineWidth', lw, 'LineStyle', '-.', 'Color', 'blue');
plot(masses, cummulative_G, 'DisplayName', "$l^2$-error="+L2_G, 'LineWidth', lw, 'Color', 'black', 'LineStyle', "--");
fontsize(gcf, 20, 'pixel');
legend(Interpreter='latex', Location='southeast');
xlim([0, max(extant_mass)]);
ylim([0, 1.1]);
xlabel("$x$", Interpreter='latex');
%%

% Use bootstap resampling to obtain confidence intervals on the maximally
% likely parameters.
bootstraps = 1000;
x_ml_full = zeros(bootstraps, 5);
LL2_full = zeros(1, bootstraps);
i = 1;
warnCount = 0;
while i < bootstraps
    lastwarn('');
    disp("Bootsrap " + i + " of " + bootstraps);
    randis = randi([1, length(extant_mass)], 1, length(extant_mass));
    func2min2 = @(theta) -LogLikelihoodSSConAb(extant_mass(randis) + theta(5), theta(1), theta(2), theta(3), theta(4), timeRecent);
    x_ml_full(i, :) = fmincon(func2min2, x_ml2_i, [0,-1,0,0,0; 0,0,0,-1,0; 0,0,-1,1,0; 0,0,-1,0,1; 0,0,0,0,-1], [0, 0, 0, -max(extant_mass), 0]);
    LL2_full(i) = -func2min2(x_ml_full(i, :));
    [warnMsg, warnId] = lastwarn;
    if isempty(warnMsg)
        i = i + 1;
    else
        warnCount = warnCount + 1;
    end
end

% Calcuate mean parameter values from bootstrapping process
% (not usually a good parameter estimate)
x_mlave = sum(x_ml_full) / bootstraps;

x_mla = sort(x_ml_full(:, 1), 'ascend');
x_mlD = sort(x_ml_full(:, 2), 'ascend');
x_mlL = sort(x_ml_full(:, 3), 'ascend');
x_mlx0 = sort(x_ml_full(:, 4), 'ascend');
x_mlmsh = sort(x_ml_full(:, 5), 'ascend');
LL2_full = sort(LL2_full, 'ascend');

% Calculate median parameter values from bootsrapping process
if round(bootstraps / 2) == (bootstraps / 2)
    med1 = bootstraps / 2;
    med2 = med1 + 1;
    x_mlmed = zeros(1, length(x_ml2));
    x_mlmed(1) = 0.5 * (x_mla(med1) + x_mla(med2));
    x_mlmed(2) = 0.5 * (x_mlD(med1) + x_mlD(med2));
    x_mlmed(3) = 0.5 * (x_mlL(med1) + x_mlL(med2));
    x_mlmed(4) = 0.5 * (x_mlx0(med1) + x_mlx0(med2));
    x_mlmed(5) = 0.5 * (x_mlmsh(med1) + x_mlmsh(med2));
else
    med = ceil(bootstraps / 2);
    x_mlmed = zeros(1, length(x_ml2));
    x_mlmed(1) = x_mla(med);
    x_mlmed(2) = x_mlD(med);
    x_mlmed(3) = x_mlL(med);
    x_mlmed(4) = x_mlx0(med);
    x_mlmed(5) = x_mlmsh(med);
end

standardGaussian = @(z) exp(-(z .* z) / 2) / sqrt(2 * pi);
confidence_area = @(ci) NumericalIntegrator(linspace(-ci, ci, 1000000), standardGaussian(linspace(-ci, ci, 1000000)));

% Number standard_deviation gives number of sigma confidence interval to calculate.
standard_deviation = 1;
l1s = floor(bootstraps * (0.5 - (confidence_area(standard_deviation) / 2)));
u1s = ceil(bootstraps * (0.5 + (confidence_area(standard_deviation) / 2)));

disp("Basic CI method:");
disp("a = " + x_ml2(1) + "(" + (x_mla(l1s)) + "," + (x_mla(u1s)) + ")");
disp("D = " + x_ml2(2) + "(" + (x_mlD(l1s)) + "," + (x_mlD(u1s)) + ")");
disp("L = " + x_ml2(3) + "(" + (x_mlL(l1s)) + "," + (x_mlL(u1s)) + ")");
disp("x0 = " + x_ml2(4) + "(" + (x_mlx0(l1s)) + "," + (x_mlx0(u1s)) + ")");
disp("m_sh = " + x_ml2(5) + "(" + (x_mlmsh(l1s)) + "," + (x_mlmsh(u1s)) + ")");
disp("LogLikelihood2 = " + LL2 + "(" + (LL2_full(l1s)) + "," + (LL2_full(u1s)) + ")");
%%

% Plot a likelihood surface in a-D space over the 2-sigma confidence range
a_m2s = -90.9;
a_p2s = -0.0739;
D_m2s = 0.248;
D_p2s = 259;

a = linspace(a_m2s, a_p2s, 100);
D = linspace(D_m2s, D_p2s, 100);
x_ml2_ired = x_ml2_i(3:5);
x_ml_ad = zeros(length(D), length(a));
k = 1;

for i = 1:length(D)
    for j = 1:length(a)
        disp(k+" of "+(length(a)*length(D)));
        func2min2 = @(theta) -LogLikelihoodSSConAb(extant_mass + theta(3), a(j), D(i), theta(1), theta(2), timeRecent);
        try
            x_ml = fmincon(func2min2, x_ml2_ired, [0,-1,0; -1,1,0; -1,0,1; 0,0,-1; 1,0,-1], [0, 0, -max(extant_mass), 0, 1.5*max(extant_mass)]);
            x_ml_ad(j, i) = -func2min2(x_ml);
        catch
            x_ml_ad(j, i) = NaN;
        end
        k = k + 1;
    end
end
%%

% Plot the likelihood surface with contours specified by levels
% ratio is used to plot the line of constant a/D - set values
% according to maximum log likelihood analysis.  Position of
% maximaum likelihood is also plotted.
levels = -11000:200:-9800;
levels(end+1:end+3) = [-9700, -9600, -9500];

figure;
contourf(D, a, x_ml_ad, levels, ShowText='on', LabelSpacing=700);
fontsize(gcf, 20, 'pixel');
set(gcf, 'Position', [600, 400, 600, 300]);
xlabel("$D$", Interpreter='latex');
ylabel("$a$", Interpreter='latex');
ratio = -0.137 / 0.406;
hold on;
plot(D, ratio * D, Color='r', LineStyle='--');
xlim([min(D), max(D)]);
ylim([min(a), max(a)]);
plot(0.406, -0.137, Marker='o', Color='r', MarkerFaceColor='r');
colorbar;
