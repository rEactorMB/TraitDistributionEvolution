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

% Take log of all mass values with minimum value mapping to 0
extant_mass = log(extant_mass / min(extant_mass));

[number, edges] = histcounts(extant_mass, 20, "Normalization", "pdf");
[number2, edges2] = histcounts(extant_mass, 50, "Normalization", "pdf");
[number3, edges3] = histcounts(extant_mass, 100, "Normalization", "pdf");
bins = edges(1:end-1) + (0.5 * (edges(2) - edges(1)));
bins2 = edges2(1:end-1) + (0.5 * (edges2(2) - edges2(1)));
bins3 = edges3(1:end-1) + (0.5 * (edges3(2) - edges3(1)));

x_range = linspace(0, max([max(edges), max(edges2), max(edges3)]), 1000);
%%

f = figure;
f.Position = [400, 100, 1000, 800];
lw = 3;
bar(bins2, number2, 'DisplayName', 'Data');
xlim([0, max(edges)]);
legend;
hold on;
[KDFx, KDFy] = KernelDensityFunction(bins, number, true);
% pause(1);
% plot(KDFx, KDFy, 'DisplayName', 'KDF', 'LineWidth', lw, 'Color', 'k');

Fabsorbing = @(x, xdata) LTConAb(xdata, x(1), x(2), x(3));
Fneumann = @(x, xdata) LTConNeu(xdata, x(1), x(2), x(3));
Frobin = @(x, xdata) LTConRob(xdata, x(1), x(2), x(3));
FabsorbingSS = @(x, xdata) SSConAb(xdata, x(1), x(2), x(3), x(4), x(5));

% pause(1);
% x_fit = lsqcurvefit(Fabsorbing, [-1, 1, 10], KDFx, KDFy, [-1000, 0, max(KDFx)]);
% % x_fit = lsqcurvefit(Fabsorbing, [-1, 1, 10], bins, Nnumber, [-1000, 0, max(edges)]);
% plot(x_range, Fabsorbing(x_fit, x_range), 'DisplayName', 'Absorbing b.c.', 'LineWidth', lw, 'LineStyle', '--', 'Color', 'r');
% disp("[a, D, L] = [" + x_fit(1) + ", " + x_fit(2) + ", " + x_fit(3) + ']');

% pause(1);
% % x_fit = lsqcurvefit(Fneumann, [-1, 1, 10], KDFx, KDFy, [-1000, 0, 0]);
% x_fit = lsqcurvefit(Fneumann, [-1, 1, 10], bins, Nnumber, [-1000, 0, 0]);
% plot(x_range, Fneumann(x_fit, x_range), 'DisplayName', 'Neumann b.c.', 'LineWidth', lw, 'LineStyle', '-.', 'Color', 'g');
% disp("[a, D, L] = [" + x_fit(1) + ", " + x_fit(2) + ", " + x_fit(3) + ']');

% pause(1);
% % x_fit = lsqcurvefit(Frobin, [-10, 0.1, 10], KDFx, KDFy, [-1000, 0, 0]);
% x_fit = lsqcurvefit(Frobin, [-1, 0.1, 10], bins, Nnumber, [-1000, 0, 0]);
% plot(x_range, Frobin(x_fit, x_range), 'DisplayName', 'Robin b.c.', 'LineWidth', lw, 'LineStyle', ':', 'Color', 'm');
% disp("[a, D, L] = [" + x_fit(1) + ", " + x_fit(2) + ", " + x_fit(3) + ']');

% pause(1);
% x_fit = lsqcurvefit(FabsorbingSS, [-10, 1, 10, 1, 20], KDFx, KDFy, [-1000, 0, max(edges), 0, 0]);
% x_fit = lsqcurvefit(FabsorbingSS, [-1, 1, 10, 1, 20], bins, number, [-1000, 0, max(edges), 0, 0]);
% plot(x_range, FabsorbingSS(x_fit, x_range), 'DisplayName', 'Absorbing b.c. - SS', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'red');
% disp("[a, D, L, x0, t] = [" + x_fit(1) + ", " + x_fit(2) + ", " + x_fit(3) + ", " + x_fit(4) + ", " + x_fit(5) + ']');

% x_fit2 = lsqcurvefit(FabsorbingSS, [-1, 1, 10, 1, 20], bins2, number2, [-1000, 0, max(edges2), 0, 0]);
% plot(x_range, FabsorbingSS(x_fit2, x_range), 'DisplayName', 'Absorbing b.c. - SS (50)', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'green');
% disp("[a, D, L, x0, t] = [" + x_fit2(1) + ", " + x_fit2(2) + ", " + x_fit2(3) + ", " + x_fit2(4) + ", " + x_fit2(5) + ']');

% x_fit3 = lsqcurvefit(FabsorbingSS, [-1, 1, 10, 1, 20], bins3, number3, [-1000, 0, max(edges3), 0, 0]);
% plot(x_range, FabsorbingSS(x_fit3, x_range), 'DisplayName', 'Absorbing b.c. - SS (100)', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'black');
% disp("[a, D, L, x0, t] = [" + x_fit3(1) + ", " + x_fit3(2) + ", " + x_fit3(3) + ", " + x_fit3(4) + ", " + x_fit3(5) + ']');

FabsorbingSSt1 = @(x, xdata) SSConAb(xdata, x(1), x(2), x(3), x(4), 1);
[x_fitt1, resnorm_t1, residuals_t1, ~, ~, ~, jac_t1] = lsqcurvefit(FabsorbingSSt1, [-10, 100, 20, 1], bins2, number2, [-1000, 0, max(edges2), 0]);
plot(x_range, FabsorbingSSt1(x_fitt1, x_range), 'DisplayName', 'Absorbing b.c. - SS (50)', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'green');
disp("[a, D, L, x0] = [" + x_fitt1(1) + ", " + x_fitt1(2) + ", " + x_fitt1(3) + ", " + x_fitt1(4) + ']');
conf_t1 = nlparci(x_fitt1, residuals_t1, 'jacobian', jac_t1);
% cov_2 = ((jac_2' * jac_2) \ ones(length(x_fit2))) * ((residuals_2 * residuals_2') / (length(number2) - length(x_fit2)));
% conf_2 = sqrt(diag(cov_2));

% [beta, R, J, CovB, MSE, ~] = nlinfit(bins2, number2, FabsorbingSS, x_fit2);
% ci = nlparci(beta, R, covariance=CovB);

%%
hold off;

xlabel("Natural logarithm of body mass");
ylabel("Density");

figure;
plot(bins2, residuals_2, LineStyle=':', Marker='x', Color="r");
yline(0, Color='black');
title("Residuals sum = " + sum(abs(residuals_2)));

%%
Gaussian = @(x, xdata) (1 / (x(2) * sqrt(2 * pi))) * exp(-((xdata - x(1)) .* (xdata - x(1))) / (2 * x(2) * x(2)));
[x_fitG, resnorm_G, residuals_G, ~, ~, ~, jac_G] = lsqcurvefit(Gaussian, [2.7, 1.5], bins2, number2);

f = figure;
f.Position = [400, 100, 1000, 800];
lw = 3;
bar(bins2, number2, 'DisplayName', 'Data');
xlim([0, max(edges)]);
legend;
hold on;
plot(x_range, Gaussian(x_fitG, x_range), DisplayName="Gaussian", LineWidth=lw, color='black');

FabsorbingSSa0 = @(x, xdata) SSConAb(xdata, 0, x(1), x(2), x(3), x(4));
[x_fita0, resnorm_a0, residuals_a0, ~, ~, ~, jac_a0] = lsqcurvefit(FabsorbingSSa0, [1, 10, 3.1, 20], bins2, number2, [0, max(edges2), 0, 0]);
disp("[D, L, x0, t] = [" + x_fita0(1) + ", " + x_fita0(2) + ", " + x_fita0(3) + ", " + x_fita0(4) + "]");

figure;
set(gcf, 'Position', [600, 400, 600, 500]);
tiledlayout(3, 1, TileSpacing='none');
nexttile([2, 1]);
bar(bins2, number2, 'DisplayName', 'Mammal mass');
xlim([0, max(edges)]);
legend(Interpreter='latex', Location='northeast');
hold on;
plot(x_range, FabsorbingSS(x_fit2, x_range), 'DisplayName', 'Dirichlet boundaries', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'red');
plot(x_range, FabsorbingSSa0(x_fita0, x_range), 'DisplayName', 'Dirichlet boundaries, $a=0$', 'LineWidth', lw, 'LineStyle', '-.', 'Color', 'blue');
plot(x_range, Gaussian(x_fitG, x_range), 'DisplayName', "No boundaries", 'LineWidth', lw, 'Color', 'black', 'LineStyle', "--");

hold off;
xticks([]);
ylabel("Density", Interpreter='latex');
fontsize(gcf, 18, 'pixel');

nexttile;
xlim([0, max(edges)]);
legend(Interpreter='latex', BackgroundAlpha=0.5, Location='northeast');
hold on;
plot(bins2, residuals_2, 'DisplayName', "sum="+resnorm_2, 'LineWidth', lw, 'LineStyle', '-', 'Color', 'red');
plot(bins2, residuals_a0, 'DisplayName', "sum="+resnorm_a0, 'LineWidth', lw, 'LineStyle', '-.', 'Color', 'blue');
plot(bins2, residuals_G, 'DisplayName', "sum="+resnorm_G, 'LineWidth', lw, 'Color', 'black', 'LineStyle', "--");
hold off;
xlabel("$\ln(\mathrm{mass}/1.8\mathrm{~g})$", Interpreter='latex');
ylabel("Residuals", Interpreter='latex');
ylim([-0.05, 0.08]);
grid on;
grid minor;
fontsize(gcf, 20, 'pixel');
%%
f = figure;
f.Position = [400, 100, 1000, 800];
lw = 3;
bar(bins2, number2, 'DisplayName', 'Data');
xlim([0, max(edges)]);
legend;
hold on;

x0 = log(2/1.8);

FabsorbingSS = @(x, xdata) SSConAb(xdata, x(1), x(2), x(3), x0, x(4));

x_range = linspace(0, max([max(edges), max(edges2), max(edges3)]), 1000);

pause(1);
x_fit = lsqcurvefit(FabsorbingSS, [-1, 1, 10, 20], bins, number, [-1000, 0, max(edges), 0]);
plot(x_range, FabsorbingSS(x_fit, x_range), 'DisplayName', 'Absorbing b.c. - SS', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'red');
disp("[a, D, L, x0, t] = [" + x_fit(1) + ", " + x_fit(2) + ", " + x_fit(3) + ", " + x0 + ", " + x_fit(4) + ']');

x_fit2 = lsqcurvefit(FabsorbingSS, [-1, 1, 10, 20], bins2, number2, [-1000, 0, max(edges2), 0]);
plot(x_range, FabsorbingSS(x_fit2, x_range), 'DisplayName', 'Absorbing b.c. - SS (50)', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'green');
disp("[a, D, L, x0, t] = [" + x_fit2(1) + ", " + x_fit2(2) + ", " + x_fit2(3) + ", " + x0 + ", " + x_fit2(4) + ']');

x_fit3 = lsqcurvefit(FabsorbingSS, [-1, 1, 10, 20], bins3, number3, [-1000, 0, max(edges3), 0]);
plot(x_range, FabsorbingSS(x_fit3, x_range), 'DisplayName', 'Absorbing b.c. - SS (100)', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'black');
disp("[a, D, L, x0, t] = [" + x_fit3(1) + ", " + x_fit3(2) + ", " + x_fit3(3) + ", " + x0 + ", " + x_fit3(4) + ']');

hold off;

xlabel("Natural logarithm of body mass");
ylabel("Density");

%%
f = figure;
f.Position = [400, 100, 1000, 800];
lw = 3;
bar(bins3, number3, 'DisplayName', 'Data');
xlim([0, max(edges)]);
legend;
hold on;

x0 = log(2);
t = 400;

FabsorbingSS = @(x, xdata) SSConAb(xdata, x(1), x(2), x(3), x0, t);

x_range = linspace(0, max([max(edges), max(edges2), max(edges3)]), 1000);

pause(1);
x_fit = lsqcurvefit(FabsorbingSS, [-1, 1, 10], bins, number, [-1000, 0, max(edges)]);
plot(x_range, FabsorbingSS(x_fit, x_range), 'DisplayName', 'Absorbing b.c. - SS', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'red');
disp("[a, D, L, x0, t] = [" + x_fit(1) + ", " + x_fit(2) + ", " + x_fit(3) + ", " + x0 + ", " + t + ']');

x_fit2 = lsqcurvefit(FabsorbingSS, [-1, 1, 10], bins2, number2, [-1000, 0, max(edges2)]);
plot(x_range, FabsorbingSS(x_fit2, x_range), 'DisplayName', 'Absorbing b.c. - SS (50)', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'green');
disp("[a, D, L, x0, t] = [" + x_fit2(1) + ", " + x_fit2(2) + ", " + x_fit2(3) + ", " + x0 + ", " + t + ']');

x_fit3 = lsqcurvefit(FabsorbingSS, [-1, 1, 10], bins3, number3, [-1000, 0, max(edges3)]);
plot(x_range, FabsorbingSS(x_fit3, x_range), 'DisplayName', 'Absorbing b.c. - SS (100)', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'black');
disp("[a, D, L, x0, t] = [" + x_fit3(1) + ", " + x_fit3(2) + ", " + x_fit3(3) + ", " + x0 + ", " + t + ']');

hold off;

xlabel("Natural logarithm of body mass");
ylabel("Density");

%%
timeRecent = 225;

x_ml2_i = [-0.1, 0.3, 20, 1, 0.01];
func2min2 = @(theta) -LogLikelihoodSSConAb(extant_mass + theta(5), theta(1), theta(2), theta(3), theta(4), timeRecent);
x_ml2 = fmincon(func2min2, x_ml2_i, [0,-1,0,0,0; 0,0,0,-1,0; 0,0,-1,1,0; 0,0,-1,0,1; 0,0,0,0,-1], [0, 0, 0, -max(extant_mass), 0]);

func2mina0 = @(theta) -LogLikelihoodSSConAb(extant_mass + theta(4), 0, theta(1), theta(2), theta(3), timeRecent);
x_mla0 = fmincon(func2mina0, [0.3, 20, 1, 0.01], [-1,0,0,0; 0,0,-1,0; 0,-1,1,0; 0,-1,0,1; 0,0,0,-1; 0,1,0,-1], [0, 0, 0, -max(extant_mass), 0, 1.5*max(extant_mass)]);

func2minG = @(theta) -LogLikelihoodGauss(extant_mass, theta(1), theta(2));
x_mlG = fminsearch(func2minG, [2.8, 4]);

disp("[a, D, L, x0, x_sh] = [" + x_ml2(1) + ", " + x_ml2(2) + ", " + x_ml2(3) + ", " + x_ml2(4) + ', ' + x_ml2(5) + ']');
disp("[D, L, x0, x_sh] = [" + x_mla0(1) + ", " + x_mla0(2) + ", " + x_mla0(3) + ", " + x_mla0(4) + ']');
disp("[mean, sigma] = [" + x_mlG(1) + ", " + x_mlG(2) + ']');

figure;
set(gcf, 'Position', [600, 400, 600, 300]);
lw = 3;
bar(bins2, number2, 'DisplayName', 'Mammal mass');
xlim([0, max(extant_mass)]);
legend(Interpreter='latex', Location='northeast');
hold on;
plot(x_range, SSConAb(x_range + x_ml2(5), x_ml2(1), x_ml2(2), x_ml2(3), x_ml2(4), timeRecent), 'DisplayName', 'Dirichlet boundaries', 'LineWidth', lw, 'LineStyle', '-', 'Color', 'red');
plot(x_range, SSConAb(x_range + x_mla0(4), 0, x_mla0(1), x_mla0(2), x_mla0(3), timeRecent), 'DisplayName', 'Dirichlet boundaries, $a=0$', 'LineWidth', lw, 'LineStyle', '-.', 'Color', 'blue');
Gaussian = @(x, xdata) (1 / (x(2) * sqrt(2 * pi))) * exp(-((xdata - x(1)) .* (xdata - x(1))) / (2 * x(2) * x(2)));
plot(x_range, Gaussian(x_mlG, x_range), 'DisplayName', "No boundaries", 'LineWidth', lw, 'Color', 'black', 'LineStyle', "--");
fontsize(gcf, 20, 'pixel');

LL2 = -func2min2(x_ml2);
LLa0 = -func2mina0(x_mla0);
LLG = -func2minG(x_mlG);
disp("LogLikelihood2 = " + LL2);
disp("LogLikelihooda0 = " + LLa0);
disp("LogLikelihoodG = " + LLG);

%%
% Bootstrap confidence intervals
bootstraps = 1000;
x_ml_full = zeros(bootstraps, 5);
LL2_full = zeros(1, bootstraps);
% for i = 1:bootstraps
%     disp("Bootsrap " + i + " of " + bootstraps);
%     randis = randi([1, length(extant_mass)], 1, length(extant_mass));
%     func2min2 = @(theta) -LogLikelihoodSSConAb(extant_mass(randis) + theta(5), theta(1), theta(2), theta(3), theta(4), timeRecent);
%     x_ml_full(i, :) = fmincon(func2min2, x_ml2_i, [0,-1,0,0,0; 0,0,0,-1,0; 0,0,-1,1,0; 0,0,-1,0,1; 0,0,0,0,-1], [0, 0, 0, -max(extant_mass), 0]);
%     LL2_full(i) = -func2min2(x_ml_full(i, :));
%     histogram(extant_mass(randis), 50, Normalization='pdf');
%     hold on;
%     plot(x_range, SSConAb(x_range + x_ml_full(i, 5), x_ml_full(i, 1), x_ml_full(i, 2), x_ml_full(i, 3), x_ml_full(i, 4), timeRecent));
%     hold off;
%     pause(0.00001);
% end

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
        % histogram(extant_mass(randis), 50, Normalization='pdf');
        % hold on;
        % plot(x_range, SSConAb(x_range + x_ml_full(i, 5), x_ml_full(i, 1), x_ml_full(i, 2), x_ml_full(i, 3), x_ml_full(i, 4), timeRecent));
        % hold off;
        % pause(0.00001);
        % cummulative_data_randi = CummulativeFunction(extant_mass(randis), extant_mass);
        % cummulative_full = CummulativeFunction(@(x) SSConAb(x + x_ml_full(i, 5), x_ml_full(i, 1), x_ml_full(i, 2), x_ml_full(i, 3), x_ml_full(i, 4), timeRecent), extant_mass, -x_ml_full(i, 5));
        % L2_full = sqrt(sum((cummulative_data_randi - cummulative_full).^2));
        % disp('l^2-norm = ' + string(L2_full));
        i = i + 1;
    else
        warnCount = warnCount + 1;
    end
end

x_mlave = sum(x_ml_full) / bootstraps;
disp("[a, D, L, x0, m_sh] = [" + x_mlave(1) + ", " + x_mlave(2) + ", " + x_mlave(3) + ", " + x_mlave(4) + ", " + x_mlave(5) + ']');

x_mla = sort(x_ml_full(:, 1), 'ascend');
x_mlD = sort(x_ml_full(:, 2), 'ascend');
x_mlL = sort(x_ml_full(:, 3), 'ascend');
x_mlx0 = sort(x_ml_full(:, 4), 'ascend');
x_mlmsh = sort(x_ml_full(:, 5), 'ascend');
LL2_full = sort(LL2_full, 'ascend');

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
disp("[a, D, L, x0, m_sh] = [" + x_mlmed(1) + ", " + x_mlmed(2) + ", " + x_mlmed(3) + ", " + x_mlmed(4) + ", " + x_mlmed(5) + ']');
%%
standardGaussian = @(z) exp(-(z .* z) / 2) / sqrt(2 * pi);
confidence_area = @(ci) NumericalIntegrator(linspace(-ci, ci, 1000000), standardGaussian(linspace(-ci, ci, 1000000)));

% Number in confidence area gives number of sigma confidence interval.
standard_deviation = 1;
l1s = floor(bootstraps * (0.5 - (confidence_area(standard_deviation) / 2)));
u1s = ceil(bootstraps * (0.5 + (confidence_area(standard_deviation) / 2)));

disp("Pivot CI method:");
disp("a = " + x_ml2(1) + "(" + (2*x_ml2(1) - x_mla(u1s)) + "," + (2*x_ml2(1) - x_mla(l1s)) + ")");
disp("D = " + x_ml2(2) + "(" + (2*x_ml2(2) - x_mlD(u1s)) + "," + (2*x_ml2(2) - x_mlD(l1s)) + ")");
disp("L = " + x_ml2(3) + "(" + (2*x_ml2(3) - x_mlL(u1s)) + "," + (2*x_ml2(3) - x_mlL(l1s)) + ")");
disp("x0 = " + x_ml2(4) + "(" + (2*x_ml2(4) - x_mlx0(u1s)) + "," + (2*x_ml2(4) - x_mlx0(l1s)) + ")");
disp("m_sh = " + x_ml2(5) + "(" + (2*x_ml2(5) - x_mlmsh(u1s)) + "," + (2*x_ml2(5) - x_mlmsh(l1s)) + ")");
disp("LogLikelihood2 = " + LL2 + "(" + (2*LL2 - LL2_full(u1s)) + "," + (2*LL2 - LL2_full(l1s)) + ")");

disp("Basic CI method:");
disp("a = " + x_ml2(1) + "(" + (x_mla(l1s)) + "," + (x_mla(u1s)) + ")");
disp("D = " + x_ml2(2) + "(" + (x_mlD(l1s)) + "," + (x_mlD(u1s)) + ")");
disp("L = " + x_ml2(3) + "(" + (x_mlL(l1s)) + "," + (x_mlL(u1s)) + ")");
disp("x0 = " + x_ml2(4) + "(" + (x_mlx0(l1s)) + "," + (x_mlx0(u1s)) + ")");
disp("m_sh = " + x_ml2(5) + "(" + (x_mlmsh(l1s)) + "," + (x_mlmsh(u1s)) + ")");
disp("LogLikelihood2 = " + LL2 + "(" + (LL2_full(l1s)) + "," + (LL2_full(u1s)) + ")");
%%
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
disp("L2_2 = " + L2_2);
disp("L2_a0 = " + L2_a0);
disp("L2_G = " + L2_G);

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

% nexttile;
% masses = linspace(0, 20, 1000);
% cummulative_data = CummulativeFunction(extant_mass, masses);
% % bar(masses, cummulative_data, 'HandleVisibility', 'off');
% hold on;
% cummulative_2 = CummulativeFunction(@(x) SSConAb(x + x_ml2(5), x_ml2(1), x_ml2(2), x_ml2(3), x_ml2(4), timeRecent), masses, -x_ml2(5));
% cummulative_a0 = CummulativeFunction(@(x) SSConAb(x + x_mla0(4), 0, x_mla0(1), x_mla0(2), x_mla0(3), timeRecent), masses, -x_mla0(4));
% cummulative_G = CummulativeFunction(@(x) Gaussian(x_mlG, x), masses, -20);
% plot(masses, cummulative_2 - cummulative_data, 'DisplayName', "$l^2=$"+L2_2, 'LineWidth', lw, 'LineStyle', '-', 'Color', 'red');
% plot(masses, cummulative_a0 - cummulative_data, 'DisplayName', "$l^2=$"+L2_a0, 'LineWidth', lw, 'LineStyle', '-.', 'Color', 'blue');
% plot(masses, cummulative_G - cummulative_data, 'DisplayName', "$l^2=$"+L2_G, 'LineWidth', lw, 'Color', 'black', 'LineStyle', "--");
% fontsize(gcf, 20, 'pixel');
% legend(Interpreter='latex', Location='southeast');
% xlim([0, max(extant_mass)]);
% ylim([-0.125, 0.1]);

%%
a_m2s = -90.9;
a_p2s = -0.0739;
D_m2s = 0.248;
D_p2s = 259;

% a = linspace(a_m2s, a_p2s, 100);
% D = linspace(D_m2s, D_p2s, 100);
a = linspace(0, -0.25, 50);
D = linspace(0.2, 1, 50);
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
a_m2s = -90.9;
a_p2s = -0.0739;
D_m2s = 0.248;
D_p2s = 259;

% a = linspace(a_m2s, a_p2s, 50);
% D = linspace(D_m2s, D_p2s, 50);
a = linspace(0, -0.25, 50);
D = linspace(0.2, 1, 50);
x_ml2_ired = x_ml2_i([3,5]);
x_ml_ad = zeros(length(D), length(a));
k = 1;

x0 = 3.52;

for i = 1:length(D)
    for j = 1:length(a)
        disp(k+" of "+(length(a)*length(D)));
        func2min2 = @(theta) -LogLikelihoodSSConAb(extant_mass + theta(2), a(j), D(i), theta(1), x0, timeRecent);
        try
            x_ml = fmincon(func2min2, x_ml2_ired, [-1,0; -1,1; 0,-1; 1,-1], [-x0, -max(extant_mass), 0, 1.5*max(extant_mass)]);
            x_ml_ad(j, i) = -func2min2(x_ml);
        catch
            x_ml_ad(j, i) = NaN;
        end
        k = k + 1;
    end
end

%%
levels = -11000:200:-9800;
% levels(end+1:end+3) = [-9700, -9600, -9500];
levels(end+1:end+11) = -9700:20:-9500;
% figure;
% contour(D, a, x_ml_ad, levels, ShowText='off');
% fontsize(gcf, 20, 'pixel');
% set(gcf, 'Position', [600, 400, 600, 300]);
% xlabel("$D$", Interpreter='latex');
% ylabel("$a$", Interpreter='latex');
% % xlim([D_m2s, D_p2s]);
% % ylim([a_m2s, a_p2s]);
% hold on;
% plot(0.406, -0.137, Marker='o', Color='r', MarkerFaceColor='r');
% colorbar;

figure;
contourf(D, a, x_ml_ad, levels, ShowText='on', LabelSpacing=700);
% contourf(D, a, x_ml_ad, ShowText='on', LabelSpacing=700);
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
% colorbar;

%%
a_m2s = -90.9;
a_p2s = -0.0739;
D_m2s = 0.248;
D_p2s = 259;

% a = linspace(a_m2s, a_p2s, 10);
% D = linspace(D_m2s, D_p2s, 10);
a = linspace(-0.2, 0, 50);
D = linspace(0.1, 10, 50);
x_ml2_ired = x_ml2_i(3:5);
x_ml_ad = zeros(length(D), length(a));
k = 1;

for i = 1:length(D)
    for j = 1:length(a)
        disp(k+" of "+(length(a)*length(D)));
        x_ml_ad(j, i) = LogLikelihoodSSConAb(extant_mass + x_ml2(5), a(j), D(i), x_ml2(3), x_ml2(4), timeRecent);
        k = k + 1;
    end
end
%%
levels = -15000:100:-9500;
figure;
contour(D, a, x_ml_ad, levels, ShowText='off');
figure;
contourf(D, a, x_ml_ad, levels, ShowText='off');

%%
for i = 1:length(D)
    for j = 1:length(a)
        if (i + j < 20) || (i + j > 80)
            x_ml_ad(j, i) = NaN;
        end
    end
end
