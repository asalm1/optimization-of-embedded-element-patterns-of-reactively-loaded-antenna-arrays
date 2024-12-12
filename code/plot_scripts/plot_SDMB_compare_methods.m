%%
% Plot SDMB results comparing different optimization methods
% ------------------------------------------------------------------------
% 06.06.2024 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%% Clear
clear
clc
close all

%% Load data
load '.\results\workspaces\SDMB_result.mat'; % saved workspace

%% Where to save results
tikzfile = '.\results\tikz_figures\SDMB_result.tex';

%% Common plot options
phicut = deg2rad(0); % cutting plane phi
theta_plt = deg2rad(linspace(-180, 180, 181)); % theta values to plot
phi_plt = phicut*ones(size(theta_plt)); % phi values to plot
pltidx = direction_indices(theta, phi, theta_plt, phi_plt); % indices of data to plot

%% Parse SDR results
g_sdr_realizable = (4*pi/377) * abs(E_sdr_realizable).^2;

%% Parse manopt optimization results
[~, idx_sort_cost_manopt] = sort(cost_manopt);

idx_manopt_mid = idx_sort_cost_manopt(ceil(end/2)); % index of "average" result
cost_manopt_mid = cost_manopt(idx_manopt_mid); % "average" cost
ropt_manopt_mid = ropt_manopt(:, idx_manopt_mid); % reflection coefficients
g_manopt_mid = (4*pi/377) * abs(E_manopt(:,:,idx_manopt_mid)).^2;

idx_manopt_worst = idx_sort_cost_manopt(end); % index of "average" result
cost_manopt_worst = cost_manopt(idx_manopt_worst); % "average" cost
ropt_manopt_worst = ropt_manopt(:, idx_manopt_worst); % reflection coefficients
g_manopt_worst = (4*pi/377) * abs(E_manopt(:,:,idx_manopt_worst)).^2;

idx_manopt_best = idx_sort_cost_manopt(1); % index of "average" result
cost_manopt_best = cost_manopt(idx_manopt_best); % "average" cost
ropt_manopt_best = ropt_manopt(:, idx_manopt_best); % reflection coefficients
g_manopt_best = (4*pi/377) * abs(E_manopt(:,:,idx_manopt_best)).^2;

%% Parse GA optimization results
[~, idx_sort_cost_ga] = sort(cost_ga);

idx_ga_mid = idx_sort_cost_ga(ceil(end/2)); % index of "average" result
cost_ga_mid = cost_ga(idx_ga_mid); % "average" cost
ropt_ga_mid = ropt_ga(:, idx_ga_mid); % reflection coefficients
g_ga_mid = (4*pi/377) * abs(E_ga(:,:,idx_ga_mid)).^2;

idx_ga_worst = idx_sort_cost_ga(end); % index of "average" result
cost_ga_worst = cost_ga(idx_ga_worst); % "average" cost
ropt_ga_worst = ropt_ga(:, idx_ga_worst); % reflection coefficients
g_ga_worst = (4*pi/377) * abs(E_ga(:,:,idx_ga_worst)).^2;

idx_ga_best = idx_sort_cost_ga(1); % index of "average" result
cost_ga_best = cost_ga(idx_ga_best); % "average" cost
ropt_ga_best = ropt_ga(:, idx_ga_best); % reflection coefficients
g_ga_best = (4*pi/377) * abs(E_ga(:,:,idx_ga_best)).^2;

%% Display cost results
cost_table = table(cost_sdr, cost_sdr_realizable, cost_manopt_best, cost_manopt_mid, cost_manopt_worst, cost_ga_best, cost_ga_mid, cost_ga_worst);

cost_table_rel = table(cost_sdr, cost_sdr_realizable / cost_sdr, ...
                            cost_manopt_best / cost_sdr, ...
                            cost_manopt_mid / cost_sdr, ...
                            cost_manopt_worst / cost_sdr, ...
                            cost_ga_best / cost_sdr, ...
                            cost_ga_mid / cost_sdr, ...
                            cost_ga_worst / cost_sdr);
disp(cost_table_rel)
disp(cost_manopt)
disp(cost_ga)

%% Display transmission line lengths for best manopt result
beta = 230*ones(NP,1);
Z_line = 48*ones(NP,1);
Z_ref = 50*ones(NP,1);
ending = "short";
len = tllength(ropt_manopt_best, beta, Z_line, Z_ref, ending);
len = real(len)*1e3; %mm

termination_table = table(angle(ropt_manopt_best)*180/pi, len);
disp(termination_table)

%% Plot scan gain engelopes
figure('Position', [1000,400,800,640]);
plot(rad2deg(theta_target), 10*log10(abs(g_sdr)), 'r-x', 'MarkerSize', 4, 'LineWidth', 2)
hold on
plot(rad2deg(theta_plt), 10*log10(g_sdr_realizable(pltidx)), 'r', 'LineWidth', 1)

plot(rad2deg(theta_plt), 10*log10(g_ga_best(pltidx)), 'g-', 'LineWidth', 1)
plot(rad2deg(theta_plt), 10*log10(g_ga_mid(pltidx)), 'g--', 'LineWidth', 1)

plot(rad2deg(theta_plt), 10*log10(g_manopt_best(pltidx)), 'b-', 'LineWidth', 2)
plot(rad2deg(theta_plt), 10*log10(g_manopt_mid(pltidx)), 'b--', 'LineWidth', 1)

plot(rad2deg(theta_target), 10*log10(g_target), 'k-o', 'MarkerSize', 3, 'LineWidth', 2 )

lgnds = {};
lgnds{1} = 'Bound';
lgnds{2} = 'SDR';

lgnds{3} = 'GA best';
lgnds{4} = 'GA';

lgnds{5} = 'Manopt best';
lgnds{6} = 'Manopt';

lgnds{7} = 'Target';

hold off
xlabel('$\theta$', 'Interpreter', 'latex')
ylabel('Realized gain (dB)', 'Interpreter','latex')
% ylim([-20, 26])
xlim([-90, 90])
ttl = sprintf('Scan gain envelope, $\\varphi = $%.0f$^{\\circ}$', rad2deg(phicut));
title(ttl, 'Interpreter','latex')

leg = legend(lgnds);
leg.Interpreter = 'latex';
leg.Location = 'south';
leg.NumColumns = 2;

grid on
set(gca, 'FontSize', 11)

matlab2tikz('filename',tikzfile, 'externalData', true, 'standalone',true);