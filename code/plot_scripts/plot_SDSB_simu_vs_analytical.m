%%
% Compare analytical Matlab-computed result curves to simulated result
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
load '.\results\workspaces\SDSB_result_scan20.mat'; % saved workspace

%% Where to save results
file_analytical = '.\results\tikz_figures\SDSB_polardata_analytical.txt';
file_simu = '.\results\tikz_figures\SDSB_polardata_simu.txt';

%% Import Simulated far-fields
simufile = '.\simulation_data\scatterer_ports_terminated\SDSB_scan20\farfield_source_(f=5)_[1].ffs'; % where to find data from
[theta_simu, phi_simu, Et_simu, Ep_simu, ff_info_simu] = import_ff(simufile, "cst");

Eco_simu = Et_simu.*sin(phi_simu) + Ep_simu.*cos(phi_simu); % co-polarized component according to Ludwig 3

%% Common plot options
phicut = deg2rad(0); % cutting plane phi
theta_plt = deg2rad(linspace(-180, 180, 181)); % theta values to plot
phi_plt = phicut*ones(size(theta_plt)); % phi values to plot
pltidx = direction_indices(theta, phi, theta_plt, phi_plt); % indices of data to plot
pltidx_simu = direction_indices(theta_simu, phi_simu, theta_plt, phi_plt); % indices of data to plot

%% Parse analytical results
[~, idx_sort_cost_manopt] = sort(cost_manopt);
idx_manopt_best = idx_sort_cost_manopt(1); % index of "average" result
cost_manopt_best = cost_manopt(idx_manopt_best); % "average" cost
ropt_manopt_best = ropt_manopt(:, idx_manopt_best); % reflection coefficients
g_manopt_best = (4*pi/377) * abs(E_manopt(:,:,idx_manopt_best)).^2;

%% Parse simulated results
g_simu = (4*pi/377) * abs(Eco_simu).^2;

%% Cut values that are too low to show
MING = -20;
g_analyt_dB = 10*log10(g_manopt_best(pltidx));
g_analyt_dB(g_analyt_dB < MING) = MING;

g_simu_dB = 10*log10(g_simu(pltidx_simu));
g_simu_dB(g_simu_dB < MING) = MING;

MAXG = round(max(g_simu_dB), 1);

%% Plot patterns in polar coordinates
figure('Position', [1000,400,800,640]);

polarplot(theta_plt, g_simu_dB, 'k', 'LineWidth', 4)
hold on
polarplot(theta_plt, g_analyt_dB, 'k--', 'LineWidth', 2)

legend({'Analytical', 'Simulated'}, 'Interpreter', 'latex');

hold off

Pax = gca;

Pax.ThetaLim = [-180,180];
Pax.RLim = [MING, MAXG+1];
Pax.RTick = [MAXG-10, MAXG-3, MAXG];
Pax.ThetaTick = -180:20:180;
% Pax.ThetaTickLabel = {'$-60^{\circ}$', '$-30^{\circ}$', '$0^{\circ}$', '$30^{\circ}$', '$60^{\circ}$'};
% % Pax.RTickLabel = {'-3', '0 dB'};
% Pax.RTickLabel = [];

Pax.ThetaZeroLocation = 'top';
Pax.ThetaDir = 'clockwise';
Pax.FontSize = 30;
Pax.GridAlpha = 0.6;
Pax.TickLabelInterpreter = 'latex';

%% Write data to text file for tikz plot
writematrix([theta_plt.', g_simu_dB], file_simu);
writematrix([theta_plt.', g_analyt_dB], file_analytical);