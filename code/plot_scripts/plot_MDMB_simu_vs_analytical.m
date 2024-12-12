%%
% MDMB, Compare analytical Matlab-computed result curves to simulated result
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
load '.\results\workspaces\MDMB_result.mat'; % saved workspace

%% Where to save results
savefolder = '.\results\tikz_figures\';

%% Import Simulated far-fields
simufolder = '.\simulation_data\scatterer_ports_terminated\MDMB\farfield\'; % where to find data from
[theta_simu, phi_simu, Et_simu, Ep_simu, ff_info_simu] = import_ff(simufolder, "cst");

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
G_manopt_best = (4*pi/377) * abs(E_manopt(:,:,idx_manopt_best)).^2;

%% Parse simulated results
G_simu = (4*pi/377) * abs(Eco_simu).^2;

%% Cut values that are too low to show
MING = -20;
G_analyt_dB = 10*log10(G_manopt_best(pltidx,:));
G_analyt_dB(G_analyt_dB < MING) = MING;

G_simu_dB = 10*log10(G_simu(pltidx_simu,:));
G_simu_dB(G_simu_dB < MING) = MING;

MAXG = round(max(G_simu_dB, [], 'all'), 1);

%% Plot patterns in polar coordinates
figure('Position', [1000,400,800,640]);

for it = 1:ND
    polarplot(theta_plt, G_simu_dB(:,it), 'LineStyle', '-', 'LineWidth', 2)
    hold on
    polarplot(theta_plt, G_analyt_dB(:,it), 'LineStyle', ':', 'LineWidth', 1)
end

% legend({'Analytical', 'Simulated'}, 'Interpreter', 'latex');

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
for it = 1:ND
    simufname = sprintf('MDMB_polardata_simu_%d.txt', it);
    anafname = sprintf('MDMB_polardata_analytical_%d.txt', it);
    writematrix([theta_plt.', G_simu_dB(:,it)], [savefolder, simufname])
    writematrix([theta_plt.', G_analyt_dB(:,it)], [savefolder, anafname])
end