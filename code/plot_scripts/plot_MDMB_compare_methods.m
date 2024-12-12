%%
% Plot MDMB results comparing different optimization methods
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
tikzfile_scan = '.\results\tikz_figures\MDMB_result_scan.tex';
tikzfile = '.\results\tikz_figures\MDMB_result.tex';
tikzdata = 'MDMB_result_data';
tikzdata_scan = 'MDMB_result_scan_data';

%% Scan directions to illustrate
theta_steer = deg2rad([-10, 0, 15]);
steeridx = direction_indices(theta, phi, theta_steer, zeros(size(theta_steer)));

%% Common plot options
phicut = deg2rad(0); % cutting plane phi
theta_plt = deg2rad(linspace(-90, 90, 361)); % theta values to plot
phi_plt = phicut*ones(size(theta_plt)); % phi values to plot
pltidx = direction_indices(theta, phi, theta_plt, phi_plt); % indices of data to plot

%% Parse SDR results
[G_sdr_realizable, g_sdr_realizable, gscan_sdr_realizable] = gain_metrics(E_sdr_realizable, steeridx, eta);

%% Parse manopt optimization results
[~, idx_sort_cost_manopt] = sort(cost_manopt);

idx_manopt_mid = idx_sort_cost_manopt(ceil(end/2)); % index of "average" result
cost_manopt_mid = cost_manopt(idx_manopt_mid); % "average" cost
ropt_manopt_mid = ropt_manopt(:, idx_manopt_mid); % reflection coefficients
[G_manopt_mid, g_manopt_mid, gscan_manopt_mid] = gain_metrics(E_manopt(:,:,idx_manopt_mid), steeridx, eta);

idx_manopt_worst = idx_sort_cost_manopt(end); % index of "average" result
cost_manopt_worst = cost_manopt(idx_manopt_worst); % "average" cost
ropt_manopt_worst = ropt_manopt(:, idx_manopt_worst); % reflection coefficients
[G_manopt_worst, g_manopt_worst, gscan_manopt_worst] = gain_metrics(E_manopt(:,:,idx_manopt_worst),steeridx, eta);

idx_manopt_best = idx_sort_cost_manopt(1); % index of "average" result
cost_manopt_best = cost_manopt(idx_manopt_best); % "average" cost
ropt_manopt_best = ropt_manopt(:, idx_manopt_best); % reflection coefficients
[G_manopt_best, g_manopt_best, gscan_manopt_best] = gain_metrics(E_manopt(:,:,idx_manopt_best),steeridx, eta);

%% Parse GA optimization results
[~, idx_sort_cost_ga] = sort(cost_ga);

idx_ga_mid = idx_sort_cost_ga(ceil(end/2)); % index of "average" result
cost_ga_mid = cost_ga(idx_ga_mid); % "average" cost
ropt_ga_mid = ropt_ga(:, idx_ga_mid); % reflection coefficients
[G_ga_mid, g_ga_mid, gscan_ga_mid] = gain_metrics(E_ga(:,:,idx_ga_mid), steeridx,eta);

idx_ga_worst = idx_sort_cost_ga(end); % index of "average" result
cost_ga_worst = cost_ga(idx_ga_worst); % "average" cost
ropt_ga_worst = ropt_ga(:, idx_ga_worst); % reflection coefficients
[G_ga_worst, g_ga_worst, gscan_ga_worst] = gain_metrics(E_ga(:,:,idx_ga_worst), steeridx,eta);

idx_ga_best = idx_sort_cost_ga(1); % index of "average" result
cost_ga_best = cost_ga(idx_ga_best); % "average" cost
ropt_ga_best = ropt_ga(:, idx_ga_best); % reflection coefficients
[G_ga_best, g_ga_best, gscan_ga_best] = gain_metrics(E_ga(:,:,idx_ga_worst),steeridx, eta);

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

%% Plot scan patterns
figure('Position', [1000,400,800,640]);
plot(rad2deg(theta_target), 10*log10(abs(g_sdr)), 'r-x', 'MarkerSize', 4, 'LineWidth', 2)
hold on
plot(rad2deg(theta_plt), 10*log10(g_manopt_best(pltidx)), 'b-', 'LineWidth', 2)

lgnds = {};
lgnds{1} = 'Bound';
lgnds{2} = 'Scan gain envelope';

clrs = {'c', 'm', 'y'};

for it = 1:length(theta_steer)
    plot(rad2deg(theta_plt), 10*log10(gscan_manopt_best(pltidx, it)), 'LineStyle', '-', 'Color', clrs{it}, 'LineWidth', 1)
    lgnds{end+1} = sprintf('Array gain, $\\theta_0 = %d^{\\circ}$', rad2deg(theta_steer(it)));
end
plot(rad2deg(theta_target), 10*log10(g_target), 'k-o', 'MarkerSize', 3, 'LineWidth', 2 )
lgnds{end+1} = 'Target';

hold off
xlabel('$\theta$', 'Interpreter', 'latex')
ylabel('Realized gain (dB)', 'Interpreter','latex')
% ylim([-20, 26])
xlim([-90, 90])
ttl = sprintf('Manopt best result and bound, $\\varphi = $%.0f$^{\\circ}$', rad2deg(phicut));
title(ttl, 'Interpreter','latex')

leg = legend(lgnds);
leg.Interpreter = 'latex';
leg.Location = 'south';
leg.NumColumns = 2;

grid on
set(gca, 'FontSize', 11)

matlab2tikz('filename',tikzfile_scan, 'externalData', true, 'standalone',true);
%% Functions

function [G, genv, gscan] = gain_metrics(E, steeridx, eta)
    G = (4*pi/eta) * abs(E).^2;
    genv = sum(G, 2);
    Ascan = conj(E(steeridx,:));
    gscan = zeros(length(E(:,1)), length(steeridx));
    for it = 1:length(steeridx)
        gscan(:,it) = (4*pi/eta) * abs( sum( Ascan(it,:) .* E, 2 ) ).^2 ./ norm(Ascan(it,:))^2;
    end
end