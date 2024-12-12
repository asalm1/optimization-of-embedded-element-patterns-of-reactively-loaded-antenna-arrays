%%
% Visualize and analyze computational time of SDR, MO, and GA algorithms.
% Find a polynomial function describing the time required to the scatterer
% port terminations using MO combined with SDR.
% ------------------------------------------------------------------------
% 20.11.2024 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%% Clear
clear
clc
close all

%% Load data
load '.\results\workspaces\time_comparison.mat'; % saved workspace

%% Where to save results
outfolder = '.\results\timecomparison';

%% Extract time data
t_sdr = zeros(NumTests, 1);
t_mo = zeros(NumTests, 1);
t_ga = zeros(NumTests, 1);
t_tot = zeros(NumTests, 1);

for it = 1:NumTests
    t_sdr(it) = info_sdr{it}.elapsed_time;
    t_mo(it) = info_manopt{it}.time;
    t_ga(it) = info_ga{it}.time;
    t_tot(it) = t_sdr(it) + t_mo(it);
end

%% Fit polynomial p(N, M) to total computational time
N = zeros(NumTests, 1); M = zeros(NumTests, 1);
for it = 1:NumTests
    N(it) = length(feds{it}); M(it) = length(pars{it});
end

[fitobject, gof] = fit([N, M], t_tot, 'poly34')   %'poly34'

%% Plot computation time with fitted polynomial and cost function values
figure
plot(testID, t_sdr, 'r', 'LineWidth', 2)
hold on
plot(testID, t_mo, 'b', 'LineWidth', 2)
plot(testID, t_tot, 'k', 'LineWidth', 2)
plot(testID, t_ga, 'g', 'LineWidth', 2)
plot(testID, fitobject(N, M), 'm:', 'LineWidth', 2)
hold off
ylim([0, Inf])
ylabel('time (s)')

yyaxis right
plot(testID, [cost_sdr{:}], 'r--', 'LineWidth', 2)
hold on
plot(testID, [cost_manopt{:}], 'b--', 'LineWidth', 2)
plot(testID, [cost_ga{:}], 'g--', 'LineWidth', 2)
hold off
ylabel('cost')
legend({'t, SDR', 't, MO', 't, SDR+MO', 't, GA', 'Fit', 'bound, SDR', 'cost, MO', 'cost, GA'})

%% Extrapolate data
N_extra = 1:20;
M_extra = 5:5:200;
[NE, ME] = meshgrid(N_extra, M_extra);
FITE = fitobject(NE,ME);

figure
surf(NE,ME, FITE)
hold on
plot(fitobject, [N, M], t_tot)
hold off
xlabel('N')
ylabel('M')
zlabel('time')
ylim([1, Inf])
xlim([1, Inf])

fprintf('With N=%d, M=%d, computation time would be %.1f days. \n\n', NE(end,end), ME(end,end), FITE(end,end) /(60*60*24))

%% Save data

% From vector to matrix form
NN = reshape(N, [length(NumM), length(NumN)]);
MM = reshape(M, [length(NumM), length(NumN)]);

NN_extra = repmat(NN(1,:), length(M_extra), 1);
MM_extra = repmat(M_extra.', 1, length(NumN));

T_SDR = reshape(t_sdr, [length(NumM), length(NumN)]);
T_MO = reshape(t_mo, [length(NumM), length(NumN)]);
T_GA = reshape(t_ga, [length(NumM), length(NumN)]);
T_TOT = reshape(t_tot, [length(NumM), length(NumN)]);
T_FIT = fitobject(NN_extra, MM_extra);

COST_SDR = reshape([cost_sdr{:}], [length(NumM), length(NumN)]);
COST_MO = reshape([cost_manopt{:}], [length(NumM), length(NumN)]);
COST_GA = reshape([cost_ga{:}], [length(NumM), length(NumN)]);

for it = 1:length(NumN)
    fname_t_sdr = sprintf('t_sdr_N%d.txt', NumN(it));
    fname_t_mo = sprintf('t_mo_N%d.txt', NumN(it));
    fname_t_ga = sprintf('t_ga_N%d.txt', NumN(it));
    fname_t_tot = sprintf('t_tot_N%d.txt', NumN(it));
    fname_t_fit = sprintf('t_fit_N%d.txt', NumN(it));

    fname_cost_sdr = sprintf('cost_sdr_N%d.txt', NumN(it));
    fname_cost_mo = sprintf('cost_mo_N%d.txt', NumN(it));
    fname_cost_ga = sprintf('cost_ga_N%d.txt', NumN(it));

    writematrix([MM(:,it), T_SDR(:,it)], fullfile(outfolder, fname_t_sdr));
    writematrix([MM(:,it), T_MO(:,it)], fullfile(outfolder, fname_t_mo));
    writematrix([MM(:,it), T_GA(:,it)], fullfile(outfolder, fname_t_ga));
    writematrix([MM(:,it), T_TOT(:,it)], fullfile(outfolder, fname_t_tot));
    writematrix([MM_extra(:,it), T_FIT(:,it)], fullfile(outfolder, fname_t_fit));

    writematrix([MM(:,it), COST_SDR(:,it)], fullfile(outfolder, fname_cost_sdr));
    writematrix([MM(:,it), COST_MO(:,it)], fullfile(outfolder, fname_cost_mo));
    writematrix([MM(:,it), COST_GA(:,it)], fullfile(outfolder, fname_cost_ga));
end

