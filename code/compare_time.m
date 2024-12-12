%%
% Compare computational time with different optimization methods and
% different size problems.
% ------------------------------------------------------------------------
% 18.11.2024 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%% Clear and initialize
clear
clc
close all

% Clear optimization variables
clear problem
clear options
clear suboptions

% Choose solver
% cvx_solver sedumi

%% Give input and output folders and define target scan sector and target gain
wspacefile = '.\results\workspaces\time_comparison.mat'; % where to save results
datadir = '.\simulation_data\all_ports_active\multidriven'; % where to find data from
rng_state_file = '.\results\rng_state.mat'; % where to find random number generator state file
sparameter_fname = "MWS-run-0001.s55p";

theta_target = deg2rad(linspace(-19.5, 19.5, 7));

target_gain_array_dB = 16; % dB, Target main beam realized gain of the array, or target element gain to the directions theta_target
target_gain_array_dB_edges = target_gain_array_dB - 6; % Target gain at the edges of theta_steer, only relevant in multi-driven cases

%% Read simulated EEPs of all ports and S-parameters
% farfields
[theta, phi, Et_all, Ep_all, ff_info] = import_ff([datadir, '\farfield\'], "cst");
Eco_all = Et_all.*sin(phi) + Ep_all.*cos(phi); % co-polarized component according to Ludwig 3

L = length(theta); % number of coordinates
freq = ff_info{1}.frequency; % frequency

% Read S-matrix of frequency freq
Sobj = sparameters(fullfile(datadir, sparameter_fname));
fidx = find( (abs(Sobj.Frequencies - freq)) == (min(abs(Sobj.Frequencies - freq))) == 1, 1); % index of closest frequency
S_all = Sobj.Parameters(:,:,fidx);
S_all = 0.5*(S_all+S_all.'); % Force symmetry

%% Define test cases
N_all = 55;                             % Number of all ports
feds_all = 1:5;                         % indices of all driven ports
pars_all = setdiff(1:N_all, feds_all);  % intides of all scatterer ports

% indices of scatterer ports in the order in which they are selected in the
% study. For example, when driven ports 1,2,3 are used, then scatterer
% ports are chosen based on order of pars_global{3}.
pars_order{1} = [51, 10, 11, 9, 12, 8, 13, 7, 14, 6, 52, 19, 20, 18, 21, 17, 22, 16, 23, 15, 53, 28, 29, 27, 30, 26, 31, 25, 32, 24, 54, 37, 38, 36, 39, 35, 40, 34, 41, 33, 55, 46, 47, 45, 48, 44, 49, 43, 50, 42];
pars_order{2} = [51, 52, 10, 19, 11, 20, 9, 18, 12, 21, 8, 17, 13, 22, 7, 16, 14, 23, 6, 15, 53, 28, 29, 27, 30, 26, 31, 25, 32, 24, 54, 37, 38, 36, 39, 35, 40, 34, 41, 33, 55, 46, 47, 45, 48, 44, 49, 43, 50, 42];
pars_order{3} = [51, 52, 53, 10, 19, 28, 11, 20, 29, 9, 18, 27, 12, 21, 30, 8, 17, 26, 13, 22, 31, 7, 16, 25, 14, 23, 32, 6, 15, 24, 54, 37, 38, 36, 39, 35, 40, 34, 41, 33, 55, 46, 47, 45, 48, 44, 49, 43, 50, 42];
pars_order{4} = [51, 52, 53, 54, 10, 19, 28, 37, 11, 20, 29, 38, 9, 18, 27, 36, 12, 21, 30, 39, 8, 17, 26, 35, 13, 22, 31, 40, 7, 16, 25, 34, 14, 23, 32, 41, 6, 15, 24, 33, 55, 46, 47, 45, 48, 44, 49, 43, 50, 42];
pars_order{5} = [51, 52, 53, 54, 55, 10, 19, 28, 37, 46, 11, 20, 29, 38, 47, 9, 18, 27, 36, 45, 12, 21, 30, 39, 48, 8, 17, 26, 35, 44, 13, 22, 31, 40, 49, 7, 16, 25, 34, 43, 14, 23, 32, 41, 50, 6, 15, 24, 33, 42];

NumN = 1:5;     % How many driven ports in tests
NumM = 1:1:50;    % How many scatterer ports in tests

NumTests = length(NumN) * length(NumM);
testID = 1 : NumTests;

S = cell(NumTests, 1);
Eco = cell(NumTests, 1);
feds = cell(NumTests, 1);
pars = cell(NumTests, 1);

kt = 1;
for it = 1:length(NumN)

    feds_global_it = 1 : NumN(it); % indices of driven ports in these tests
    pars_order_it = pars_order{NumN(it)}; % selection order for scatterer ports

    for jt = 1:length(NumM)
        pars_global_it_jt = pars_order_it( 1 : NumM(jt) ); % indices of scatterer ports

        % Terminate unused ports to short circuits. First find indices of
        % used and unused ports
        used = [feds_global_it, pars_global_it_jt];
        others = setdiff(1:N_all, used);

        % Diagonal reflection coefficient matrix of unused ports'
        % terminations (short circuits)
        R_shorts = -1 * eye(length(others));

        % Resulting E and S matrices of used ports after shorting the
        % unused ports
        Eco{kt} = (Eco_all(:, used).' + S_all(used, others) * ( ( inv(R_shorts) - S_all(others, others) ) \ (Eco_all(:, others).') )).';
        S{kt} = S_all(used, used) + S_all(used, others) * ( ( inv(R_shorts) - S_all(others, others) ) \ S_all(others, used) );


        % Define indices of driven and scatterer pports related to the new
        % E and S matrices
        feds{kt} = 1:NumN(it);
        pars{kt} = setdiff(1:length(used), feds{kt});

        kt = kt+1;
    end
end

%% Manifold optimization settings
manopt_method = "MMF"; % default MMF (magnitude minimax fitting)
options.subsolver = @rlbfgs;
options.theta_rho = 10; % penalty coefficient's increment factor
options.theta_eps = 0.9; % accuracy tolerance's increment factor
options.theta_sigma = 0.8; % sigma-indicator's increment factor
options.rho0 = 0.1; % starting penalty coefficient
options.eps0 = 0.9; % starting accuracy tolerance
options.d_min = 1e-4; % minimum step size
options.eps_min = 1e-4; % minimum accuracy tolerance
options.lambda_min_magnitude = 1e-6; % minimum mulptiplier value
options.lambda_max_magnitude = 1e6; % max multiplier value
options.lambda0_magnitude = 1; % starting multiplier value
options.verbosity = 1; % 0, 1, or 2. use this also in SDR.
options.maxiter = 1000; % Maximum number of iterations
warning('off', 'manopt:getHessian:approx');

% manifold optimization sub-solver settings
suboptions.maxiter = 50; % max number of iterations
suboptions.miniter = 2; % min number of sub solver iterations
suboptions.memory = Inf; % maximum memory usage
suboptions.minstepsize = 1e-6;
suboptions.verbosity = 1; % 0, 1, or 2

%% Load random number generator state or generate new and save
if exist(rng_state_file)
    load(rng_state_file, "rng_state");
    rng(rng_state)
else
    rng_state = rng("default");
    save(rng_state_file, "rng_state");
end


%% Define target
phi_target = zeros(size(theta_target)); % phi coordinates

% Indices of target directions
scanidx = direction_indices(theta, phi, theta_target, phi_target);

%% Run optimizations
info_sdr = cell(length(testID), 1);
info_manopt = cell(length(testID), 1);
info_ga = cell(length(testID), 1);
cost_sdr = cell(length(testID), 1);
cost_manopt = cell(length(testID), 1);
cost_ga = cell(length(testID), 1);

for it = 1:length(testID)
    fprintf('Run %d / %d, N = %d, M = %d \n\n', it, length(testID), length(feds{it}), length(pars{it}));
    [info_sdr{it}, info_manopt{it}, info_ga{it}, cost_sdr{it}, cost_manopt{it}, cost_ga{it}] = optimize_terminations(theta, phi, Eco{it}, S{it}, feds{it}, pars{it}, scanidx, target_gain_array_dB, target_gain_array_dB_edges, options, suboptions);
end

%% Save results
save(wspacefile, "cost_ga", "cost_manopt", "cost_sdr", "info_ga", "info_manopt", "info_sdr", "feds", "pars", "testID", "NumTests", "NumN", "NumM");

%% Function for optimizing loads
function [info_sdr, info_manopt, info_ga, cost_sdr, cost_manopt, cost_ga] = optimize_terminations(theta, phi, Eco, S, feds, pars, scanidx, target_gain_array_dB, target_gain_array_dB_edges, options, suboptions)
    eta = 377; % wave impedance

    ND = length(feds); NP = length(pars);
    theta_target = theta(scanidx); phi_target = phi(scanidx);

    
    % Transform to linear scale and divide by number of driven elements
    target_gain_element = eta/(4*pi*ND) * 10^(target_gain_array_dB/10);
    target_gain_element_edges = eta/(4*pi*ND) * 10^(target_gain_array_dB_edges/10);
    
    
    % |\tilde{e}_{nl}|^2 and |\tilde{e}_{nl}|
    EEP_target_squared = target_gain_element * ones(length(theta_target), ND);
    if length(theta_target) > 1
        EEP_target_squared([1, end], :) = target_gain_element_edges;
    end
    
    
    % Optimization data in correct form
    C = Eco(scanidx,feds).';
    SDP = S(feds, pars);
    SPP = S(pars, pars);
    Q = Eco(scanidx,pars).';

    % Start SDR optimization
    fprintf('\n Running SDR optimization \n \n')
    [a, A, EEP_squared_sdr, info_sdr] = optimize_scateq_SDR(C, SDP, SPP, Q, EEP_target_squared.', options);
    
    cost_sdr = max(abs( EEP_squared_sdr - EEP_target_squared ), [], 'all');

    % take first NP terms of a, and compute the correcponding reflection
    % coefficients
    bar_ropt_sdr = a ./ ( (kron(eye(ND), SPP))*a + reshape(SDP.', ND*NP, 1)); % non-realizable, repeated R
    ropt_sdr = bar_ropt_sdr(1:NP); % non-realizable
    ropt_sdr_realizable = ropt_sdr ./ abs(ropt_sdr); % reactive reflection coefficients


    fprintf('\n Running Manopt \n \n')   
    [~, cost_manopt, ~, info_manopt] = optimize_scateq_manopt(C, SDP, SPP, Q, ropt_sdr_realizable, EEP_target_squared.', options, "MMF", suboptions);

    % Genetic algorithm
    % objective function
    gafun = @(psi)  max( abs( abs( C + SDP * ( ( diag(1./(exp(1j*psi))) - SPP ) \ Q) ).^2 - EEP_target_squared.'), [], 'all' );
    rng(1)
    fprintf('\n Running GA \n \n')
    gastart = tic;
    psi_opt_ga = ga(gafun, NP, [], [], [], [], 0, 2*pi);
    cost_ga = gafun(psi_opt_ga);
    info_ga.time = toc(gastart);

    if options.verbosity > 0
        if info_ga.time > 60
            elapsed_string_ga = seconds(info_ga.time);
            elapsed_string_ga.Format = 'hh:mm:ss';
        else
            elapsed_string_ga = string(info_ga.time);
        end
        
        fprintf('Total elapsed time in SDR optimization: %s \n \n', elapsed_string_ga);
    end

end