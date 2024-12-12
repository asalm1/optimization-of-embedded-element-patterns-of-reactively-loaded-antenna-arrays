%%
% Tune Manopt to converge fast.
% ------------------------------------------------------------------------
% 19.11.2024 Albert Salmi, Department of Electronics and Nanoengineering,
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

%% Give input and output folders and define target scan sector and target gain
datadir = '.\simulation_data\all_ports_active\multidriven'; % where to find data from
rng_state_file = '.\results\rng_state.mat'; % where to find random number generator state file
sparameter_fname = "MWS-run-0001.s55p";

theta_target = deg2rad(linspace(-19.5, 19.5, 7));

target_gain_array_dB = 16; % dB, Target main beam realized gain of the array, or target element gain to the directions theta_target
target_gain_array_dB_edges = target_gain_array_dB - 6; % Target gain at the edges of theta_steer, only relevant in multi-driven cases

%% Define test cases
N = 55;
feds = 1:5;
pars = setdiff(1:N, feds);

%% Read simulated EEPs of all ports and S-parameters
% farfields
[theta, phi, Et, Ep, ff_info] = import_ff([datadir, '\farfield\'], "cst");

L = length(theta); % number of coordinates
freq = ff_info{1}.frequency; % frequency

% Read S-matrix of frequency freq
Sobj = sparameters(fullfile(datadir, sparameter_fname));
fidx = find( (abs(Sobj.Frequencies - freq)) == (min(abs(Sobj.Frequencies - freq))) == 1, 1); % index of closest frequency
S = Sobj.Parameters(:,:,fidx);
S = 0.5*(S+S.'); % Force symmetry

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


%% Define target
phi_target = zeros(size(theta_target)); % phi coordinates

% Indices of target directions
scanidx = direction_indices(theta, phi, theta_target, phi_target);

Eco = Et.*sin(phi) + Ep.*cos(phi); % co-polarized component according to Ludwig 3

%% Run optimizations
[info_manopt, cost_manopt] = optimize_terminations(theta, phi, Eco, S, feds, pars, scanidx, target_gain_array_dB, target_gain_array_dB_edges, options, suboptions);

disp(info_manopt.time)
disp(cost_manopt)

%% Function for optimizing loads
function [info_manopt,cost_manopt] = optimize_terminations(theta, phi, Eco, S, feds, pars, scanidx, target_gain_array_dB, target_gain_array_dB_edges, options, suboptions)
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

    r0 = exp(1j*2*pi*rand(NP,1));

    fprintf('\n Running Manopt \n \n')   
    [~, cost_manopt, ~, info_manopt] = optimize_scateq_manopt(C, SDP, SPP, Q, r0, EEP_target_squared.', options, "MMF", suboptions);

end