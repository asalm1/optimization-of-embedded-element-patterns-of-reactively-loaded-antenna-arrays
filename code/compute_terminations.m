%%
% Compute reactive termiations for scatterer ports. Save .mat workspace.
% ------------------------------------------------------------------------
% 06.06.2024 Albert Salmi, Department of Electronics and Nanoengineering,
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
cvx_solver sedumi

eta = 377; % wave impedance

%% Most important configurations, modify these ones
wspacefile = '.\results\workspaces\MDMB_result.mat'; % where to save results
datadir = '.\simulation_data\all_ports_active\multidriven'; % where to find data from
rng_state_file = '.\results\rng_state.mat'; % where to find random number generator state file
sparameter_fname = "MWS-run-0001.s55p";

feds = 1:5; % indices of driven elements
% feds = 1;

% theta_target = deg2rad(20); % Target directions in radians
theta_target = deg2rad(linspace(-19.5, 19.5, 7));

target_gain_array_dB = 16; % dB, Target main beam realized gain of the array, or target element gain to the directions theta_target
target_gain_array_dB_edges = target_gain_array_dB - 6; % Target gain at the edges of theta_steer, only relevant in multi-driven cases

evs = 10; % Number optimization runs with manopt and ga

%% Read simulated EEPs of all ports and S-parameters
% farfields
[theta, phi, Et, Ep, ff_info] = import_ff([datadir, '\farfield\'], "cst");

N = length(ff_info); % total number of elements
L = length(theta); % number of coordinates
pars = setdiff(1:N, feds); % indices of passive scatterer elements
ND = length(feds); % number of driven elements
NP = length(pars); % number of passive ports
freq = ff_info{1}.frequency; % frequency

% Read S-matrix of frequency freq
Sobj = sparameters(fullfile(datadir, sparameter_fname));
fidx = find( (abs(Sobj.Frequencies - freq)) == (min(abs(Sobj.Frequencies - freq))) == 1, 1); % index of closest frequency
S = Sobj.Parameters(:,:,fidx);
S = 0.5*(S+S.'); % Force symmetry


%% Manifold optimization settings
manopt_method = "MMF"; % default MMF (magnitude minimax fitting)
options.subsolver = @rlbfgs;
options.theta_rho = 3.3; % penalty coefficient's increment factor
options.theta_eps = 0.8; % accuracy tolerance's increment factor
options.theta_sigma = 0.8; % sigma-indicator's increment factor
options.rho0 = 1; % starting penalty coefficient
options.eps0 = 1e-3; % starting accuracy tolerance
options.d_min = 1e-8; % minimum step size
options.eps_min = 1e-6; % minimum accuracy tolerance
options.lambda_min_magnitude = 1e-4; % minimum mulptiplier value
options.lambda_max_magnitude = 1e6; % max multiplier value
options.lambda0_magnitude = 0.1; % starting multiplier value
options.verbosity = 1; % 0, 1, or 2. use this also in SDR.
options.maxiter = 1000; % Maximum number of iterations
warning('off', 'manopt:getHessian:approx');

% manifold optimization sub-solver settings
suboptions.maxiter = 200; % max number of iterations
suboptions.miniter = 10; % min number of sub solver iterations
suboptions.memory = Inf; % maximum memory usage
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

% Transform to linear scale and divide by number of driven elements
target_gain_element = eta/(4*pi*ND) * 10^(target_gain_array_dB/10);
target_gain_element_edges = eta/(4*pi*ND) * 10^(target_gain_array_dB_edges/10);

% Indices of target directions
scanidx = direction_indices(theta, phi, theta_target, phi_target);

% |\tilde{e}_{nl}|^2 and |\tilde{e}_{nl}|
EEP_target_squared = target_gain_element * ones(length(theta_target), ND);
if length(theta_target) > 1
    EEP_target_squared([1, end], :) = target_gain_element_edges;
end
EEP_target = sqrt(EEP_target_squared);

% Array target mainbeam realized gain in linear scale
g_target = 4*pi/eta * sum(EEP_target_squared, 2);

%% Formulate optimization data
Eco = Et.*sin(phi) + Ep.*cos(phi); % co-polarized component according to Ludwig 3

% Optimization data in correct form
C = Eco(scanidx,feds).';
SDP = S(feds, pars);
SPP = S(pars, pars);
Q = Eco(scanidx,pars).';

% Function handle to compute the new co-pol EEPs
Eco_new = @(r) (Eco(:,feds).' + SDP * ((diag(1./r) - SPP ) \ Eco(:,pars).')).' ;

%% Start SDR optimization
fprintf('\n Running SDR optimization \n \n')
[a, A, EEP_squared_sdr, info_sdr] = optimize_scateq_SDR(C, SDP, SPP, Q, EEP_target_squared.', options);

%% Extract realizable reflection coefficients from SDR result
G_sdr = (4*pi/377)*EEP_squared_sdr; % bound, element-wise
g_sdr = sum(G_sdr, 2); % bound, total scan gain

% take first NP terms of a, and compute the correcponding reflection
% coefficients
bar_ropt_sdr = a ./ ( (kron(eye(ND), SPP))*a + reshape(SDP.', ND*NP, 1)); % non-realizable, repeated R
ropt_sdr = bar_ropt_sdr(1:NP); % non-realizable
ropt_sdr_realizable = ropt_sdr ./ abs(ropt_sdr); % reactive reflection coefficients

E_sdr_realizable = Eco_new(ropt_sdr_realizable); % EEP

cost_sdr_realizable = max(abs( abs(E_sdr_realizable(scanidx,:)).^2 - EEP_target_squared ), [], 'all');

cost_sdr = max(abs( EEP_squared_sdr - EEP_target_squared ), [], 'all');

%% Minimum-power SDR optimization (Corcoles' method) if number of target directions is 1
if length(theta_target) == 1
    fprintf('\n Running Minimum-power SDR optimization \n \n')
    [a_cor, A_cor, cost_cor_raw, elapsed_time_cor] = optimize_Corcoles(S, Eco.', feds, pars, scanidx, [], [], options);
    
    % Extract realizable reflection coefficients from Minimum power SDR result
    ropt_cor_all = a_cor ./ (S*a_cor);
    ropt_cor = ropt_cor_all(pars);
    ropt_cor_realizable = ropt_cor ./ abs(ropt_cor);
    
    E_cor_realizable = Eco_new(ropt_cor_realizable);
    cost_cor = max(abs( abs(E_cor_realizable(scanidx,:)).^2 - EEP_target_squared ), [], 'all');
end

%% Manopt optimization
% initialize
ropt_manopt = zeros(NP, evs); % reflection coefficients
cost_manopt = zeros(evs, 1); % final cost function values
info_manopt = cell(evs, 1); % optimizer info
E_manopt = zeros(L, ND, evs); % EEPs

% Run manopt evs times
for it = 1:evs
    fprintf('\n Running Manopt (%d) \n \n', it)

    rng(it)

    % At run 1, use SDR result as initial guess, run 2: Min-power SDR, run
    % 3... use random 
    if it == 1
        r0 = ropt_sdr_realizable; % SDR result as initial guess
    elseif it == 2 && length(theta_target) == 1
        r0 = ropt_cor_realizable; % MP-SDR as initial guess
    else 
        r0 = exp(1j*2*pi*rand(NP,1)); % random
    end
   
    [result_manopt, cost_manopt(it), ~, info_manopt{it}] = optimize_scateq_manopt(C, SDP, SPP, Q, r0, EEP_target_squared.', options, manopt_method, suboptions);
    if manopt_method == "MMF"
        ropt_manopt(:,it) = result_manopt.B;
    else % if using MLSF formulation, the result_manopt is not a struct
        ropt_manopt(:,it) = result_manopt;
    end
    E_manopt(:,:,it) = Eco_new(ropt_manopt(:,it)); % EEP
end


%% Genetic algorithm
% objective function
gafun = @(psi)  max( abs( abs( C + SDP * ( ( diag(1./(exp(1j*psi))) - SPP ) \ Q) ).^2 - EEP_target_squared.'), [], 'all' );

psi_opt_ga = zeros(NP, evs); % phases of reflection coefficients
cost_ga = zeros(evs,1); % cost values
E_ga = zeros(L, ND, evs); % EEPs
for it = 1:evs
    rng(it)
    fprintf('\n Running GA (%d) \n \n', it)
    psi_opt_ga(:,it) = ga(gafun, NP, [], [], [], [], 0, 2*pi);
    cost_ga(it) = gafun(psi_opt_ga(:,it));
    E_ga(:,:,it) = Eco_new( exp(1j*psi_opt_ga(:,it)) );
end

ropt_ga = exp(1j*psi_opt_ga); % reflection coefficients

%% Save results
save(wspacefile);

