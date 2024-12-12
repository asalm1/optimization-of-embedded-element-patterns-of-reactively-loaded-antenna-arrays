function [ropt, cost, elapsed_time, info] = optimize_scateq_manopt(C, SDP, SPP, Q, r0, target, options, method, suboptions)
%
% Optimize the scattering equation. Provide target envelope optionally.
% Otherwise, minimize or maximize the scattering equation by setting target
% = +- Inf. This function uses SDR and returns the raw SDR output.
%-------------------------------------------------------------------------
% INPUT  C       : (ND,:)    driven port initial data
%        SDP     : (ND, NP) coupling between driven and scatterer ports
%        SPP     : (NP, NP) ceoupling between scatterer ports
%        Q       : (NP, :)  scatterer port data
%        r0      : (NP, 1) initial guess
%        target  : (ND, :)  target for driven ports (default +Inf)
%        options : struct   optimization settings
%        method  : string, "MMF" or "MLSF" (default)
%        suboptions : struct, options for sub-solver (MMF only)
%
% OUTPUT  ropt          : (NP, 1) optimal reflection coefficients or struct
%         cost          : (1,1) last objective value
%         elapsed time  : string
%         info          : struct, optimizer output at iterations
% ------------------------------------------------------------------------
% 06.06.2023 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    C           (:,:)
    SDP         (:,:)
    SPP         (:,:)
    Q           (:,:)
    r0          (:,1)
    target      (:,:) = Inf
    options     (:,:) = []
    method      (:,:) = "MLSF"
    suboptions  (:,:) = []
end

tstart = tic;

% Set default options if not set
if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
end
if ~isfield(options, 'tolgradnorm')
    options.tolgradnorm = 1e-8;
end
if ~isfield(options, 'verbosity')
    options.verbosity = 1;
end
if ~isfield(options, 'solver')
    options.solver = @(prob, x0, opts) trustregions(prob, x0, opts);
end

info = {}; % init info struct

if target == Inf % maximize frobenius norm of the scattering equation

    [ropt, cost, info] = minimize_scateq_fronorm_manopt(C, SDP, SPP, Q, r0, options, 1);

elseif target == -Inf % minimize frobenius norm of the scattering equation

    [ropt, cost, info] = minimize_scateq_fronorm_manopt(C, SDP, SPP, Q, r0, options, 0);

else % minimize distance to the target

    if method == "MLSF"
        [ropt, cost] = optimize_scateq_MLSF_manopt(C, SDP, SPP, Q, r0, target, options);
    else
        [ropt, cost, info] = optimize_scateq_MMF_manopt(C, SDP, SPP, Q, r0, target, options, suboptions);
    end

end

elapsed_time = toc(tstart);
if options.verbosity > 0
    if elapsed_time > 60
        elapsed_string = seconds(elapsed_time);
        elapsed_string.Format = 'hh:mm:ss';
    else
        elapsed_string = string(elapsed_time);
    end
    
    fprintf('Total elapsed time in Manifold optimization: %s \n \n', elapsed_string);
end

end


function [x_opt, x_cost, info] = minimize_scateq_fronorm_manopt(C, SDP, SPP, Q, x0, options, maximize)
%
% Minimize f(x) = || C + SDP (diag(x)^{-1} - SPP)^{-1} Q ||_F^2 using
% manifold optimization
%-------------------------------------------------------------------------
% INPUT  C        : (:,:) driven port data
%        SDP      : (:,:) coupling from scatterers to driven
%        SPP      : (:,:) coupling between the scatterers
%        Q        : (:,:) scatterer's data
%        x0       : (:,1) initial guess, default ines
%        options  : struct, for manopt algorithm
%        maximize : bool, maximize instead of minimization
%
% OUTPUT  x_opt   : (:,1) optimized vector
%         x_cost  : (1,1) final cost function value
%         info    : struct, manopt info
% ------------------------------------------------------------------------
% 10.10.2023 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    C       (:,:)
    SDP     (:,:)
    SPP     (:,:)
    Q       (:,:)
    x0      (:,1) = ones(length(SPP(:,1)),1)
    options (:,:) = []
    maximize (1,1) = 0
end

function [cost, rgrad] = scateq_costgrad(C, P, S, Q, x, maximize)
% Compute f(x) = || C + P (diag(x)^{-1} - S)^{-1} Q ||_F^2 and it's
% Riemannian gradient.
    M = diag(1./x) - S;
    MinvQ = M \ Q; % = M^{-1} * Q
    PMinv = ( (M.') \ (P.') ).'; % = P * M^{-1}
    alpha = -1 ./ (imag(x) - 1j*real(x)).^2;
    G = C + P * MinvQ;
    cost = (1-2*maximize) * norm(G , "fro")^2;
    egrad = (1-2*maximize) *  2*conj( alpha .* diag( MinvQ * G' * PMinv ));
    rgrad = egrad - real( egrad .* conj(x) ) .* x; % Projection to the manifold
end

options.statsfun = statsfunhelper('current_point', @(x) x);

clear problem % clear old problem

problem.M = complexcirclefactory(length(x0));

problem.costgrad = @(r) scateq_costgrad(C, SDP, SPP, Q, r, maximize);

[x_opt, x_cost, info] = options.solver(problem, x0, options);

end


function [ropt, cost] = optimize_scateq_MLSF_manopt(C, SDP, SPP, Q, r0, target, options)
%
% Optimize scattering equation using magnitude least squares fittinf (MLSF)
% formulation
%-------------------------------------------------------------------------
% INPUT  C       : (ND,:)    driven port initial data
%        SDP     : (ND, NP) coupling between driven and scatterer ports
%        SPP     : (NP, NP) ceoupling between scatterer ports
%        Q       : (NP, :)  scatterer port data
%        r0      : (NP, 1) initial guess
%        target  : (ND, :)  target for driven ports (default +Inf)
%        options : struct
%
% OUTPUT  x_opt   : (:,1) optimized vector
%         x_cost  : (1,1) final cost function value
%         info    : struct, manopt info
% ------------------------------------------------------------------------
% 06.06.2023 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    C           (:,:)
    SDP         (:,:)
    SPP         (:,:)
    Q           (:,:)
    r0          (:,1)
    target      (:,:)
    options     (:,:) = []
end

if ~isfield(options, 'COST_DIFF_TOL')
    options.COST_DIFF_TOL = 0.01;
end

if ~isfield(options, 'MAX_PHASE_ITER')
    options.MAX_PHASE_ITER = 100;
end

ropt = r0;

cost_last = 0;
cost_diff = 2*options.COST_DIFF_TOL;
kit = 1;
while (cost_diff >= options.COST_DIFF_TOL) && (kit <= options.MAX_PHASE_ITER)
    % X-step
    X = exp(1j*angle( C + SDP*( ((diag(1./ropt)) - SPP) \ Q ) ));

    % a-step
    Y = C - target .* X;
    [ropt, ~, ~] = minimize_scateq_fronorm_manopt(Y, SDP, SPP, Q, ropt, options, 0);

    cost = norm(abs(SDP*((diag(1./ropt) - SPP)\Q) + C) - target, 'fro')^2;

    % Cost diff
    cost_diff = abs(cost - cost_last) / abs(cost);
    if options.verbosity > 0
        fprintf('Iteration %d, Manopt cost difference = %f \n \n', kit, cost_diff);
    end
    
    cost_last = cost;
    kit = kit+1;
end

end


function [x_opt, cost_opt, info] = optimize_scateq_MMF_manopt(C, SDP, SPP, Q, r0, target, options, suboptions)
%
% Optimize the scattering equation using the magnitude minimax fitting
% (MMF) formulation and optimization on Riemannian manifold. Use Riemannian
% Augmented Lagrangiand Method (RALM)
%-------------------------------------------------------------------------
% INPUT  C      : (ND,:)    driven port initial data
%        SDP    : (ND, NP) coupling between driven and scatterer ports
%        SPP    : (NP, NP) ceoupling between scatterer ports
%        Q      : (NP, :)  scatterer port data
%        r0     : (NP, 1) initial guess
%        target : (ND, :)  target for driven ports (default +Inf)
%        options : struct   RALM outer solver optimization settings
%        suboptions : struct, RALM sub-solver settings
%
% OUTPUT  x_opt          : struct optimal reflection coefficients
%         cost_opt          : (1,1) last objective value
%         info          : struct, optimizer output at iterations
% ------------------------------------------------------------------------
% 06.06.2024 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    C   (:,:)
    SDP     (:,:)
    SPP     (:,:)
    Q       (:,:)
    r0    (:,1)
    target  (:,:) = Inf
    options     (:,:) = []
    suboptions (:,:) = []
end

function [egrad] = scateq_egrad(C, SDP, SPP, Q, x)
    % Compute the Euclidean gradient of
    % f(x) = || C + P (diag(x)^{-1} - S)^{-1} Q ||_F^2
    M = diag(1./x) - SPP;
    MinvQ = M \ Q; % = M^{-1} * Q
    PMinv = ( (M.') \ (SDP.') ).'; % = P * M^{-1}
    alpha = -1 ./ (imag(x) - 1j*real(x)).^2;
    G = C + SDP * MinvQ;
    egrad = 2*conj( alpha .* diag( MinvQ * G' * PMinv ));
end

clear problem

NP = length(r0);
ND = length(C(:,1));
L = length(C(1,:));

% Since the manifold is a product manifold, cost function, gradients,
% optimization variables etc must be structs with fields A and B. The field
% A refers to the "first" manifold (R) and B refers to the "second"
% manifold (C^NP)

cost = @(x) x.A; % cost function with respect to struct variable x
egrad.A = @(x) 1; % euclidean gradient of the cost function with respect to struct variable's A-part
egrad.B = @(x) zeros(NP,1);  % euclidean gradient of the cost function with respect to struct variable's B-part

% Initialize equality constraints
eq_const = {}; % cell array of inequality constraint function handles
eq_const_egrad.A = {}; % cell array of inequality constraint function's gradient function handles (A-part)
eq_const_egrad.B = {}; % cell array of inequality constraint function's gradient function handles (B-part)

% Initialize inequality constraints, separately for low and upper bounds
% because -t <= g(x) <= t induce two sort of inequality constraints. The
ineq_const_low = cell(ND, L);
ineq_const_low_egrad.A = cell(ND, L);
ineq_const_low_egrad.B = cell(ND, L);

ineq_const_hi = cell(ND, L);
ineq_const_hi_egrad.A = cell(ND, L);
ineq_const_hi_egrad.B = cell(ND, L);

for jt = 1:ND

    for it = 1:L

        ineq_const_low{jt, it} = @(x) - abs( C(jt, it) + SDP(jt, :) * (( diag(1./(x.B)) - SPP ) \ Q(:,it))  )^2 + target(jt, it) - x.A;
        ineq_const_hi{jt, it} = @(x) abs(C(jt, it) + (SDP(jt, :) * (( diag(1./(x.B)) - SPP ) \ Q(:,it)) ))^2  - target(jt, it) - x.A;
    
        ineq_const_low_egrad.A{jt, it} = @(x) -1;
        ineq_const_low_egrad.B{jt, it} = @(x) -scateq_egrad(C(jt,it), SDP(jt,:), SPP, Q(:,it), x.B);
    
        ineq_const_hi_egrad.A{jt, it} = @(x) -1;
        ineq_const_hi_egrad.B{jt, it} = @(x) scateq_egrad(C(jt,it), SDP(jt,:), SPP, Q(:,it), x.B);
    end

end

% Reshape to one-column cell arrays and combine the low and hi inequality
% constraints
ineq_const = [reshape(ineq_const_low, [], 1); reshape(ineq_const_hi, [], 1)  ];
ineq_const_egrad.A = [reshape(ineq_const_low_egrad.A, [], 1); reshape(ineq_const_hi_egrad.A, [], 1)];
ineq_const_egrad.B = [reshape(ineq_const_low_egrad.B, [], 1); reshape(ineq_const_hi_egrad.B, [], 1)];

% Manifold construction
manifold_struct.A = euclideanfactory(1);
manifold_struct.B = complexcirclefactory(NP);
problem.M = productmanifold(manifold_struct);

problem.cost = cost;
problem.ineq = ineq_const;
problem.eq = eq_const;

problem.egrad_cost = egrad;
problem.egrad_ineq = ineq_const_egrad;
problem.egrad_eq = eq_const_egrad;

% Initial guess
x0.B = r0;
EEP0 = abs( C + SDP * (( diag(1./(x0.B)) - SPP ) \ Q)  ).^2;
x0.A = max(abs(target-EEP0), [], 'all');

lambda0 = options.lambda0_magnitude * ones(length(problem.ineq), 1);
gamma0 = [];

% Solve using RALM
[x_opt, cost_opt, info] = ralm(problem, x0, lambda0, gamma0, options, suboptions);

end