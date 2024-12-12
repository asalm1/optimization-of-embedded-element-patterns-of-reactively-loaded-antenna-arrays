function [x_opt, cost_opt, info] = ralm(problem, x0, lambda0, gamma0, options, suboptions)
%
% Riemannian Augmented Lagrangian Method. Based on 
% C. Liu and N. Boumal, “Simple algorithms for optimization on riemannian 
% manifolds with constraints,” Applied Mathematics & Optimization,
% vol. 82, pp. 949–981, 2020.
%-------------------------------------------------------------------------
% INPUT  problem        : struct, problem data
%        x0             : struct or array, initial quess
%        lambda0        : (:,1), intial guess of inequality constraint
%                                multipliers
%        gamma0         : (:,1), intial guess of equality constraint
%                                multipliers
%
% OUTPUT  x_opt         : struct or array, resulting point
%         cost_opt      : (1,1) final objective function value
%         info          : struct, optimizer output at iterations
% ------------------------------------------------------------------------
% 06.06.2024 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    problem     struct
    x0          (:,:)   = []
    lambda0     (:,1)   = []
    gamma0      (:,1)   = []
    options     struct  = []
    suboptions  struct  = []
end

elapsed_time = tic;

% Local defaults for the outer solver
localdefaults.maxtime = Inf;
localdefaults.miniter = 0;
localdefaults.maxiter = 1000;
localdefaults.eps_min = 1e-6;
localdefaults.eps0 = 1e-3;
localdefaults.rho0 = 1;
localdefaults.theta_eps = (localdefaults.eps_min/localdefaults.eps0)^(1/30);
localdefaults.theta_rho = 3.3;
localdefaults.theta_sigma = 0.8;

localdefaults.lambda_min_magnitude = 0;
localdefaults.lambda_min = localdefaults.lambda_min_magnitude * ones(length(problem.ineq), 1);
localdefaults.lambda_max_magnitude = 1e6;
localdefaults.lambda_max = localdefaults.lambda_max_magnitude * ones(length(problem.ineq), 1);
localdefaults.lambda0_magnitude = 10;

localdefaults.gamma_min_magnitude = 0;
localdefaults.gamma_min = localdefaults.gamma_min_magnitude * ones(length(problem.eq), 1);
localdefaults.gamma_max_magnitude = 1e6;
localdefaults.gamma_max = localdefaults.gamma_max_magnitude * ones(length(problem.eq), 1);
localdefaults.gamma0_magnitude = 10;

localdefaults.d_min = 1e-10;
localdefaults.subsolver = @rlbfgs;

% Local defaults for subsolver
subdefaults.maxtime = Inf;
subdefaults.miniter = 10;
subdefaults.maxiter = 200;
subdefaults.verbosity = 1;

% Merge outer solver options
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults, options);

% Merge subsolver options
subdefaults = mergeOptions(getGlobalDefaults(), subdefaults);
if ~exist('suboptions', 'var') || isempty(suboptions)
    suboptions = struct();
end
suboptions = mergeOptions(subdefaults, suboptions);

% Init primal and dual variables if empty
if isempty(x0)
    x0 = problem.M.rand();
end
if isempty(lambda0)
    lambda0 = options.lambda0_magnitude * ones(length(problem.ineq), 1);
end
if isempty(gamma0)
    gamma0 = options.gamma0_magnitude * ones(length(problem.eq), 1);
end

% Initialize data storage
x = cell(1,options.maxiter);
lambda = cell(1,options.maxiter);
gamma = cell(1,options.maxiter);
sigma = cell(1,options.maxiter);
eps = cell(1,options.maxiter);
rho = cell(1,options.maxiter);

% Append intial guess to data storage
x{1} = x0;
lambda{1} = lambda0;
gamma{1} = gamma0;
sigma{1} = update_sigma(problem.ineq, x0, lambda0, options.rho0);
rho{1} = options.rho0;
eps{1} = options.eps0;

% Start optimization
k = 1; % iteration number
while true
   
    % Liu2020 step 4
    x{k+1} = update_x(x{k}, rho{k}, lambda{k}, gamma{k}, eps{k}, problem, suboptions, options.subsolver);

    % Liu2020 step 8
    gamma{k+1} = update_gamma(options.gamma_min, options.gamma_max,...
                                gamma{k}, rho{k}, problem.eq, x{k+1});

    lambda{k+1} = update_lambda(options.lambda_min, options.lambda_max,...
                                lambda{k}, rho{k}, problem.ineq, x{k+1});

    % Liu2020 step 10
    sigma{k+1} = update_sigma(problem.ineq, x{k+1}, lambda{k}, rho{k});

    % Liu2020 step 11
    eps{k+1} = update_epsilon(options.eps_min, options.theta_eps, eps{k});

    % Liu2020 step 12-16
    rho{k+1} = update_rho(k, problem.eq, x{k+1}, sigma{k+1}, ...
                            options.theta_sigma, options.theta_rho, x{k}, sigma{k}, rho{k});

    % Liu2020 step 5-7
    if stop_now(k, x, eps, toc(elapsed_time), options, problem.M)
        break;
    end

    k = k+1;
end

% Remove zeros
info.x = x(1:k+1);
info.time = toc(elapsed_time);
info.gamma = gamma(1:k+1);
info.lambda = lambda(1:k+1);
info.sigma = sigma(1:k+1);
info.eps = eps(1:k+1);
info.rho = rho(1:k+1);

x_opt = x{k+1};
cost_opt = problem.cost(x_opt);

end

% Update functions
function xk1 = update_x(xk, rhok, lambdak, gammak, epsk, problem, suboptions, solver)
% Liu2020 step 4. Minimize the augmented Lagrangian function with respect
% to x on the given manifold.

    function L = augmented_lagrangian_fun(cost, ineq, eq, rho, lambda, gamma, x)
    % Augmented Lagrangian function L_{rho_k}(x, lambda^k, gamma^k)
    
        sum_eq = 0;
        for jt = 1:length(eq)
            sum_eq = sum_eq + ( eq{jt}(x) + gamma(jt)/rho )^2;
        end
        
        sum_ineq = 0;
        for it = 1:length(ineq)
            sum_ineq = sum_ineq + max([0, lambda(it)/rho + ineq{it}(x)])^2;
        end
        
        L = cost(x) + (rho/2) * (sum_eq + sum_ineq);
    
    end

    function grad_L = augmented_lagrangian_egrad(ineq, eq, egrad_cost, egrad_ineq, egrad_eq, rho, lambda, gamma, x)
    % Euclidean gradient of the Augmented Lagrangian function, for product
    % manifolds (struct variables)

        function grad_L_nostruct = augmented_lagrangian_egrad_nostruct(ineq, eq, egrad_cost, egrad_ineq, egrad_eq, rho, lambda, gamma, x)
        % Euclidean gradient of the Augmented Lagrangian function, for simple
        % manifolds
         
            sum_eq = 0;
            if ~isempty(eq)
                sum_eq = (eq{1}(x) + gamma(1)/rho) * egrad_eq{1}(x);
                for jt = 2:length(eq)
                    sum_eq = sum_eq + (eq{jt}(x) + gamma(jt)/rho) * egrad_eq{jt}(x);
                end
            end
            
            sum_ineq = 0;
            if ~isempty(ineq)
                sum_ineq = max([0,lambda(1)/rho+ineq{1}(x)])*egrad_ineq{1}(x);
                for jt = 2:length(ineq) 
                    sum_ineq = sum_ineq + max([0,lambda(jt)/rho+ineq{jt}(x)])*egrad_ineq{jt}(x);
                end
            end
            
            grad_L_nostruct = egrad_cost(x) + rho * (sum_eq + sum_ineq);
            
        end

        if isstruct(x)
            fnames = fieldnames(x);
            for it = 1:length(fnames)
                grad_L.(fnames{it}) = augmented_lagrangian_egrad_nostruct(ineq, eq, ...
                                        egrad_cost.(fnames{it}), egrad_ineq.(fnames{it}), egrad_eq.(fnames{it}),...
                                        rho, lambda, gamma, x);
            end
        else
            grad_L = augmented_lagrangian_egrad_nostruct(ineq, eq, egrad_cost, egrad_ineq, egrad_eq, rho, lambda, gamma, x);
        end
        
    end

    subproblem.cost = @(x) augmented_lagrangian_fun(problem.cost, problem.ineq, problem.eq,...
                                                     rhok, lambdak, gammak, x);

    subproblem.egrad = @(x) augmented_lagrangian_egrad(problem.ineq, problem.eq, ...
                                                       problem.egrad_cost, problem.egrad_ineq, problem.egrad_eq,...
                                                       rhok, lambdak, gammak, x);

    subproblem.M = problem.M;

    suboptions.tolgradnorm = epsk;
    [xk1, ~, ~] = solver(subproblem, xk, suboptions);
end

function gammak1 = update_gamma(gamma_min, gamma_max, gammak, rhok, h, xk1)
% Liu2020 step 8. Update multipliers of equality constraints
    gammak1 = zeros(size(gammak));
    for it = 1:length(gammak1)
        gammak1(it) = clip(gamma_min(it), gamma_max(it), ...
                        gammak(it) + rhok*h{it}(xk1) );
    end
end

function lambdak1 = update_lambda(lambda_min, lambda_max, lambdak, rhok, g, xk1)
% Liu2020 step 9. Update multipliers of inequality constraints
    lambdak1 = zeros(size(lambdak));
    for it = 1:length(lambdak1)
        lambdak1(it) = clip(lambda_min(it), lambda_max(it), ...
                        lambdak(it) + rhok*g{it}(xk1) );
    end
end

function sigmak1 = update_sigma(g, xk1, lambdak, rhok)
% Liu2020 step 10. Update sigma-coefficients
    sigmak1 = zeros(length(g), 1);

    for it = 1:length(g)
        sigmak1(it) = max([g{it}(xk1), -lambdak(it)/rhok]);
    end
end

function epsk1 = update_epsilon(eps_min, theta_eps, epsk)
% Liu2020 step 11. Update accuracy tolerance
    epsk1 = max([eps_min, theta_eps*epsk]);
end

function rhok1 = update_rho(k, h, xk1, sigmak1, theta_sigma, theta_rho, xk, sigmak, rhok)
% Liu2020 step 8. Update penalty parameter

    % Maximum of equality constraint functions
    if isempty(h)
        max_abs_h_k1 = 0;
        max_abs_h_k = 0;
    else
        abs_h_k1 = zeros(size(h));
        abs_h_k = zeros(size(h));
        for jt = 1:length(h)
            abs_h_k1(jt) = abs(h{jt}(xk1));
            abs_h_k(jt) = abs(h{jt}(xk));
        end
    
        max_abs_h_k1 = max(abs_h_k1);
        max_abs_h_k = max(abs_h_k);
    end
    
    % Maximum of sigmas
    if isempty(sigmak)
        max_abs_sigma_k = 0;
        max_abs_sigma_k1 = 0;
    else
        max_abs_sigma_k = max(abs(sigmak));
        max_abs_sigma_k1 = max(abs(sigmak1));
    end
    
    max_abs_k = max([max_abs_h_k, max_abs_sigma_k]);
    max_abs_k1 = max([max_abs_h_k1, max_abs_sigma_k1]);
    
    if (k==1) || ( max_abs_k1 <= theta_sigma * max_abs_k )
        
        rhok1 = rhok;
    
    else
        rhok1 = theta_rho*rhok;
    end

end

function s = stop_now(k, x, eps, t, options, M)
% Terminate iteration if following conditions: 

    if k < options.miniter
        s = false;
        return
    end

    if (M.dist(x{k+1}, x{k}) < options.d_min) && (eps{k} <= options.eps_min)

        s = true;
        return

    end

    if k >= options.maxiter
        s = true;
        return
    end

    if t >= options.maxtime
        s = true;
        return
    end

    s = false;
    return
end

function c = clip(a,b,x)
    c = max(a, min(b,x));
end