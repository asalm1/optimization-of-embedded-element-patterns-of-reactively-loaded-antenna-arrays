function [aopt, Aopt, result, info] = optimize_scateq_SDR(C, SDP, SPP, Q, target, options)
%
% Optimize the scattering equation using semidefinite relaxation. 
% Provide target data optionally.
% Otherwise, minimize or maximize the frobenius norm of the scattering 
% equation by setting target = +- Inf. If target is given, 
% minimize maximum distance to the target. Target can be defined
% element-wise or as a sum of elements' values (EEPs or S-parameters)
% This function uses SDR and returns the raw SDR output. Require CVX being
% installed.
%-------------------------------------------------------------------------
% INPUT  C      : (ND,:)   driven port initial data
%        SDP    : (ND, NP) coupling between driven and scatterer ports
%        SPP    : (NP, NP) coupling between scatterer ports
%        Q      : (NP, :)  scatterer port data
%        target : (ND, :)  target for driven ports (default +Inf)
%        options : struct  optimization settings
%
% OUTPUT  aopt          : (ND*NP, 1) raw SDR solution, incident waves vector
%         Aopt          : (ND*NP, ND*NP) raw SDR solution, incident waves matrix
%         result        : (L, ND) raw SDR result
%         info          : struct, (cost, status, tightness)
% ------------------------------------------------------------------------
% 05.08.2024 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    C       (:,:)
    SDP     (:,:)
    SPP     (:,:)
    Q       (:,:)
    target  (:,:) = Inf
    options (:,:) = []
end

tstart = tic; % start timing

% problem dimensions
NP = length(SPP(1,:));
ND = length(C(:,1));
L = length(C(1,:));

% Set verbosity if not defined
if ~isfield(options, 'verbosity')
    options.verbosity = 1;
end

% Vectorized optimization data
bar_c = reshape(C.', ND*L, 1);  % ED or SDD
bar_Q = kron(eye(ND), Q.');    % EP or SPD
bar_p = reshape(SDP.', ND*NP, 1);
bar_SPP = kron(eye(ND), SPP);

% Begin CVX
if options.verbosity > 0
    cvx_begin sdp
else
    cvx_begin sdp quiet
end
variable A(ND*NP, ND*NP) hermitian complex % matrix variable
variable a(ND*NP) complex % vector variable

if target == Inf % maximize fro norm of scateq
    F0 = - bar_Q' * bar_Q;
    F0 = 0.5 * (F0 + F0'); % Force Hermitean symmetry

    f0 = - bar_c' * bar_Q; % row
    g0 = - bar_c' * bar_c;

    minimize( trace(F0*A) + 2*real(f0*a) + g0 ) % frobenius norm in quadratic form
    subject to % Constraints outside of if-else, coming later
      
elseif target == -Inf % minimize fro norm of scateq
    F0 = bar_Q' * bar_Q;
    F0 = 0.5 * (F0 + F0'); % Force Hermitean symmetry

    f0 = bar_c' * bar_Q; % row
    g0 = bar_c' * bar_c;

    minimize( trace(F0*A) + 2*real(f0*a) + g0 ) 
    subject to % Constraints outside of if-else

else % minimize max distance to target
    bar_target = reshape(target.', ND*L, 1); % ND x L -> ND*L x 1

    variable t(1)

    minimize(t)

    subject to

        if all(size(target) == size(C)) % Target is defined for each matrix element

            for it = 1:ND*L
                F0 = (bar_Q(it,:))' * bar_Q(it,:);
                F0 = 0.5 * (F0 + F0'); % Force Hermitean symmetry
            
                f0 = (bar_c(it))' * bar_Q(it,:); % row
                g0 = abs(bar_c(it))^2;

                % add constraint
                -t <= trace(F0*A) + 2*real(f0*a) + g0 - bar_target(it) <= t;
            end

        else % target is definited for the sum of the element results

            for n = 1:L
                u = zeros(ND,ND*L);
                for it = 1:ND
                    u(it, (it-1)*L + n) = 1;
                end

                U = u*u';

                F0 = bar_Q' * U * bar_Q;
                F0 = 0.5 * (F0 + F0'); % Force Hermitean symmetry
            
                f0 = bar_c' * U * bar_Q; % row
                g0 = bar_c' * U * bar_c;

                % add constraint
                -t <= trace(F0*A) + 2*real(f0*a) + g0 - bar_target(n) <= t;
            end

        end

end

% Default constraints for passivity and common terminations
for n = 1:ND*NP

    m = n - NP*(ceil(n/NP)-1);
    un = zeros(ND*NP,1);
    um = zeros(ND*NP,1);
    un(n) = 1;
    um(m) = 1;
    
    % Fn = bar_SPP' * um*(un.') * bar_SPP - um*(un.');
    % f1n = bar_p' * um*(un.') * bar_SPP; % row
    % f2n = bar_SPP' * um*(un.') * bar_p; % col
    % gn = bar_p' * um*(un.') * bar_p;

    Fn = (bar_SPP(m,:))' * bar_SPP(n,:) - um*(un.');
    f1n = (bar_p(m))' * bar_SPP(n,:); % row
    f2n = (bar_SPP(m,:))' * bar_p(n); % col
    gn = (bar_p(m))' * bar_p(n);
    
    % Hn = bar_SPP' * un*(un.') * bar_SPP - un*(un.');
    % hn = bar_p' * un*(un.') * bar_SPP; % row
    % mn = bar_p' * un*(un.') * bar_p;

    Hn = (bar_SPP(n,:))' * (bar_SPP(n,:)) - un*(un.');
    hn = (bar_p(n))' * (bar_SPP(n,:)); % row
    mn = abs(bar_p(n))^2;

    trace(Fn*A) + f1n*a + a'*f2n + gn == 0;
    trace(Hn*A) + 2*real(hn*a) + mn == 0;
end

[A, a; a', 1] >= 0;

cvx_end

Aopt = A;
aopt = a;

info.cost = cvx_optval;
info.status = cvx_status;
info.tightness = norm(Aopt - aopt*aopt', 'fro') / (aopt'*aopt);

% Raw SDR result
if all(size(target) == size(C)) % Target is defined for each matrix element
    result = zeros(ND*L, 1);
    for it = 1:ND*L
        F0 = (bar_Q(it,:))' * bar_Q(it,:);
        F0 = 0.5 * (F0 + F0'); % Force Hermitean symmetry
        
        f0 = (bar_c(it))' * bar_Q(it,:); % row
        g0 = abs(bar_c(it))^2;

        % add constraint
        result(it) = trace(F0*Aopt) + 2*real(f0*aopt) + g0;
    end

    result = reshape(result, L, ND);

else % total gain target
    result = zeros(L, 1);

    for n = 1:L
        u = zeros(ND,ND*L);
        for it = 1:ND
            u(it, (it-1)*L + n) = 1;
        end

        U = u*u';

        F0 = bar_Q' * U * bar_Q;
        F0 = 0.5 * (F0 + F0'); % Force Hermitean symmetry
    
        f0 = bar_c' * U * bar_Q; % row
        g0 = bar_c' * U * bar_c;

        result(n) = trace(F0*Aopt) + 2*real(f0*aopt) + g0;
    end

end

if options.verbosity > 1
    % Check the constraints
    rep_const = ones(ND*NP,1);
    cmc_const = ones(ND*NP,1);
    for n = 1:ND*NP
        m = n - NP*(ceil(n/NP)-1);
        un = zeros(ND*NP,1);
        um = zeros(ND*NP,1);
        un(n) = 1;
        um(m) = 1;
        
        Fn = bar_SPP' * um*(un.') * bar_SPP - um*(un.');
        f1n = bar_p' * um*(un.') * bar_SPP; % row
        f2n = bar_SPP' * um*(un.') * bar_p; % col
        gn = bar_p' * um*(un.') * bar_p;

        % Fn = (bar_SPP(m,:))' * bar_SPP(n,:) - um*(un.');
        % f1n = (bar_p(m))' * bar_SPP(n,:); % row
        % f2n = (bar_SPP(m,:))' * bar_p(n); % col
        % gn = (bar_p(m))' * bar_p(n);
        
        Hn = bar_SPP' * un*(un.') * bar_SPP - un*(un.');
        hn = bar_p' * un*(un.') * bar_SPP; % row
        mn = bar_p' * un*(un.') * bar_p;

        % Hn = (bar_SPP(n,:))' * (bar_SPP(n,:)) - un*(un.');
        % hn = (bar_p(n))' * (bar_SPP(n,:)); % row
        % mn = abs(bar_p(n))^2;
    
        rep_const(n) = trace(Fn*A) + f1n*a + a'*f2n + gn;
        cmc_const(n) = trace(Hn*A) + 2*real(hn*a) + mn;
    end

    disp('Constraint function values (should be zeros): \n')
    disp([rep_const, cmc_const])
end

% Show elapsed time
info.elapsed_time = toc(tstart);
if options.verbosity > 0
    if info.elapsed_time > 60
        elapsed_string = seconds(info.elapsed_time);
        elapsed_string.Format = 'hh:mm:ss';
    else
        elapsed_string = string(info.elapsed_time);
    end
    
    fprintf('Total elapsed time in SDR optimization: %s \n \n', elapsed_string);
end

end
