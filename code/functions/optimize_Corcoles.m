function [aopt, Aopt, cost, elapsed_time] = optimize_Corcoles(S, E, feds, pars, scanidx, outidx, SLL, options)
%
% Optimize incident waves of the scatterer ports using Corcoles
% formulation: https://doi.org/10.1109/TAP.2015.2478487. This optimizes
% both driven port waves and scatterer port waves. Uses semidefinite
% relaxation
%-------------------------------------------------------------------------
% INPUT  S          : (N, N) S-parameters
%        E          : (N, L) Embedded element patterns
%        feds, pars : (:,1) indices of driven and scatterer ports
%        scanidx    : (:,1) indices of E referring to scan directions
%        outidx     : (:,1) indices of E referring to out directions
%        SLL        : (:,1) maximum sidelobe levels (linear scale)
%        options    : struct, optimization options
%        
% OUTPUT  a             : (NP, 1) optimal incident waves to scatterer ports
%         A             : (NP, NP) optimizer variable, A ~ a a^T
%         cost          : (1,1) last objective value
%         elapsed time  : string
%         info          : struct, optimizer output at iterations
% ------------------------------------------------------------------------
% 07.12.2023 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    S       (:,:)
    E       (:,:)
    feds    (:,1)
    pars    (:,1)
    scanidx (:,1)
    outidx  (:,1) = []
    SLL     (:,1) = []
    options     (:,:) = []
end

tstart = tic;

if ~isfield(options, 'verbosity')
    options.verbosity = 1;
end

NP = length(pars); ND = length(feds); N = NP+ND; 
Lout = length(outidx);

if length(SLL) ~= length(outidx)
    SLL = SLL(1) .* ones(size(outidx));
end

% Normalize EEPs
Enorm = zeros(size(E));
for it = 1:N
    Enorm = E(it, :) ./ max(abs(E(it, :)));
end
% Enorm = E;

if options.verbosity > 0
    cvx_begin sdp
else
    cvx_begin sdp quiet
end

    variable A(N, N) hermitian complex
    minimize( trace( (eye(N) - S'*S) * A ) );

    subject to
        trace(conj(E(:,scanidx)) * (E(:,scanidx).') * A) == 1; % Squared

        % Sidelobe constraints
        for m = 1:Lout
            Fm = conj(Enorm(:,outidx(m))) * Enorm(:,outidx(m)).';
            real(trace( Fm * A)) - SLL(m) <= 0;
        end

        % Reactivity constraints
        for l = 1:NP
            ul = zeros(N,1);
            ul(pars(l)) = 1;
            Fl = ul*ul.' - S(pars(l), :)' * S(pars(l), :);
            trace( Fl * A) == 0;
        end

        A >= 0;

cvx_end

% Data extraction
Aopt = A;
[evecs, evals] = eig(Aopt);
[maxeval, maxidx] = max(diag(evals));
aopt = sqrt(maxeval) * evecs(:, maxidx);

cost = cvx_optval;

elapsed_time = toc(tstart);
if options.verbosity > 0
    if elapsed_time > 60
        elapsed_string = seconds(elapsed_time);
        elapsed_string.Format = 'hh:mm:ss';
    else
        elapsed_string = string(elapsed_time);
    end
    
    fprintf('Total elapsed time in Corcoles optimization: %s \n \n', elapsed_string);
end

end