function obs_idcs = direction_indices(theta, phi, obs_theta, obs_phi)
%
% Find indices that describe the desired points on measurement sphere
%-------------------------------------------------------------------------
% INPUT  theta, phi     : (:,1) Whole measurement sphere (eg. from CST)
%        obs_theta, obs_phi   : (L,1) Measurement points to be observed
%
% OUTPUT obs_idcs : (L,1) Indices that points to observation points on
%                         original theta and phi
% ------------------------------------------------------------------------
% 15.08.2022 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    theta   (:,1)
    phi     (:,1)
    obs_theta   (:,1)
    obs_phi     (:,1)
end

if ( length(theta) ~= length(phi) ) || ( length(obs_theta) ~= length(obs_phi) )
    error('The sizes of theta and phi vectors must be equal.')
end

L = length(obs_theta);

% Wrap arbitrary spherical coordinates to (0..pi) and (0..2*pi)
[thetaS, phiS] = wrap2sphere(obs_theta, obs_phi);

% Find correct indices
obs_idcs = zeros(L,1);
for lit = 1:L
    obs_idcs(lit) = find( (abs(thetaS(lit) - theta) == min(abs(thetaS(lit) - theta))) .*...
                 ( abs(phiS(lit) - phi) == min(abs(phiS(lit) - phi) ) ), 1 );
end

end