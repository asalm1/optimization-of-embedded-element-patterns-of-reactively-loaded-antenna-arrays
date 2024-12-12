function [w_theta, w_phi] = wrap2sphere(theta, phi)
%
% Wrap fuckful theta and phi coordinates (-inf...inf) to normal theta
% (0...pi) and phi (0...2*pi)
%-------------------------------------------------------------------------
% INPUT  theta,     (:,1), fuckful theta vector
%        phi,       (:,1), fuckful phi vector
%
% OUTPUT m_theta,     (:,:), kosher theta
%        m_phi,       (:,:), kosher phi
% ------------------------------------------------------------------------
% 09.02.2023 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    theta   (:,1)
    phi     (:,1)
end

% Transform to cartesian
x = sin(theta).*cos(phi);
y = sin(theta).*sin(phi);
z = cos(theta);

% Back to spherical
back_theta = atan2(sqrt(x.^2 + y.^2) , z); % [-pi, pi]
back_phi = atan2(y , x); % [-pi, pi]

% Find indices
theta_neg = (back_theta < 0) .* (back_phi >= 0);
phi_neg = (back_phi < 0) .* (back_theta >= 0);
both_neg = (back_phi < 0) .* (back_theta < 0);
both_pos = (back_phi >= 0) .* (back_theta >= 0);

% Check that all dots found
sumidx = theta_neg + phi_neg + both_neg + both_pos;
prodidx = theta_neg .* phi_neg .* both_neg .* both_pos;

if any(sumidx==0) || any(prodidx==1)
    error("Not all indices found")
end

idx_theta_neg = find(theta_neg==1);
idx_phi_neg = find(phi_neg==1);
idx_both_neg = find(both_neg==1);
idx_both_pos = find(both_pos==1);

% Initialize
w_theta = -66*ones(size(theta));
w_phi = -66*ones(size(phi));

% Both positive
w_theta(idx_both_pos) = back_theta(idx_both_pos);
w_phi(idx_both_pos) = back_phi(idx_both_pos);

% Both negative
w_theta(idx_both_neg) = -back_theta(idx_both_neg);
w_phi(idx_both_neg) = back_phi(idx_both_neg) + pi;

% Phi negative
w_theta(idx_phi_neg) = back_theta(idx_phi_neg);
w_phi(idx_phi_neg) = back_phi(idx_phi_neg) + 2*pi;

% Theta negative
w_theta(idx_theta_neg) = back_theta(idx_theta_neg);
w_phi(idx_theta_neg) = back_phi(idx_theta_neg) + 2*pi;

if any(w_theta==-66) || any(w_phi==-66)
    error('Nyt meni joku päin helvettiä...')
end

end