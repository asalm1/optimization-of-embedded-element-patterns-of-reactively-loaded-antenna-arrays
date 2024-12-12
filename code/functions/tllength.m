function len = tllength(r, beta, Z_line, Z_ref, ending)
%
% Calculate length of transmission line for realizing a given reflection
% coefficient
%-------------------------------------------------------------------------
% INPUT  rho    : (N,1) diagonal matrix of reflection coefficients
%        beta   : (N,N) Diagonal matrix of propagation constants
%        ZC     : (N,N) Tline characteristic impedance
%        Z0     : (N,N) Port impedance
%        l1     : (N,N) Length from port to end, initially
%        endint : string, open or short
%
% OUTPUT len : (N,N) Diegonal matrix of slot line lengths
% ------------------------------------------------------------------------
% 04.08.2022 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%
arguments
    r             (:,1)
    beta            (:,1)
    Z_line          (:,1)
    Z_ref           (:,1)
    ending          (1,1) string
end

if ending == "short"
    len = atan( (Z_ref./(1j*Z_line)) .* (1+r) ./ (1-r) );
elseif ending == "open"
    len = 1j*acoth((Z_ref.*(r+1)) ./ (Z_line.*(r-1)));
else
    error("ending must be open or short")
end


len = len ./ beta;

% len(len < 0) = len(len<0) + (pi./beta(len<0));


end