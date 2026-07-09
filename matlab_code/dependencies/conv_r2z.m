function z = conv_r2z(r)
% conv_r2z  Fisher r-to-z transform (elementwise).
%   z = conv_r2z(r) returns the Fisher variance-stabilizing transform
%   z = atanh(r) = 0.5*log((1+r)./(1-r)) of a correlation coefficient (or
%   array of them). Operates elementwise; NaNs pass through unchanged.
%   Inputs of exactly +/-1 map to +/-Inf, so callers should NaN out
%   self-correlations (diagonal 1s) beforehand if averaging. Inverse of conv_z2r.
z = atanh(r);
end
