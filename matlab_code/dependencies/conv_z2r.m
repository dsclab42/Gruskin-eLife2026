function r = conv_z2r(z)
% conv_z2r  Fisher z-to-r transform (elementwise).
%   r = conv_z2r(z) returns r = tanh(z), the inverse of the Fisher r-to-z
%   transform (conv_r2z). Operates elementwise; NaNs pass through unchanged.
%   Typically used to convert a mean-of-z back to a correlation for reporting
%   or visualization.
r = tanh(z);
end
