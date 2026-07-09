% Custom function to generate the blue-white-red colormap
function c = bluewhitered(m)
if nargin < 1, m = size(get(gcf,'colormap'),1); end
bottom = [0 0 1]; top = [1 0 0]; middle = [1 1 1];
% Interpolate between colors
interp = @(top, bottom, m) [linspace(bottom(1), top(1), m)', linspace(bottom(2), top(2), m)', linspace(bottom(3), top(3), m)'];
c = [interp(middle, bottom, ceil(m/2)); interp(top, middle, floor(m/2))];
end