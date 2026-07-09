% Custom function to generate the white-red colormap
function c = whitered(m)
    if nargin < 1, m = size(get(gcf,'colormap'),1); end
    bottom = [1 1 1]; top = [1 0 0]; middle = [.5 0 0];
    % Interpolate between colors
    interp = @(top, bottom, m) [linspace(bottom(1), top(1), m)', linspace(bottom(2), top(2), m)', linspace(bottom(3), top(3), m)'];
    c = [interp(middle, bottom, ceil(m/2)); interp(top, middle, floor(m/2))];
end