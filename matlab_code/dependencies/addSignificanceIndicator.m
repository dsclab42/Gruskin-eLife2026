% Define function to add significance indicators
function addSignificanceIndicator(barGroup1, barGroup2, barIndex1, barIndex2, pValue, maxY)
significanceThreshold = .05;
if pValue < significanceThreshold
    x1 = barGroup1.XEndPoints(barIndex1);
    x2 = barGroup2.XEndPoints(barIndex2);
    y = maxY; % Use maximum Y value for placement
    plot([x1, x2], [y, y], 'k', 'LineWidth', 1); % Line
    plot(mean([x1, x2]), y * 1.05, 'k*', 'LineWidth', 1); % Asterisk
end
end