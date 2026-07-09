function outliers = find_outliers(data)
Q1  = quantile(data, 0.25);
Q3  = quantile(data, 0.75);
IQR = Q3 - Q1;

lower_fence = Q1 - 1.5 * IQR;
upper_fence = Q3 + 1.5 * IQR;
outliers = data < lower_fence | data > upper_fence;
end
