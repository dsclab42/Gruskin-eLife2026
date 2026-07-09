function [p_value, jack_se1,jack_se2] = compareHeritability(h2_1, se_1, h2_2, se_2)
n_fam    = length(se_1);
A        = (n_fam-1)/n_fam;
B1       = sum((se_1 - sum(se_1/n_fam)).^2);
h2_var1  = A*B1;
jack_se1 = sqrt(h2_var1);

B2       = sum((se_2 - sum(se_2/n_fam)).^2);
h2_var2  = A*B2;
jack_se2 = sqrt(h2_var2);

% Calculate the Z-score
z_score = (h2_1 - h2_2) / sqrt(jack_se1^2 + jack_se2^2);

% Calculate the p-value for the two-tailed test
p_value = 2 * (1 - normcdf(abs(z_score)));
end
