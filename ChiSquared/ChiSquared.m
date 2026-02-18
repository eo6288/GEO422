N = 3; % determines "skewness" of the resulting chi-squared
sigma = 2;
mu = 5;


%number of trials (scalar values of (N-1)* var_est / var to realize
M = 10000;
Y = []; %place to collect the data 

for i = 1:M
    X = mu + sigma * rand(N, 1);
    variance = (N-1) * var(X, 1) / (sigma ^2);
    Y(i) = variance;
end


[ah, ha, H] = krijetem(subnum(1,1));
axes(ah(1))
histogram(Y)
title("Distribution of scaled ratio of variances")