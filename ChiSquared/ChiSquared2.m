clear; clf;
set(groot, "defaultTextInterpreter", "latex")

X =  random("Exponential", 3, [1, 1000]); %1.5 = 1/lambda (it is the mean)
K = 100; % number of bins
p = 1; % only 1 parameter to account for, plus 1
p_value = 0.01;

[frequencies, edges] = histcounts(X, K);

N = length(X);
expected = [];
for i = 1:K
    left = edges(i);
    right = edges(i+1);
    expected(i) = N * ( ...
        cdf("Exponential", right, 3) - ... %lambda still 10
        cdf("Exponential", left, 3));
end

valid = expected > 5; %avoids a divide by zero(or near zero)
chi2 = sum((frequencies(valid) - expected(valid)).^2 ./ (expected(valid)));

df = K - p - 1;
critical = chi2inv(1 - p_value, df);
if chi2 > critical
    disp("Reject null hypothesis")
else
    disp("Fail to reject null")
end


[ah, ha, H] = krijetem(subnum(1,3));

axes(ah(1))
bar(frequencies)
ylim([0,300])
title("sampled distribution")

axes(ah(2))
bar(expected)
ylim([0,300])
title("expected distribution")

axes(ah(3))
bar(abs(frequencies - expected))
ylim([0,300])
title("|difference|")
xpos = 0.6 * K;            
ypos = 0.9 * max(ylim);
text(xpos, ypos, ...
     sprintf("$\\chi_{%i - %i}^2 = %.3f$", K, p+1, chi2))

