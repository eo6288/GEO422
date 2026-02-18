rng default   % reproducibility

M = 500; %each distribution has 500 samples
X = cell(1,10); %place to store all the distributions

X{1} = random("Exponential", 1, [1, M]); 
X{2} = random("Uniform", 0, 1, [1, M]);  
X{3} = random('Binomial', 5, 0.5, [1, M]); 
X{4} = random("Poisson", 5, [1, M]);
X{5} = random("Exponential", 1, [1, M]);  
X{6} = random("Noncentral Chi-square", 3, 2, [1, M]);
X{7} = random("Lognormal", 0, 0.25, [1, M]); 
X{8} = random("Uniform", -0.5, 0.5, [1, M]);
X{9} = random("Beta", 2, 5, [1, M]);       
X{10} = random("Rayleigh", 1, [1, M]);  

clf;
figure
tiledlayout(5,2)

S = zeros(1, M); %this will store the running sum of all the distributions
S_means = 0;
S_vars = 0;


for k = 1:10
    S_means = S_means + k * mean(X{k});
    S_vars = S_vars + k^2 * std(X{k})^2;

    S = S + k * X{k}; % progressive sum, I will use the k value arbitrarily as a coefficient
    
    % normalize?
    %S_norm = (S - mean(S)) / std(S);
    
    nexttile
    histogram(S, 40)
    %ylim([0 1]) %use this line if normalizing
    ylim([0 50])
    xlim([0 180])
    title(sprintf('%d distributions', k))
    ylabel('PDF')
end

disp("totalled mean of each distribution: " + S_means)
disp("mean of plot 10: " + mean(S))

disp("totalled variance of each distribution: "+ S_vars)
disp("variance of plot 10: " + std(S)^2)