
rng default   % reproducibility

M = 500;

X = cell(1,10);

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


goldenrod = [0.85 0.65 0.13];

figure
tiledlayout(5,2)

S = zeros(1, M);

for k = 1:10
    S = S + X{k}; % progressive sum
    
    % normalize (CLT form)
    S_norm = (S - mean(S)) / std(S);
    
    nexttile
    histogram(S_norm, 40, ...
        'Normalization','pdf', ...
        'FaceColor', goldenrod, ...
        'EdgeColor','none')
    ylim([0 1])
    title(sprintf('Sum of %d distributions', k))
    %xlabel('Normalized sum')
    ylabel('PDF')
end

disp(mean(X{1}))