dx = 0.01; % spacing
x0 = -0.5:0.01:0.5; %grid

% % Uniform convolved 5 times
% p0 = ones(size(x0));
% p0 = p0 / trapz(x0, p0); %force distribution to sum to 1

% Exponential distributions
p0 = exp(-x0);
p0 = p0 / trapz(x0, p0); %force distribution to sum to 1

mu0 = trapz(x0, x0 .* p0);
var0 = trapz(x0, (x0 - mu0).^2 .* p0); 
disp(mu0)
disp(var0)

p = p0;
x = x0;

for n = 2:10
    p = conv(p, p0) * dx;
    x = linspace(x(1)+x0(1), x(end)+x0(end), length(p));
end

p = p / trapz(x, p); %normalize again

mu = trapz(x, x .* p); %should be n*mu0
var = trapz(x, (x - mu).^2 .* p); %should be n*var0

disp(mu)
disp(var)

[ah, ha, ~] = krijetem(subnum(1, 2));

axes(ah(1))
plot(p0)
ylim([])
title("Starting Exponential Distribution")

axes(ah(2))
plot(p)
title("Distribution after 10 Convolutions")