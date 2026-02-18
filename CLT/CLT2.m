clear; clf;
set(groot, "defaultTextInterpreter", "latex")

dx = 0.1;
x0 = -5:dx:5;
% 
% %create a uniform distribution
% p0 = ones(size(x0));
% X = p0 / trapz(x0,p0);

%create an exponential distribution
p0 = exp(-x0);
X = p0 / trapz(x0, p0); %force distribution to sum to 1

% make the x axis grid points
x2 = 2*x0(1):dx:2*x0(end);
x3 = 3*x0(1):dx:3*x0(end);
x4 = 4*x0(1):dx:4*x0(end);
x5 = 5*x0(1):dx:5*x0(end);
x6 = 6*x0(1):dx:6*x0(end);

% successive convolutions
f1 = X;
f2 = conv(X,X) * dx;
f2 = f2 / trapz(x2,f2);

f3 = conv(f2,X) * dx;
f3 = f3 / trapz(x3,f3);

f4 = conv(f3,X) * dx;
f4 = f4 / trapz(x4,f4);

f5 = conv(f4,X) * dx;
f5 = f5 / trapz(x5,f5);

f6 = conv(f5,X) * dx;
f6 = f6 / trapz(x6,f6);


% calculate moments to be able to plug into normpdf()
mu_2 = trapz(x2, x2 .* f2);
sigma_2 = sqrt(trapz(x2, (x2 - mu_2).^2 .* f2));

mu_3 = trapz(x3, x3 .* f3);
sigma_3 = sqrt(trapz(x3, (x3 - mu_3).^2 .* f3));

mu_4 = trapz(x4, x4 .* f4);
sigma_4 = sqrt(trapz(x4, (x4 - mu_4).^2 .* f4));

mu_5 = trapz(x5, x5 .* f5);
sigma_5 = sqrt(trapz(x5, (x5 - mu_5).^2 .* f5));

mu_6 = trapz(x6, x6 .* f6);
sigma_6 = sqrt(trapz(x6, (x6 - mu_6).^2 .* f6));

% gaussians to plot 
g2 = normpdf(x2, mu_2, sigma_2);
g3 = normpdf(x3, mu_3, sigma_3);
g4 = normpdf(x4, mu_4, sigma_4);
g5 = normpdf(x5, mu_5, sigma_5);
g6 = normpdf(x6, mu_6, sigma_6);



[ah, ha, H] = krijetem(subnum(3, 2)); 

axes(ah(1)) 
hold on 
plot(x0, f1) 
title("Exponential Distribution") 

axes(ah(2)) 
hold on 
plot(x2, f2) 
plot(x2, g2) 
legend("convolution", "gaussian", "Location", "north")
title("x1 = conv(x, x)")
hold off 

axes(ah(3)) 
hold on 
plot(x3, f3) 
plot(x3, g3) 
title("x2 = conv(x1, x)")
hold off 

axes(ah(4)) 
hold on 
plot(x4, f4) 
plot(x4, g4) 
title("x3 = conv(x2, x)")
hold off 

axes(ah(5)) 
hold on 
plot(x5, f5) 
plot(x5, g5) 
xlim([-25 25])
title("x4 = conv(x3, x)")
hold off 

axes(ah(6)) 
hold on 
plot(x6, f6) 
plot(x6, g6) 
title("x5 = conv(x4, x)")
hold off 



