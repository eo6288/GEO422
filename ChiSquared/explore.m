A = 2;
T = 2; % takes 2 seconds to go all the way around, frequency, we only see 1/2 of the 
% plot per unit time
t = 0:.0625:1;
phase = 2 * pi * t / T;

exponential = A * exp(1i * phase);
sinusoidal = A * cos(phase) + A * 1i * sin(phase);

[ah, ha, H] = krijetem(subnum(1, 2));
axes(ah(1))
plot(exponential, 'o')
title("exponential")
xlim([-2.5 2.5])
ylim([-2.5 2.5])

axes(ah(2))
plot(sinusoidal, 'o')
title("sinusoidal")
xlim([-2.5 2.5])
ylim([-2.5 2.5])