clear; clc;

f0 = 5; %frequency in signal
T  = 5; % the signal will run for 5 seconds

%Nyquist frequency is 1/2*5 = 0.1
delta_x = [0.05 0.15 0.3]; % above, near, below Nyquist

for k = 1:length(delta_x)

    dt = delta_x(k);
    t = 0:dt:T;
    fs = 1/dt;

    x = sin(2*pi*f0*t);

    figure
    subplot(2,1,1)
    plot(t,x)
    title(['Time series, dt = ', num2str(dt)])

    subplot(2,1,2)
    %estimate frequency content of signal
    periodogram(x,[],[],fs)
end

figure;
t = 0:0.001:T;
plot(t,sin(2*pi*f0*t))
title("true signal (frequency of 5)")