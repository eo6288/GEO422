clear; clf;

N = 16;
unit_time = 0.0566;
t = 0 : unit_time : (N - 1)*unit_time; %seconds
m_true = [10; 500; -981];  % made up truth

G = [ones(N,1) t' 0.5*t'.^2]; % design matrix
y_true = G * m_true;

% noise grows over time
sigma_t = 2 + 2 * (0:N-1)';
Sigma_base = toeplitz(0.8.^(0:N-1));
Sigma = diag(sigma_t) * Sigma_base * diag(sigma_t);

%this is a weighted least squares
plot(t, y_true)
hold on
for i = 1:100
y_noisy = mvnrnd(y_true, Sigma)';
plot(t, y_noisy, 'o', "MarkerFaceColor", [0.7 0.7 0.7], "MarkerEdgeColor", [0.7 0.7 0.7])
end
hold on
title("100 Realizations (full, time-dependent covariance matrix)")
hh=errorxy(t,y_true,[],sqrt(diag(Sigma))); %one standard deviation
legend(hh, "1 standard deviation")
hold off

figure;
imagesc(Sigma);
title("Covariance matrix")
colorbar

m_est = G \ y_noisy;

residual = y_noisy - G * m_est;
G_g = inv(G'*inv(Sigma)*G) * G' * inv(Sigma);
Cov_m = G_g * Sigma * G_g';

% %plot residuals
% figure;
% hold on
% plot(t, residual)
% hold off
