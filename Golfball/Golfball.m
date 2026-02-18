% Uncomment the following if you want to collect new data points

% %get data
% URL = 'https://geoweb.princeton.edu/people/simons/GOLFBALL/';
% data = zeros(24,2); 
% 
% figure;
% for i = 31:54
%     I = imread([URL sprintf('%08d.jpg', i)]);
% 
%     clf
%     imshow(I)
% 
%     [y,t] = ginput(1);
%     data(i-30,:) = [y t];
% end
% 

%%% data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [
  596.8100  450.6845
  591.2995  443.5996
  556.6624  382.1974
  542.4926  334.9649
  483.4520  188.5443
  459.0486  143.6734
  416.5394   79.9096
  393.7103   58.6550
  329.9465   35.8260
  310.2663   32.6771
  266.1827   61.0166
  227.6095  106.6747
  192.1851  175.1617
  176.4410  217.6710
  %181.1642  212.9477
  135.5062  359.3684
  117.4004  416.0474 ];

y = data(:, 2);
y = y(end:-1:1); %reverese so that the ball starts at origin at t=0
y = (y*-1)+max(y) + 1; % the plus one is so that there are no zeros in the Cov_d 
%figure out y dimension conversion
height_ball = 43.2174 - 23.6522; % obtained through ginput
% height_ball = 19.5652, so every 19.5652 pixels = 4 cm
length_scale = 4 / height_ball;
y = y* length_scale;

% there are 25 fps and 16 photos, so we will say each photo is 1/25th of a
% second
unit_time = 0.0566;
t = 0 : unit_time : (length(y) - 1)*unit_time; %seconds


% Setting Up Inversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(y);
G = zeros(N, 3); %also called the design matrix
G(:, 1) = 1;
G(:, 2) = t;
G(:, 3) = .5*t.^2;

G_g = inv(G' * G) * G';
m_est = G_g * y;
disp(m_est)
%disp(m_est(3) /100) %this will be gravity in m/s^2

% Uncertainty Quant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
residual = y - (G * m_est); % you see this in the data misfit

sigma_2 = ( 1 / ( N - 3) ) * (residual' * residual); %sigma^2
Cov_d = sigma_2 * eye(N); %noise covariance, previously tried E(d*d^T), there is only one realization
Cov_m = G_g * Cov_d * (G_g)'; %since I assumed Cov_d has constant diagonals, I could have written Cov_m = sigma_2 inv(G' * G)

m_est2 = inv(G' * inv(Cov_d) * G) * G' * inv(Cov_d) * y;
disp(m_est2)

intervals = {};
% Confidence intervals
for i = 1:length(m_est)
    lower = m_est(i) - sqrt(sigma_2) * 1.96;
    upper = m_est(i) + sqrt(sigma_2) * 1.96;
    intervals{i} = [lower upper];
end

% X^2 test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
deg_free = N - 3; % number of data points - number of parameters
chi2_obs = (residual' * residual) / sigma_2;
p_value = 1 - chi2cdf(chi2_obs, deg_free);

% Visualizing C.L. Ellipses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1=sqrt(chi2inv(0.4,2));
s2=sqrt(chi2inv(0.68,2));
s3=sqrt(chi2inv(0.95,2));

Cov_va = Cov_m(2:3, 2:3); %covariance between velocity and accelaration
[V, D] = eig(Cov_va);
a = sqrt(D(1, 1));
b = sqrt(D(2, 2));
theta = linspace(0,2*pi,100);
ellipse1 = s1*V * [a*cos(theta); b*sin(theta)]; %these ellipses live in parameter space
ellipse2 = s2*V * [a*cos(theta); b*sin(theta)];
ellipse3 = s3*V * [a*cos(theta); b*sin(theta)];

%Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ah = gca;
scatter(t, y, 'filled')
xlabel("time (s)")
ylabel("height (cm)")
title("Projectile motion of a Golfball")

v0 = m_est(2);
a0 = m_est(3);

figure
ah2 = gca;
plot(v0 + ellipse1(1,:), a0 + ellipse1(2,:), 'b', 'LineWidth', 1.5)
hold on
plot(v0 + ellipse2(1,:), a0 + ellipse2(2,:), 'r', 'LineWidth', 1.5)
plot(v0 + ellipse3(1,:), a0 + ellipse3(2,:), 'g', 'LineWidth', 1.5)
xlabel("Velocity $(cm/s)$")
ylabel("Acceleration $(cm/s^2)$", Interpreter="latex")
legend( ...
    sprintf("CL = 0.4"), ...
    sprintf("CL = 0.68"), ...
    sprintf("CL = 0.95"));
title(['$\mathcal{X}^2$ Uncertainty Bound']);
hold off




