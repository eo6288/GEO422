function varargout=Geiger(xver)
% [ah1,ah2,t1,t2,tt]=GEIGER(xver)
%
% solves the nonlinear problem F(m) = d where m = (x, y, z, t) 
% and the data (d) is the arrival times. Also computes some
% basic uncertainty quantification
%
% INPUTS:
%   xver   1 Writes the plot to PDF
%
% OUTPUTS:
%   ah1     The axis handles to the first figure
%   ah2     The axis handles to the second figure
%   t1      The title of the first figure
%   t2      The title of the first figure
%   tt      The overarching title

defval('xver',0)


% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stationlocations = [ % station locations (x, y, z)
   10.6886   21.7918   37.8341
   -9.0488   17.5562  -80.2807
  152.0287  -99.5501   58.5911
  -26.4936   74.7840  -47.0340
   -6.0372  133.1682  -23.8966
    4.2637   67.1250   38.3537
  -56.6421  -18.8118   98.1316
 -121.0475   34.1007 -165.2877
 -196.4180 -136.5620  -59.1821
 -114.9036  189.1655   32.0711
   19.9958  -56.7188  174.0942
   40.2587   -4.6160  120.7359
  -54.7807 -130.8797   19.0181
  -12.8942  129.1239   19.4332
  265.7312   58.1433   -6.5751
   83.1159  -75.9345   76.1278
 -154.8477  242.9923   -6.0048
  168.2754  -99.1848 -100.3387
  216.2341  182.8079  -54.1840
 -112.9423   28.0698   -9.6501
    0.2106  121.1739 -267.9798
  -38.4603   85.6691 -136.2375
 -195.2332  -42.0421  -51.0887
  -53.9187   47.6592   16.3262
   31.1105   31.7220 -131.2957
  -52.9779  -92.4970 -145.9214
  -43.2357 -148.4846  157.9665 ];

noisyarrivaltimes = [ % Noisy arrival times 
   55.4982
   67.7140
   57.6544
   69.8631
   77.8886
   60.6434
   45.4901
   80.5806
   58.4471
   83.6240
   41.5813
   52.0467
   35.3534
   71.8921
   91.2715
   48.3367
   96.5975
   73.9880
  100.5424
   62.5134
  106.6134
   85.0471
   65.2416
   57.5509
   76.2200
   68.3029
   26.4729 ];

N = size(stationlocations, 1); 

% The velocity of the wave through a homogenous medium (km/s)
% Homogenous assumption gives us simple, straight-line ray paths
v = 5.6150;


% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = mean(stationlocations(:, 1)); % initial guesses m = (x, y, z, t)
y = mean(stationlocations(:, 2));
z = mean(stationlocations(:, 3));
t = mean(noisyarrivaltimes);
m = [x; y; z; t];
m0 = m;

maxIter = 20;
tol = 1e-4;

% Keep track of model misfits at each step
model_misfits = [];
model_misfits_t = [];

% Keep track of data misfits at each step
data_misfits = [];

for iter = 1:maxIter

    % Empty Jacobian G matrix
    G = zeros(N,4);

    % Empty vector of guessed arrival times
    d_pred = zeros(N,1);

    for i = 1:N
        xi = stationlocations(i,1);
        yi = stationlocations(i,2);
        zi = stationlocations(i,3);

        Di = sqrt((xi - x)^2 + (yi - y)^2 + (zi - z)^2);

        d_pred(i) = t + Di / v;

        % Jacobian
        G(i, 1) = (x - xi) / (v * Di); % Partial derivative with respect to x
        G(i, 2) = (y - yi) / (v * Di); 
        G(i, 3) = (z - zi) / (v * Di); 
        G(i, 4) = 1;
    end

    % for misfit
    data_residuals = noisyarrivaltimes - d_pred; %delta d

    % least square inversion
    dm = G \ data_residuals; %delta d = G * delta m
    m = m + dm; % Update the model parameters

    x = m(1);
    y = m(2);
    z = m(3);
    t = m(4);

    % Model misfit
    model_residuals = sqrt((m(1)-m0(1))^2 + (m(2)-m0(2))^2 + (m(3)-m0(3))^2);
    model_residuals_t = sqrt((m(1)-m0(1))^2 + (m(2)-m0(2))^2 + (m(3)-m0(3))^2) + (v*(m(4)-m0(4))^2); % space-time

    model_misfits(iter) = model_residuals;
    model_misfits_t(iter) = model_residuals_t;

    m0 = m; 

    % Data misfit
    data_misfits(iter) = norm(data_residuals);

    if data_misfits(iter) < tol
        break; % Convergence check
    end
 
end
 

% Uncertainty Quantification

% data covariance (noise level)
sigma_d = std(data_residuals); 

% model covariance matrix as a mapping from the data covariance using the
% linear operator which you have implicitly done by right-dividing
Cm = sigma_d^2 * inv(G' * G); 

% Scale for the ellipse
s1=sqrt(chi2inv(0.4,2));
s2=sqrt(chi2inv(0.68,2));
s3=sqrt(chi2inv(0.95,2));


Cxy = Cm(1:2,1:2);
% You're getting the ellipticity
[V,D] = eig(Cxy);
% semi-axis proportional length
a = sqrt(D(1,1));   
b = sqrt(D(2,2));
theta = linspace(0,2*pi,200);
% Finally apply the scaling that that maps it to the chosen confidence
% level
ellipse1 = s1*V * [a*cos(theta); b*sin(theta)];
ellipse2 = s2*V * [a*cos(theta); b*sin(theta)];
ellipse3 = s3*V * [a*cos(theta); b*sin(theta)];

% Monte Carlo
num_MC = 1000;
MC_models = zeros(4, num_MC);
for i = 1:num_MC
    % Generate random "data" from being just noise
    % and you attempt to map them to model perturbation
    dm = G \ [sigma_d * randn(N,1)];
    % Prior line is the step that maps data to models linearly
    % which is why the data covariance maps to the model covariance 
    % quadratically, which is nothing but lin 168 since inv(G'*G)*G' is the
    % linear map, and inv(G'*G)*G'*G*inv(G'G) is the quadratic and if it
    % all works well G'*G*inv(G'G) cancels and then the cloud you make
    % below should be well summarized by the ellipse AT THE RIGHT CL
    % which means you can use INPOLYGON to check that it does.
    % Produce a new model
    MC_models(:,i) = m + dm;
end

counts1 = sum(inpolygon(MC_models(1,:), MC_models(2,:), x + ellipse1(1,:), y + ellipse1(2,:)));
counts2 = sum(inpolygon(MC_models(1,:), MC_models(2,:), x + ellipse2(1,:), y + ellipse2(2,:)));
counts3 = sum(inpolygon(MC_models(1,:), MC_models(2,:), x + ellipse3(1,:), y + ellipse3(2,:)));

% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(groot, 'defaultTextInterpreter','latex') % All in latex font

station_x = stationlocations(:, 1);
station_y = stationlocations(:, 2);
station_z = stationlocations(:, 3);

fh(1)=figure(1);
clf
[ah1, ~, ~] = krijetem(subnum(2,2));

axes(ah1(1))
s1 = scatter(x, y, "green", "filled", "diamond");
hold on
s2 = scatter(-48.8540, -156.7437, "red", "o");
legend([s1, s2], 'Predition', "Truth")
scatter(station_x, station_y, "black", "filled", "^", "DisplayName", "Station Location")
xlim([-250 250])
ylim([-250 250])
xlabel("X Coordinate")
ylabel("Y Coordinate")
t1(1) = title("X-Y plane");
hold off

axes(ah1(2))
scatter(station_x, station_z, "black", "filled", "^")
hold on
scatter(x, z, "green", "filled", "diamond")
scatter(-48.8540, 124.1183, "red", "o")
ylim([-250 250])
xlim([-250 250])
xlabel("X Coordinate")
ylabel("Z Coordinate")
t1(1) = title("X-Z plane");
hold off

axes(ah1(3))
scatter(station_y, station_z, "black", "filled", "^")
hold on
scatter(y, z, "green", "filled", "diamond")
scatter(-156.7437,  124.1183, "red", "o")
ylim([-250 250])
xlim([-250 250])
xlabel("Y Coordinate")
ylabel("Z Coordinate")
t1(1) = title("Y-Z plane");
hold off

axes(ah1(4))
plot(data_misfits, "Marker","o", "MarkerFaceColor","blue")
xlabel("Iteration Number")
ylabel("Data Misfit");
%yscale('log')
ylim([0 max(data_misfits)+20])
hold off

tt = supertit(ah1(1:2),...
    sprintf(['Predicted Earthquake = [%.4f, %.4f, %.4f] with ' ...
    't = %.2f s \n True Earthquake =  [-48.8540 -156.7437 124.1183]' ...
    ' with t = 18.5044 s'], m(1), m(2), m(3), m(4)));

% Plotting the uncertainty appraisal
fh(2)=figure(2);
clf
[ah2, ~, ~] = krijetem(subnum(2,2));

axes(ah2(1))
hist(data_residuals, 15)
xlabel('Bins')
ylabel('$t_{noisy} - t_{pred}$', 'FontSize', 14)
xlim([-4, 4])
ylim([0,6])
t2(1) = title('Data Residuals');
text(1, 5.5, ['data $\sigma =$', num2str(sigma_d)])

axes(ah2(2))
plot(model_misfits, "Marker","o", "MarkerFaceColor","blue")
hold on
plot(model_misfits_t,  "Marker","^", "MarkerFaceColor","r")
yscale('log')
ylabel('Model Misfit (log-scale)')
xlabel('Iteration Number')
t2(2) = title('$L_2$ Norm of Model Misfit Variables');
legend('X, Y, Z', 'X, Y, Z, t')
hold off

axes(ah2(3))
plot(x + ellipse1(1,:), y + ellipse1(2,:), 'b', 'LineWidth', 1.5)
hold on
plot(x + ellipse2(1,:), y + ellipse2(2,:), 'r', 'LineWidth', 1.5)
plot(x + ellipse3(1,:), y + ellipse3(2,:), 'g', 'LineWidth', 1.5)
scatter(-48.8540, -156.7437, "red", "o")
xlabel("X Coordinate")
ylabel("Y Coordinate ")
xlim([-70, -30])
ylim([-190, -130])
legend( ...
    sprintf("CL = 0.4, %d points", counts1), ...
    sprintf("CL = 0.68, %d points", counts2), ...
    sprintf("CL = 0.95, %d points", counts3), ...
    "Truth")
t2(3) = title(['$\mathcal{X}^2$ Uncertainty Bound']);
hold off

axes(ah2(4))
scatter(MC_models(1,:), MC_models(2,:), 'filled', 'MarkerFaceAlpha', 0.2)
hold on
scatter(m(1), m(2), 'r', 'diamond')
xlabel("X Coordinate")
ylabel("Y Coordinate")
xlim([-70, -30])
ylim([-190, -130])
legend("Simulated", "Best X,Y model guesses")
t2(4) = title('Monte Carlo Method (1000 Simulations)');
hold off

% COSMETICS
axes(ah1(1))
box on
axes(ah1(2))
box on
axes(ah1(3))
box on
axes(ah1(4))
box on
serre(ah1([1 3]),1/3,'down') %adjust for title
serre(ah1([2 4]),1/3,'down')
movev(tt,0.2)

axes(ah2(1))
box on
axes(ah2(2))
box on
axes(ah2(3))
box on
axes(ah2(4))
box on

if xver==1
    tic
    %  MAKE THE PLOT
    figure(1)
    exportgraphics(fh(1),sprintf('%s_%i.pdf',mfilename,1))
    figure(2)
    exportgraphics(fh(2),sprintf('%s_%i.pdf',mfilename,2))
    disp(sprintf('Writing the PDF took %4.2f s',toc))
end

% OPTIONAL OUTPUTS
varns={ah1,ah2,t1,t2,tt};
varargout=varns(1:nargout);
