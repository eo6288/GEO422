clear; close all;

N = 128; %number of grid points in each direction
[kx, ky] = meshgrid((-N/2):(N/2-1), (-N/2):(N/2-1)); %makes coordinate grid with 0 in the center
kr = sqrt(kx.^2 + ky.^2); %radial wavenumber, also dimension 128x128

% Fourier coefficients (i.e. just generate a random complex number)
% randn(N, N) = generate 128 x128 random points, each one pulled from N(0,1)
Z = randn(N,N) + 1i*randn(N,N);

%choose a specific wavenumber to dominate the ring
spike_loc = 40; % preferred radial wavenumber location
width = 1;      % thickness of ring

% basically, the intuition for this is, if the radial wavenumber is close 
% in value to spike_loc, it is going to value an exponent of zero (or
% near zero), so that's going to be a 1. 
% If you slice the domain along the center line and look at one quadrant
% it's a Gaussian hill centered at the spike with the thickness (variance)
% of your choice
S_ring = exp(-(kr - spike_loc).^2/(2*width^2));

% %have monotonically decreasing background
% background = 5*exp(-kr/40);
% S_ring = background + exp(-(kr - spike_loc).^2/(2*width^2));

% we have to do the following to taper off to zero in some places and
% emphasize the spike
F_ring = Z .* sqrt(S_ring); %random coefficients weighted by the ring

%bring to space
field_ring = real(ifft2(ifftshift(F_ring)));

[ah,ha,H]=krijetem(subnum(1,2));

axes(ah(1))
imagesc(kx(1,:), ky(:,1), S_ring)
axis image
colorbar
title('Ring spectrum')
xlabel("kx")
ylabel("ky")

axes(ah(2))
imagesc(field_ring)
axis image
colorbar
title('Realization from ring spectrum')
xlabel("spatial index in x")
ylabel("spatial index in y")
