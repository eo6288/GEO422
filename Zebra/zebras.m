function zebras(wt)
% zebras(wt)
%
% INPUT
%
% wt      Threshold color
%
% Loads the image of a zebra and plots histogram of colors
%
% Last modified by erinoneil on 01/13/2026

% PRELIMINARIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
if length(wt) == 1
    wt = [wt wt wt]; % Ensure wt is a 3-element vector for RGB thresholds
end

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = imread(fullfile('/Users','eo6288','Documents','MATLAB','ZebraPractice','zebra.jpeg'));
cX = whos('X'); %this checks the class and saves it

% CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
red = X(:, :, 1);
green = X(:, :, 2);
blue = X(:, :, 3);
colx = {'red', 'green', 'blue'};

% Use the input parameter wt
condr = red<=wt(1);
condg = green<=wt(2);
condb = blue<=wt(3);

red(condr)=0;
green(condg)=0;
blue(condb)=255;

X(:,:,1) = red;
X(:,:,2) = green;
X(:,:,3) = blue;

% Define histogram parameters
histbins = linspace(0, double(intmax(cX.class)), 15);
histxlim = [histbins(1) histbins(end)];

% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(1) = subplot(211);
image(X); axis image;
xl(1) = xlabel('column (pixels)');
yl(1) = ylabel('rows (pixels)');
t(1) = title(sprintf('This image is of size %ix%i', size(X,1), size(X,2)));
ah(1).XTick = unique([1 [1:floor(size(X,2)/100)]*100 size(X,2)]);
ah(1).YTick = unique([1 [1:floor(size(X,1)/100)]*100 size(X,1)]);

% Plot histograms
ah(2) = subplot(234);
h(1) = histogram(red, histbins);
ah(3) = subplot(235);
h(2) = histogram(green, histbins);
ah(4) = subplot(236);
h(3) = histogram(blue, histbins);

%Cosmetic
set(ah(2:4), 'xlim', histxlim)
set(ah(2:4), 'ylim', [0 1.1*max(cat(2, h(:).Values))])
set(ah(2:4), 'xtick', ...
    unique([histxlim histxlim(1) + [1:3]*round(histxlim(2)/4)]))
for index=1:length(colx)
    h(index).FaceColor = colx{index};
    axes(ah(index+1))
    xl(index+1)=xlabel(colx{index});
end

print(gcf, sprintf('%s_%i_%i_%i', mfilename,wt(1), wt(2), wt(3)), '-dpdf', '-bestfit');










