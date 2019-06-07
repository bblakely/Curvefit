function [fitresult, gof] = MosaicAlbFit(LAIMo, SNOWMo, AlbindexMo, COUNTSMo)
%CREATEFIT1(LAIMO,SNOWMO,ALBINDEXMO)
%  Create a fit.
%
%  Data for 'MosaicAlbFit' fit:
%      X Input : LAIMo
%      Y Input : SNOWMo
%      Z Output: AlbindexMo
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 02-Dec-2015 13:32:07


%% Fit: 'MosaicAlbFit'.
[xData, yData, zData, weights] = prepareSurfaceData( LAIMo, SNOWMo, AlbindexMo, COUNTSMo );

% Set up fittype and options.
ft = fittype( 'lowess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';
opts.Span = 0.15;
opts.weights = weights;

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% %Plot fit with data.
% figure( 'Name', 'MosaicAlbFit' );
% h = plot( fitresult, [xData, yData], zData );
% zlim([0,0.7]);
% legend( h, 'MosaicAlbFit', 'AlbindexMo vs. LAIMo, SNOWMo', 'Location', 'NorthEast' );
% % Label axes
% xlabel LAIMo
% ylabel SNOWMo
% zlabel AlbindexMo
% grid on
% view( 40.0, 20.0 );
% 
% %Saving r2 
global r2Alb
r2Alb=gof;

