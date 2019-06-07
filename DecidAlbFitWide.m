function [fitresult, gof] = DecidAlbFitWide(LAID, SNOWD, AlbindexD,COUNTSD)
%CREATEFIT1(LAID,SNOWD,ALBINDEXD,COUNTSD)
%  Create a fit.
%
%  Data for 'DecidAlbFitWide' fit:
%      X Input : LAID
%      Y Input : SNOWD
%      Z Output: AlbindexD
%      Weights : COUNTSD
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 16-Jun-2016 13:20:13


%% Fit: 'DecidAlbFitWide'.
[xData, yData, zData, weights] = prepareSurfaceData( LAID, SNOWD, AlbindexD, COUNTSD);

% Set up fittype and options.
ft = fittype( 'lowess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';
opts.Span = 0.25;
opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'DecidAlbFitWide' );
% h = plot( fitresult, [xData, yData], zData );
% zlim([0,0.6]);
% legend( h, 'DecidAlbFitWide', 'AlbindexD vs. LAID, SNOWD with COUNTSD', 'Location', 'NorthEast' );
% % Label axes
% xlabel LAID
% ylabel SNOWD
% zlabel AlbindexD
% grid on
% view( 33.5, 32.0 );
% 

%Saving r2 (terrible practice, I know)
global r2Alb
r2Alb=gof;

