function [fitresult, gof] = EvgreenAlbFitWide(LAIE, SNOWE, AlbindexE, COUNTSE)
%CREATEFIT1(LAIE,SNOWE,ALBINDEXE,COUNTSE)
%  Create a fit.
%
%  Data for 'EvgreenAlbFitWide' fit:
%      X Input : LAIE
%      Y Input : SNOWE
%      Z Output: AlbindexE
%      Weights : COUNTSE
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 01-Dec-2015 18:20:52


%% Fit: 'EvgreenAlbFitWide'.
[xData, yData, zData, weights] = prepareSurfaceData( LAIE, SNOWE, AlbindexE, COUNTSE );

% Set up fittype and options.
ft = fittype( 'lowess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';
opts.Span = 0.5;
opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'EvgreenAlbFitWide' );
% h = plot( fitresult, [xData, yData], zData );
% zlim([0,0.6]);
% legend( h, 'EvgreenAlbFitWide', 'AlbindexE vs. LAIE, SNOWE with COUNTSE', 'Location', 'NorthEast' );
% % Label axes
% xlabel LAIE
% ylabel SNOWE
% zlabel AlbindexE
% grid on
% view( 29.5, 30.0 );
%Saving R2
global r2Alb
r2Alb=gof;
