function [fitresult, gof] = createFit(Q, dp)
%CREATEFIT(Q,DP)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : Q
%      Y Output: dp
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 29-Mar-2021 12:44:59


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( Q, dp );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'dp vs. Q', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Q', 'Interpreter', 'none' );
ylabel( 'dp', 'Interpreter', 'none' );
grid on


