function createfigure(X1, Y1)
%CREATEFIGURE(X1, Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 07-Jan-2021 12:44:56

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot(X1,Y1,'LineWidth',1);

% Create ylabel
ylabel({'Number of people in UCI'});

% Create xlabel
xlabel({'date'});

% Create title
title({'People in UCI due to Covid-19 in Basque country'});

box(axes1,'on');
grid(axes1,'on');
