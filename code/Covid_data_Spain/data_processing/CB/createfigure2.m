function createfigure2(X1, Y1)
%CREATEFIGURE(X1, Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 08-Jan-2021 15:31:11

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot(X1,Y1,'LineWidth',1);

% Create ylabel
ylabel({'number of positives/number of test'});

% Create xlabel
xlabel({'date'});

% Create title
title({'Positivity ratio in Cantabria'});

box(axes1,'on');
grid(axes1,'on');
saveas(figure1, "figure2_CB.png");
saveas(figure1, "figure2_CB", "epsc");
