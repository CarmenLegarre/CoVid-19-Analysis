% Date: 25/10/2020
% Author: Carmen Legarreta
% 
% Script which has been created to solve the system of differential
% equations when R>1 and R<1. 

clear all; close all;

%% Constant definition 
global alpha gamma la mu p1 p2 rho sigma ba be bs

R = 'R_unstable.mat';
c = load(R);

alpha = c.alpha; gamma = c.gamma; la = c.la; mu = c.mu; p1 = c.p1;
p2 = c.p2; rho = c.rho; sigma = c.sigma; ba = c.ba; be = c.be; 
bs = c.bs;

spectral = be * la * (1 - p1) / mu / (mu + sigma) +...
    la * sigma * (p2 * ba + (1 - p2) * bs) *...
    ((1 - p1) / (mu + sigma) + p1 * (1 - p2) / ...
    ( mu + sigma + p2 * gamma))/ mu / (mu+alpha*(1-p2)+gamma);


%% Solve the system
ts = 0.01;
[t, y] = ode45(@vdp1, [0:ts:200], [1e7; 0; 0; 100; 0]);


%% Plot the results
if strcmp(R, 'R_unstable.mat')
    createfigure_unstable(t, [y(:, 1) y(:, 5)], y(:, 2:4))
else 
    createfigure_stable(t, y(:, 1), y(:, 2:5))
end


%% System definition in a function

function dxdt = vdp1(t, x)

    global alpha gamma la mu p1 p2 rho sigma ba be bs
    
    S = x(1); Q = x(2); E = x(3); I = x(4); R = x(5);
    
    dxdt = zeros(5, 1);
    
    dxdt(1) = la + R * rho - be * S * E - S * I * ...
        (p2 * ba + (1 - p2)*bs) - mu * S;
    dxdt(2) = (be * S * E + S * I * (p2 * ba + (1 - p2)*bs)) * ...
        p1 - Q * (mu + sigma + gamma * p2);
    dxdt(3) = (be * S * E + S * I * (p2 * ba + (1 - p2)*bs)) * ...
        (1-p1) - E * (sigma + mu);
    dxdt(4) = sigma * (E + Q * (1 - p2)) - I * (mu + alpha * ...
        (1 - p2) + gamma);
    dxdt(5) = Q * (sigma + gamma) * p2 + gamma * I - ...
        (rho + mu) * R;
    
end