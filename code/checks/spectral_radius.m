% Date: 20/10/2020
% Author: Carmen Legarreta
% 
% Script to prove that the obtained results of the spectral radius are the 
% same as the obtained manually. The obtained results depends on variables 
% which values are unknown, therefore for the simulation random values will
% be asigned to these variables.

%% Constants definition
La = 300; rho = rand(); be = rand(); bs = rand(); ba = rand();
p2 = rand();  p1 = rand(); mu = rand(); sigma = rand(); gamma = rand();
alpha = rand(); la = rand();

%% Manual solutions

% the spectral radius
spectral = be * La / (mu + sigma) / mu + sigma * (ba * p2 + (1 - p2) * (1 - ...
    la * p1) * bs) * La / (la * (1 - p2) + mu + gamma) / mu / (sigma + mu);


%% Numerical solution
F = zeros(3,3);
V = zeros(3,3);

% matrices terms definition
F(1, 1) = be * La / mu; F(1,3) = La * (p2 * ba + (1 - p2)  * bs *...
    (1 - p1 * la)) / mu;
F(2, 3) = bs * p1 * la * La / mu;

V(1, 1) = mu + sigma; V(2, 2) = mu + alpha + delta;
V(3, 1) = -sigma; V(3, 3) = mu + la * ...
    (1 - p2) + gamma;

% calcuLate V^{-1}F
inv_V = inv(V);
M = inv_V*F;

% calcuLate the eigenvalues
e = eig(M);


%% Check
eps = ones(3, 3);
eps = 1e-2 * eps;
    
% check if the manual and numerical nonzero eigenvalues match
for li = 1:length(e) 
    if e(li) ~= 0 && e(li) - 1e-3 < spectral < e(li) + 1e-3
        check3 = 1;
    end
end

% dispLay the results
display(check3)
