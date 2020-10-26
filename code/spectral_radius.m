% Date: 20/10/2020
% Author: Carmen Legarreta
% 
% Script to prove that the obtained results of the spectral radius are the 
% same as the obtained manually. The obtained results depends on variables 
% which values are unknown, therefore for the simulation random values will
% be asigned to these variables.

%% Constants definition
la = 1000; rho = rand(); be = rand(); bs = rand(); ba = rand();
p2 = rand();  p1 = rand(); mu = rand(); sigma = rand(); gamma = rand();
alpha = rand(); 

%% Manual solutions

% the inverse of V
manual_inv_V = zeros(3, 3);
manual_inv_V(1, 1) = 1/(mu + sigma); 
manual_inv_V(2, 2) = 1/(mu + sigma + p2*gamma);
manual_inv_V(3, 1) = sigma / (sigma + mu) / (mu + alpha * (1 - p2) + gamma);
manual_inv_V(3, 2) = sigma*(1 - p2)/ (mu + sigma + gamma * p2) / ...
    (mu + alpha * (1 - p2) + gamma);
manual_inv_V(3, 3) = 1 / (mu + alpha * (1 - p2) + gamma);

% the product of V^{-1} and F
manual_M = zeros(3, 3);
manual_M(1, 1) = be * la * (1 - p1) / (mu + sigma) / mu;
manual_M(1, 3) = la * (p2 * ba + (1 - p2) * bs) * (1 - p1) / mu / (mu + sigma);
manual_M(2, 1) = be * la * p1 / (mu + sigma + p2 * gamma) / mu;
manual_M(2, 3) = la * (ba * p2 + (1 - p2) * bs) * p1 / mu / ...
    (mu + sigma + p2 * gamma);
manual_M(3, 1) = be * sigma * la * ((1 - p1) / (sigma + mu) + ...
    p1 * (1 - p2) / (mu + sigma + p2 * gamma)) / ...
    (mu + alpha * (1- p2) + gamma) / mu;
manual_M(3, 3) = la * (p2 * ba + (1 - p2) * bs) * sigma * ((1 - p1) / ...
    (sigma + mu) + p1 * (1 - p2) / (mu + sigma + p2 * gamma)) / ... 
    (mu + alpha * (1 - p2) + gamma) / mu;

% the spectral radius
spectral = be * la * (1 - p1) / mu / (mu + sigma) +...
    la * sigma * (p2 * ba + (1 - p2) * bs) *...
    ((1 - p1) / (mu + sigma) + p1 * (1 - p2) / ...
    ( mu + sigma + p2 * gamma))/ mu / (mu+alpha*(1-p2)+gamma);


%% Numerical solution
F = zeros(3,3);
V = zeros(3,3);

% matrices terms definition
F(1, 1) = be * la *(1-p1)/mu; F(1,3) = la * (p2 * ba + (1 - p2) * bs) *...
    (1 - p1) / mu;
F(2, 1) = be * p1 * la / mu; F(2, 3) = la * (p2 * ba + (1 - p2) * bs) *...
    p1 / mu;

V(1, 1) = mu + sigma; V(2, 2) = mu + sigma + gamma * p2;
V(3, 1) = -sigma; V(3, 2) = -sigma * (1 - p2); V(3, 3) = mu + alpha * ...
    (1 - p2) + gamma;

% calculate V^{-1}F
inv_V = inv(V);
M = inv_V*F;

% calculate the eigenvalues
e = eig(M);


%% Check
eps = ones(3, 3);
eps = 1e-2 * eps;

% check if the manual and numerical V^{-1} is the shame
if inv_V - eps < manual_inv_V
    if manual_inv_V < inv_V + eps
        check1 = 1;
    end
end
    
% check if the manual and numerical V^{-1}F is the shame
if M - eps < manual_M 
    if manual_M < M + eps
        check2 = 1;
    end
end
    
% check if the manual and numerical nonzero eigenvalues match
for li = 1:length(e) 
    if e(li) ~= 0 && e(li) - 1e-3 < spectral < e(li) + 1e-3
        check3 = 1;
    end
end

% display the results
display(check1 * check2 * check3)
