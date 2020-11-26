% Date: 13/10/2020
% Author: Carmen Legarreta
% 
% Create an script to prove that the obtained results from the system in
% the endemic equilibium case are the same as the obtained manually. The
% system equations depende of variables which values are unknown, for the
% simulation random values will be asigned to these variables.

%% Constants definition
la = 1000; rho = rand(); be = rand(); bs = rand(); ba = rand();
p2 = rand();  p1 = rand(); p3 = rand();  mu = rand(); sigma = rand(); 
gamma = rand(); alpha = rand(); delta = rand();

%% Manual solution
e = 100;

s = 1.0 / mu * (rho* (gamma - p3 * delta) / (rho + mu) / ...
    (mu + (1 - p2) * alpha + gamma - p3*delta) * ((sigma + p1 * p3 * ...
    delta * (sigma + mu) / (1 - p1) / (mu + (1 - p2) * alpha + p3 * ...
    delta)) + p1 * (sigma + mu) * delta * (1 - p3) / (gamma - p3 * delta))- ...
    (sigma + mu)/(1 - p1)) * e + la/mu ;

q = (sigma + mu)*p1*e/(mu + (1-p2) * alpha + p3 * delta)/(1-p1);

i = (sigma + p1 * p3 * delta * (sigma + mu) / (1 - p1) / ...
    (mu + (1 - p2) * alpha + p3 *delta)) * e / (mu + (1 - p2) * alpha + ...
    gamma - p3 * delta); 

r = ((gamma - p3 * delta)*(sigma + p1 * p3 * delta * (sigma + mu) / (1 - p1) / ...
    (mu + (1 - p2) * alpha + p3 *delta)) + p1 * (sigma + mu) * delta * ...
    (1 - p3)) * e / (rho + mu) / (mu + (1 - p2) * alpha + gamma - p3 * delta);


%% Numerical solution
system_a = zeros (5, 5);

% first array
system_a(1, 3) = (mu + (1 - p2) * alpha + p3 * delta)*(1-p1);
system_a(1, 2) = -(sigma + mu)*p1;

% second array
system_a(2, 1) = -mu;
system_a(2, 2) = -(sigma + mu);
system_a(2, 3) = -(mu + (1 - p2) * alpha + p3 * delta);
system_a(2, 5) =  rho;

% third array
system_a(3, 2) = sigma;
system_a(3, 3) = delta * p3;
system_a(3, 4) = -(mu + alpha*(1 - p2) + gamma - p3 * delta);

% fourth array
system_a(4, 3) = delta * (1 - p3);
system_a(4, 4) = (gamma - p3 * delta);
system_a(4, 5) = -(mu + rho);

% fifth array
system_a(5, 2) = 1;

% system
system_b = [0, -la, 0, 0, e]';

% solution
solution = linsolve(system_a, system_b);


%% Verification
eps = 1.0e-3;
solution2 = [s, e, q, i, r];
v = 1;

for i = 1:length(solution)
    v = v * solution(i) - eps < solution2(i) < solution(i) +eps;
end



