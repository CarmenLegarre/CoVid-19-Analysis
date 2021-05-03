% Date: 13/10/2020
% Author: Carmen Legarreta
% 
% Create an script to prove that the obtained results from the system in
% the endemic equilibium case are the same as the obtained manually. The
% system equations depende of variables which values are unknown, for the
% simulation random values will be asigned to these variables.

%% Constants definition
La = 1000; rho = rand(); be = rand(); bs = rand(); ba = rand();
p2 = rand();  p1 = rand(); p3 = rand();  mu = rand(); sigma = rand(); 
gamma = rand(); alpha = rand(); delta = rand(); la = rand();

%% Manual solution
b1 = (la*(1-p2) + mu + gamma);
b2 = (la*(1-p2)/b1 + (sigma + mu) * (1-p2) * p1 * la * bs)/(b1 * be + ...
    (ba * p2 + (bs * (1-p2) * (1 - p1 * la))) * sigma) / (mu + delta + alpha);

e = 100;

s = La/mu + sigma * (b2 * ((rho * delta)/(rho + mu) - mu - alpha ...
    - delta) + ((rho * gamma)/(rho + mu) + la * (1 - p2))/b1 - ...
    (sigma + mu)/sigma) * e / mu ;

q = sigma * b2 * e;

i = sigma * e / b1; 

r = sigma * (gamma / b1 + delta * b2) * e/ (rho + mu);


%% Numerical solution
% x = S E Q I R

system_a = zeros (5, 5);

% first array
system_a(1, 2) = -sigma * b1 * b2 * (mu + alpha + delta);
system_a(1, 3) =  b1 * (mu + alpha + delta);

% second array
system_a(2, 1) = -mu;
system_a(2, 2) = -(sigma + mu);
system_a(2, 3) = -(mu + alpha + delta);
system_a(2, 4) = la * (1 - p2);
system_a(2, 5) =  rho;

% third array
system_a(3, 2) = sigma;
system_a(3, 4) = -(la * (1 - p2) + mu + gamma);

% fourth array
system_a(4, 3) = delta;
system_a(4, 4) = gamma;
system_a(4, 5) = -(mu + rho);

% fifth array
system_a(5, 2) = 1;

% system
system_b = [0, -La, 0, 0, e]';

% solution
solution = linsolve(system_a, system_b);


%% Verification
eps = 1.0e-3;
solution2 = [s, e, q, i, r];
v = 1;

for i = 1:length(solution)
    v = v * solution(i) - eps < solution2(i) < solution(i) +eps;
end



