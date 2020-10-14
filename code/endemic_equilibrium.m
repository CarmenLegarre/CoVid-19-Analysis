% Date: 13/10/2020
% Author: Carmen Legarreta
% 
% Create an script to prove that the obtained results from the system in
% the endemic equilibium case are the same as the obtained manually. The
% system equations depende of variables which values are unknown, for the
% simulation random values will be asigned to these variables.

%% Constants definition
la = 1000; rho = rand(); be = rand(); bs = rand(); ba = rand();
p2 = rand();  p1 = rand(); mu = rand(); sigma = rand(); gamma = rand();
alpha = rand(); 

%% Manual solution
e = 100;
s = 1/mu * (rho*e/(mu + rho)*((gamma/(mu+alpha*(1-p2)+gamma)*(1+(sigma+mu)*(1-p2)*p1/(mu+sigma+gamma*p2)/(1-p1))*sigma) + (sigma + mu)*(sigma + gamma)*p1*p2/(mu+sigma+gamma*p2)/(1-p1)) - (sigma + mu)*e - p1*(sigma+ mu)*e/(1-p1) + la);
q = (sigma + mu)*p1*e/(sigma + mu + gamma*p2)/(1-p1);
i = (1 + (sigma + mu)*(1 - p2)*p1/(sigma + mu + gamma*p2)/(1 - p1))*e*sigma/(mu + alpha*(1 - p2) + gamma); 
r = e/(mu + rho)*((gamma/(mu+alpha*(1-p2)+gamma)*(1+(sigma+mu)*(1-p2)*p1/(mu+sigma+gamma*p2)/(1-p1))*sigma) + (sigma + mu)*(sigma + gamma)*p1*p2/(mu+sigma+gamma*p2)/(1-p1));



%% Numerical solution
system_a = zeros (5, 5);

% first array
system_a(1, 3) = (mu + sigma + gamma*p2)*(1-p1);
system_a(1, 2) = -(sigma + mu)*p1;

% second array 
system_a(2, 1) = -mu;
system_a(2, 2) = -(sigma + mu);
system_a(2, 3) = -(mu + sigma + gamma*p2);
system_a(2, 5) =  rho;

% third array
system_a(3, 2) = sigma;
system_a(3, 3) = sigma*(1 - p2);
system_a(3, 4) = -(mu + alpha*(1 - p2) + gamma);

% fourth array
system_a(4, 3) = (sigma + gamma)*p2;
system_a(4, 4) = gamma;
system_a(4, 5) = -(mu + rho);

% fifth array
system_a(5, 2) = 1;

% system
system_b = [0, -la, 0, 0, 100]';

% solution
solution = linsolve(system_a, system_b);



