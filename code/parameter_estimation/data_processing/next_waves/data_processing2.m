% Update: 03/05/2020
% Author: Carmen Legarreta
% Institution: University of Basque Country (UPV/EHU)
% 
% On this script epidemiological data from the Autonomous Community of
% Cantabria is gathered and is used to estimate the people in quarantine Q 
% on the first wave. It is also analysed possible delays between estimated
% Q and the reported deaths d(k) and recovered people r(k). Moreover, the
% result from the linear regression of the data sets {Q, r(k)} and {Q, dp(k)}
% is displayed.


%% Gather data

% include paths to call the function generate_data() 
addpath("/home/carmen/Repositories/CoVid/CoVid-19-Analysis/code/Covid_data_Spain/data_processing/CB");
addpath("/home/carmen/Repositories/CoVid/CoVid-19-Analysis/code/parameter_estimation");

% obtain data from cantabria
[date_cb, test_cb, tp_cb, ph_cb, pi_cb, dp_cb, r_cb] = generate_data();

% apply movemean technique over data. The window width is 14
w = 14;
tp_cb = movmean(tp_cb, w); ph_cb = movmean(ph_cb, w);
pi_cb = movmean(pi_cb, w); dp_cb = movmean(dp_cb, w);
r_cb = movmean(r_cb, w); test_cb = movmean(test_cb, w);

% compact data
data1 = [tp_cb, dp_cb, r_cb, test_cb, ph_cb];


%% Following waves
% time range
ts = 140; tf = 352;

% call function optimal_Q() to find the optimal delays
[data, tau1, tau2, tau3] = optimal_Q(data1, ts, tf);

% extract data
Q = data(:, 1); tp = data(:, 2); test = data(:, 5);
r = data(:, 3); dp = data(:, 4);

% intervals where it is assumed that p1 and k are constant
t_cb = [1 21;25 48;50 63;63 90;91 106;106 116;116 130;132 150;156 174;174 188;188 196];


%% Create figures and fit
% Q vs. r
[fitresult_r, gof_r] = createFit2r(Q, r);

coeff = coeffvalues(fitresult_r);
delta = coeff(1);

error = confint(fitresult_r);
error_delta = abs(diff(error(:, 1)));

display(coeff); display(error); display(gof_r.rsquare)

% Q vs. dp
[fitresult_dp, gof_dp] = createFit2dp(Q, dp);
coeff = coeffvalues(fitresult_dp);
alpha = coeff(1);

error = confint(fitresult_dp);
error_alpha = abs(diff(error(:, 1)));

display(coeff); display(error); display(gof_dp.rsquare)


%% deviation standard of Q
% obtain data from Cantabria
[date_cb, test_cb, tp_cb, ph_cb, pi_cb, dp_cb, r_cb] = generate_data();

% calculate people in quarantine with the previously calculated delays and
% without prefiltering data
[Q1, tp2, r2, dp2, test2] = quarantine([tp_cb dp_cb r_cb test_cb],ts, tf, tau1, tau3);

% first derivative approximation
Q2 = Q1(2:end)-Q1(1:end-1);
% standard deviation of the first derivative mean value
st2 = movstd(Q2, w)./sqrt(w);
% filter first derivative
Q2 = movmean(Q2, w);

% second derivative approximation
Q3 = Q1(3:end)-2*Q1(2:end-1)+Q1(1:end-2);
% standard deviation of the second derivative mean value
st3 = movstd(Q3, w)./sqrt(w);
% filter second derivative
Q3 = movmean(Q3, w); 

% standard deviation of the mean values 
est1 = movstd(tp2, w)./sqrt(w); est2 = movstd(dp2, w)./sqrt(w);
est3 = movstd(r2, w)./sqrt(w);
% standard deviation of Q
st = sqrt(est1.^2+est2.^2+est3.^2);


%% Functions
function [data2,tau1, tau2, tau3] = optimal_Q(data1, ts, tf)
% this function calculates the people in quarantine
%
% Inputs
%   -data1: array which holds the following information
%       data(:, 1): number of postive tests per day
%       data(:, 2): number of death people per day
%       data(:, 3): number of recovered people per day
%       data(:, 4): number of performed test per day
%   -ts: when to start gathering data
%   -tf: when to stop gathering data
%
% Outputs
%   -data2: array which holds the following information
%       data(:, 1): number of people in quarantine
%       data(:, 2): number of postive tests per day
%       data(:, 3): number of death people per day
%       data(:, 4): number of recovered people per day
%       data(:, 5): number of performed test per day
%   -tau1: delay of the number of recovered people with respect the number
%   of deaths
%   -tau2: delay between the number of people in quarantine and the number
%   of recovered people (or death people)
%   -tau3: delay of the number of positive test with respect the number of
%   deaths

% Variable to maximize
Z_max = 0;

% data break down
tp_cb = data1(:, 1); dp_cb = data1(:, 2);
r_cb = data1(:, 3); test_cb = data1(:, 4);

% delay between dp and r
dd = finddelay(dp_cb(ts:tf), r_cb(ts:tf));
for d = dd
tau1 = d;

% try different delays between dp and tp
l =15;

% only positive values of tau3 are considered because for tau3>0 Z is 
% increased only a little, which is insignificant since at the end
% some points are deleted.
Tau3 = [-l:-1,0];
for j = 1:length(Tau3)
    
    tau3 = Tau3(j);
    [Q1, tp2, r2, dp2, test2] = quarantine([tp_cb dp_cb r_cb test_cb],ts, tf, tau1, tau3);
    [val1 index1] = max(r2);
    [val3 index2] = max(Q1);
    tau2 = index1-index2;
    tau2 = round(tau2);
    
    % apply delay 
    if  tau2 > 0
        Q1 = Q1(1:end-tau2); dp2 = dp2(tau2+1:end); r2 = r2(tau2+1:end);
        tp2 = tp2(tau2+1:end); test2 = test2(tau2+1:end);
    elseif  tau2 < 0
        Q1 = Q1(abs(tau2)+1:end); dp2 = dp2(1:end-abs(tau2)); r2 = r2(1:end-abs(tau2));
        tp2 = tp2(1:end-abs(tau2)); test2 = test2(1:end-abs(tau2));  
    end
    

    % apply linear regression 
    [fitobject,gof1]=fit(Q1,dp2,'poly1');
    [fitobject,gof12]=fit(Q1,r2,'poly1');
    
    % linear regression residual (cost function)
    Z = (1*gof12.rsquare+1*gof1.rsquare)/2;
    
    % condition to optimize the residual
    if Z > 1*Z_max 
        Z_max = Z;
        t1 = tau1;
        t2 = tau2;
        t3 = tau3;
        data2 = [Q1, tp2, r2, dp2, test2];
    end
end
end
tau1 = t1;
tau2 = t2;
tau3 = t3;

end

function [Q1, tp, r, dp, test] = quarantine(data, ts, tf, tau1, tau2)
% this function calculates the people in quarantine
%
% Inputs
%   -data: array wich holds the following information
%       data(:, 1): number of postive tests per day
%       data(:, 2): number of death people per day
%       data(:, 3): number of recovered people per day
%       data(:, 4): number of performed test per day
%   -tau1: delay of r with respect to dp
%   -tau2: delay of tp and test with respect to dp and r
%
% Outputs
%   -Q1: number of people in quarantine
%   -tp, r, dp, test: number of positive test, recovered people, death 
%       people, and performed test per day, respectively,
%       with the delays especified at the input.


tp = data(ts:tf, 1); dp = data(ts:tf, 2); r = data(ts:tf, 3); 
test = data(ts:tf, 4);


if tau1>0
    tp = tp(1:end-tau1); dp = dp(1:end-tau1); r = r(tau1+1:end);
    test = test(1:end-tau1);
elseif tau1<0
    tau1 = -tau1;
    tp = tp(tau1+1:end); dp = dp(tau1+1:end); r = r(1:end-tau1);
    test = test(tau1+1:end);
end


if tau2>0
    tp = tp(tau2+1:end); dp = dp(1:end-tau2); r = r(1:end-tau2);
    test = test(tau2+1:end);
elseif tau2<0
    tau2 = -tau2;
    tp = tp(1:end-tau2); dp = dp(tau2+1:end); r = r(tau2+1:end);
    test = test(tau2+1:end);
end


Q = tp-dp-r;

Q1 = zeros(length(Q), 1);
for jj = 2:length(Q)
    Q1(jj) = sum(Q(jj:-1:1));
end

if min(Q1) < 0
    Q1 = Q1-min(Q1)*ones(length(Q1), 1);
end

end