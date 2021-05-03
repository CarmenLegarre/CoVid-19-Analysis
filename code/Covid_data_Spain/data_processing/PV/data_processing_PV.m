clear all;
table1 = readtable('../../data/PV/situacion-epidemiologica.xlsx');
dimension = length(table1.KORONABIRUSARENEGOERAEPIDEMILOGIKOAEUSKADIN_COVID_19__SITUACI_N(:, 1));
table2 = readtable("../../data/PV/01.csv");
table3 = readtable("../../data/PV/04.csv");

%% dates
dates = table1.KORONABIRUSARENEGOERAEPIDEMILOGIKOAEUSKADIN_COVID_19__SITUACI_N(2:dimension);
dates = cell2mat(dates);
dt = datetime(dates, 'InputFormat', 'yyyy/MM/dd');

%% total test
% cumulative
tt = table1.Var2(2:dimension);
tt = str2double(tt);

% daily
tt_new = zeros(dimension-1, 1);
tt_new(1, 1)= tt(1, 1);
for i = 2:dimension-1
    tt_new(i, 1) = tt(i, 1)-tt(i-1, 1);
end

%% Positives

% Araba
tp_A = table1.Var10(2:dimension);
tp_A = str2double(tp_A);
% Bizkaia
tp_B = table1.Var11(2:dimension);
tp_B = str2double(tp_B);
% Gipuzkoa
tp_G = table1.Var12(2:dimension);
tp_G = str2double(tp_G);
% total
tp_t = table1.Var9(2:dimension);
tp_t = str2double(tp_t);
% cumulative
tp = zeros(dimension-1, 1);
for i = 1:dimension-1
    tp(i, 1) = tp_A(i, 1) + tp_B(i, 1) + tp_G(i, 1);
end

% daily
tp_new = nan(dimension-1, 1);
tp_new(1, 1) = tp(1, 1);
for i = 2:dimension-1
    tp_new(i, 1) = tp(i, 1)-tp(i-1, 1);
end



%% people in hospital (daily)
% ph gathers data coming from different files.
ph = table1.Var18(2:dimension);
ph = str2double(ph);
ph1 = table2.Var2(2:end, 1);
ph1 = str2double(ph1);
dates1 = table2.Evoluci_nN_DeIngresadosEnPlanta(2:end, 1);
dates1 = cell2mat(dates1);
dt1 = datetime(dates1, 'InputFormat', 'yyyy/MM/dd');
dt2 = dt;
ph2 = NaN(length(dt2), 1);
j = 0;
for i = 1:length(ph2)
    if dt1(1) > dt(i)
        j = j+1;
    end
    if dt2(i) <= dt1(end) && dt2(i)>=dt1(1)
        ph2(i) = ph1(i-j);
    else
        ph2(i) = ph(i);
    end
end

%% people in UCI (daily)
pi = table1.Var19(2:dimension);
pi = str2double(pi);
pi1 = table3.Var2(3:end, 1);
pi1 = str2double(pi1);
dates1 = table3.Evoluci_nN_DeIngresadosEnUCI(3:end, 1);
dates1 = cell2mat(dates1);
dt1 = datetime(dates1, 'InputFormat', 'yyyy/MM/dd');
pi2 = NaN(length(dt2), 1);
j = 0;
for i = 1:length(ph2)
    if dt1(1) > dt(i)
        j = j+1;
    end
    if dt2(i) <= dt1(end) && dt2(i)>=dt1(1)
        pi2(i) = pi1(i-j);
    else
        pi2(i) = pi(i);
    end
end


%% deaths (acumulative)
dp = table1.Var16(2:dimension);
dp = str2double(dp);
dp_daily = zeros(dimension-1, 1);
dp_daily(1)= dp(1);
dp_daily_s = nan(dimension-1, 1);
dp_daily_s (1) = dp(1);

j= 1;
for i = 2:dimension-1
    if isnan(dp(i))
        dp_daily(i) = 0;
    else 
        dp_daily(i) = (dp(i) - dp(j))/(i-j);
        dp_daily_s(i) = (dp(i)-dp(j))/(i-j);
        % linear interpolation
        if j+1 < i
            m = (dp_daily(i) - dp_daily(j)) / (i-j);
            n = dp_daily(i) - m*i;
            for ii = i-1:-1:j
                dp_daily(ii) = (dp(i) - dp(j))*ii/(i-j) + n;
                dp_daily(ii) = m*ii+n;
            end
        end
        j = i;
    end
end



%% Figure
createfigure1(dt, tp_new);
createfigure2(dt, tp_new./tt_new);
createfigure3(dt, pi);
createfigure4(dt2, ph2);
createfigure5(dt, dp);
createfigure6(dt, [ph2(:), pi2(:), dp_daily(:)]);
close all