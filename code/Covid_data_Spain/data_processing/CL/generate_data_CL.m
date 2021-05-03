function [date, tests, tp, ph, pi, dp, r] = generate_data_CL()
table1 = readtable("../../data/CL/pruebas-realizados-coronavirus.csv");
table2 = readtable("../../data/CL/situacion-de-hospitalizados-por-coronavirus-en-castilla-y-leon.csv");

date1 = table1.Fecha(:);
date2 = table2.fecha(:);

index2 = [1];

ii =1;
for i = 2:length(date2)
    if date2(index2(ii)) ~= date2(i)
        index2 = [index2, i];
        ii = ii + 1;
    end
end

reduced_dates2 = NaT(length(index2)-1, 1);

ph = zeros(length(index2)-1, 1);
pi = zeros(length(index2)-1, 1);
dp1 = zeros(length(index2)-1, 1);
dp = zeros(length(index2)-1, 1);
n_ph = zeros(length(index2)-1, 1);
n_uci = zeros(length(index2)-1, 1);

total_h = table2.hospitalizados_planta(:);
total_uci = table2.hospitalizados_uci(:);
total_deaths = table2.fallecimientos(:);
total_r = table2.altas(:);
total_n_ph = table2.nuevos_hospitalizados_planta(:);
total_n_uci = table2.nuevos_hospitalizados_uci(:);
for i = 1:length(index2)-1
    reduced_dates2(i) = date2(index2(i));
    ph(i) = sum(total_h(index2(i):index2(i+1)-1));
    pi(i) = sum(total_uci(index2(i):index2(i+1)-1));
    dp1(i) = sum(total_deaths(index2(i):index2(i+1)-1));
    r(i) = sum(total_r(index2(i):index2(i+1)-1));
    n_ph(i) = sum(total_n_ph(index2(i):index2(i+1)-1));
    n_uci(i) = sum(total_n_uci(index2(i):index2(i+1)-1));
end

dp1 = dp1(length(index2)-1:-1:1);
dp(1) = dp1(1);

r = r(end:-1:1);
r(2:end) = r(2:end)'-r(1:end-1)';

for i = 2:length(index2)-1
    dp(i) = dp1(i)-dp1(i-1);
end

date1 = date1(end:-1:1);
reduced_dates2 = reduced_dates2(end:-1:1);
ph = ph(end:-1:1);
pi = pi(end:-1:1);
n_ph = n_ph(end:-1:1);
n_uci = n_uci(end:-1:1);


test_rapido = table1.test_rapidos_total(end:-1:1);
test_antigeno = table1.test_antigenos_total(end:-1:1);
test_pcr = table1.pcr_total(end:-1:1);
tests = zeros(length(date1), 1);

positive_rapido = table1.test_rapidos_positivos(end:-1:1);
positive_antigeno = table1.test_antigenos_positivos(end:-1:1);
positive_pcr = table1.pcr_positivos(end:-1:1);
tp = zeros(length(date1), 1);

for i = 1:length(date1)
    if isnan(test_rapido(i))
        test_rapido(i) = 0;
        positive_rapido(i) = 0;
    end
    if isnan(test_antigeno(i))
        test_antigeno(i) = 0;
        positive_antigeno(i) = 0;
    end
    if isnan(test_pcr(i))
        test_pcr(i) = 0;
        positive_pcr(i) = 0;
    end
    if i == 1
        tests(i) = test_rapido(i) + test_antigeno(i) + test_pcr(i); 
        tp(i) = positive_antigeno(i) + positive_rapido(i) + positive_pcr(i);
    else
        tests(i) = test_rapido(i) + test_antigeno(i) + test_pcr(i) - ...
            test_rapido(i-1) - test_antigeno(i-1) - test_pcr(i-1);
        tp(i) = positive_antigeno(i) + positive_rapido(i) + ... 
            positive_pcr(i) - positive_antigeno(i -1) - ...
            positive_rapido(i-1) - positive_pcr(i-1);
    end
end

date = date1(1):date1(end);
test_n = zeros(length(date), 1);
tp_n = zeros(length(date), 1);

for i = 1:length(date)
    n = true;
    for j = 1:length(date1) 
        if date(i) == date1(j)
            n = false;
            test_n(i) = tests(j);
            tp_n(i) = tp(j);
            k = j;
        end
    end
    if n
        test_n(i) = interp1([date1(k) date1(k+1)], [tests(k) tests(k+1)], date(i));
        tp_n(i) = interp1([date1(k) date1(k+1)], [tp(k) tp(k+1)], date(i));
    end
end

date1 = date;

date = reduced_dates2(1):reduced_dates2(end);
ph_n = zeros(length(date), 1);
pi_n = zeros(length(date), 1);
dp_n = zeros(length(date), 1);
r_n = zeros(length(date), 1);

for i = 1:length(date)
    n = true;
    for j = 1:length(reduced_dates2)
        if date(i) == reduced_dates2(j)
            n = false;
            ph_n(i) = ph(j);
            pi_n(i) = pi(j);
            dp_n(i) = dp(j);
            r_n(i) = r(j);
            k = j;
            
        end
    end
    if n
        ph_n(i) = interp1([reduced_dates2(k) reduced_dates2(k+1)], [ph(k) ph(k+1)], date(i));
        tp_n(i) = interp1([reduced_dates2(k) reduced_dates2(k+1)], [pi(k) pi(k+1)], date(i));
        test_n(i) = interp1([reduced_dates2(k) reduced_dates2(k+1)], [dp(k) dp(k+1)], date(i));
        r_n(i) = interp1([reduced_dates2(k) reduced_dates2(k+1)], [r(k) r(k+1)], date(i));
    end
end

reduced_dates2 = date;

date = reduced_dates2(1:end-1);
tests = test_n(7:end);tests(isnan(tests)) = 0;
tp = tp_n(7:end); tp(isnan(tp)) = 0;
ph = ph_n(1:end-1);ph(isnan(ph)) = 0;
pi = pi_n(1:end-1);pi(isnan(pi)) = 0;
dp = dp_n(1:end-1);dp(isnan(dp)) = 0;
r = r_n(1:end-1);r(isnan(r)) = 0;

end 

