function Q = quarantine_CL()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
table = readtable("../../data/CL/situacion-epidemiologica-coronavirus-en-castilla-y-leon.csv");
date = table.fecha(end):table.fecha(1);
tp = table.casos_confirmados(end:-1:1);
dp = table.fallecimientos(end:-1:1);
r = table.altas(end:-1:1);

tp_n = zeros(length(tp(1:9:end)), 1); dp_n = zeros(length(dp(1:9:end)), 1);
r_n = zeros(length(r(1:9:end)), 1);

j = 1;
for i = 1:9:length(tp)
    tp_n(j) = sum(tp(i:i+8));
    dp_n(j) = sum(dp(i:i+8));
    r_n(j) = sum(r(i:i+8));
    j = j + 1;
end

j1 = 0;
j2 = 0;
n = true;
for i = 1:length(tp_n)
    if isnan(tp_n(i)) && n
        n = false;
        j1 = i-1;
    elseif ~isnan(tp_n(i)) && ~n
        n = true;
        j2 = i;
        tp_n(j1+1:j2-1) = interp1([j1 j2], [tp_n(j1) tp_n(j2)], j1+1:j2-1);
    end
    
end

for i = 1:length(dp_n)
    if isnan(dp_n(i)) && n
        n = false;
        j1 = i-1;
    elseif ~isnan(dp_n(i)) && ~n
        n = true;
        j2 = i;
        dp_n(j1+1:j2-1) = interp1([j1 j2], [dp_n(j1) dp_n(j2)], j1+1:j2-1);
    end
    
end

for i = 1:length(r_n)
    if isnan(r_n(i)) && n
        n = false;
        j1 = i-1;
    elseif ~isnan(r_n(i)) && ~n
        n = true;
        j2 = i;
        r_n(j1+1:j2-1) = interp1([j1 j2], [r_n(j1) r_n(j2)], j1+1:j2-1);
    end
    
end

tp_n(2:end) = movmean(tp_n(2:end)-tp_n(1:end-1), 7);
dp_n(2:end) = movmean(dp_n(2:end)-dp_n(1:end-1),  7);
r_n(2:end) = movmean(r_n(2:end)-r_n(1:end-1), 7);

Q = tp_n(1:end-3) - dp_n(4:end) - r_n(1:end-3);
Q(2:end) = Q(2:end)-Q(1:end-1);

end

