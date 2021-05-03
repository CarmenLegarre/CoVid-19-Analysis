function [reduced_dates, test, tp, ph, pi, dp, r] = generate_data()
table = readtable("/home/carmen/Repositories/CoVid/CoVid-19-Analysis/code/Covid_data_Spain/data/CB/COVID19_historico.csv");
reduced_dates = table.FECHA(:);
dimension = length(reduced_dates);
test_cumulative = table.TOTALTEST(:);
ph = table.TOTALHOSPITALIZADOS(:);
tp_cumulative = table.TOTALCASOS(:);
dp_cumulative = table.FALLECIDOS(:);
pi = table.HOSPITALIZADOSUCI(:);
r_cumulative = table.RECUPERADOS(:);

test = zeros(dimension, 1);
tp = zeros(dimension, 1);
dp = zeros(dimension, 1);

test(1) =  test_cumulative(1);
tp(1) = tp_cumulative(1);
dp(1) = dp_cumulative(1);
for i = 2:dimension
   test(i) = test_cumulative(i)-test_cumulative(i-1);
   tp(i) = tp_cumulative(i)-tp_cumulative(i-1);
   dp(i) = dp_cumulative(i)-dp_cumulative(i-1);  
end

r = r_cumulative(:) - [0 r_cumulative(1:end-1)']';

end

