%% import data from the folder

% table1 contains the numer reported cases with respect time
table1 = readtable('AN/cs_export.xls');
table1 = table1(5::, :)

% table2 contains  the number of reporterd, hospitalized and deaths with
% respect to time.
table2 = readtable('casos_hosp_uci_def_sexo_edad_provres.csv');