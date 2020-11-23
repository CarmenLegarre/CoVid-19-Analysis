% Date: 01/11/2020
% Author: Carmen Legarreta
%
% Script to process data from the excel file called
% situacion-epidemiologica.xlsx.
% data to depict; infected, uci, death. 

%% Load the data
data = load('data.mat');
data = data.situacionepidemiologica;
