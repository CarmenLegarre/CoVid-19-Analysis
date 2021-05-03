These are the steps to be followed to estimate the number of people in quarantine in the FIRST WAVE, the values of the parameters and the validation process.

1- Open Matlab

2- Open three scripts: CoVid/CoVid-19-Analysis/code/parameters_estimation/data_processing/first_wave/data_processing1.m
                       CoVid/CoVid-19-Analysis/code/parameters_estimation/parameter_estimation.m
                       CoVid/CoVid-19-Analysis/code/validation/validation.m
                      
In case you want to evaluate the FOLLOWING WAVES instead of opening data_processing1.m you MUST open the following script: 
CoVid/CoVid-19-Analysis/code/parameters_estimation/data_processing/next_waves/data_processing2.m

3- Run data_processing1.m (or data_processing2.m if you want to analyse the following waves). 3 Figures will appear: the first one will display the number of people 
in quarantine Q, the second and third present the linear regression between Q and the number of recovered people and death people, respectively. Moreover the bash 
will show the values of delta and lambda (model's parameters).

4- Without clearing the workspace run data_estimation.m. Then, the variables lambda_v, k_v and p1_v will contains the values of the parametes lambda, k and p1, 
respectively for the periods in which k and p1 mantain approximately constant.

5- Without clearing the workspace run validation.m. Then, 3 Figures will be displayed: the first one will show the computed Q with the estimated parameters and it 
will be compared with the previously estimated Q. The second and third will compare the computed number of deaths and recovered with the real data, respectively. 
