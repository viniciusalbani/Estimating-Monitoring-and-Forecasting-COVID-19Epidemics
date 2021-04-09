# Estimating-Monitoring-and-Forecasting-COVID-19Epidemics

All the scripts listed here are written in MATLAB and are part of a project that estimates susceptible-exposed-infected-recovered-like (SEIR-like) models using COVID-19 datasets found in https://www1.nyc.gov/site/doh/covid/covid-19-data.page. The model accounts for different levels of disease severity and deaths. An updated version of the manuscript draft can be fond in https://arxiv.org/abs/2011.08664. The final version of this article shall be published in Scientific Reports.

The project members are Vinicius V.L. Albani, Roberto M. Velho, and Jorge P. Zubelli.

#Contents

There are three main folders: "With_Age", "With_AgeAndGender" and "With_AgeAndSpace". The scripts in these folder access the daily reports of COVID-19 infections, hospitalizations and deaths during the period 01-Mar-2020 to 22-Aug-2020 (https://www1.nyc.gov/site/doh/covid/covid-19-data.page). 

"With_Age": The scripts in this folder estimate the dynamics of a SEIR-like epidemiological model, accounting for different age rages.
"With_AgeAndGender": The scripts in this folder estimate the dynamics of a SEIR-like epidemiological model, accounting for different age rages and sex.
"With_AgeAndSpace": The scripts in this folder estimate the dynamics of a SEIR-like epidemiological model, accounting for different age rages in the five New York City boroughs.

Inside each folder there is a main file, named as "mySEIR20...", that reads the data and calls the optimization functions.
Objective functions are named as "ObjFun_...". 
Time time interpolation of the time-dependent rates of hospitalization and death are evaluated using, respectively, the functions "factorWorse" and "factorDeath" in the folders  "With_Age" and "With_AgeAndGender". The same rates are evaluated by the functions "hospA" and "deathA" in the folder "With_AgeAndSpace". The time interpolation of the beta parameter (transmission parameter) is evaluated by "betaA" in the folder "With_AgeAndSpace".
The SEIR-like model terms are evaluated using the functions "seir_death_age_beta3" and "seir_death_age_beta_b3". To implement the scripts for the proposed SEIR-like model, we consider as an preliminary example the scripts in https://cs.uwaterloo.ca/~paforsyt/SEIR.html.
How to use this repository
The time-dependent effective reproduction number is evaluated by the function "basic_reproduction_rate_beta2".

#The main objective of these scripts is to estimate the parameters of a SEIR-like model, with time-dependent parameters, accounting for different levels of disease severity and death, as well as age range, sex, and spatial distribution of the disease.


Contact: Prof. Vinicius V.L. Albani: v.albani@ufsc.br.
    
