clear all; clc; close all; format long e; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% We estimate the parameters of a SEIR-type epidemiological model by
%%%% using a maximum a posteriori estimator. All the estimation procedures
%%%% are carried out by LSQNONLIN, although the likelihood function (or
%%%% data misfit) is log-Poisson. The model parameters are estimated from
%%%% the daily records of infections.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup for the ODE solver and the least-square function
tol = 1.e-6;  % ode solver tolerance
options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats',...
                                                         'off','Refine',1);
options2 = optimset('MaxFunEvals',15000,'MaxIter',10000,'TolFun',...
                                                       1e-30,'TolX',1e-30);

options3 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experimental Data

%%%% Daily Cases, Hospitalizations and Deaths
%%%% NYC, USA:

DATA = importdata('boroughs-case-hosp-death_20200822.csv');
data = DATA.data;

%%% We shall delete the last 10 days.

t_actual = 0:size(data,1);

%%%% Smoothing the data - averaging every seven consecutive days (2 times):
data2 = data;
for jj=4:size(data,1)-3
for ii = 1:size(data,2) 
data2(jj,ii) = mean(data(jj-3:jj+3,ii));
end
end

t_span = datetime(2020,2,29) + caldays(0:length(t_actual)-1);

%%%% Total Population:
N = 8622698;      % NYC, USA

%%%% Population proportion on each age range:
Proportion = [0.0762	0.0471	0.0319	0.0605	0.0132
0.1192	0.0621 0.0848	0.0993	0.0182
0.0703	0.0404 0.0458	0.0728	0.0155
0.0233	0.0119 0.0170	0.0233	0.0053
0.0181	0.0094 0.0137	0.0181	0.0036]; 
Proportion = reshape(Proportion,1,size(Proportion,1)*size(Proportion,2));

PropInfections = [385.17	566.79	273.09	411	480.57
2251.41	3477.97	1521.88	3044.96	3667.37
3637.23	5499.35	2467.89	4237.76	4253.65
4111.77	5639.42	2589.84	3882.3	3227.94
4883.21	6340.88	3601.47	4596.62	4165.09];
for jj =1:size(PropInfections,2)
PropInfections(:,jj) = PropInfections(:,jj)/sum(PropInfections(:,jj));
end
factor = [0.328074788	0.298685783	0.202229299	0.304609871	0.300665457
0.872047535	0.941124909	0.835845154	0.939431041	1.066737476
1.117558363	1.160525989	1.163044813	1.103399675	1.090323839
1.419594595	1.277592716	1.328426218	1.190752327	0.881998796
1.623443538	1.463131027	1.892795883	1.391144625	1.165194346];
PropInfections = PropInfections.*factor;
PropInfections = reshape(PropInfections,1,...
                            size(PropInfections,1)*size(PropInfections,2));

PropHosp = [0.093499555	0.092	0.107086614	0.071657754	0.054325956
0.108806974	0.119631584	0.102085911	0.110778902	0.06057086
0.231459942	0.236023517	0.260318801	0.241377976	0.154589372
0.438600666	0.47783765	0.482585752	0.467384339	0.354266212
0.59865115	0.60651289	0.615541459	0.629186935	0.504927976];

factor = [3.387096774	3.228070175	1.837837838	3.435897436	5.4
4.0171875	4.09929078	3.620689655	3.991712707	6.517857143
3.877630553	4.051376147	3.589519651	3.931904161	6.496240602
4.47873633	4.2953125	3.916488223	4.206459054	5.40625
4.894054054	4.712349398	4.775043937	4.780359029	6.9375];
PropHosp = PropHosp.*factor;
PropHosp = reshape(PropHosp,1,size(PropHosp,1)*size(PropHosp,2));

PropICU = ones(size(PropHosp));

PropDeath = [0.019047619	0.027173913	0.014705882	0.014925373	0
0.080124465	0.07482699	0.054545455	0.087889273	0.084931507
0.257286432	0.216485507	0.158556367	0.233483002	0.21412037
0.389582203	0.370316479	0.314379442	0.403345215	0.393063584
0.600839408	0.577181208	0.550607287	0.617848465	0.717717718];
factor = [2	2.5	2	2	1
2.675324675	3.035087719	3	2.760869565	2.583333333
2.65010352	3.025316456	3.128	2.696296296	2.569444444
2.639705882	3.020771513	3.108108108	2.729128015	2.582278481
2.733668342	3.124567474	3.071868583	3.063526835	2.569892473];
PropDeath = PropDeath.*factor;
PropDeath = reshape(PropDeath,1,size(PropDeath,1)*size(PropDeath,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial parameters

I_M0 = 7;      % Potential initial infected and mild at t=0
E_0 = 0;       % initial exposed
I_H0 = 0;      % initial infected and hospitalized at t=0
I_I0 = 0;      % initial infected and in ICU at t=0
R_0 = 0;       % initial recovered 
D_0 = 0;       % initial dead

%  params is a structure used to pass parameters to the
%   ODE solver

S_0 = N-(I_M0+I_I0+I_H0+R_0+D_0+E_0);    % Suceptible pop.,  excluding initial infected 
params.N =   N;  % N = total population

NumberOfAgeClasses = 5;  % total age ranges
NumberOfPlaces = 5;

aux = NumberOfAgeClasses*NumberOfPlaces;
yinit(1:aux) = S_0*Proportion;
yinit(aux+1:2*aux) = E_0*Proportion;
yinit(2*aux+1:3*aux) = I_M0*Proportion;
yinit(3*aux+1:4*aux) = I_H0*Proportion;
yinit(4*aux+1:5*aux) = I_I0*Proportion;
yinit(5*aux+1:6*aux) = R_0*Proportion;
yinit(6*aux+1:7*aux) = D_0*Proportion;
yinit = yinit/N;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Model Parameters

params.sigma = 1./5.1;   % inverse of mean incubation time
params.NumberOfAgeClasses = NumberOfAgeClasses;
params.NumberOfPlaces = NumberOfPlaces;
%--------------------------------------------------------------------------
% Mean time until recovery
nu_M = 1./14;
nu_H = 1./12;
nu_I = 1./9;

%--------------------------------------------------------------------------
% Mean time until death
mu_M = zeros;
mu_H = zeros;
mu_I = 1./7;

%--------------------------------------------------------------------------
% Mean time until passing to another infectious Class
gamma_M = 1./1.2;
gamma_H = 1./3.5;

%--------------------------------------------------------------------------
% Proportion of individuals that will recover
p_M  = (1-PropHosp);
params.p_M = p_M; % In Mild conditions

p_H = 1-PropICU;
params.p_H = p_H; % Hospitalized individuals

p_I = 1-PropDeath;
params.p_I = p_I; % In ICU

%--------------------------------------------------------------------------
% Proportion of individuals that will die
q_M = zeros(1,NumberOfAgeClasses*NumberOfPlaces);
params.q_M = q_M;  % In Mild conditions

q_H = zeros(1,NumberOfAgeClasses*NumberOfPlaces);
params.q_H = q_H; % Hospitalized individuals

%--------------------------------------------------------------------------
%%% RATES
%%% Recovery Rate

Recovery_M = ones-PropHosp;
Recovery_H = (nu_H+mu_H+gamma_H)*p_H;
Recovery_I = ones-PropDeath;

params.Recovery_M = Recovery_M; % in Mild conditions
params.Recovery_H = ones-PropICU;%Recovery_H; % Hospitalized individuals
params.Recovery_I = Recovery_I; % in ICU individuals

%%% Getting Worse Rate

GetWorse_M = PropHosp;
GetWorse_H = PropICU;

params.GetWorse_M = GetWorse_M; % Mild conditions
params.GetWorse_H = GetWorse_H; % Hospitalized individuals

%%% Death Rate

Death_M = (nu_M+mu_M+gamma_M)*q_M;
Death_H = (nu_H+mu_H+gamma_H)*q_H;
Death_I = PropDeath; 

params.Death_M = Death_M; % in Mild conditions
params.Death_H = Death_H; % Hospitalized individuals
params.Death_I = Death_I; % in ICU individuals

%%%%%%%%
Hosp = zeros(length(t_actual),NumberOfPlaces);
Death = zeros(length(t_actual),NumberOfPlaces);
for jj = 1:NumberOfPlaces
Hosp(3:end,jj) = min(0.4,data2(2:end,(jj-1)*3+2)./data2(1:end-1,(jj-1)*3+1));
Death(3:end,jj) = min(1,data2(2:end,(jj-1)*3+3)./data2(1:end-1,(jj-1)*3+2));
end
Death(Hosp==0)=zeros;
params.factorWorse = @(t)hospA(Hosp,t_actual,t,NumberOfPlaces);
params.factorDeath = @(t)deathA(Death,t_actual,t,NumberOfPlaces);

%--------------------------------------------------------------------------
% A priori value for the transmission parameter
R_zero = 1.4*2.8/0.1782;
gamma = 1/18;

%--------------------------------------------------------------------------
% Transmission parameters
beta = 2.2911*R_zero*gamma;
bb = [0.5,0.25,0.20,0.15];
beta_M = 0.5*diag(PropInfections);

for ii=1:NumberOfPlaces
    aux = (ii-1)*NumberOfAgeClasses;
for jj = 1:NumberOfAgeClasses-1
beta_M(aux+jj,aux+jj+1:ii*NumberOfAgeClasses) = ...
                                   PropInfections(aux+jj)*bb(1:end-jj+1);
end
beta_M(aux+1:ii*NumberOfAgeClasses,ii*NumberOfAgeClasses+1:end) = ...
    min(PropInfections(aux+1:ii*NumberOfAgeClasses))*min(bb);
end
beta_M = beta_M + beta_M';
params.beta_M = beta_M; % In Mild conditions
params.beta_H = beta_M; % Hospitalized individuals
params.beta_I = beta_M; % In ICU

params.a = ones;
params.b = 0.1;
params.c = 0.01;
paramsOld = params;
yinitOld = yinit;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating basic parameters from day 1 until 20:
Number = NumberOfAgeClasses;
NP = NumberOfPlaces;
day = 20;
LB1 = [ones(1,NP)/N,1E-3*ones(1,NP)];
UB1 = [50*ones(1,NP)/N,50*ones(1,NP)];
unknowns10=[ones(1,NP)/N,beta*ones(1,NP)];
                          
priors1 = unknowns10;

unknowns20 = [];
for jj = 1:NumberOfPlaces
unknowns20 = [unknowns20,0.5,0.25,0.20,0.15];
end
priors2 = unknowns20;
LB2 = zeros(size(unknowns20));
UB2 = ones(size(unknowns20));

unknowns0 = [unknowns10,unknowns20]; % initial parameters
LB = [LB1,LB2]; % lower bounds
UB = [UB1,UB2]; % upper bounds

priors = [priors1,priors2];

%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:

OF = @(unknowns)ObjFun_InitialPopBetaM7(t_actual(1:day),params,...
data2(1:day-1,:),options,priors,yinit,Proportion,PropInfections,unknowns);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);

I_M = unknowns(1:NP);
N = params.N;
for jj = 1:NP
yinit(2*NP*Number+(jj-1)*Number+1:2*NP*Number+jj*Number) = ...
                       I_M(jj)*PropInfections((jj-1)*Number+1:jj*Number);
yinit((jj-1)*Number+1:jj*Number) = Proportion((jj-1)*Number+1:jj*Number)-...
                       I_M(jj)*PropInfections((jj-1)*Number+1:jj*Number);
end
a = unknowns(NP+1:2*NP);
params.a = a;
beta_M = 0.5*diag(PropInfections);
for ii=1:NP
aux = (ii-1)*Number;
bb = unknowns(2*NP+(ii-1)*(Number-1)+1:2*NP+ii*(Number-1));
for jj = 1:Number-1
beta_M(aux+jj,aux+jj+1:ii*Number) = PropInfections(aux+jj)*bb(1:end-jj+1);
end
beta_M(aux+1:ii*Number,ii*Number+1:end) = ...
                 min(PropInfections(aux+1:ii*Number))*min(bb);
beta_M(aux+1:ii*Number,:) = a(ii)*beta_M(aux+1:ii*Number,:);
end
beta_M = beta_M + beta_M';
params.beta_M = beta_M; % In Mild conditions
params.beta_H = beta_M; % Hospitalized individuals
params.beta_I = beta_M; % In ICU

PARAMS1 = unknowns;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 1
% %%%% until 19:

BETA = zeros(length(t_actual),NP);
BETA(1,:) = ones(1,NP); 
unknowns0 = ones(1,NP);
priors = unknowns0;
yinit2 = yinit;
yb = zeros(length(t_actual),7*Number*NP);
yb(1,:) = yinit;
for jj =1:day-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta5(t_actual2,params,unknowns2,...
                    data2(jj,:),options,priors,BETA(jj,:),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-3*ones(1,NP),10*ones(1,NP),options3);
unknowns0 = unknowns;
priors = unknowns;
BETA(jj+1,:) = unknowns;
beta2 = @(t)betaA(BETA(jj:jj+1,:),t_actual2,t,NP);
[~,y2] = ode45(@(t,y)seir_death_age_beta3(t,y,params,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
disp(num2str(unknowns))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Total Cases for each day
sigma=params.sigma;
NewInfections = zeros(day,NP);
for jj = 1:NP
NewInfections(:,jj) = ...
     sigma*sum(yb(1:day,NP*Number+(jj-1)*Number+1:NP*Number+jj*Number),2)*N;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating basic parameters from day 21 until 80:


LB1 = 1E-3*ones(1,NumberOfPlaces);
UB1 = 50*ones(1,NumberOfPlaces);
unknowns10 = beta*ones(1,NumberOfPlaces);
priors1 = unknowns10;
unknowns20 = [];
for jj = 1:NP
unknowns20 = [unknowns20,0.5,0.25,0.20,0.15];
end
priors2 = unknowns20;
LB2 = zeros(size(unknowns20));
UB2 = ones(size(unknowns20));

unknowns0 = [unknowns10,unknowns20]; % initial parameters
LB = [LB1,LB2]; % lower bounds
UB = [UB1,UB2]; % upper bounds
priors = [priors1,priors2];
yinit2 = y2(end,:);

%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:

OF = @(unknowns)ObjFun_BetaM7(t_actual(day:end),params,...
           data2(day:end,:),options,priors,yinit2,unknowns,PropInfections);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
a = unknowns(1:NP);
params.a = a;
b = params.b;
c = params.c;

beta_M = 0.5*diag(PropInfections);
for ii=1:NP
aux = (ii-1)*Number;
bb = unknowns(NP+(ii-1)*(Number-1)+1:NP+ii*(Number-1));
for jj = 1:Number-1
beta_M(aux+jj,aux+jj+1:ii*Number) = PropInfections(aux+jj)*bb(1:end-jj+1);
end
beta_M(aux+1:ii*Number,ii*Number+1:end) = ...
                             min(PropInfections(aux+1:ii*Number))*min(bb);
beta_M(aux+1:ii*Number,:) = a(ii)*beta_M(aux+1:ii*Number,:);
end
beta_M = beta_M + beta_M';
params.beta_M = beta_M; % In Mild conditions
params.beta_H = beta_M; % Hospitalized individuals
params.beta_I = beta_M; % In ICU

PARAMS2 = unknowns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 21
%%%% until 80:

unknowns0 = unknowns(1:NP);
priors = unknowns0;
for jj = day:length(t_actual)-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns)ObjFun_beta5(t_actual2,params,unknowns,...
                  data2(jj,:),options,priors,BETA(jj,:),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-3*ones(1,NP),10*ones(1,NP),options3);
unknowns0 = unknowns;
priors = unknowns;
BETA(jj+1,:) = unknowns;
beta2 = @(t)betaA(BETA(jj:jj+1,:),t_actual2,t,NP);
[~,y2] = ode45(@(t,y)seir_death_age_beta3(t,y,params,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
disp(num2str(unknowns))
end


%%% Final Number of Cases for each Age Range
%%% Total Number of Deaths for each day

NewInfections = zeros(size(yb,1),NP);
NewDeaths = zeros(size(yb(:,1),1),NP);
NewHospB = zeros(size(yb(:,1),1),NP*Number);
NewHosp = zeros(size(yb(:,1),1),NP);
Deaths = zeros(length(yb(:,1)),NP*Number);

TotalHosp = zeros(NP,Number);
TotalInfections = zeros(NP,Number);
TotalDeaths = zeros(NP,Number);

for ii = 1:NP
NewInfections(:,ii) = ...
     sigma*sum(yb(:,NP*Number+(ii-1)*Number+1:NP*Number+ii*Number),2)*N;
 NewCasesB = sigma*sum(yb(:,NP*Number+(ii-1)*Number+1:NP*Number+ii*Number))'*N;
for jj=1:Number
NewDeaths(:,ii) = NewDeaths(:,ii) + Death_M((ii-1)*NP+jj)*yb(:,2*NP*Number+(ii-1)*NP+jj)...
+ Death_H((ii-1)*NP+jj)*yb(:,3*NP*Number+(ii-1)*NP+jj)...
+ Death_I((ii-1)*NP+jj)*(Death(:,ii).*yb(:,4*NP*Number+(ii-1)*NP+jj));
NewHosp(:,ii) = NewHosp(:,ii) + GetWorse_M((ii-1)*NP+jj)*(Hosp(:,ii).*yb(:,2*NP*Number+(ii-1)*NP+jj));
NewHospB(:,(ii-1)*NP+jj) = GetWorse_M((ii-1)*NP+jj)*(Hosp(:,ii).*yb(:,2*NP*Number+(ii-1)*NP+jj));
Deaths(:,(ii-1)*NP+jj) = Death_M((ii-1)*NP+jj)*yb(:,2*NP*Number+(ii-1)*NP+jj)...
+Death_H((ii-1)*NP+jj)*yb(:,3*NP*NumberOfAgeClasses+(ii-1)*NP+jj)...
+Death_I((ii-1)*NP+jj)*(Death(:,ii).*yb(:,4*NP*NumberOfAgeClasses+(ii-1)*NP+jj));
end
TotalHosp(ii,:) = sum(NewHospB(:,(ii-1)*NP+1:ii*NP))*N;
TotalInfections(ii,:) = ...
      sum(sigma*yb(:,Number*NP+(ii-1)*NP+1:NP*Number+ii*NP))*N;
TotalDeaths(ii,:) = sum(Deaths(:,(ii-1)*NP+1:ii*NP))*N;

disp(['  Infections   Hosp.    Deaths  Borough #',num2str(ii)])
disp(num2str(round([TotalInfections(ii,:)',TotalHosp(ii,:)',TotalDeaths(ii,:)'])))


figure
hold on
title(['New Infections - Borough #',num2str(ii)])
plot(t_actual(2:end),data2(:,(ii-1)*3+1),'--r')
plot(t_actual,NewInfections(:,ii),'b')
hold off

figure
hold on
title(['Hospitalizations - Borough #',num2str(ii)])
plot(t_actual(2:end),data2(:,(ii-1)*3+2),'--r')
plot(t_actual,NewHosp(:,ii)*N,'b')
hold off

figure
hold on
title(['Deaths - Borough #',num2str(ii)])
plot(t_actual(2:end),data2(:,(ii-1)*3+3),'--r')
plot(t_actual,NewDeaths(:,ii)*N,'b')
hold off
end

elapsedtime = toc;
disp(['Elapsed time: ',num2str(elapsedtime/3600),' hours.'])
