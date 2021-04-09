clear all; clc; close all; format long e; tic
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
options2 = optimset('MaxFunEvals',10000,'MaxIter',7000,'TolFun',...
                                                       1e-30,'TolX',1e-30);
options3 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experimental Data

%%%% Daily Cases, Hospitalizations and Deaths
%%%% NYC, USA:

DATA = importdata('case-hosp-death_20200822.csv');
data = DATA.data;
ndays = 0; %%% Keeping the last 10 days for forecast.
data = data(1:end-ndays,:);

%%% We shall delete the last 10 days.

t_actual = 0:size(data,1);

%%%% Smoothing the data - averaging every 7e consecutive days:
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
Proportion = [22.8,38.4,12.8+11.7,8,4.4+2]/100; % Queens-NY, USA
Proportion = [(4510159/8622698)*Proportion,(4112539/8622698)*Proportion]; 
PropInfections = [419.94,2564.42,3945.68,3893.62,4689.98];
PropInfections = PropInfections/max(PropInfections);
factor = [0.873840222,0.883821544,1.005266419,1.182288047,1.364884278];
PropInfections = PropInfections.*factor;
factor = [0.951331497,0.959185259,0.961367542,1.082044759,1.202503164];
PropInfections = PropInfections.*factor;
PropInfections = [0.961360722*0.9362*PropInfections,1.034274427*1.0682*PropInfections];

PropHosp = [0.08591,0.10759,0.23369,0.45678,0.60687];
factor = [4.019354839,3.919477234,4.028474288,4.342996856,4.970625798];
PropHosp = factor.*PropHosp;
factor = [0.235360786,0.244169556,0.23681515,0.246066771,0.237879135];
PropHosp = factor.*PropHosp;
PropHosp = [0.838085254*PropHosp,1.177093331*PropHosp];
PropICU = ones(size(PropHosp));
PropDeath = [0.019261637,0.07830483,0.225129233,0.37855879,0.597250771];
factor = [3.0,3.124463519,2.883783784,3.106930693,3.084936961];
PropDeath = factor.*PropDeath;
factor = [1.5,1.491803279,1.621580547,1.50768738,1.522515147];
PropDeath = factor.*PropDeath;
PropDeath = [0.7681*PropDeath,0.8576*1.1254*PropDeath];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial parameters

I_M0 = 7;      % Potential initial infected and mild at t=0
E_0 = 0;       % initial exposed
I_H0 = 0;      % initial infected and hospitalized at t=0
I_I0 = 0;      % initial infected and in ICU at t=0
R_0 = 0;       % initial recovered 
D_0 = 0;       % initial dead

S_0 = N-(I_M0+I_I0+I_H0+R_0+D_0+E_0);    % Suceptible pop.,  excluding initial infected 
params.N =   N;  % N = total population

NumberOfAgeClasses = 10;  % total age ranges

yinit(1:NumberOfAgeClasses) = S_0*Proportion;
yinit(NumberOfAgeClasses+1:2*NumberOfAgeClasses) = E_0*Proportion;
yinit(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = I_M0*Proportion;
yinit(3*NumberOfAgeClasses+1:4*NumberOfAgeClasses) = I_H0*Proportion;
yinit(4*NumberOfAgeClasses+1:5*NumberOfAgeClasses) = I_I0*Proportion;
yinit(5*NumberOfAgeClasses+1:6*NumberOfAgeClasses) = R_0*Proportion;
yinit(6*NumberOfAgeClasses+1:7*NumberOfAgeClasses) = D_0*Proportion;
yinit = yinit/N;

NumberOfAgeClasses = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Model Parameters

params.sigma = 1./5.1;   % inverse of mean incubation time
params.NumberOfAgeClasses = NumberOfAgeClasses;

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

NumberOfAgeClasses = 10;
%--------------------------------------------------------------------------
% Proportion of individuals that will die
q_M = zeros(1,NumberOfAgeClasses);
params.q_M = q_M;  % In Mild conditions

q_H = zeros(1,NumberOfAgeClasses);
params.q_H = q_H; % Hospitalized individuals

%--------------------------------------------------------------------------
%%% RATES
%%% Recovery Rate

Recovery_M = ones-PropHosp;
Recovery_H = (nu_H+mu_H+gamma_H)*p_H;
Recovery_I = ones-PropDeath;

params.Recovery_M = Recovery_M; % in Mild conditions
params.Recovery_H = ones-PropICU; % Hospitalized individuals
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

%--------------------------------------------------------------------------
%%%% Corrections
Hosp = ones(length(t_actual),1);
Hosp(3:end) = min(1,data2(2:end,2)./data2(1:end-1,1));
Death = ones(length(t_actual),1);
Death(3:end) = data2(2:end,3)./data2(1:end-1,2);
params.factorWorse = @(t)factorWorse(t,t_actual,Hosp);
params.factorDeath = @(t)factorDeath(t,t_actual,Death);

%--------------------------------------------------------------------------
% A priori value for the transmission parameter
R_zero = 1.4*2.8/0.1782;
gamma = 1/18;

%--------------------------------------------------------------------------
% Transmission parameters
beta = 2.2911*R_zero*gamma;
bb = [0.5,0.25,0.20,0.15];
aux1 = 0.5*diag(PropInfections(1:NumberOfAgeClasses/2));
aux2 = 0.5*diag(PropInfections(NumberOfAgeClasses/2+1:end));
aux3 = aux1+aux2;
for jj = 1:NumberOfAgeClasses/2-1
aux1(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
aux2(jj,jj+1:end) = PropInfections(jj+NumberOfAgeClasses/2)*bb(1:end-jj+1);
aux3(jj,jj+1:end) = 0.5*(aux1(jj,jj+1:end)+aux2(jj,jj+1:end));
end
aux3 = aux3 + aux3';
beta_M = [aux1,aux3;zeros(size(aux2)),aux2];
beta_M = beta_M + beta_M';
params.beta_M = beta_M;      % In Mild conditions
params.beta_H = beta_M;  % Hospitalized individuals
params.beta_I = beta_M; % In ICU

params.a = ones;
params.b = 0.1;
params.c = 0.01;

paramsOld = params;
yinitOld = yinit;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating basic parameters from day 1 until 20:

day = 20;

LB1 = [1/N,1E-3];
UB1 = [50/N,10];
unknowns10=[1/N,beta];
priors1 = unknowns10;

unknowns20 = [0.5,0.25,0.20,0.15];
priors2 = unknowns20;
LB2 = zeros(size(unknowns20));
UB2 = ones(size(unknowns20));

unknowns0 = [unknowns10,unknowns20]; % initial parameters
LB = [LB1,LB2]; % lower bounds
UB = [UB1,UB2]; % upper bounds

priors = [priors1,priors2];

%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:

OF = @(unknowns)ObjFun_InitialPopBetaM5(t_actual(1:day),params,...
data2(1:day-1,:),options,priors,yinit,Proportion,PropInfections,unknowns);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);

I_M = unknowns(1);
yinit(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = I_M*PropInfections/sum(PropInfections);
yinit(1:NumberOfAgeClasses) = Proportion-I_M*PropInfections/sum(PropInfections);
params.a = unknowns(2);

bb = unknowns(3:end);
aux1 = 0.5*diag(PropInfections(1:NumberOfAgeClasses/2));
aux2 = 0.5*diag(PropInfections(NumberOfAgeClasses/2+1:end));
aux3 = aux1+aux2;
for jj = 1:NumberOfAgeClasses/2-1
aux1(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
aux2(jj,jj+1:end) = PropInfections(jj+NumberOfAgeClasses/2)*bb(1:end-jj+1);
aux3(jj,jj+1:end) = 0.5*(aux1(jj,jj+1:end)+aux2(jj,jj+1:end));
end
aux3 = aux3 + aux3';
beta_M = [aux1,aux3;zeros(size(aux2)),aux2];
beta_M = beta_M + beta_M';
params.beta_M = params.a*beta_M;
params.beta_H = params.a*beta_M;
params.beta_I = params.a*beta_M;

PARAMS1 = unknowns;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 1
%%%% until 19:

BETA = zeros(length(t_actual),1);
BETA(1) = unknowns(2); 
unknowns0 = unknowns(2);
priors = unknowns0;
yinit2 = yinit;
yb = zeros(length(t_actual),7*NumberOfAgeClasses);
yb(1,:) = yinit;
for jj =1:day-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta5(t_actual2,params,unknowns2,...
                    data2(jj,:),options,priors,BETA(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-3,50,options3);
unknowns0 = unknowns;
priors = unknowns;
BETA(jj+1) = unknowns;
beta2 = @(t)interp1(t_actual2,BETA(jj:jj+1),t);
[~,y2] = ode45(@(t,y)seir_death_age_beta3(t,y,params,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
disp(num2str(unknowns))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating basic parameters from day 21 until 80:

LB1 = 1E-3;
UB1 = 10;
unknowns10 = beta;
priors1 = unknowns10;

unknowns20 = [0.5,0.25,0.20,0.15];
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

OF = @(unknowns)ObjFun_BetaM5(t_actual(day:end),params,...
           data2(day:end,:),options,priors,yinit2,unknowns,PropInfections);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
params.a = unknowns(1);

bb = unknowns(2:end);
aux1 = 0.5*diag(PropInfections(1:NumberOfAgeClasses/2));
aux2 = 0.5*diag(PropInfections(NumberOfAgeClasses/2+1:end));
aux3 = aux1+aux2;
for jj = 1:NumberOfAgeClasses/2-1
aux1(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
aux2(jj,jj+1:end) = PropInfections(jj+NumberOfAgeClasses/2)*bb(1:end-jj+1);
aux3(jj,jj+1:end) = 0.5*(aux1(jj,jj+1:end)+aux2(jj,jj+1:end));
end
aux3 = aux3 + aux3';
beta_M = [aux1,aux3;zeros(size(aux2)),aux2];
beta_M = beta_M + beta_M';
params.beta_M = beta_M;
params.beta_H = beta_M;
params.beta_I = beta_M;

PARAMS2 = unknowns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 21
%%%% until 80:

unknowns0 = unknowns(1);
priors = unknowns0;
for jj = day:length(t_actual)-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns)ObjFun_beta5(t_actual2,params,unknowns,...
                  data2(jj,:),options,priors,BETA(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-3,50,options2);
unknowns0 = unknowns;
priors = unknowns;
BETA(jj+1,:) = unknowns;
beta2 = @(t)interp1(t_actual2,BETA(jj:jj+1),t);
[~,y2] = ode45(@(t,y)seir_death_age_beta3(t,y,params,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
disp(num2str(unknowns))
end

%%% Total Cases for each day
sigma=params.sigma;
NewCases = sigma*sum(yb(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses),2);
   
%%% Final Number of Cases for each Age Range
NewCasesB = sigma*sum(yb(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses)'*N);

%%% Total Number of Deaths for each day
NewDeaths = zeros(size(yb(:,1)));
NewHospB = zeros(size(yb(:,1),1),NumberOfAgeClasses);
NewHosp = zeros(size(yb(:,1)));
Deaths = zeros(length(yb(:,1)),NumberOfAgeClasses);
for jj=1:NumberOfAgeClasses
NewDeaths = NewDeaths + Death_M(jj)*yb(:,2*NumberOfAgeClasses+jj)...
+ Death_H(jj)*yb(:,3*NumberOfAgeClasses+jj)...
+ Death_I(jj)*yb(:,4*NumberOfAgeClasses+jj);
NewHosp = NewHosp + GetWorse_M(jj)*yb(:,2*NumberOfAgeClasses+jj);
NewHospB(:,jj) = GetWorse_M(jj)*yb(:,2*NumberOfAgeClasses+jj);
Deaths(:,jj) = Death_M(jj)*yb(:,2*NumberOfAgeClasses+jj)...
+Death_H(jj)*yb(:,3*NumberOfAgeClasses+jj)...
+Death_I(jj)*yb(:,4*NumberOfAgeClasses+jj);
end

%%% Total Hospitalizations for each Age Range
TotalHosp = sum(NewHospB)'*N; 
TotalInfections = ...
      sum(params.sigma*yb(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses))'*N;
TotalDeaths = sum(Deaths)'*N;

disp('  Infections   Hosp.    Deaths  (By Gender')
disp(num2str(round([TotalInfections(1:5),TotalHosp(1:5),TotalDeaths(1:5)])))
disp(num2str(round([sum(TotalInfections(1:5)),sum(TotalHosp(1:5)),sum(TotalDeaths(1:5))])))
disp(num2str(round([TotalInfections(6:end),TotalHosp(6:end),TotalDeaths(6:end)])))
disp(num2str(round([sum(TotalInfections(6:end)),sum(TotalHosp(6:end)),sum(TotalDeaths(6:end))])))

disp('  Infections   Hosp.    Deaths  (Total)')
disp(num2str(round([TotalInfections(1:5)+TotalInfections(6:end),TotalHosp(1:5)+TotalHosp(6:end),TotalDeaths(1:5)+TotalDeaths(6:end)])))
disp(num2str(round([sum(TotalInfections),sum(TotalHosp),sum(TotalDeaths)])))
%% Estimating Bootstrap Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data; tic;
%%% Bootstrap Sampling
NSamples = 200;
trajectories = zeros(length(t_actual)-1,NSamples);
for jj = 1:NSamples
trajectories(:,jj) = poissrnd(NewCases(2:end)*N);
end

%%% Estimating Bootstrap Parameters
PARAMS1Boot = zeros(NSamples,6);
BETABoot = zeros(NSamples,length(BETA));
yinit2B = zeros(NSamples,7*NumberOfAgeClasses);
PARAMS2Boot = zeros(NSamples,5);
AUX = zeros(size(BETABoot));
t_actualB = t_actual(day:end);

TotalHospBoot = zeros(NumberOfAgeClasses,NSamples); 
TotalInfectionsBoot = zeros(NumberOfAgeClasses,NSamples);
TotalDeathsBoot = zeros(NumberOfAgeClasses,NSamples);

NewCasesBoot = zeros(length(t_actual),NSamples);
NewDeathsBoot = zeros(length(t_actual),NSamples);
NewHospBoot = zeros(length(t_actual),NSamples);


parfor ll = 1:NSamples
% %%% Estimating the transmission constant parameters (M,H,I), the initial
% %%% infecve population (I_M0) and the transmission matrix:
params2=paramsOld;
yinitB = yinitOld;
t_actual2 = t_actual(1:day);
LB1 = [1/N,1E-3];
UB1 = [50/N,10];
unknowns10=[1/N,beta];
priors1 = unknowns10;

unknowns20 = [0.5,0.25,0.20,0.15];
priors2 = unknowns20;
LB2 = zeros(size(unknowns20));
UB2 = ones(size(unknowns20));

unknowns0 = [unknowns10,unknowns20]; % initial parameters
LB = [LB1,LB2]; % lower bounds
UB = [UB1,UB2]; % upper bounds

priors = [priors1,priors2];

OF = @(unknowns)ObjFun_InitialPopBetaM5(t_actual2,params2,...
trajectories(1:day-1,ll),options,priors,yinitB,Proportion,PropInfections,unknowns);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
I_M = unknowns(1);
yinitB(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = I_M*PropInfections/sum(PropInfections);
yinitB(1:NumberOfAgeClasses) = Proportion-I_M*PropInfections/sum(PropInfections);
params2.a = unknowns(2);
bb = unknowns(3:end);
aux1 = 0.5*diag(PropInfections(1:NumberOfAgeClasses/2));
aux2 = 0.5*diag(PropInfections(NumberOfAgeClasses/2+1:NumberOfAgeClasses));
aux3 = aux1+aux2;
for jj = 1:NumberOfAgeClasses/2-1
aux1(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
aux2(jj,jj+1:end) = PropInfections(jj+NumberOfAgeClasses/2)*bb(1:end-jj+1);
aux3(jj,jj+1:end) = 0.5*(aux1(jj,jj+1:end)+aux2(jj,jj+1:end));
end
aux3 = aux3 + aux3';
beta_M = [aux1,aux3;zeros(size(aux2)),aux2];
beta_M = beta_M + beta_M';
params2.beta_M = params2.a.*beta_M;
params2.beta_H = params2.a.*beta_M;
params2.beta_I = params2.a.*beta_M;

PARAMS1Boot(ll,:) = unknowns;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 1
%%%% until 19:

Beta = zeros(length(t_actual),1);
Beta(1) = unknowns(2); 
unknowns0 = unknowns(2);
priors = unknowns0;
yinit2 = yinitB;
yb2 = zeros(length(t_actual),7*NumberOfAgeClasses);
yb2(1,:) = yinit;
for jj =1:day-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta5(t_actual2,params2,unknowns2,...
                    trajectories(jj,ll),options,priors,Beta(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-3,50,options3);
unknowns0 = unknowns;
priors = unknowns;
Beta(jj+1) = unknowns;
beta2 = @(t)interp1(t_actual2,Beta(jj:jj+1),t);
[~,y2B] = ode45(@(t,y)seir_death_age_beta3(t,y,params2,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2B(end,:);
yb2(jj+1,:) = yinit2;
disp(num2str(unknowns))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LB1 = 1E-3;
UB1 = 10;
unknowns10 = beta;
priors1 = unknowns10;
unknowns20B = [0.5,0.25,0.20,0.15];
priors2 = unknowns20B;
LB2 = zeros(size(unknowns20B));
UB2 = ones(size(unknowns20B));
unknowns0 = [unknowns10,unknowns20B]; % initial parameters
LB = [LB1,LB2]; % lower bounds
UB = [UB1,UB2]; % upper bounds

priors = [priors1,priors2];
OF = @(unknowns)ObjFun_BetaM5(t_actualB,params2,...
trajectories(day:end,ll),options,priors,yinit2,unknowns,PropInfections);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
params2.a = unknowns(1);
bb = unknowns(2:end);
aux1 = 0.5*diag(PropInfections(1:NumberOfAgeClasses/2));
aux2 = 0.5*diag(PropInfections(NumberOfAgeClasses/2+1:NumberOfAgeClasses));
aux3 = aux1+aux2;
for jj = 1:NumberOfAgeClasses/2-1
aux1(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
aux2(jj,jj+1:end) = PropInfections(jj+NumberOfAgeClasses/2)*bb(1:end-jj+1);
aux3(jj,jj+1:end) = 0.5*(aux1(jj,jj+1:end)+aux2(jj,jj+1:end));
end
aux3 = aux3 + aux3';
beta_M = [aux1,aux3;zeros(size(aux2)),aux2];
beta_M = beta_M + beta_M';
params2.beta_M = params2.a.*beta_M;
params2.beta_H = params2.a.*beta_M;
params2.beta_I = params2.a.*beta_M;
PARAMS2Boot(ll,:) = unknowns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating a time-dependent transmission parameter (beta) from day 1
%%% until 19:

unknowns0 = Beta(day);
priors = unknowns0;
yb2(day,:) = yinit2;
for jj =day:length(t_actual)-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta5(t_actual2,params2,unknowns2,...
                    trajectories(jj,ll),options,priors,Beta(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-3,50,options3);
unknowns0 = unknowns;
priors = unknowns;
Beta(jj+1) = unknowns;
beta2 = @(t)interp1(t_actual2,Beta(jj:jj+1),t);
[~,y2] = ode45(@(t,y)seir_death_age_beta3(t,y,params2,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb2(jj+1,:) = yinit2;
disp(num2str(unknowns))
end
BETABoot(ll,:) = Beta;



NewCasesBoot(:,ll) = ...
             sigma*sum(yb2(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses),2);

factor = zeros(length(t_actual),1);
factorD = factor;
factorW = params2.factorWorse;
factorDeath = params2.factorDeath;
for jj = 1:length(t_actual)
factor(jj) = factorW(t_actual(jj));
factorD(jj) = factorDeath(t_actual(jj));
end          

%%% Total Number of Deaths for each day
NewHospBb = zeros(size(yb2(:,1),1),NumberOfAgeClasses);
Deathsb = zeros(length(yb2(:,1)),NumberOfAgeClasses);
for jj=1:NumberOfAgeClasses
NewDeathsBoot(:,ll) = NewDeathsBoot(:,ll)...
+ Death_M(jj)*yb2(:,2*NumberOfAgeClasses+jj)...
+ Death_H(jj)*yb2(:,3*NumberOfAgeClasses+jj)...
+ Death_I(jj)*factorD.*yb2(:,4*NumberOfAgeClasses+jj);
NewHospBoot(:,ll) = NewHospBoot(:,ll)...
+ GetWorse_M(jj)*factor.*yb2(:,2*NumberOfAgeClasses+jj);
NewHospBb(:,jj) = GetWorse_M(jj)*factor.*yb2(:,2*NumberOfAgeClasses+jj);
Deathsb(:,jj) = Death_M(jj)*yb2(:,2*NumberOfAgeClasses+jj)...
+Death_H(jj)*yb2(:,3*NumberOfAgeClasses+jj)...
+Death_I(jj)*factorD.*yb2(:,4*NumberOfAgeClasses+jj);
end

%%% Total Hospitalizations for each Age Range
TotalHospBoot(:,ll) = sum(NewHospBb)'*N; 
TotalInfectionsBoot(:,ll) = ...
      sum(sigma*yb2(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses))'*N;
TotalDeathsBoot(:,ll) = sum(Deathsb)'*N;

end
disp(['Elapsed time: ',num2str(toc),' secons'])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CI Evaluation
%%% Beta
aux = sort(BETABoot);
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95BETA = [min(aux);max(aux)];

%%% New Cases:
aux = sort(NewCasesBoot');
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95NewCases = [min(aux);max(aux)];

%%% New Deaths:
aux = sort(NewDeathsBoot');
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95NewDeaths = [min(aux);max(aux)]*N;

%%% New Hospitalizations:
aux = sort(NewHospBoot');
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95NewHosp = [min(aux);max(aux)]*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux = sort(PARAMS1Boot);
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95PARAMS1 = [min(aux);max(aux)];

aux = sort(PARAMS2Boot);
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95PARAMS2 = [min(aux);max(aux)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting Results:

figure
hold on
box on
title('Daily New Infections')
h1=area(t_span,CI95NewCases(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewCases(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data2(:,1),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
plot(t_span,N*NewCases,'r','LineWidth',2)
plot([t_span(day),t_span(day)],[0,1.3*max(data(:,1))],'--k','LineWidth',2)
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Reported','Estimated')
ylabel('Number of Individuals')
xlim([t_span(1),t_span(end)])
ylim([0,1.3*max(data(:,1))])
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
hold off


figure
hold on
box on
title('Daily New Deaths')
h1=area(t_span,CI95NewDeaths(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewDeaths(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data2(:,3),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
plot(t_span,N*NewDeaths,'r','LineWidth',2)
plot([t_span(day),t_span(day)],[0,1.3*max(data(:,3))],'--k','LineWidth',2)
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Reported','Estimated')
ylabel('Number of Individuals')
xlim([t_span(1),t_span(end)])
ylim([0,1.3*max(data(:,3))])
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
hold off

%%%%%%%
figure
hold on
box on
title('Daily New Hospitalizations')
h1=area(t_span,CI95NewHosp(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewHosp(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data2(:,2),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t_span,N*NewHosp,'r','LineWidth',2)
plot([t_span(day),t_span(day)],[0,1.3*max(data(:,2))],'--k','LineWidth',2)
legend('Reported','Estimated')
ylabel('Number of Individuals')
xlim([t_span(1),t_span(end)])
ylim([0,1.3*max(data(:,2))])
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
hold off
