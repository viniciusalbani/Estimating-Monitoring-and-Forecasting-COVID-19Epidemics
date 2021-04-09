function f = ObjFun_beta5(t_actual,params,unknowns,data,options,priors,beta,yinit)
Number = params.NumberOfAgeClasses;
NP = params.NumberOfPlaces;
tspan = [t_actual(1),t_actual(end)];
beta = @(t)betaA([beta;unknowns],t_actual,t,NP);

N = params.N;

sigma = params.sigma;

[~,y] = ode45(@(t,y)seir_death_age_beta3(t,y, params,beta),tspan,yinit,options);

f =zeros(1,NP);
for jj = 1:NP
NewInfections = ...
    sigma*sum(y(end,NP*Number+(jj-1)*Number+1:NP*Number+jj*Number),2)*N;

aux = data(end,(jj-1)*3+1);
Stirling = 0.5*log(2*pi*aux) + aux.*log(aux) - aux;
f(jj) = aux.*log(NewInfections) - NewInfections - Stirling;
end
f = [f,1E-1*(unknowns-priors)];
f(isnan(f))=zeros;