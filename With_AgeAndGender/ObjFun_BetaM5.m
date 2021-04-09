function f = ObjFun_BetaM5(t_actual,params,data,options,priors,...
                                             yinit,unknowns,PropInfections)
Number = 2*params.NumberOfAgeClasses;
params.a = unknowns(1);

bb = unknowns(2:end);
aux1 = 0.5*diag(PropInfections(1:Number/2));
aux2 = 0.5*diag(PropInfections(Number/2+1:end));
aux3 = aux1+aux2;
for jj = 1:Number/2-1
aux1(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
aux2(jj,jj+1:end) = PropInfections(jj+Number/2)*bb(1:end-jj+1);
aux3(jj,jj+1:end) = 0.5*(aux1(jj,jj+1:end)+aux2(jj,jj+1:end));
end
aux3 = aux3 + aux3';
beta_M = [aux1,aux3;zeros(size(aux2)),aux2];
beta_M = beta_M + beta_M';
params.beta_M = beta_M;
params.beta_H = beta_M;
params.beta_I = beta_M;

tspan = [t_actual(1),t_actual(end)];

N = params.N;


sigma = params.sigma;

[t,y]=ode45(@(t,y)seir_death_age_beta_b3(t,y, params),tspan,yinit,options);
NewInfections = sigma*sum(y(:,Number+1:2*Number),2)*N;
NewInfections = interp1(t,NewInfections,t_actual(2:end)');

% log-Poisson Misfit of Likelihood
Stirling = 0.5*log(2*pi*data(:,1)) + data(:,1).*log(data(:,1)) - data(:,1);
f = data(:,1).*log(NewInfections) - NewInfections - Stirling;
f = [f;1E-10*(unknowns-priors)'];