function f = ObjFun_BetaM7(t_actual,params,data,options,priors,...
                                             yinit,unknowns,PropInfections)

NP = params.NumberOfPlaces;
Number = params.NumberOfAgeClasses;
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
params.beta_M = beta_M;   % In Mild conditions
params.beta_H = b*beta_M; % Hospitalized individuals
params.beta_I = c*beta_M; % In ICU

%%%%%

tspan = [t_actual(1),t_actual(end)];
N = params.N;
sigma = params.sigma;

[t,y]=ode45(@(t,y)seir_death_age_beta_b3(t,y, params),tspan,yinit,options);
NewInfections = zeros(length(t_actual)-1,NP);
f =[];
for jj = 1:NP
aux = sigma*sum(y(:,NP*Number+(jj-1)*Number+1:NP*Number+jj*Number),2)*N;
NewInfections(:,jj) = interp1(t,aux,t_actual(2:end)');

aux = data(:,(jj-1)*3+1);
Stirling = 0.5*log(2*pi*aux) + aux.*log(aux) - aux;
f = [f;aux.*log(NewInfections(:,jj)) - NewInfections(:,jj) - Stirling];
end
f = [f;1E-10*(unknowns-priors)'];
f(isnan(f))=zeros;