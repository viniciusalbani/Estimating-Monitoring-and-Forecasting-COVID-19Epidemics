function R0 = basic_reproduction_rate_beta2(S,params,beta,t)

a = params.b;
b = params.c;

factor = params.factorWorse;
factorD = params.factorDeath;

Number = params.NumberOfAgeClasses;
beta_M = beta*params.beta_M;
beta_H = beta*a*params.beta_H;
beta_I = beta*b*params.beta_I;
sigma = params.sigma;
gamma_M = factor(t).*params.GetWorse_M;
gamma_H = params.GetWorse_H;
nu_M = ones-gamma_M;
nu_H = params.Recovery_H;
mu_M = params.Death_M;
mu_H = params.Death_H;
mu_I = factorD(t).*params.Death_I;
nu_I = ones-mu_I;

f = zeros(4*Number);
for jj = 1:Number
for ii =1:Number
f(jj,(ii-1)*4+1:ii*4) = [0,beta_M(jj,ii)*S(jj),beta_H(jj,ii)*S(jj),beta_I(jj,ii)*S(jj)];
end
end
v = zeros(4*Number);
for jj = 1:Number
aux = diag([sigma,nu_M(jj) + mu_M(jj) + gamma_M(jj),nu_H(jj) + mu_H(jj)...
                                       + gamma_H(jj),nu_I(jj) + mu_I(jj)]);
aux = aux - [zeros(1,4);diag([sigma,gamma_M(jj),gamma_H(jj)]),zeros(3,1)];
v((jj-1)*4+1:jj*4,(jj-1)*4+1:jj*4) = aux;
end

aux = f/v;
eigen = eig(aux);
R0 = max(eigen);
