function  yprime = seir_death_age_beta3(t,y,params,beta)
factor = params.factorWorse;
factorD = params.factorDeath;

a = params.b;
b = params.c;

Number = 2*params.NumberOfAgeClasses;
beta_M = beta(t)*params.beta_M;
beta_H = beta(t)*a*params.beta_H;
beta_I = beta(t)*b*params.beta_I;
sigma = params.sigma;
gamma_M = factor(t).*params.GetWorse_M;
gamma_H = params.GetWorse_H;
nu_M = ones-gamma_M;
nu_H = params.Recovery_H;
mu_M = params.Death_M;
mu_H = params.Death_H;
mu_I = factorD(t).*params.Death_I;
nu_I = ones-mu_I;
S = y(1:Number);
E = y(Number+1:2*Number);
I_M = y(2*Number+1:3*Number);
I_H = y(3*Number+1:4*Number);
I_I = y(4*Number+1:5*Number);

yprime = zeros(2*7,1);
for jj = 1:Number
yprime(jj) = -S(jj)*(beta_M(jj,:)*I_M +beta_H(jj,:)*I_H +beta_I(jj,:)*I_I);
yprime(Number+jj) = S(jj)*(beta_M(jj,:)*I_M...
                  + beta_H(jj,:)*I_H + beta_I(jj,:)*I_I) - sigma*E(jj);
yprime(2*Number+jj) = sigma*E(jj)-(nu_M(jj)+mu_M(jj)+gamma_M(jj))*I_M(jj);
yprime(3*Number+jj) = ...
    gamma_M(jj)*I_M(jj) - (nu_H(jj)+mu_H(jj)+gamma_H(jj))*I_H(jj);
yprime(4*Number+jj) = ...
    gamma_H(jj)*I_H(jj)-(nu_I(jj) + mu_I(jj))*I_I(jj);
yprime(5*Number+jj) = nu_M(jj)*I_M(jj)+nu_H(jj)*I_H(jj)+nu_I(jj)*I_I(jj);
yprime(6*Number+jj) = mu_M(jj)*I_M(jj)+mu_H(jj)*I_H(jj)+mu_I(jj)*I_I(jj);
end
