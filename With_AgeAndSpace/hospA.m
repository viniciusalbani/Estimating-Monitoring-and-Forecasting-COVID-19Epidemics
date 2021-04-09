function f = hospA(Hosp,TA,t,NP)
f = zeros(length(t),NP);
for jj = 1:NP
f(:,jj) = interp1(TA,Hosp(:,jj),t);
end