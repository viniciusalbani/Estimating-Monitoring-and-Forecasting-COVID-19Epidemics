function f = deathA(Death,TA,t,NP)
f = zeros(length(t),NP);
for jj = 1:NP
f(:,jj) = interp1(TA,Death(:,jj),t);
end