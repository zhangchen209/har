function cv = cmax_finite(alpha,K,h,Nsim)

%%parameters


% simu finite sample
T = 10000;
cutoff = 8000;
v = zeros(T,Nsim);

e = normrnd(0,1,[T,Nsim]);
v(1,:)  = e(1,:);

for t=2:T
    v(t,:) = (1-h/T)*v(t-1,:) + e(t,:);
end

y = sum(v(T-cutoff+1:T,:));
v = v(T-cutoff+1:T,:);

F = F_infty(v,y,K);

cv = quantile(F,1-alpha);