clear, clc

T = 2000;
nh = 500;
alpha = 0.05;
Nsim = 2000;
v = zeros(T,Nsim);
seed = 872013;
mu = 0;
sig = 1;
v(1,:)  = normrnd(mu,sig,[1,Nsim]);
e = normrnd(0,1,[T,Nsim]);


cv = zeros(3,5);
step_h = T/nh;
h = 0;
for n_h = 1:5
    rng(seed+n_h,'twister');
    for t=2:T
        v(t,:) = (1-h/T)*v(t-1,:) + e(t,:);
    end
    F = zeros(3,Nsim);
    k = 1;
    for K = [6,12,24]
        F(k,:) = F_infty(v,K);
        cv(k,n_h) = quantile(F(k,:),1-alpha);
        k = k+1;
    end
    h = h + step_h;
end


f1 = figure
hold on
plot(cv(1,:));
plot(cv(2,:));
plot(cv(3,:));



