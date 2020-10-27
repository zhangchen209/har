function cv = ppAlgorithm(alpha,beta0,X,AR_MaxOrder)


T = size(X,1);
% v = zeros(T,Nsim);
Beta_bar = [beta0;0;0];
delta = 0.05;
% mu = 0;
% sig = 1;
% v(1,:)  = normrnd(mu,sig,[1,Nsim]);
% e = normrnd(0,1,[T,Nsim]);
% 
% for t=2:T
%     v(t,:) = 0.9*v(t-1,:) + e(t,:);
% end


M0 = 500;
M1 = 10;
M2 = 2;
N0 = 200;
N1 = 1000;
N2 = 5000;
rho = (1-delta)*(rand(M0,AR_MaxOrder)*2-ones(M0,AR_MaxOrder));
rhoInit0 = [rho, zeros(M0,1)];
Z0 = normrnd(0,1,[T,N0]);
for n_rho = 1:M0
    rhoInit0(n_rho,AR_MaxOrder+1) = ppCval(alpha,rho(n_rho,1:AR_MaxOrder)',Beta_bar,X,Z0);
end
[~,sortIndex] = sort(rhoInit0(:,AR_MaxOrder+1),1,'descend');
rhoInit0 = rhoInit0(sortIndex,1:AR_MaxOrder);


rhoInit1 = [rhoInit0(1:M1,1:AR_MaxOrder), zeros(M1,1)];
Z1 = normrnd(0,1,[T,N1]);
Finv1 = @(r) -ppCval(alpha,r,Beta_bar,X,Z1);
options = optimset('Display','iter','Maxiter',5);
for n_rho = 1:M1
    [rho_min,fval] = fminsearch(Finv1,rhoInit0(n_rho,:)',options);
    rhoInit1(n_rho,:) = [rho_min',-fval];
end
[~,sortIndex] = sort(rhoInit1(:,AR_MaxOrder+1),1,'descend');
rhoInit1 = rhoInit1(sortIndex,1:AR_MaxOrder);


rhoInit2 = [rhoInit1(1:M2,1:AR_MaxOrder), zeros(M2,1)];
Z2 = normrnd(0,1,[T,N2]);
Finv2 = @(r) -ppCval(alpha,r,Beta_bar,X,Z2);
for n_rho = 1:M2
    [rho_min,fval] = fminsearch(Finv2,rhoInit1(n_rho,:)',options);
    rhoInit2(n_rho,:) = [rho_min',-fval];
end

cv = max(rhoInit2(:,AR_MaxOrder+1));







