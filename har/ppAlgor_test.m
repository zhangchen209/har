function cv = ppAlgor_test(alpha,beta0,X)

T = size(X,1);
% v = zeros(T,Nsim);
Beta_bar = [beta0;0;0];
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
N1 = 500;
N2 = 2500;
rho = rand(M0,1)*2-ones(M0,1);
rhoInit0 = [rho, zeros(M0,1)];
Z0 = normrnd(0,1,[T,N0]);
for n_rho = 1:M0
    rhoInit0(n_rho,2) = ppCval(alpha,rho(n_rho),Beta_bar,X,Z0);
end
[~,sortIndex] = sort(rhoInit0(:,2),1,'descend');
rhoInit0 = rhoInit0(sortIndex,2);

cv = max(rhoInit0);









