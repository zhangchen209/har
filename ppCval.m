function [cv,Tw] = ppCval(alpha,rho,beta0,X,Z)

% PP18 algorithm
T = size(X,1);
p = size(X,2)-2;
AR_MaxOrder = length(rho);
rho = [rho;zeros(T-1-AR_MaxOrder,1)];
Mn = T/10; %bandwidth equals to 1/10
Nsim = size(Z,2);
Tw = zeros(Nsim,1);




mu0 = X*beta0;

%- compute Sigma via Patial Autocovariance-----------------
delta = 0.05;


sig = zeros(AR_MaxOrder,1); % autocovariance function
sig(1) = rho(1);

% Sigma = eye(T);
% Sigma(2,1) = gamma(1);
% Sigma(1,2) = Sigma(2,1)';

% R2 = 1 + delta;
% sig(2) = sig(1)^2/R2+rho(2)*(1-sig(1)^2/R2);
for k=2:T-1
    r1 = sig(1:k-1);
    r3 = flip(r1);
%     R2 = [R2,flip(sig(1:k-2));flip(sig(1:k-2)'),1+delta];
    R2 = eye(k-1)*(1+delta);
    for i=1:k-1
        R2(i,i+1:k-1) = sig(1:k-i-1)';
        R2(i+1:k-1,i) = R2(i,i+1:k-1)';
    end
    D = (1-r1'*R2^-1*r1)*(1-r3'*R2^-1*r3);
    sig(k) = r1'*R2^-1*r3+rho(k)*D;
end
Sigma = [R2,r3;r3',1+delta];
Sigma = [Sigma,flip(sig);flip(sig'),1+delta];

%-------------------------------
y = mu0*ones(1,Nsim)+chol(Sigma)*Z;

J = (1:T-1)';
w = (ones(T-1,1) - J/Mn).*(abs(J/Mn)<1);

beta = (X'*X)^-1*X'*y;

Psi = zeros(p+2,p+2);
for sim = 1:Nsim
    v = X.*((y(:,sim)-X*beta(:,sim))*ones(1,size(beta,1)));
    for j=1:T-1
        Psi = Psi +T^-1*w(j)*(v(j+1:T,:)'*v(1:T-j,:));
    end
    Psi = Psi + Psi' + T^-1*(v'*v);
    Tw(sim) = (det(Psi)~=0)*(beta(1:end-2,sim)-...
        beta0)'*(T*[ones(1,p),0,0]*(X'*X)^-1*...
        Psi*(X'*X)^-1*[ones(1,p),0,0]')^-1*(beta(1:end-2,sim)-beta0);
end

cv = quantile(Tw,1-alpha);






