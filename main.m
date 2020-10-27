clear, clc

%% dgp
Tsim = 1000;
T = 50;
I = 5000; % replications
Nsim = 5000; % number of simulations for calculating critival values
alpha = 0.05;
beta0=0.1;
rho = 0.9;
mu = 0;
sig = 1;
K =24;
p = length(mu);
AR_MaxOrder = T-1; % parameter space: max order of AR process

seed = 82019;

v = zeros(T,1);

reject_cmax = zeros(I,1);
reject_pp = zeros(I,1);
power_cmax = zeros(21,1);
power_pp = zeros(21,1);
cv_pp = zeros(I,1);
Tw = zeros(I,1);

phi = zeros(T,K);
for t=1:T  % basis funcitons
    for j=1:K/2
        phi(t,2*j-1) = sqrt(2)*cos(2*j*pi*t/T);
        phi(t,2*j) = sqrt(2)*sin(2*j*pi*t/T);
    end
end

% cv_cmax = cmax(alpha,K,10,Nsim,5000);
b = -0.5:0.05:0.5;
b = b';
B = length(b);
for n_b = 13
    bb = b(n_b);
    
    parfor (i=1:100,16)
        i
        rng(seed+i,'twister');
        e = single(normrnd(0,1,[Tsim,1]));
        u = zeros(Tsim,1);
        u(1)  = e(1);
        X = single(normrnd(mu,sig,[T,1]));
        %X = ones(T,1);
        for t = 2:Tsim
            u(t) = rho*u(t-1) + e(t);
        end
        u = u(Tsim+1-T:end);
        Y = bb*X + u;

        %% cmax test
%         vhat = X.*(Y-X*beta0);
%         F = (sum(vhat,1)).^2.*(K^-1*(sum((phi'*vhat).^2,1)).^-1);
%         reject_cmax(i) = F>cv_cmax;

        %% pp test
        Xtilde = [X,ones(T,1),(-1).^(1:T)'];
        cv_pp(i) = ppAlgorithm(alpha,beta0,Xtilde,AR_MaxOrder);
        J = (1:T-1)';
        Mn = T/10;
        w = (ones(T-1,1) - J/Mn).*(abs(J/Mn)<1);

        Beta = (Xtilde'*Xtilde)^-1*Xtilde'*Y;
        
        Psi = zeros(p+2,p+2);
        vbar_pp = Xtilde.*((Y-Xtilde*Beta)*ones(1,size(Beta,1)));
        for j=1:T-1
            Psi = Psi +T^-1*w(j)*(vbar_pp(j+1:T,:)'*vbar_pp(1:T-j,:));
        end
        Psi = Psi + Psi' + T^-1*(vbar_pp'*vbar_pp);
        Tw(i) = (det(Psi)~=0)*(Beta(1:end-2)-beta0)'*(T*...
            [ones(1,p),0,0]*(Xtilde'*Xtilde)^-1*Psi*...
            (Xtilde'*Xtilde)^-1*[ones(1,p),0,0]')^-1*...
            (Beta(1:end-2)-beta0);
        reject_pp(i) = Tw(i)>cv_pp(i);
    end

    power_cmax(n_b) = mean(reject_cmax);
    power_pp(n_b) = mean(reject_pp);
end

save cvpp_r90