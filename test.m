clear,clc


T = 50;
delta = 0.05;
rho = (1-delta)*(rand(1,T-1)*2-ones(1,T-1));

sig = zeros(T-1,1); % autocovariance function
sig(1) = rho(1);

% Sigma = eye(T);
% Sigma(2,1) = gamma(1);
% Sigma(1,2) = Sigma(2,1)';

% R2 = 1 + delta;
% gamma(2) = gamma(1)^2/R2+rho(2)*(1-gamma(1)^2/R2);
for k=2:T-1
    r1 = sig(1:k-1);
    r3 = flip(r1);
%     R2 = [R2,flip(gamma(1:k-2));flip(gamma(1:k-2)'),1+delta];
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


