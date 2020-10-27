
gamma = zeros(T-1,1);
gamma(1) = rho(1);


for k=2:T
    r1 = gamma(1:k-1);
    r3 = flip(r1);
    R2 = eye(k-1);
    for i=1:k-1
        R2(i,i+1:k-1) = gamma(1:k-2)';
        R2(i+1:k-1,i) = R2(i,i+1:k-1)';
    end
    D = sqrt((1-r1'*R2^-1*r1)*(1-r3'*R2^-1*r3));
    gamma(k) = r1'*R2^-1*r3+rho(k)*D;
end
