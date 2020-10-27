function F = F_infty(v,y,K,dt)


T = size(v,1);
phi = zeros(T,K);

% basis funcitons
for t=1:T
    for j=1:K/2
        phi(t,2*j-1) = sqrt(2)*cos(2*j*pi*t/T);
        phi(t,2*j) = sqrt(2)*sin(2*j*pi*t/T);
    end
end

F = (y.^2)*K.*((sum((phi'*v*dt).^2,1)).^-1);


