function cv = cmax(alpha,K,h,Nsim)

%%parameters
t_start = 0;          %simulation start time
t_end = 1;          %simuation end time
dt = 0.005;            %time step
% h = 0.1;            %relaxation time
c = 1;                 %diffusion constant
x0 = 0;               %initial value for stochastic variable x
%mu = 0;               %mean of stochatic process x

% start_dist = -2.0;    %start of OU pdf 
% end_dist = 2.0;       %end of OU pdf

%time
T = t_start:dt:t_end;
%compute x and y
i = 1;
x = zeros(length(t_start+dt:dt:t_end)+1,Nsim);
x(1,:) = x0;
y = zeros(1,Nsim);
for t=t_start+dt:dt:t_end
    i = i+1;
   r1 = randn(1,Nsim,'single');
   %r2 = randn(1,Nsim,'single');
   x(i,:) = exp(-h*dt)*x(i-1,:) + sqrt((c/h*0.5)*(1-exp(-2*h*dt)))*r1;
%    y = y + x(i-1,:)*h*(1-exp(-dt*h))+sqrt((c*h^3*(dt/h-2*...
%        (1-exp(-dt/h))+0.5*(1-exp(-2*dt/h))))-((0.5*c*h^2)*...
%        (1-exp(-dt/h))^2)^2/((c*h/2)*(1-exp(-2*dt/h))))*r2+...
%        ((0.5*c*h^2)*(1-exp(-dt/h))^2)/(sqrt((c*h/2)*...
%        (1-(exp(-dt/h))^2)))*r1;
end

y = dt*ones(1,size(x,1)-1)*x(2:end,:);
%pdf for OU process
% k = 0; j = start_dist:dt:end_dist;
% for l=start_dist:dt:end_dist
%     k = k + 1;
%     p(k) = sqrt((1/tau)/(pi*c))*exp(-(1/tau)*(l-mu)^2/(c)); 
% end


%% simu finite sample
% v = zeros(2000,Nsim);
% 
% e = normrnd(0,1,[T,Nsim]);
% v(1,:)  = e(1,:);
% 
% for t=2:T
%     v(t,:) = (1-h/T)*v(t-1,:) + e(t,:);
% end


F = F_infty(x,y,K,dt);

cv = quantile(F,1-alpha);



