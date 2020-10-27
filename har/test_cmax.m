clear, clc


alpha = 0.05;
h = 10;
Nsim = 5000;
K = 24;

%%parameters
t_start = 0;          %simulation start time
t_end = 1;          %simuation end time
dt = 0.0002;            %time step
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
for t=t_start+dt:dt:t_end
    i = i+1;
   r1 = randn(1,Nsim,'single');
   %r2 = randn(1,Nsim,'single');
   x(i,:) = exp(-h*dt)*x(i-1,:) + sqrt((c/h*0.5)*(1-exp(-2*h*dt)))*r1;
end

y = dt*ones(1,size(x,1)-1)*x(2:end,:);


F = F_infty(x,y,K);

cv = quantile(F,1-alpha);



