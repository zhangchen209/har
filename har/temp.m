th = 1;
mu = 1.2;
sig = 0.3;
dt = 1e-2;
t = 0:dt:2;             % Time vector
x = zeros(1,length(t)); % Allocate output vector, set initial condition
rng(1);                 % Set random seed
for i = 1:length(t)-1
    x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);