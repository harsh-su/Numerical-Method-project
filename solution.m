 clc, clear;
n = 5;  % eta = n
t = [0 n];
p1 = 0.1;
p2 = 1;
error = 1;
tol = 1e-3;
iterations = 1;

while error > tol
    if iterations == 1
        x = Blasius_RK4([0, 0, p1], t);
        m1 = x(end,2);
        %-------------%-------------%-----------%-----------%---------%----------%
        x = Blasius_RK4([0, 0, p2], t);
        m2 = x(end,2);
    else
        p2  = (p2 -p1)/(m2-m1)*(1 - m1) + p1;
        x = Blasius_RK4([0, 0, p2], t);
        m2 = x(end,2);
        error = abs(1-m2);
    end
    iterations = iterations + 1;
end



function x = Blasius_RK4(y_init, tSpan)

dt              = 1e-3;                             % timestep remains uniform
N               = 1 + (tSpan(2) - tSpan(1))/dt;     % number of time instances where we will find our values
num_steps       = N - 1;
x               = zeros(N, 4);

% we need to initiate the vector y, and also vector t
t               = zeros(N, 1);                              % initiating the time vectors
y               = zeros(N, length(y_init));                 % define the y vector so that it is equal to the number of initial equations at hand

for j1 = 1:length(y_init)
    y(1, j1)      = y_init(j1);
end

for i = 1:num_steps
    t(i+1)      = t(i) + dt;
    k1          = blasius(t(i), y(i,:));                      % define k1; implementation of RK-4 method
    k2          = blasius(t(i) + dt/2, 0.5*k1*dt + y(i,:));   % define k2
    k3          = blasius(t(i) + dt/2, 0.5*k2*dt + y(i,:));   % define k3
    k4          = blasius(t(i) + dt, k3*dt + y(i,:));         % define k4
    y(i+1, :)   = y(i, :) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
%     x(i, 4)     = -(1/2)*y(1)*y(3) ;
end


x(:, 1:3)               = y;
x(:, 4)                 = -(1/2).*y(:, 1).*y(:, 3);

%hold on; plot(t,x(:,1), '-k',t, x(:,2), '--r', t, x(:,3), '-.b');

figure(1), plot(t, y(:, 1), 'r'); hold on;
xlabel('\eta')
ylabel('f')
figure(2), plot(t, y(:, 2), 'g'); hold on;
xlabel('\eta')
ylabel("f'")
figure(3), plot(t, y(:, 3), 'b'); hold on;
xlabel('\eta')
ylabel('f"')

figure(4), plot(t, x(:, 4), 'black'); hold on;
xlabel('\eta')
ylabel("f'''")

end

% Define the blasius equation

function xdot = blasius(t,x)
    xdot        = zeros(size(x));               % initiating xdot

    xdot(1)     = x(2);                         % dx1/dt = x2
    xdot(2)     = x(3);                         % dx2/dt = x3
    xdot(3)     = -(1/2)*x(1)*x(3);             % dx3/dt = -(1/2)*x1*x3
end
