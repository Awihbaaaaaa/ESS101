% PSS 11 exercise 7.4 (Euler). Written by Robert Hult.
clc; clear; close all

lambda = -2;

% Simulation parameters:
tFinal  = 2;
dt = 0.1;
x_0 = 1;

% Defining the function:
syms x
x_dot = lambda*x;
f = matlabFunction(x_dot,'Vars',{x});

% The exact Solution:
tExact = 0:0.1:tFinal;
xExact = exp(lambda*tExact);

% Explicit Euler Method:
N = tFinal / dt;

%tEuler = zeros(nEuler,1);
xRK1 = zeros(N,1);
xRK1(1) = 1;

for j = 2:N+1
    %tEuler(j) = tEuler(j-1) + dtEuler; % time increment
    xRK1(j) = xRK1(j-1) + dt*f(xRK1(j-1));
end

% Runge Kutta Method of order 2:
a=1;
b1 = 0;
b2 = 1;
c = 0.5;

%tRK2 = zeros(nRK2,1);
xRK2 = zeros(N,1);
xRK2(1) = 1;

for j = 2:N+1
    %tRK2(j) = tRK2(j-1) + dtRK2;
    K1 = f(xRK2(j-1));
    K2 = f(xRK2(j-1)+(a*dt*K1));
    xRK2(j) = xRK2(j-1)+(dt*b1*K1)+(dt*b2*K2);
end

% Runge Kutta Method of order 4:
a1 = 1/2; a2 =0.5; a3 = 1;
b1 = 1/6; b2 = 1/3; b3 = 1/3; b4 = 1/6;
c1 = 0.5; c2 = 0.5; c3 = 1;

% tRK2 = zeros(nRK2,1);
xRK4 = zeros(N,1);
xRK4(1)  = 1;

for j = 2:N+1
    %tRK2(j) = tRK2(j-1) + dtRK2;
    K1 = f(xRK4(j-1));
    K2 = f(xRK4(j-1)+(a1*dt*K1));
    K3 = f(xRK4(j-1)+(a2*dt*K2));
    K4 = f(xRK4(j-1)+(a3*dt*K3));
    xRK4(j) = xRK4(j-1)+(dt*b1*K1)+(dt*b2*K2)+(dt*b3*K3)+(dt*b4*K4);
end

%% Compare the results:
plot(tExact,xExact  ,'marker','.','markersize',20)
hold on
plot(tExact,xRK1  ,'marker','.','markersize',10)
hold on
plot(tExact,xRK2      ,'marker','.','markersize',10)
hold on
plot(tExact,xRK4    ,'marker','.','markersize',10)
hold on
grid on
xlabel('t');
ylabel('x');
legend('Exact','Euler Method','RK2','RK4')

%% Errors comparision
N_table = floor(logspace(1,3,50));

for index = 1:length(N_table)
    N = N_table(index);
    dt(index) = tFinal / N_table(index);
    xRK1 = zeros(N,1)
    xRK2 = zeros(N,1)
    xRK4 = zeros(N,1)
    xRK1(1) = 1; 
    xRK2(1) = 1;
    xRK4(1) = 1;
    
    for j = 2:N+1
        % Euler step
        xRK1(j) = xRK1(j-1) + dt(index)*f(xRK1(j-1));
        
        % RK2 step
        K1 = f(xRK2(j-1));
        K2 = f(xRK2(j-1)+(a*dt(index)*K1));
        xRK2(j) = xRK2(j-1)+(dt(index)*b1*K1)+(dt(index)*b2*K2);
        
        % RK4 step
        K1 = f(xRK4(j-1));
        K2 = f(xRK4(j-1)+(a1*dt(index)*K1));
        K3 = f(xRK4(j-1)+(a2*dt(index)*K2));
        K4 = f(xRK4(j-1)+(a3*dt(index)*K3));
        xRK4(j) = xRK4(j-1)+(dt(index)*b1*K1)+(dt(index)*b2*K2)+(dt(index)*b3*K3)+(dt(index)*b4*K4);

    end
    
    % Compute the global error
    err(1,index) = norm(xRK1(end)-xExact(end),inf);
    err(2,index) = norm(xRK2(end)-xExact(end),inf);
    err(3,index) = norm(xRK4(end)-xExact(end),inf);
end

figure(1);clf; 
loglog(dt,err,'marker','.','markersize',15,'linestyle','none')
grid on
set(gca,'XDir','reverse')
xlabel('stepsize')
ylabel('error')
legend('Euler', 'RK2','RK4')

%% Stability check
close all; clc
lambda = (-0.5:-0.5:-200);
dt = 0.1;
N = tFinal/dt;
err = zeros(3,length(lambda));

for index = 1:length(lambda)
    xdot = lambda(index)*x;
    f = matlabFunction(xdot,'var',{x});
    xRK1 = zeros(N,1);
    xRK2 = zeros(N,1);
    xRK4 = zeros(N,1);
    xRK1(1) = 1; 
    xRK2(1) = 1;
    xRK4(1) = 1;
    
    for j = 2:N+1
        % Euler step
        xRK1(j) = xRK1(j-1) + dt*f(xRK1(j-1));
        
        % RK2 step
        K1 = f(xRK2(j-1));
        K2 = f(xRK2(j-1)+(a*dt*K1));
        xRK2(j) = xRK2(j-1)+(dt*b1*K1)+(dt*b2*K2);
        
        % RK4 step
        K1 = f(xRK4(j-1));
        K2 = f(xRK4(j-1)+(a1*dt*K1));
        K3 = f(xRK4(j-1)+(a2*dt*K2));
        K4 = f(xRK4(j-1)+(a3*dt*K3));
        xRK4(j) = xRK4(j-1)+(dt*b1*K1)+(dt*b2*K2)+(dt*b3*K3)+(dt*b4*K4);

    end
    
    % Compute the global error
    err(1,index) = norm(xRK1(end)-xExact(end),inf);
    err(2,index) = norm(xRK2(end)-xExact(end),inf);
    err(3,index) = norm(xRK4(end)-xExact(end),inf);
    
    if abs(lambda(index)*0.1)<=1
        RK1_unstable=lambda(index);
    end
    
    if abs(lambda(index)*0.1+((lambda(index)*0.1)^2)/2)<=1
        RK2_unstable=lambda(index);
    end
    
    if abs(lambda(index)*0.1+((lambda(index)*0.1)^2/2)+((lambda(index)*0.1)^3)/6+((lambda(index)*0.1)^4)/24)<=1
        RK4_unstable=lambda(index);
    end
end
RK1_unstable
RK2_unstable
RK4_unstable

figure(1);clf; 
plot(lambda,err)
grid on
set(gca,'XDir')
xlabel('stepsize')
ylabel('error')
legend('Euler', 'RK2','RK4')

