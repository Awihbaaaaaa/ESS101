%% 
clear all
close all
clc
load input.mat;
load output.mat;
% In this exercise we will work with the ARX define as 

% We will work with simulated data, hence let's generate the data:
% (For a real SysID problem, data are not simulated but they are real, there is no model generating them, 
% they come from the true system)

N = length(y); % number of data we want to generate

uest = u(1:N/2); % Estimation data
yest = y(1:N/2);
uval = u(N/2+1:end); % Validation data
yval = y(N/2+1:end);

% ARX 1
% y(t) = -a1 y(t-1) -a2 y(t-2) + b0 u(t) + e(t)
PHI = zeros(3,N/2); % The matrix has 3 rows (3 parameters), N/2 columns (time instants)
PHI(:,1) = [ uest(1);0;0];
PHI(:,2) = [ uest(2);yest(1);0];

for i=3:N/2
    PHI(:,i) = [uest(i) ; yest(i-1) ;yest(i-2)] ;  
end

th = (PHI*PHI')\PHI*yest; 

% th should be be a vector with three elements, estimated value for a and
% estimate value for b. Lets pick this three elements separetely:
b1hat = th(1);
a1hat = th(2);
a2hat = th(3);

% Prediction 
%  y(t) = -a1 y(t-1) -a2 y(t-2) + b0 u(t) + e(t)
% Lets do prediction first, Write your code for 1-step-ahead predictor
ypred1 = zeros(N/2,1); % that's the vector where we will store the predicted output
% Simulation
ysim1 = zeros(N/2,1); % that's the vector where we will store the simulated output
%ypred1(1) = [uest(1) 0 0]*th
for i=3:N/2
    ypred1(i) = a1hat*yval(i-1) +a2hat*yval(i-2) + b1hat*uval(i);
    ysim1(i) = a1hat*ysim1(i-1) + a2hat*ysim1(i-2) + b1hat*uval(i);
end

predERROR1 = yval-ypred1;
predRMSE1  = rms(predERROR1)

% compare with real data and compute RMSE
simERROR1 = yval-ysim1;
simRMSE1  = rms(simERROR1)
%% 
%       clc
% y(t) = -a1 y(t-1) -a2 y(t-2) + b0 u(t) + b1 u(t-1) + e(t)
PHI = zeros(4,N/2); % The matrix has 4 rows (4 parameters), N/2 columns (time instants)
PHI(:,1) = [ uest(1);0;0;0];
PHI(:,2) = [ uest(2);uest(1);yest(1);0];

for i=3:N/2
    PHI(:,i) = [uest(i); uest(i-1) ; yest(i-1) ;yest(i-2)] ;  
end

th = (PHI*PHI')\PHI*yest; 

% th should be be a vector with three elements, estimated value for a and
% estimate value for b. Lets pick this three elements separetely:
b0hat = th(1);
b1hat = th(2);
a1hat = th(3);
a2hat = th(4);

% fprintf('====================================\n')
% fprintf('a1hat = %.6f \na2hat = %.6f \nb0hat = %.6f \nb1hat = %.6f \n',a1hat,a2hat,b0hat,b1hat);
% fprintf('====================================\n')

% Prediction 
%  % y(t) = -a1 y(t-1) -a2 y(t-2) + b0 u(t) + b1 u(t-1) + e(t)
% Lets do prediction first, Write your code for 1-step-ahead predictor
ypred2 = zeros(N/2,1); % that's the vector where we will store the predicted output

% Simulation
ysim2 = zeros(N/2,1); % that's the vector where we will store the simulated output

for i=3:N/2
    ypred2(i) = a1hat*yval(i-1) +a2hat*yval(i-2) + b0hat*uval(i) + b1hat*uval(i-1);
    ysim2(i) =  a1hat*ysim2(i-1) +a2hat*ysim2(i-2) + b0hat*uval(i) + b1hat*uval(i-1);
end

predERROR2 = yval-ypred2;
predRMSE2  = rms(predERROR2)
simERROR2 = yval - ysim2;
simRMS2 = rms(simERROR2)

%% 
%clc
close all
% y(t) = -a1 y(t-1) -a2 y(t-2) -a3 y(t-3)+ b1 u(t-1) + e(t)
PHI = zeros(4,N/2); % The matrix has 4 rows (4 parameters), N/2 columns (time instants)
PHI(:,1) = [ 0;0;0;0];
PHI(:,2) = [ uest(1);yest(1);0;0];
PHI(:,3) = [ uest(2);yest(2);yest(1);0];

for i=4:N/2
    PHI(:,i) = [uest(i-1) ; yest(i-1) ; yest(i-2) ; yest(i-3)] ;  
end

th = (PHI*PHI')\PHI*yest; 

% th should be be a vector with three elements, estimated value for a and
% estimate value for b. Lets pick this three elements separetely:
b1hat = th(1);
a1hat = th(2);
a2hat = th(3);
a3hat = th(4);

% % fprintf('====================================\n')
% % fprintf('a1hat = %.6f \na2hat = %.6f \nb1hat = %.6f \na3hat = %.6f \n',a1hat,a2hat,b1hat,a3hat);
% % fprintf('====================================\n')

% Prediction 
% y(t) = -a1 y(t-1) -a2 y(t-2) -a3 y(t-3)+ b1 u(t-1) + e(t)
% Lets do prediction first, Write your code for 1-step-ahead predictor
ypred3 = zeros(N/2,1); % that's the vector where we will store the predicted output
ysim3 = zeros(N/2,1);
for i=4:N/2
    ypred3(i) = a1hat*yval(i-1) +a2hat*yval(i-2) + a3hat*yval(i-3) + b1hat*uval(i-1);
    ysim3(i) =  a1hat*ysim3(i-1) +a2hat*ysim3(i-2) + a3hat*ysim3(i-3) + b1hat*uval(i-1);
end

predERROR3 = yval-ypred3;
predRMSE3  = rms(predERROR3)
simERROR3 = yval-ysim3;
simRMS3 = rms(simERROR3)

%% Comparing
figure (1)
subplot(3,1,1)
plot(yval)
hold on
plot(ypred1)
legend('DATA','Model prediction 1')
title('Output')
xlabel('Samples')
ylabel('output')
grid on

subplot(3,1,2)
plot(yval)
hold on
plot(ypred2)
legend('DATA','Model prediction 2')
title('Output')
xlabel('Samples')
ylabel('output')
grid on

subplot(3,1,3)
plot(yval)
hold on
plot(ypred3)
legend('DATA','Model prediction 3')
title('Output')
xlabel('Samples')
ylabel('output')
grid on 

figure (2)
subplot(3,1,1)
plot(yval)
hold on
plot(ysim1)
legend('DATA','Model simulation 1')
title('Output')
xlabel('Samples')
ylabel('output')
grid on

subplot(3,1,2)
plot(yval)
hold on
plot(ysim2)
legend('DATA','Model simulation 2')
title('Output')
xlabel('Samples')
ylabel('output')
grid on

subplot(3,1,3)
plot(yval)
hold on
plot(ysim3)
legend('DATA','Model simulation 3')
title('Output')
xlabel('Samples')
ylabel('output')
grid on 