clc;clear;
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Lets do simulation now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the data samples and identified parameters:
load pss8_data_part2.mat

N = length(y);
NN = N/2;

% lets redifine the data variables so it will be easy to switch from
% validation to estimation sets.
yn = yval;
un = uval;
% yn = yest; 
% un = uest;

ysim = zeros(NN,1); % that's the vector where we will store the simulated output
% Write the code for simulation
for i=2:NN 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Write your code here: %%%%%%%%%%%%%%
    % use the variables un, ahat and bhat.
    ysim(i) = ahat*ysim(i-1) + bhat*un(i);
end

% compare with real data and compute RMSE
simERROR = yn-ysim;
simRMSE  = rms(simERROR)

% plot DATA vs MODEL prediction
figure
subplot(2,1,1)
plot(yn)
hold on
plot(ysim)
legend('DATA','Model simulation')
title('Output')
xlabel('Samples')
ylabel('output')
subplot(2,1,2)
plot(simERROR)
legend('Simulation error')
xlabel('Samples')
ylabel('error')

