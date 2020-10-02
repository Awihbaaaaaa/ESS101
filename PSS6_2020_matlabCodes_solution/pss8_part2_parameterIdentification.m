clc;clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Identification using Least Squares formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the data samples:
load pss8_data_part1.mat

% Now we do some identification!
% First of all lets split the data in estimation and validation sets
% (half and half)
N = length(y);  % number of data
uest = u(1:N/2);
yest = y(1:N/2);
uval = u(N/2+1:end);
yval = y(N/2+1:end);

% As model for identification we will use the 1-step-ahead predictor of
% ARX, which, as we know, is linear in the parameters.
% Hence, first let's write this predictor as a linear regression, in the
% form: y(t) = phi(t).'*theta
% where theta is the parameter vector and phi(t) the regressors vector.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% With pen and paper, write the predictor in that form %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now that we have the predictor in the y(t) = phi(t).'*theta form, we can
% derive the vector form as well: Y = PHI.'*theta, where, now Y is a vector and
% PHI is matrix.

% Let's build the matrix PHI:
PHI = zeros(2,N/2);           % The matrix has 2 rows (2 parameters), N/2 columns (time instants)
PHI(:,1) = [ uest(1) ; 0 ];
for i=2:N/2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Write your code here: %%%%%%%%%%%%
    PHI(:,i) = [uest(i) ; yest(i-1)] ;  
end

%%
% This vector form of the predictor is useful for deriving the least
% squares estimate of theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% With pen and paper, write the LS estimate of theta in %%%%%%%%%%%%
%%%%%% terms of PHI and y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write the expression in code as well, so we get our estimate!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Write your code here: %%%%%%%%%%%%
th = (PHI*PHI')\PHI*yest; 

% th should be be a vector with two elements, estimated value for a and
% estimate value for b. Lets pick this two elements separetely:
bhat = th(1);
ahat = th(2);


% print the results in the command window:
fprintf('====================================\n')
fprintf('ahat = %.6f \nbhat = %.6f \n',ahat,bhat);
fprintf('====================================\n')

save pss8_data_part2.mat  u  y  uest yest uval yval ahat  bhat