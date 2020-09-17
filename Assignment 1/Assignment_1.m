%% Task 1.a
clc
clear all
syms x y z phi theta l dtheta dphi m2 m1 g u1 u2 u3 
Z = [0 0 1];
p1 = [x y z].';% what this is????

% Masses positions
%p1 = sym('p_',[3,1],'real');
p2 = [p1(1)+l*cos(theta)*sin(phi)
    p1(2)+l*sin(theta)*sin(phi)
    p1(3)-l*cos(phi)];

% dp1 = [1;1;1] what is this ????? 
dp1 = sym('dp1_',[3,1],'real');% [dp1;dp2;dp3]
dq = [dp1; dtheta; dphi];
q = [p1;theta;phi];

dp2dq = jacobian(p2,q)*dq;

% Kinetic energies
T1 = 0.5*m1*dp1.'*dp1;
T2 = 0.5*m2*(dp2dq).'*dp2dq;
T = T1+T2;

% Potential energies
V1 = m1*g*Z*p1;
V2 = m2*g*Z*p2; 
V = V1+V2;

% Lagrange function
L = T-V;

% gradient(L)_dq
dLdq = jacobian(L,dq).';

% d/dt * gradient(L)_dq
ddLdq= jacobian(dLdq,q)*dq;

% grad(L)_q
dLq = jacobian(L,q).';

% What the forces are ? 3x1 ot 5x1
u = [u1;u2;u3;0;0];


% Euler - Lagrange model
M_v = simplify(u + dLq -ddLdq)
%M_v_latex = latex(M_v)

%% Task 1.b
clc
clear all
syms z phi theta l dtheta dphi m2 m1 g u1 u2 u3 ddl dl
Z = [0 0 1]';

% Masses positions
p1 = sym('p_1',[3,1],'real');
p2 = sym('p_2',[3,1],'real');

% dp1 derivates of masses positions 
dp1 = sym('dp1_',[3,1],'real');
dp2 = sym('dp2_',[3,1],'real');

% dp1 secund derivates of masses positions
ddp1 = sym('ddp1_',[3,1],'real');
ddp2 = sym('ddp2_',[3,1],'real');
 
q = [p1; p2];
dq = [dp1;dp2];
ddq = [ddp1; ddp2];


u = [u1;u2;u3;0;0;0];


% Constraint
e = p1-p2;
C = 0.5*(e.'*e-(l.^2));

% Kinetic energies
T1 = 0.5*m1*dp1'*dp1;
T2 = 0.5*m2*dp2'*dp2;
T = T1+T2;

% Potential energies
V1 = m1*g*Z.'*p1;
V2 = m2*g*Z.'*p2;
V = V1+V2;

% Lagrange function
L = simplify(T-V-z*C)

% Model equations
% Derivatives of constraint
% Derivative of C wrt q
dcq = jacobian(C,q);

% Derivative of C wrt q 
dcdq = dcq*ddq;

% Derivative of C wrt dq
ddcdq = jacobian(dcdq,q)*dq;

% Derivative of C wrt L
dcdL = diff(C,l)*ddl;

% Derivative of C wrt dl
ddcdl = diff(dcdL,l)*dl;

% Derivative of C wrt q,l
dcdql= dcdq + ddcdq + dcdL + ddcdl;

% Derivatives of L 
% Derivatives of L wrt dq
dLdq = jacobian(L,dq);

% Why r we getting zeros?
ddLdq = jacobian(dLdq,q)'*dq;

% Derivative of L wrt q
dLq = jacobian(L,q).';

% Maybe 
Mv2 = simplify(u+dLq - ddLdq)

% c(q,l)=dcdql
%Mv2_latex = latex(Mv2)

%% 2.a
clc
w = jacobian(dLdq,dq);

a = dcq';

S = [w a
    a' 0];

c = [u-ddLdq+jacobian(T,q).'-jacobian(V,q).'
    -ddcdq]

%2.b
ddqz = simplify(inv(S)*c)





