clear;clc;close all;
%% Q1
% define H and V that are correspondent to the polyhedron for H and v representation. Plot the polyhedron and
% explain the difference in the report.
A = [0 1; -1 0; -1 -1; 1 1];
b = [0; 0; 1; 1];

P1 = Polyhedron(A, b);
figure(1)
plot(P1)

H = P1.H;
V = P1.V; % corner points
P2 = Polyhedron('V', V);
figure(2)
plot(P2)

%% Q2
% define Q and P, find the sum and the difference and plot the results.
% Sketch the plots in the report as well
clc

A1 = [0 1; 1 0; 0 -1; -1 0];
b1 = [2; 2; 2; 2];
A2 = [-1 -1; 1 1; 1 -1; -1 1];
b2 = [1; 1; 1; 1];

P = Polyhedron(A1, b1);
Q = Polyhedron(A2, b2);

sum = plus(P,Q)
diff = minus(P,Q)

figure(1)
plot(sum)

figure(2)
plot(diff)

%% Q3
% write a code that shows S is invariant and explain your approach in the
% report
Ain = [1 0 -1 0 1 1 -1 -1; 0 1 0 -1 1 -1 1 -1]';
bin = [1 1 1 1 1.5 1.5 1.5 1.5]';
S = Polyhedron(Ain,bin);

A = [0.8 0.4; -0.4 0.8];
sys = LTISystem('A',A);
XReach = sys.reachableSet('X', S, 'direction', 'forward', 'N', 1);
figure(1)
plot(S, 'color', 'r', XReach, 'color', 'b', 'alpha', 0.5)
legend('State set', 'Reachable set')

%% Q4
% Fill in the Reach_X function and Plot S and its one step reachable set.
% Note that you are not supposed to change the inputs and outputs of the
% function.
% Q3 has to be run before

A = [0.8 0.4; -0.4 0.8];
B = [0; 1];

% For verification
sys = LTISystem('A', A, 'B', B);
sys.u.min = -1;
sys.u.max = 1;

% Feasable inputs 
U = Polyhedron('lb', -1, 'ub', 1);

XReach = Reach_14(A,B,S,U);
XReach2 = sys.reachableSet('X', S, 'direction', 'forward', 'N', 1);

figure(1)
subplot(1,2,1);
plot(S, 'color', 'r', XReach, 'color', 'b', 'alpha', 0.5);
legend('State set', 'Reacable set')
title('User implementation')

subplot(1,2,2);
plot(S, 'color', 'r', XReach2, 'color', 'b', 'alpha', 0.5);
legend('State set', 'Reacable set')
title('MPT implementation')




%% Q5
% Fill in the Pre_X function and Plot S and its Pre set.
% Note that you are not supposed to change the inputs and outputs of the
% function.
% Q4 has to be run before

XPre = Pre_14(A,B,S,U);
XPre2 = sys.reachableSet('X', S, 'direction', 'backward', 'N', 1);

figure(2)
subplot(1,2,1);
plot(S, 'color', 'r', XPre, 'color', 'b', 'alpha', 0.5);
legend('State set', 'Pre set')
title('User implementation')

subplot(1,2,2);
plot(S, 'color', 'r', XPre2, 'color', 'b', 'alpha', 0.5);
legend('State set', 'Pre set')
title('MPT implementation')



%% Q6
clc

disp("====== 6a ======");
A = [0.9 0.4; -0.4 0.9];
B = [0; 1];
Pf = zeros(size(A));
Q = eye(2);
R = 1;

x_ub = [3; 3];
u_ub = 0.1;

x0 = [2; 0];


% part 1: Fill in the function ShorterstN_X.m and use it to find the shortest N
% that is feasible. Note that you are not supposed to change the inputs and outputs of the
% function.

for N = 1:100
    [Z, exitflag] = ShortestN_14(A, B, N, Q, R, Pf,x_ub,u_ub,x0);
    if exitflag == 1
        break
    end
end
N

%%
% part 2: Fill in the function RHCXf_X.m and use it to check the
% feasibility. Note that you are not supposed to change the inputs and outputs of the
% function.
disp("====== 6b ======");

N = 2;

sys = LTISystem('A', A, 'B', B);
sys.u.min = -u_ub;
sys.u.max = u_ub;
sys.x.min = -x_ub;
sys.x.max = x_ub;



Xf = sys.invariantSet();
[Z, exitflag] = RHCXf_14(A, B, N, Q, R, Pf, x_ub, u_ub, Xf, x0)

fprintf("Feasable: %d\n", exitflag);



%%
% part 3: Plot the feasible sets for the initial condition in part 1 and 2
% and plot those sets. Answer to the rest of the question in the report.

Xf26 = Polyhedron('A', [1 0; -1 0; 0 1; 0 -1], 'b', [0; 0; 0; 0]);
Xf2 = sys.invariantSet();

XPre26 = sys.reachableSet('X', Xf26, 'direction', 'backward', 'N', 26);
XPre2 = sys.reachableSet('X', Xf2, 'direction', 'backward', 'N', 2);

figure(1)
plot(XPre26, 'color', 'r', XPre2, 'color', 'r', 'alpha', 0.5);
legend('X_{f-26}', 'X_{f-2}')




























