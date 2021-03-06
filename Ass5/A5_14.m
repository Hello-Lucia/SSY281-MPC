clc;clear;close all;
%% Q1
% you are supposed to find Matrix S such that V(x(k)) is a Lyapunov
% function that V is decreasing (except at the origin) in the report;
%In addition, define matrix S here as well.
A = [0.5, 1; -0.1, 0.2];
S = idare(A, zeros(2,1), eye(2))
S2 = dlyap(A', eye(2))

%% Q2
% answer in the report;


%% Q3 
clc
disp("==== Question 3a ====");
% part a:
% Define N (the shortest horizon that ...) here as N3. You can use DP_XX.m
% that you have writen in previous assignments. Do note that when I run
% your code, I should be able to see the result and just writing a number
% would not be enough. Mention this N in your report as well.
Q = eye(2);
R = 1;
Pf = Q;
A = [1 1; -1 5];
B = [0; 1];



N3 = 1;
while 1
    [K, P] = DP_14(A, B, N3, Q, R, Pf);
    poles = eig(A + B*K);
    
    if ~any(poles(:) >= 1)
        fprintf('Solution found with N = %x\n', N3);
        P
        K
        break
    end
    
    if(N3 == 50)
        warning('Solution not found after 50 iterations')
        break
    end
    
    N3 = N3 + 1;
end
% part b:
% explain in the report

% part c:
% Fill in the function Pf_X.m in which X is your group number. Motivate
% the concept behind your code in the report.
fprintf("\n==== Question 3c ====\n");
Pf = Pf_14(A,B,Q,R);
K = -inv(R + B'*Pf*B)*B'*Pf*A;
% Double check
[K, P] = DP_14(A, B, 1, Q, R, Pf);
poles = eig(A + B*K);
if any(poles(:) >= 1)
    warning('Found Pf does NOT stabalize with N = 1\n');
else
    fprintf('Found Pf does stablize with N = 1');
    
end
Pf
K


% part d:
% you can use trial and error to find this R; just provide the R that works
% and check the stability by checking the eigenvalues of the closed loop
% system with this new R; define it in the code as Rnew
fprintf("\n==== Question 3d ====\n");
Pf = Q;
Rnew = 0.1337;
N = 1;
[K, P] = DP_14(A, B, N, Q, Rnew, Pf);
poles = eig(A + B*K);

if any(poles(:) >= 1)
    warning('Found R does NOT stabalize with R = %d\n', Rnew);
else
    fprintf('Found R does stablize with R = %d', Rnew);
    Rnew
end

%% Q4 
% write your proof in the report

%% Q5 a
% answer to the question in the report. Do note that you can verify your
% answer by checking it numerically in Matlab (this is just for your own
% and you may not provide any code regarding this)
clc
syms R real positive
syms lambda
K = -inv(R+B'*Pf*B)*B'*Pf*A
test = det(A + B*K - lambda*eye(2))
sol = solve(test, lambda) 

% plug sol into wolfram and find answer

%complex conjugate solution
R_sol = [-16/20 + sqrt(16/20+1/5); -16/20 - sqrt(16/20+1/5)]
%% b 
clc

P = Q + A'*Pf*A - A'*Pf*B*inv(R + B'*Pf*B)*B'*Pf*A;
K = -inv(R+B'*P*B)*B'*P*A;
eigs = det(A + B*K - lambda*eye(2))
sol = solve(eigs, lambda) 
% plug sol into wolfram and find answer

%% Q6
% answer in the report



%% do not comment this! when you run this script, the printed values should be the correct answers!
%S
%N3
%Rnew



