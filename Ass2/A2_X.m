clear;clc;
%% Note that when you run this file, N2, N3, Pf3, and N5 should be defined in Workspace and have the correct values.
%% Q1: Fill the DP_X.m function using the dynamic programic approach.

%% Q2: Find the shortest N that stabilizes the system using DP_X.m ; define N2 that gives the shortest N, i.e. N2=min(N) subject to the system stability
clc
A=[1.0025 0.1001;0.05 1.0025];
B=[0.005;0.1001];
C=[1 0];
D=0;

Q=[5 0;0 1];
Pf=[5 0;0 1];
R=0.5;

% Find shortest N
N2 = 1;

while 1
    [K, P] = DP_X(A, B, N2, Q, R, Pf);
    poles = eig(A + B*K);
    
    if ~any(poles(:) >= 1)
        fprintf('Solution found with N = %x\n', N2);
        P
        K
        break
    end
    
    if(N2 == 50)
        disp('Solution not found after 50 iterations')
        break
    end
    
    N2 = N2 + 1;
end

%% Q3: Define Pf3 as the solution to the Riccati Equation; define N3 that gives the shortest N, i.e. N3=min(N) subject to the system stability
clc
N3 = 1; % Initial
Pf3 = dare(A, B, Q, R) % from Q2

while 1
    [K, P] = DP_X(A, B, N3, Q, R, Pf3);
    poles = eig(A + B*K);
    
    if ~any(poles(:) >= 1)
        fprintf('Solution found with N = %x\n', N3);
        P
        K
        break
    end
    
    if(N3 == 50)
        disp('Solution not found after 50 iterations')
        break
    end
    N3 = N3 + 1;
end

Pf3 = P;


%% Q4: Fill the BS_XX.m function using the batch solution approach.
%% Q5: Find the shortest N that stabilizes the system using BS_XX.m; define N5 that gives the shortest N, i.e. N5=min(N) subject to the system stability
clc
N5 = 1;
while 1
    [K, P] = BS_X(A, B, N5, Q, R, Pf);
    poles = eig(A + B*K);
    
    if ~any(poles(:) >= 1)
        fprintf('Solution found with N = %x\n', N5);
        P
        K
        break
    end
    
    if(N5 == 50)
        disp('Solution not found after 50 iterations')
        break
    end
    N5 = N5 + 1;
end

%N5=the shortest N!

%% Q6: Use BS_XX.m or DP_XX.m and simulate the system for 20 steps; plot the inputs and the states for these four cases.
tf = 20;
h = 0.1;
time = 1:(tf/h);
x0 = [1;0];

x1(:,1) = x0;
x2(:,1) = x0;
x3(:,1) = x0;
x4(:,1) = x0;

R=0.5;N=5;
[K1,P1] = DP_X(A,B,N,Q,R,Pf);

R=0.5;N=15;
[K2,P2] = DP_X(A,B,N,Q,R,Pf);

R=0.05;N=5;
[K3,P3] = DP_X(A,B,N,Q,R,Pf);

R=0.05;N=15;
[K4,P4] = DP_X(A,B,N,Q,R,Pf);

for k = time
    u1(k) = K1*x1(:,k);
    u2(k) = K2*x2(:,k);
    u3(k) = K3*x3(:,k);
    u4(k) = K4*x4(:,k);
    
    
    x1(:,k+1) = (A+B*K1)*x1(:,k);
    x2(:,k+1) = (A+B*K2)*x2(:,k);
    x3(:,k+1) = (A+B*K3)*x3(:,k);
    x4(:,k+1) = (A+B*K4)*x4(:,k);
    
    y1(:,k) = C*x1(:,k);
    y2(:,k) = C*x2(:,k);
    y3(:,k) = C*x3(:,k);
    y4(:,k) = C*x4(:,k);
end

fig = figure(1);
fig.Position = [200 100 1000 500];

subplot(1,2,1)
plot(time*h,[y1; y2; y3; y4;])
title('States')
xlabel('Time [s]')
ylabel('Angle [deg]')
legend('R=0.5     N=5','R=0.5     N=15', 'R=0.05   N=5', 'R=0.05   N=15')


subplot(1,2,2)
plot(time*h, [u1(:) u2(:) u3(:) u4(:)])
title('Inputs')
xlabel('Time [s]')
ylabel('Torque [Nm]')
legend('u_1', 'u_2', 'u_3', 'u_4')



%% Q7: Fill the CRHC1_X.m function
%% Q8: Fill the CRHC2_X.m function
%% Q9: Solve Q6 using CRHC1_X.m or CRHC2_X.m considering the given constraints for 100 sample times
clc
clf
A=[1.0025 0.1001;0.05 1.0025];
B=[0.005;0.1001];
C=[1 0];
D=0;

Q=[5 0;0 1];
Pf=[5 0;0 1];

tf = 10;
h = 0.1;
time = 1:(tf/h);

F1 = []; 
G1 = [];
h1 = [];

fig = figure(1);
fig.Position = [200 100 1000 500];
for s = 1:4
    if s == 1
        R = 0.5; N = 5; 
    elseif s == 2
        R=0.5; N=15;
    elseif s == 3
        R=0.05; N=5;
    else
        R=0.05; N=15;
    end
    F2 = kron(eye(N), [0 1;0 -1]);
    F2 = [F2; zeros(size(F2))];
    G2 = kron(eye(N), [1; -1]);
    G2 = [zeros(size(G2)); G2];
    h2 = [0.5*ones(2*N,1); 0.7*ones(2*N,1)];

    x1 = [];
    x2 = [];
    u = [];
    x0 = [1; 0];
    for k = time
        [Z, VN] = CRHC1_X(A,B,N,Q,R,Pf,F1,G1,h1,F2,G2,h2,x0);
        x0 = [Z(1); Z(2)];
        x1(k) = Z(1);
        x2(k) = Z(2);
        u(k) = Z(2*N+1);
    end
    subplot(1,2,1)
    hold on
    plot(time*h, rad2deg(x1(:)), 'LineWidth', 1);
    %plot(time*h, x2(:), '--');
    subplot(1,2,2)
    hold on
    plot(time*h, u(:), 'LineWidth',1);
end
subplot(1,2,1)
yl = ylim;
ylim([yl(1)*1.1 yl(2)*1.1])
hold off
title('Output ( State x_1 = \theta )')
xlabel('Time [s]')
ylabel('Angle [deg]')
legend('R=0.5     N=5','R=0.5     N=15', 'R=0.05   N=5', 'R=0.05   N=15')

subplot(1,2,2)
hold off
title('Input (\tau)')
xlabel('Time [s]')
ylabel('Torque [Nm]')
legend('u_1', 'u_2', 'u_3', 'u_4')


%% Unit Test
clc
A = 1; B = 1;
Q = 1; R = 1; Pf = 1;
N = 2;

F2 = zeros(4,2);
G2 = [1 0; -1 0; 0 1; 0 -1];
h2 = [1; 1; 1; 1];

F1 = [];
G1 = [];
h1 = [];

x0 = [-5/3];

[Z_b, VN_b] = CRHC2_X(A,B,N,Q,R,Pf,F1,G1,h1,F2,G2,h2,x0)
[Z_dp, VN_dp] = CRHC1_X(A,B,N,Q,R,Pf,F1,G1,h1,F2,G2,h2,x0)

    



