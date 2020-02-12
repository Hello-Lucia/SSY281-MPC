clear;
close all;
clc;
%% Question 1
A=diag([0.5 0.6 0.5 0.6]);
B=[diag([0.5 0.4]);diag([0.25 0.6])];
C=[1 1 0 0;0 0 1 1];
z_sp=[1;-1];

controllable = rank(ctrb(A,B)) == length(A)
rank(ctrb(A,B))
% calculate the steady state (xs,us) and define them as follows; Do not change
% the names of variables; Report these values in the report

A_s = [eye(length(A)) - A, -B; C, zeros(size(C,1),2)];
b =  [0; 0; 0; 0; z_sp];

x = A_s \ b;

xs1 = x(1:4)
us1 = x(5:end)
%% Question 2
A=diag([0.5 0.6 0.5 0.6]);
B=[0.5;0;0.25;0];
C=[1 1 0 0;0 0 1 1];
z_sp=[1;-1];
% calculate the steady state (xs,us) and define them as follows; Do not change
% the names of variables; Report these values in the report


syms u real; 
n = length(A);
e = sqrt(z_sp - C*inv(eye(n)-A)*B*u);
Vn = e'*e;

% Find minimum analitically
derv = diff(Vn, u);
u2s = solve(derv);

% Find minimum using quadprog
H = 2 * 5/4;
f = -1;

us2 = quadprog(H,f)
xs2 = inv(eye(n)-A)*B*us2

%xs2 =
%us2 = 
%% Question 3
clc
A=diag([0.5 0.6 0.5 0.6]);
B=[diag([0.5 0.4]);diag([0.25 0.6])];
C=[1 1 0 0];
z_sp=1;

% calculate the steady state (xs,us) and define them as follows; Do not change
% the names of variables; Report these values in the report
syms u [2 1] real;
H = eye(2);

Aeq = C*inv(eye(n) - A)*B;
beq = z_sp;

u = quadprog(H,[],[],[],Aeq,beq);

xs3 = inv(eye(n) - A)*B*u
us3 = u
%% second part
tf=50;                
%==========================================================================
% Process model
%==========================================================================

h = 1; %sampling time in minutes

A = [ 0.2681   -0.00338   -0.00728;
      9.703    0.3279   -25.44;
         0         0       1   ];
B = [ -0.00537  0.1655;
       1.297   97.91 ;
       0       -6.637];
C = [ 1 0 0;
      0 1 0;
      0 0 1];
Bp = [-0.1175;
      69.74;
       6.637 ];
   
n = size(A,1); % n is the dimension of the state
m = size(B,2); % m is the dimension of the control signal
p = size(C,1); % p is the dimension of the measured output

d=0.01*[zeros(1*tf/5,1);ones(4*tf/5,1)]; %unmeasured disturbance trajectory

x0 = [0.01;1;0.1]; % initial condition of system's state

%==========================================================================
% Observer model
%==========================================================================
%% choose one case, i.e. a, b, or c and then write the code for that case! for the other ones
% you just need to change the example case!
example = 'c';
switch example
    case 'a'
        nd = 2;
        Bd = zeros(n,nd);
        Cd = [1 0;0 0; 0 1]; 
    case 'b'
        nd=3;
        Bd = zeros(n,nd); 
        Cd = [1 0 0;0 0 1;0 1 0];
    case 'c'
        nd=3;
        %Bd = [0 0 0.1655;0 0 97.91; 0 0 -6.637]; 
        Bd = [zeros(3,2) Bp];
        Cd = [1 0 0;0 0 0;0 1 0];
end

%% Question 4
%Augment the model with constant disturbances; check the detectability for case "a", "b", and "c" and
%report the detectability of each one in the report

% define Ae, Be, and Ce which are the matrices for the augmented system. No
% need to report these in the report

Ae = [A Bd; zeros(nd,n) eye(nd)];

Be = [B; zeros(nd,m)];

Ce = [C Cd];

fprintf("============== System %s ====================\n", example)
if rank(obsv(A, C)) == n %The original system is observable (only detectable is needed)
    disp("The original system is observable")
    fprintf("The augmented system is ")
    detectable = (rank([eye(n)-A -Bd; C Cd]) == n + nd)
    if detectable
        fprintf("detectable\n");
    else
        fprintf(" NOT detectable\n");
    end
else
    disp("The original system is NOT observable")
end


%% Question 5
% Calculate Kalman filter gain and name it Le; no need to report this in
% the report
Q = eye(n+nd);
R = eye(p);

if detectable
   [P, ~, ~] = idare(Ae', Ce', Q, R, [], []);
    Le = P*Ce'*inv(Ce*P*Ce' + R) 
end


%% Question 6
%Target Selector
H = [1 0 0;0 0 1]; 

% Matrices for steady stae target calculation
Ass = [eye(n)-A, -B; H*C zeros(size(H*C,1),m)];
bss = [Bd; -H*Cd];

Mss = inv(Ass) * bss
%note that Mss is defined by [xs;us]=Mss*d and you will use it later on; no need to report this in
% the report

%% Question 7
%==========================================================================
% Setup MPC controller
%==========================================================================
clc
%sp=zeros(tf,3);         % setpoint trajectory

N=10;                   % prediction horizon
M=3;                    % control horizon

Q = diag([1 0.001 1]);  % state penalty
Pf = Q;                 % terminal state penalty
R = 0.01*eye(m);         % control penalty
%==========================================================================
% Simulation
%==========================================================================
% n is the dimension of the state
% m is the dimension of the control signal
% p is the dimension of the measured output
% nd is the dimensions of the disturbances
% Simulation
xe_hat = [zeros(n,1); zeros(nd,1)]; % Augmented state estimate
clear x;
x(:,1) = x0;


    for k = 1:tf
        
        %=============================
        % Calculate steady state target
        %=============================
        d_hat = xe_hat(end-nd+1:end,k);
        xs_us = Mss*d_hat;
        xs = xs_us(1:n);
        us = xs_us(end-m+1:end);
        
        %=============================
        % Solve the QP
        %=============================
        x_hat = xe_hat(1:n, k);
        dx = x_hat - xs;
        
        [du, ~] = CRHC(A,B,N,M,Q,R,Pf,[],[],[],[],[],[],dx);
        u(:, k) = du(1:m) + us;
        
        %=============================
        % Update the observer state
        %=============================
        y(:, k) = C*x(:, k);
        
        %Correction
        xe_hat_corr = xe_hat(:, k) + Le*(y(:,k) - Ce*xe_hat(:,k));
        %Prediction
        xe_hat(:, k+1) = Ae*xe_hat_corr + Be*u(:, k);
        
        %=============================
        % Update the process state
        %=============================
        x(:,k+1) = A*x(:,k) + B*u(:,k) + Bp*d(k);
        
        %=============================        
        % Store current variables in log 
        %=============================
        % Not needed since vector indexing is used instead
        
   end % simulation loop
 
        %%
%==========================================================================
% Plot results
%==========================================================================
clf
% plot the states, the state estimations, and the input and report them
% in the report
figure(1)

% State plot + state estimation
subplot(2,1,1)
plot(x', 'Linewidth', 1)
hold on
plot(xe_hat(1:n,:)', '--', 'Linewidth', 1)
title('States & state estimates')
xlabel('Sample [k]')
ylabel('State value')
legend('c', 'T', 'h', '$\hat{c}$', '$\hat{T}$', '$\hat{h}$','Interpreter','latex','FontSize',14)

% Inputs
subplot(2,1,2)
plot(u', 'Linewidth', 1)
title('Inputs')
xlabel('Sample [k]')
ylabel('Input value')
legend('$u_1$ - Coolant temp.', '$u_2$ - Flow rate','Interpreter','latex','FontSize',14)

% Disturbance plot + disturbance estimation



   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
