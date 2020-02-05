clear;clc;close all;
%%You are supposed to fill this template. Remeber not to change variable
%%names! When you run this file, you should see all these variables in the Workspace with
%%the right values! 
%% Discrete Model
Ac=[0 1;0.5 0];
Bc=[0;1];
Cc=[1 0];
Dc=0;
h=0.1;

%Use given parameters and find A, B, and C in Question #1
sys = ss(Ac, Bc, Cc, Dc);
sys_d = c2d(sys, h);

A = sys_d.A;
B = sys_d.B;
C = sys_d.C; 
D = 0;

%% Delayed model
% define Aa, Ba, Ca, Da according to Question #2

Aa = [expm(Ac*h) expm(Ac*(h-0.5*h)) * integral(@(t) expm(Ac*t)*Bc, 0, 0.5*h, 'ArrayValued', true); 0 0 0];
Ba = [integral(@(t) expm(Ac*t)*Bc, 0, h-0.5*h, 'ArrayValued', true); 1];
Ca = [1 0 0]; 
Da = 0;

sys_delay = ss(Aa, Ba, Ca, Da, h);

%% Controllability and Observability
%find the rank of controllablity and observability matrices according to
%Question #3


% the following parameters should be the ranks of the observability and 
% the controllability matrices, accordingly
sys2_c = rank(ctrb(sys))
sys2_o = rank(obsv(sys))

sys3_c = rank(ctrb(sys_d))
sys3_o = rank(obsv(sys_d))

sys4_c = rank(ctrb(sys_delay))
sys4_o = rank(obsv(sys_delay))

%% Non observable system
const = 1;
C_nonobsv = [const const*sqrt(2)];

sys_non_obsv = ss(A, B, C_nonobsv, D);
sys_non_obsv_o = rank(obsv(sys_non_obsv))

%% Controller design
lambda1 = -4 + 1i*6;
lambda2 = -4 - 1i*6;

%calculate desired poles for the discrete time system (3) and define them as p1 and
%p2 as it is asked in Question #6
p1 = exp(lambda1 * h);
p2 = exp(lambda2 * h);
p = [p1 p2];

%define the feedback gain for the discrete time system (3) as K1
K1 = place(A, B, p);
sys_d_cl = ss(A-B*K1, B, C, D, h);


%define the feedback gain for the delayed discrete time system (4) as K2
K = [K1, 0];
P = [p 0];
K2 = place(Aa, Ba, P);
sys_delay_cl = ss(Aa-Ba*K, Ba, Ca, Da, h);
sys_delay_cl_improved = ss(Aa-Ba*K2, Ba, Ca, Da, h);


%plot the step response of the systems in one figure. Your figure should
%have labels and legend.
figure(1)
T_final = 5;
[y1, t1] = step(sys_d_cl, T_final);
[y2, t2] = step(sys_delay_cl, T_final);
[y3, t3] = step(sys_delay_cl_improved, T_final);
plot(t1, y1, t2, y2, t3, y3)
title('Step response')
legend('Discrete','Discrete + delay', 'Discrete + delay, improved')
xlabel('t')
ylabel('Theta [\theta]')




%% Steady State
ys = pi/6;
% plot the system output as explained in Question #7. Your figure should
%have labels and legend.

x = sym('x', [3 1], 'real');
u = sym('u', 'real');

y = ys;

A_improved = sys_delay_cl_improved.A;
B_improved = sys_delay_cl_improved.B

steady_state = solve([
    x == A_improved*x + B_improved*u
    y == Ca*x
]);
double([steady_state.x1; steady_state.x2; steady_state.x3])
double(steady_state.u)
u_k = ( (1-Aa(1,1))*pi*(1-Aa(2,2)) - Aa(1,2)*Aa(2,1)*pi) / (6*Aa(1,2)*(Aa(2,3)*Ba(3) + Ba(2) + (1-Aa(2,2))*(Aa(1,3)*Ba(3)+Ba(1)) ))


clear x
x(:,1) = [0; 0; 0];

for k = 1:length(t1)
    x(:, k+1) = (Aa-Ba*K2)*x(:,k) + Ba*steady_state.u;
    y(k) = Ca*x(:,k);
end

figure(2)
title('Steady state with feed back controller')
plot(t1, y, t1, ys*ones(length(t1)), '--r')
xlabel('Time [s]')
ylabel('Angle [deg]')
legend('\theta - Measured angle', '\theta = \pi/6')


%% Disturbance
Bd = [0;1;0];
% define Ae, Be, Ce, and De as asked in Question #8
Ae = [Aa Bd; 0 0 0 1];
Be = [Ba; 0];
Ce = [Ca 0];
De = 0;


sys_aug = ss(Ae, Be, Ce, De, h)
eig(sys_aug)

%define the rank of the controllability and observability matrices as

sys6_c = rank(ctrb(sys_aug))
sys6_o = rank(obsv(sys_aug))


%% Feedback gain
%define the controller gain as K3 according to Question #9
p = [p1 p2 0 1]; % Can't move since we can't control the state. 
K3 = place(Ae, Be, p)

sys_aug_cl = ss(Ae-Be*K3, Be, Ce, De, h);

step(sys_aug_cl)




%% Observer
%define L as the observer gain according to Question #10


p = [0.1 0.2 0.3 0.4];
L = place(Ae', Ce', p)'


