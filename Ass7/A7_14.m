close all;clear;clc;
%% You can use any command in MPT in this assignment

%% Question 1
% write the code to plot X0. In the report, give Pf,  the X0, and your
% motivation for choosing this Pf.
close all; clc

A = [1.2 1; 0 1];
B = [0; 1];
    
u_max = 1;
x_max = [15; 15];

Xlim = Polyhedron('lb', -x_max, 'ub', x_max);
Ulim = Polyhedron('lb', -u_max, 'ub', u_max);

Q = eye(2);             % Weight for states
R = 100;                % Weight for control signals
N = 3;

% Terminal state is the origin
Xf = Polyhedron('A', [1 0; -1 0; 0 1; 0 -1], 'b', [0; 0; 0; 0]);
Pf =  Pf_14(A,B,Q,R);   % Terminal state weight.

sys = LTISystem('A', A, 'B', B);
sys.u.min = -u_max;
sys.u.max = u_max;
sys.x.min = -x_max;
sys.x.max = x_max;

sys.x.penalty = QuadFunction( Q );
sys.u.penalty = QuadFunction( R );
sys.x.with('terminalPenalty');
sys.x.terminalPenalty = QuadFunction( Pf );


XPre = sys.reachableSet('X', Xf, 'direction', 'backward', 'N', N);
figure(1)
plot(XPre, 'color', 'r', 'alpha', 0.5)
legend('X_0')


fig2 = figure(2);
fig2.Units = 'normalized';
fig2.Position = [0.05 0.1 0.8 0.7];

% Without terminal set
ctrl = MPCController(sys, N);
ectrl = ctrl.toExplicit();

subplot(1,2,1)
ectrl.partition.plot();
hold on
Xlim.plot('wire', true, 'linestyle', '-.', 'linewidth', 2);
title('Without terminal set')
fprintf('\n<a href="matlab:ectrl.clicksim()">%s</a>\n', "Simulate without terminal state")

% With terminal set
sys.x.with('terminalSet');
sys.x.terminalSet = Xf;
ctrl_new = MPCController(sys, N);
ectrl_new = ctrl_new.toExplicit();

subplot(1,2,2)
ectrl_new.partition.plot();
hold on
Xlim.plot('wire', true, 'linestyle', '-.', 'linewidth', 2);
title('With terminal set X_0')
fprintf('\n<a href="matlab:ectrl_new.clicksim()">%s</a>\n', "Simulate with terminal state")




%% Question 2
% write the code to plot the requested figures. Provide the figures in your
% report and explain your observation about the error between the state and
% state prediction as N increases.
close all;
x0 = [4; -2.6];

N = [10; 15; 20];

for i = 1:length(N)
    ctrl = MPCController(sys, N(i));
    
    % Prediction (open loop)
    [~, ~, pred ] = ctrl.evaluate(x0);
    
    % Simulate (closed loop)
    sim = ctrl.simulate(x0, N(i));
    
    
    fig = figure(i);
    fig.Units = 'normalized';
    fig.Position = [(i-1)*0.33 0.1 0.15 0.2];
    
    plot(pred.X(1, :), pred.X(2,:), '.', 'color', 'b');
    hold on
    plot(sim.X(1, :), sim.X(2, :), '.', 'color', 'r');
    
    grid on
    title("RH controller with N = " + num2str(N(i)))
    xlabel('x_1')
    ylabel('x_2')
    legend('Prediction', 'Simulation')
    
end

%% Question 3
% no code is needed. Answer in the report



%% Question 4
% write a code that calculates the figures and costs. Provide the figures
% and costs in the report. for costs, provide a table in the report that
% provides all costs for all different methods in the question (4 methods,
% each with three different costs as defined in A7 assignment). If you what
% to use some functions in the code, you can write them in different matlab
% files and submit them with the rest of your files