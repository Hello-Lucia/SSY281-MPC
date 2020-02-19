clear; clc;

H = eye(4);
Ale = [
       -1  0  0  0;
        1  0  0  0;
        0  1  0  0;
        0 -1  0  0;
        0  0  1  0;
        0  0 -1  0;
        0  0  0  1;
        0  0  0 -1;
    ];
ble = [
        -2.5; %1000; % Check if thre is almost no ineq const in x1
        1000;%5;
        1;
        1;
        2;
        2;
        2;
        2;
    ];

Aeq = [1 0 -1 0; -0.5 1 0 -1];
beq = [1; 0];

[z, ~, ~, ~, L_mult] = quadprog(H, [], Ale, ble, Aeq, beq)
lambda = L_mult.eqlin;
mu = L_mult.ineqlin;

syms x [2 1] real;
syms u [2 1] real;
X = [x; u];


f = 1/2* X'*H*X;
f_grad = jacobian(f, X);

g = Ale * X - ble;
g_grad = jacobian(g, X);
h = Aeq * X - beq;
h_grad = jacobian(h, X);
%%
n = 4;

mu_n = size(mu,1);
lambda_n = size(lambda,1);
tol = 1e-8;

KKT = [
        abs(f_grad' + g_grad'*mu + h_grad'*lambda) < tol*ones(n,1);
        mu >= zeros(mu_n, 1);
        g <= zeros(mu_n, 1);
        abs(h) < tol*ones(lambda_n, 1);
        mu .* g < tol*ones(mu_n,1);
    ];
%%
% Substitute solution and check logical KKT conditions
KKT = logical(subs(KKT, X, z))











