clc

syms x1 x2 mu lambda real

L = x1^2 + 2*x2^2 + mu*(x1-2) + lambda*(x1-2*x2-1);
dLdx1 = diff(L,x1);
dLdx2 = diff(L,x2);

x1_sol = solve(dLdx1, x1)
x2_sol = solve(dLdx2, x2)

q = subs(L, [x1 x2], [x1_sol x2_sol]);
q = simplify(q);
symmer(q)

Q = matlabFunction(q, 'vars',[mu lambda]);

[X,Y] = meshgrid(0:0.5:10, -10:0.5:10);
surf(X,Y,Q(X,Y))

% QP of the Langrange
H = -2 * [-1/4 -1/4; -1/4 -3/4];
f = -[-2; -1];

Ale = -[1 0];
ble = 0;

z = quadprog(H, f, Ale, ble);
x_q = [subs([x1_sol; x2_sol], [mu lambda], [z(1) z(2)])];
x_q = double(x_q)



% QP of the original problem
Aeq = [1 -2];
beq = 1;

Ale = [1 0];
ble = 2;

H = 2 * [1 0; 0 2];

x_l = quadprog(H, [], Ale, ble, Aeq, beq)

