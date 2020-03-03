function [Z,exitflag]=ShortestN_14(A,B,N,Q,R,Pf,x_ub,u_ub,x0)
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% N is the length of the horizon
%% Q, R, and Pf are the gains in the cost function
%% x_ub is the upper bound for absolute value of x elements
%% u_ub is the upper bound for absolute value of u elements
%% x0 is the initial condition
%% Z is the vector of optimal variables 
%% exitflag shows if the quadprog have a solution or not; it is one of quadprog outputs

Gamma = kron(eye(N),B);
Omega = A;
for i=1:N-1
    Gamma = Gamma + kron(diag(ones(N-i,1),-i),A^i*B);
    Omega = [Omega; A^(i+1)];
end

Q_bar = blkdiag( kron(eye(N-1),Q), Pf );
R_bar = kron(eye(N),R);



H = 2 * (R_bar + Gamma'*Q_bar*Gamma);
f = ( (Omega*x0)'*Q_bar*Gamma + (Q_bar*Omega*x0)'*Gamma)';


n = size(A,1);
m = size(B,2);

% Terminal state constraint (end up in origin)
F1 = [zeros(n,2*(N-1)), eye(n)];
G1 = zeros(n,N);
h1 = [0; 0];

Aeq = F1*Gamma + G1;
beq = h1 - F1*Omega*x0;


% State and control constraints
F2 = [kron(eye(N), [1 0; -1 0; 0 1; 0 -1]); zeros(N*2, n*N)];
G2 = [zeros(N*n*2, m); kron(eye(N),[1;-1])];
h2 = [repmat(x_ub, N*n*2, 1); repmat(u_ub, N*2, 1)];

Ale = G2 + F2*Gamma;
ble = h2 - F2*Omega*x0;



options = optimset('Display', 'none');
[Z,~, exitflag] = quadprog(H, f, Ale, ble, Aeq, beq, [], [], [], options);

end
