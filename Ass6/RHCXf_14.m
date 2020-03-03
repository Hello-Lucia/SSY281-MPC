function [Z,exitflag] = RHCXf_14(A,B,N,Q,R,Pf,x_ub,u_ub,Xf,x0)
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% Z is the vector of optimal variables 
%% x0 is the initial condition
%% x_ub is the upper bound for absolute value of x elements
%% u_ub is the upper bound for absolute value of u elements
%% exitflag shows if the quadprog have a solution or not; it is one of quadprog outputs
%% Xf is a polytope that shows the terminal set



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


neXf = size(Xf.Ae, 1);
    
F1 = [zeros(neXf, (N-1)*n) Xf.Ae];
G1 = zeros(neXf, N*m);
h1 = Xf.be;

Aeq = F1*Gamma + G1;
beq = h1 - F1*Omega*x0;





nXf = size(Xf.A, 1);
F2 = [
         eye(N*n);
        -eye(N*n);
        zeros(N*m, N*n);
        zeros(N*m, N*n);
        
        [zeros(nXf, (N-1)*n) Xf.A];
    ];
G2 = [
    zeros(N*n, N*m);
    zeros(N*n, N*m);
     eye(N*m);
    -eye(N*m);

    zeros(nXf, N*m);
];
h2 = [
    kron(ones(2*N, 1), x_ub);
    kron(ones(2*N, 1), u_ub);

    Xf.b
];

Ale = G2 + F2*Gamma;
ble = h2 - F2*Omega*x0;

options = optimset('Display', 'none');
[Z,~, exitflag] = quadprog(H, f, Ale, ble, Aeq, beq, [], [], [], options);



end
