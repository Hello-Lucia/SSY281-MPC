function [Z,VN]=CRHC1_X(A,B,N,Q,R,Pf,F1,G1,h1,F2,G2,h2,x0)
%% do not chane the inputs and outputs!
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% Z is the vector of optimal variables and VN is the cost function 
%% F1, G1, h1, F2, G2, h2 are constraint matrices
%% Be aware of the F1, F2, G1, G2, h1, and h2! 
%% x0 is the initial condition
Q_bar = kron(eye(N),Q);
n = length(Q);
Q_bar(end-n+1:end, end-n+1:end) = Pf;
R_bar = kron(eye(N),R);


H = 2 * [Q_bar zeros(size(Q_bar,1), size(R_bar,2)); zeros(size(R_bar,1), size(Q_bar,2)) R_bar];

%k1 = kron(-eye(N-1),eye(2));
k1 = -eye(2*N);
k_small = kron(eye(N-1),A);
k1 = k1 + [zeros(n,n*N); k_small zeros(n*(N-1),n)];

k2 = kron(eye(N), B);

q = [-A*x0; zeros(n*(N-1),1)];
Aeq = [k1 k2; F1 G1];
beq = [q; h1];

Ale = [F2 G2];
ble = h2;

options = optimset('Display', 'off');
[Z, VN] = quadprog(H, [], Ale, ble, Aeq, beq, [], [], [], options); 
VN = VN + x0'*Q*x0;


end
