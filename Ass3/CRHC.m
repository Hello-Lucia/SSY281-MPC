function [Z,VN] = CRHC(A,B,N,M,Q,R,Pf,F1,G1,h1,F2,G2,h2,x0)
%% do not chane the inputs and outputs!
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the prediction horizon
%% M is the length of the control horizon
%% Z is the vector of optimal variables and VN is the cost function 
%% F1, G1, h1, F2, G2, h2 are constraint matrices
%% Be aware of the F1, F2, G1, G2, h1, and h2! 
%% x0 is the initial condition


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

if ~isempty(F1)
    Aeq = [F1*Gamma + G1];
    beq = [h1 - F1*Omega*x0];
else
    Aeq = [];
    beq = [];
end


m = size(B,2);
if (M > N)
    error('Control horizon (M) is not allowed to be longer than prediction horizon (N).')
elseif (M < N)
    Aeq = [
        Aeq;
        zeros((N-M)*m,m*(M-1)), repmat(-eye(m),N-M,1), kron(eye(N-M), eye(m))
    ];
    beq = [
        beq;
        zeros((N-M)*m, 1)
    ];
end

if ~isempty(F2)
    Ale = G2 + F2*Gamma;
    ble = h2 - F2*Omega*x0;
else
    Ale = [];
    ble = [];
end


options = optimset('Display', 'none');
[Z, VN] = quadprog(H, f, Ale, ble, Aeq, beq, [], [], [], options);
VN = VN + x0'*Q*x0 + (Omega*x0)'*Q_bar*(Omega*x0);




end
