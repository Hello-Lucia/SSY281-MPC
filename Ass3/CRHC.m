function [Z,VN]=CRHC2_X(A,B,N,M,Q,R,Pf,F1,G1,h1,F2,G2,h2,x0)
%% do not chane the inputs and outputs!
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% Z is the vector of optimal variables and VN is the cost function 
%% F1, G1, h1, F2, G2, h2 are constraint matrices
%% Be aware of the F1, F2, G1, G2, h1, and h2! 
%% x0 is the initial condition


Omega = {};
Omega{1} = A;
for i = 2:N
    Omega{i} = Omega{i-1}*A;
end
Omega = vertcat(Omega{:});

Gamma = cell(N,N);
zeroCell = zeros(length(A),1); 

for i = 1 : N % Row 
    Gamma(i,i) = {B};
    for k = 1 : i % Column
        if(k>1)
            Gamma(i,k-1) = {A^(i-k+1)*B};
            Gamma(k-1,i) = {zeroCell};
        end
    end
end
Gamma = cell2mat(Gamma);


Q_bar = kron(eye(N),Q);
n = length(Q);
Q_bar(end-n+1:end, end-n+1:end) = Pf;
R_bar = kron(eye(N),R);


H = 2 * (R_bar + Gamma'*Q_bar*Gamma);
f = ( (Omega*x0)'*Q_bar*Gamma + (Q_bar*Omega*x0)'*Gamma)';

if ~~sempty(F1)
    Aeq = [Aeq; F1*Gamma + G1];
    beq = [beq; h1 - F1*Omega*x0];
else
    Aeq = F1*Gamma + G1;
    beq = h1 - F1*Omega*x0;
end


Ale = G2 + F2*Gamma;
ble = h2 - F2*Omega*x0;

options = optimset('Display', 'off');
[Z, VN] = quadprog(H, f, Ale, ble, Aeq, beq, [], [], [], options);
VN = VN + x0'*Q*x0 + (Omega*x0)'*Q_bar*(Omega*x0);




end
