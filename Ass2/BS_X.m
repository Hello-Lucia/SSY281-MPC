function [K0,P0]=BS_X(A,B,N,Q,R,Pf)
%% do not chane the inputs and outputs!
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% K0 is the controller gain when u(0)=K0x
%% P0 describes the final cost as VN=(x0^T)*P0*x0 

Omega = {};
Omega{1} = A;
for i = 2:N
    Omega{i} = Omega{i-1}*A;
end
Omega = vertcat(Omega{:});

Gamma = cell(N,N);
zeroCell = [0 ; 0]; 

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


Q_bar = [kron(eye(N-1),Q), zeros(size(Q,2)*(N-1),size(Q,2)); zeros(size(Pf,1),size(Pf,2)*(N-1)), Pf];
%Q_bar = kron(eye(N),Q);
%Q_bar(N-2:N,N-2:N) = Pf;
R_bar = kron(eye(N),R);

u0 = -inv(Gamma'*Q_bar*Gamma + R_bar)*Gamma'*Q_bar*Omega;

K0 = u0(1,:);
P0 = Q + Omega'*Q_bar*Omega - Omega'*Q_bar*Gamma*-u0;

end
