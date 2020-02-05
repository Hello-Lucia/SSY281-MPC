function [K0,P0]=DP_X(A,B,N,Q,R,Pf)
%% do not chane the inputs and outputs!
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% K0 is the controller gain when u(0)=K0x(0)
%% P0 describes the final cost as VN=(x0^T)*P0*x0 

% [P, L, G] = dare(A, B, Q, R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P - The stabalizng solution
% L - closed loop eigenvalues
% G - gain matrix

K = {};
P = {};
P{N+1} = Pf; %End cost


for k = N:-1:1
    K{k} = -inv(R + B'*P{k+1}*B)*B'*P{k+1}*A;
    P{k} = Q + A'*P{k+1}*A - A'*P{k+1}*B*-K{k} ; 
    %[X{k}, L{k}, G{k}] = dare(A,B,Q,R);
end

K0 = K{1}; 
P0 = P{1};



end

