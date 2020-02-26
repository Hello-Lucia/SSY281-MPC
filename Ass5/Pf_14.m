function Pf =Pf_14(A,B,Q,R)
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, and R are the stage cost gains in the cost function
%% Pf is the terminal cost gain

Pf = idare(A,B,Q,R);


end

