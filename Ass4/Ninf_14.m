function [x,fval]=Ninf_XX(A,b)
%% This function finds x to minimize the infinity norm
%% The inputs are Matrix a and vector b
%% The outputs are the optimal x and the norm value (fval)
%% Do not change the inputs and outputs!
%% use linprog to solve the problem

n = length(A);

f = [
        zeros(n,1);
        1
    ];


Ale = [
        A  -ones(n,1); 
       -A  -ones(n,1)
    ];
ble = [
        b; 
        -b
    ];

z = linprog(f, Ale, ble);

x = z(1:n);
fval = norm(A*x - b, Inf); %Infinity norm
end