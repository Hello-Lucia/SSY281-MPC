function [x,fval]=N1_XX(A,b)
%% This function finds x to minimize the first norm
%% The inputs are Matrix a and vector b
%% The outputs are the optimal x and the norm value (fval)
%% Do not change the inputs and outputs!
%% use linprog to solve the problem
n = length(A);

f = [
        zeros(n,1);
        ones(n,1)
    ];


Ale = [
        A  -eye(n); 
       -A  -eye(n)
    ];
ble = [
        b; 
        -b
    ];

z = linprog(f, Ale, ble);

x = z(1:n);
fval = norm(A*x - b, 1); % 1-Norm

end

