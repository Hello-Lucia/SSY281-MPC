H = 2 * [1 0; 0 1];
f = [0 0];

Ale = [1 0];
ble = 0;

Aeq = [1 1];
beq = 1;

Z = quadprog(H,f,Ale,ble,Aeq,beq)