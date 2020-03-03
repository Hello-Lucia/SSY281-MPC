function P=Reach_14(A,B,S,U)
% A and B are the system matrices x^+=Ax+Bu
% S is the polytope for set S
% U is the polytope for feasible inputs
% P is the polytope Reach(S)

 P = A*S + B*U;

end