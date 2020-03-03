function P=Pre_14(A,B,S,U)
% A and B are the system matrices x^+=Ax+Bu
% S is the polytope for set S
% U is the polytope for feasible inputs
% P is the polytope Pre(S)

AA = [S.A*A S.A*B; zeros(size(U.A,1),size(S.A*A,2)) U.A];
bb = [S.b; U.b];

P = Polyhedron(AA, bb);

P = P.projection(1:size(A));
end