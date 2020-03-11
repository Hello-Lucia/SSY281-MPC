A = [
    1.2 1
    0   1
];
B = [0; 1];

x_bound = [15; 15];
u_bound = 1;

N = 3;

Q = eye(2);
R = 100;

Pf = Pf_14(A, B, Q, R);

Xf = Polyhedron('A', [1 0; -1 0; 0 1; 0 -1], 'b', [0; 0; 0; 0]);


% Question 1
% write the code to plot X0. In the report, give Pf, the X0, and your
% motivation for choosing this Pf.
sys = LTISystem('A', A, 'B', B);

sys.x.min = -x_bound;
sys.x.max = x_bound;
sys.x.penalty = QuadFunction(Q);

sys.x.with('terminalPenalty');
sys.x.terminalPenalty = QuadFunction(Pf);

sys.x.with('terminalSet');
sys.x.terminalSet = Xf;

sys.u.min = -u_bound;
sys.u.max = u_bound;
sys.u.penalty = QuadFunction(R);

 ctrl_new = MPCController(sys, N);
ectrl_new = ctrl_new.toExplicit();