%{
A = [-1 -1; -1 1; 1 1; 1 -1];
b = ones(4,1);

P = Polyhedron(A, b);
%nexttile
%P.plot
R = ExamplePoly.randHrep('d', 2);
A = R.A
b = R.b

P1 = Polyhedron('ub', [1;1], 'lb', [-1; -1]);
W = [1 1; 1 -1];
P2 = W*P1;
%nexttile
%P2.plot

P3 = Polyhedron([0 0; 0 0; 0.5 0; 0 0.5; ], ones(4,1));
P3 = P3 + [1;1];
nexttile
P3.plot

%nexttile
%P4 = ReLU.reach_exact(P2)
%P4.plot

%}
U2 = eye(2);%ones(2,2);
U1 = [0.5 0; 0 0.5];
L1 = [0 0; 0 0];
W1 = [1 1; 1 -1];
W2 = [1 1; 1 -1];
W3 = [1 1; 0 1];
U = [1 0; 0 1];
L = [-1 0; 0 -1];
b = ones(2,1)


L = L1*W1*L;
U = U1*W1*U +U2;

L = W2*L
U = W2*U
lb = [0;0];
ub = [2;2];
P = Polyhedron('A',[L;U], 'b',ones(4,1));
P.plot