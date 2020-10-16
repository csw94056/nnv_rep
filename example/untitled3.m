P = Polyhedron('lb', [-1;-1], 'ub', [1;1]);
S = Star(P);
W1 = [1 1; 1 -1];
W2 = [1 1; 1 -1];
W3 = [1 1; 0 1];

b1 = [0; 0];
b2 = [0; 0];
b3 = [1; 0];


S = S.affineMap(W1, b1);
[lb, ub] = S.getRanges
S = ReLU.reach_approx(S, 'abs_domain');
[lb, ub] = S.getRanges
S = S.affineMap(W2, b2);
[lb, ub] = S.getRanges
S = ReLU.reach_approx(S, 'abs_domain');
[lb, ub] = S.getRanges
S = S.affineMap(W3, b3);
[lb, ub] = S.getRanges

S.plot
% S = ReLU.reach_approx(S, 'abs_domain');
