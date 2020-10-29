close all;
% clear;
% clc;

P = Polyhedron('lb', [-2;-2], 'ub', [2;2]);

W1 = [1 1; 1 -1];
W2 = [1 1; 1 -1];
W3 = [1 1; 0 1];
W4 = [1 -1];

b1 = [0; 0];
b2 = [0; 0];
b3 = [1; 0];

S = Star(P);
S = S.affineMap(W1, b1);
[lb, ub] = S.getRanges;
nexttile
S.plot
S = ReLU.reach_approx(S, 'abs_domain');
[lb, ub] = S.getRanges;
nexttile
S.plot
S = S.affineMap(W2, b2);
[lb, ub] = S.getRanges;
nexttile
S.plot
S = ReLU.reach_approx(S, 'abs_domain');
[lb, ub] = S.getRanges;
nexttile
S.plot
S = S.affineMap(W3, b3);
[lb, ub] = S.getRanges;
nexttile
S.plot
S = ReLU.reach_approx(S, 'abs_domain');
nexttile
S.plot
% nexttile
% S.plot
disp('exact')
S = Star(P);
S = S.affineMap(W1, b1);
[lb, ub] = S.getRanges;

S = ReLU.reach_exact(S);
[lb, ub] = S.getRanges;

S = S.affineMap(W2, b2);
[lb, ub] = S.getRanges;

S = ReLU.reach_exact(S);
[lb, ub] = S.getRanges;

S = S.affineMap(W3, b3);
[lb, ub] = S.getRanges;
%S = ReLU.reach_exact(S);
%  nexttile
%  S.plot
% S = ReLU.reach_approx(S, 'abs_domain');
