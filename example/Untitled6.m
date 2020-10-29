close all;
clear;
clc;

P = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);
W = [1 1; 1 -1];
P1 = W*P;
S = Star(P1);
m = size(S.get_V, 2);
C1 = [eye(m); -eye(m)]
d1 = [2*ones(m,1); ones(m,1)]
C = [S.C; C1];
d = [S.d,; d1];
S1 = Star(S.V, C, d);


figure;
S.plot('y');
S1.plot('r');
