close all;
clear;
clc;

P = ExamplePoly.randHrep('d',2);
S = Star(P);

nN = [2 3 3 2];
[W, b] = ReLU.rand_Layers(nN);

P = Polyhedron('lb', [0; 0], 'ub', [2; 2]);
U = Star(P);




%[R, ct] = ReLU.FNN_reach_BFS(S, W, b, 'approx')


[safe, R, cntr_ex] = ReLU.verify_FNN_BFS(S, W, b, 'exact', U)



figure;
nexttile;
plot(S, 'r');
hold on;
plot(U, 'g');
hold on;
plot(R, 'b');
nexttile;
plot(S, 'r');
hold on;
plot(cntr_ex, 'y');

%{
nexttile;
C = cntr_ex.C(1:4,:);
d = cntr_ex.d(1:4,:);
H = Star(cntr_ex.V, C, d);
plot(H);
%}
