close all;
clear;
clc;


% P = ExamplePoly.randHrep('d',2);
P = Polyhedron('lb',[-1;-1],'ub',[1;1]);
I_z = Zonotope(P);
I_rlxp = RlxPoly(P, inf);
I_rlxs = RlxStar(P, inf);
I_s = Star(P);

nN = [2 3 3 2];
[W, b] =  ReLU.rand_Layers(nN);
% W{1} = [1 1; 1 -1];
% W{2} = [1 1; 1 -1];
% W{3} = [1 1; 1 0];
% 
% b{1} = [0; 0];
% b{2} = [0; 0];
% b{3} = [1; 0];

Z = ReLU.FNN_reach_BFS(I_z, W, b, 'approx');
RlxP = ReLU.FNN_reach_BFS(I_rlxp, W, b, 'approx');
RlxS = ReLU.FNN_reach_BFS(I_rlxs, W, b, 'approx');
S = ReLU.FNN_reach_BFS(I_s, W, b, 'approx');
Se = ReLU.FNN_reach_BFS(I_s, W, b, 'exact');

figure('Name','All methods');

% ReLU.plot_set_hold(RlxP, 'Abs-Domain', 'r');
ReLU.plot_set_hold(RlxS, 'Star with Abs-Domain bounds','b');
ReLU.plot_set_hold(Z, 'Zonotope approx', 'c');
ReLU.plot_set_hold(S, ' ','y');
ReLU.plot_set_hold(Se, ' ','m');
title('All methods');
% legend('Zonotope', 'Star wiht Abs-Domain', 'Star approx', 'Star exact');
legend('Star w/ Abs-Dom bounds','Zonotope','Star approx', 'Star exact');
