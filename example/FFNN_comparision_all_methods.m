% FFNN comparsion for all methods: Polyhedron, Star (exact,
% over-approximation), Zonotope
close all;
clear;
clc;
 
P = ExamplePoly.randHrep('d',2);
%P = Polyhedron('lb',[-1;-1],'ub',[1;1]);
S = Star(P);
Z = Zonotope(P);

nN = [2 3 3 2];
[W, b] =  ReLU.rand_Layers(nN);
R_pe = ReLU.relu_FNN_reach_BFS(P, W, b, 'exact');
R_se = ReLU.relu_FNN_reach_BFS(S, W, b, 'exact');
R_sa = ReLU.relu_FNN_reach_BFS(S, W, b, 'approx');
R_z = ReLU.relu_FNN_reach_BFS(Z, W, b, 'approx');
R_sza = ReLU.relu_FNN_reach_BFS(S, W, b, 'zono');
R_sad = ReLU.relu_FNN_reach_BFS(S, W, b, 'abs_domain');

figure('Name','All methods');
ReLU.plot_set_hold(P, 'Input set','r');
ReLU.plot_set_hold(R_z, 'Zonotope approx', 'c');
%ReLU.plot_set_hold(R_sza, 'Star zonotope approx', 'm');
ReLU.plot_set_hold(R_sad, 'Star abs-domain approx', 'g');
ReLU.plot_set_hold(R_sa, 'Star approx','b');
%ReLU.plot_set_hold(R_pe, 'Polyhedron exact','m');
ReLU.plot_set_hold(R_se, 'Star exact','y');
title('All methods');
%legend('Input set', 'Zonotope approx', 'Star zono approx', 'Star abs-domain approx','Star approx','Star exact');
legend('Input set', 'Zonotope approx', 'Star abs-domain approx','Star approx','Star exact');

figure('Name','Seperate methods');
ReLU.plot_set(P, 'Input set');
ReLU.plot_set(R_pe, 'Polyhedron exact');
ReLU.plot_set(R_se, 'Star exact');
ReLU.plot_set(R_sa, 'Star approx');
ReLU.plot_set(R_z, 'Zonotope approx');
%ReLU.plot_set(R_sza, 'Star zonotope approx');
ReLU.plot_set(R_sad, 'Star abs-domain approx');


