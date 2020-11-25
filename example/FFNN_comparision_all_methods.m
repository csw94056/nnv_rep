% FFNN comparsion for all methods: Polyhedron, Star (exact,
% over-approximation), Zonotope
close all;
clear;
clc;

orange = [0.91 0.41 0.17];
dark_green = [0, 0.5, 0];
gray = [0.5 0.5 0.5];

P = ExamplePoly.randHrep('d',2);
%P = Polyhedron('lb',[-1;-1],'ub',[1;1]);
S = Star(P);
Z = Zonotope(P);

figure('Name','Input set');
ReLU.plot_set(P, 'Input set','r');

nN = [2 3 3 2];
[W, b] =  ReLU.rand_Layers(nN);
R_pe = ReLU.FNN_reach_BFS(P, W, b, 'exact');
R_pa = ReLU.FNN_reach_BFS(P, W, b, 'approx');
R_se = ReLU.FNN_reach_BFS(S, W, b, 'exact');
R_sa = ReLU.FNN_reach_BFS(S, W, b, 'approx');
R_z = ReLU.FNN_reach_BFS(Z, W, b, 'approx');
R_sza = ReLU.FNN_reach_BFS(S, W, b, 'zono');
R_sad = ReLU.FNN_reach_BFS(S, W, b, 'abs_domain');
R_stan = Sigmoid.FNN_reach_BFS(S, W, b, 'tansig');
R_slog = Sigmoid.FNN_reach_BFS(S, W, b, 'logsig');


figure('Name','All methods');
% ReLU.plot_set_hold(R_pa, 'Polyhedron approx','m');
% ReLU.plot_set_hold(R_stan, 'Star tansig approx', gray);
% ReLU.plot_set_hold(R_slog, 'Star logsig approx', dark_green);
ReLU.plot_set_hold(R_z, 'Zonotope approx', 'c');
%ReLU.plot_set_hold(R_sza, 'Star zonotope approx', 'm');
ReLU.plot_set_hold(R_sad, 'Star abs-domain approx', 'g');
ReLU.plot_set_hold(R_sa, 'Star approx','b');
%ReLU.plot_set_hold(R_pe, 'Polyhedron exact','m');
ReLU.plot_set_hold(R_se, 'Star exact','y');
title('All methods');
%legend('Input set', 'Zonotope approx', 'Star zono approx', 'Star abs-domain approx','Star approx','Star exact');
% legend('Input set', 'Polyhedron apporx', 'Star tansig', 'Star logsig', 'Zonotope approx', 'Star abs-domain approx','Star approx','Star exact');
legend('Zonotope', 'Abs-domain', 'Star approx', 'Star exact');


% figure('Name','Seperate methods');
% ReLU.plot_set(P, 'Input set');
% ReLU.plot_set(R_pe, 'Polyhedron exact');
% ReLU.plot_set(R_se, 'Star exact');
% ReLU.plot_set(R_sa, 'Star approx');
% ReLU.plot_set(R_z, 'Zonotope approx');
% %ReLU.plot_set(R_sza, 'Star zonotope approx');
% ReLU.plot_set(R_sad, 'Star abs-domain approx');
% ReLU.plot_set(R_slog, 'Star logsig approx');
% ReLU.plot_set(R_stan, 'Star tansig approx');
