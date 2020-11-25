close all;
clear; 
clc; 

% P = ExamplePoly.randHrep('d',2);
P = Polyhedron('lb',[-1;-1],'ub',[1;1]);
S = Star(P);
W = [1 -1; 1 0];
b = [1; 0];
S = S.affineMap(W,b);
ReLU.plot_set(S);