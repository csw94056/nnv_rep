close all;
clear;
clc;

W{1} = [1 -1; 1 -1];
W{2} = [1 1; 1 -1];
W{3} = [1 1; 1 -1];
b{1} = [0; 0]; 
b{2} = [0; 0];
b{3} = [0; 0];


P = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);
S = Star(P);

for i = 1:steps
    S = S.affineMap(W{i}, b{i});
    nexttile;
    plot(S, 'y');
end