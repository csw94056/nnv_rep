 %abstract domain
close all;
clear;
clc;

W{1} = [1 1; 1 -1];
W{2} = [1 1; 1 -1];
W{3} = [1 1; 1 0];

b{1} = [0; 0];
b{2} = [0; 0];
b{3} = [1; 0];

% W{1} = [1 1; 1 -1];
% W{2} = [1 1; 1 -1];
% W{3} = [1 1; 1 0];
% W{4} = [1 1; 1 -1];
% W{5} = [1 1; 1 0];

% b{1} = [2; 5];
% b{2} = [-5; 1];
% b{3} = [0; -2];
% b{4} = [1; 5];
% b{5} = [3; 3];

% b{1} = [0; 1];
% b{2} = [0; 0];
% b{3} = [0; 0];
% b{4} = [0; 0];
% b{5} = [1; 0];


steps = 3;
for i = 1:steps
%     W{i} = rand(2,2);
%     b{i} = rand(2,1);
   
%     W{i} = randi([-5 5], 2,2);
%     b{i} = randi([-5 5], 2,1);
    
%     W{i} = randi([0 5], 2,2);
%     b{i} = randi([0 5], 2,1);
end

% W{1} = rand(3,2);
% b{1} = rand(3,1);
% W{2} = rand(3,3);
% b{2} = rand(3,1);
% W{3} = rand(2,3);
% b{3} = rand(2,1);


P = Polyhedron('lb', [-1;-1], 'ub', [1;1]);
S = Star(P);
R = RelaxedPoly(P);


nexttile;
plot(S, 'y');
hold on;
plot(R, 'r');
title('input');


for i = 1:steps
    S = S.affineMap(W{i}, b{i});
    R = R.affineMap(W{i}, b{i});
    
    nexttile;
    plot(R, 'r');
    hold on;
    plot(S, 'y');
%     hold on;
%     plot(R, 'g');
    title('affine map');

    S = ReLU.reach_approx(S, 'abs_domain');
    R = ReLU.reach_approx(R);
    
    nexttile;
    plot(R, 'r');
    hold on;
    plot(S, 'y');
%     hold on;
%     plot(R, 'g');
    title('ReLU');
end