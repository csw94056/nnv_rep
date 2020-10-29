 %abstract domain
close all;
clear;
clc;

W{1} = [1 1; 1 -1];
W{2} = [1 1; 1 -1];
W{3} = [1 1; 1 0];
W{4} = [1 1; 1 -1];
W{5} = [1 1; 1 0];

% b1 = [0; 0];
% b2 = [0; 0];
% b3 = [0; 0];
% b4 = [0; 0];
% b5 = [1; 0];

% W1 = [3 1; 1 -5];
% W2 = [1 4; 1 -1];
% W3 = [1 1; 0 1];
% W4 = [1 1; 1 -1];
% W5 = [1 1; 0 1];
% 
% b1 = [0; 2];
% b2 = [0; 0];
% b3 = [1; 0];
% b4 = [1; 0];
% b5 = [1; 0];

steps = 5;
for i = 1:5;
%     W{i} = rand(2,2);
%     b{i} = rand(2,1);
    
%     W{i} = randi([-5 5], 2,2);
    b{i} = randi([-5 5], 2,1);
    
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




for i = 1:4
    S = S.affineMap(W{i}, b{i});
    R = R.affineMap(W{i}, b{i});
    
    nexttile;
    plot(R, 'r');
    hold on;
    plot(S, 'y');
    title('affine map');

    S = ReLU.reach_approx(S, 'abs_domain');
    R = ReLU.reach_approx(R);
    
    nexttile;
    plot(R, 'r');
    hold on;
    plot(S, 'y');
    title('ReLU');
end
