close all;
clear;
clc;

% W{1} = [1 1; 1 -1];
% W{2} = [1 1; 1 -1];
% W{3} = [1 1; 1 0];
% 
% b{1} = [0; 0];
% b{2} = [0; 0];
% b{3} = [1; 0];

steps = 5;
for i = 1:steps
%     W{i} = rand(2,2);
%     b{i} = rand(2,1);

%     W{i} = randi([-5 5], 2,2);
%     b{i} = randi([-5 5], 2,1);
    
    W{i} = randi([0 5], 2,2);
    b{i} = randi([0 5], 2,1);
end
% 
% W{1} = rand(3,2);
% b{1} = rand(3,1);
% W{2} = rand(3,3);
% b{2} = rand(3,1);
% W{3} = rand(2,3);
% b{3} = rand(2,1);

%good example
% W{1} = [-5 -5; -2 -2];
% b{1} = [2; -3];
% W{2} = [2 5; 3 4];
% b{2} = [-2; 4];
% W{3} = [4 -2; 4 3];
% b{3} = [-2; 4];
% W{4} = [1 -5; -3 -1];
% b{4} = [0; -5];
% W{5} = [4 1; 3 5];
% b{5} = [2; 0];


P = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);
S = Star(P);
R = RlxPoly(P, inf);
RlxS= RlxStar(P, inf);
Se = Star(P);

nexttile;
plot(R, 'r');
hold on;
plot(S, 'y');    
hold on;
plot(RlxS, 'g');
hold on;
plot(Se, 'c');
title('input');
legend('RlxPoly', 'Abs-dom Star','RlxStar','Exact Approx Star');

for i = 1:steps
    S = S.affineMap(W{i}, b{i});
    R = R.affineMap(W{i}, b{i});
    RlxS = RlxS.affineMap(W{i}, b{i});
    Se = Se.affineMap(W{i}, b{i});
    nexttile;
    plot(R, 'r');
    hold on;
    plot(S, 'y');    
    hold on;
    plot(RlxS, 'g');
    hold on;
    plot(Se, 'c');
    title('affine map');

    S = ReLU.reach_approx(S, 'abs_domain');
    R = ReLU.reach_approx(R);
    RlxS = ReLU.reach_approx(RlxS);
    Se = ReLU.reach_approx(Se, 'approx');
    nexttile;
    plot(R, 'r');
    hold on;
    plot(S, 'y');
    hold on;
    plot(RlxS, 'g');
    hold on;
    plot(Se, 'c');
    title('ReLU');
end