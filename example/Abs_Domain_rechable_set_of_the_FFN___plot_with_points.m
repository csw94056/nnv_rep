close all;
clear;
clc;

% steps = 3;
% for i = 1:steps
%     W{i} = rand(2,2);
%     b{i} = rand(2,1);
% 
%     W{i} = randi([-5 5], 2,2);
%     b{i} = randi([-5 5], 2,1);
%     
%     W{i} = -1 + 2*rand(2,2);
%     b{i} = -1 + 2*rand(2,1);
%     
%     W{i} = -10 + 20*rand(2,2);
%     b{i} = -10 + 20*rand(2,1);
%     
%     W{i} = randi([0 5], 2,2);
%     b{i} = randi([0 5], 2,1);
% end

W{1} = [2 0; 1 -1; 1 1];
W{2} = [-1 0 1; 1 -1 0];
b{1} = [0.5; -1; -0.5];
b{2} = [-0.5; 0.5];

% W{1} = [1 1; 1 -1];
% W{2} = [1 1; 1 -1];
% W{3} = [1 1; 1 0];
% b{1} = [0; 0];
% b{2} = [0; 0];
% b{3} = [1; 0];

% W{1} = [1 1; 1 -1];
% W{2} = [1 0; 2 1];
% b{1} = [0; 0];
% b{2} = [0; 1.5];

lb = [-1; -1];
ub = [1; 1];
P = Polyhedron('lb', lb, 'ub', ub);
p{1} = (lb + ub) * 0.5;
p{2} = lb;
p{3} = ub;
p{4} = [ub(1); lb(1)];
p{5} = [lb(1); ub(2)];
RlxS = RlxStar(P, inf);

figure('Name','Sung Star with Abs-Dom bounds');
nexttile;
plot(RlxS, 'c');
hold on;
for i = 1:length(p)
    plot(p{i}(1),p{i}(2), '*', 'color', 'black');
end
hold off;
title('Input Set');


steps = length(W);
for i = 1:steps
    nexttile;
    RlxS = ReLU.layerReach(RlxS, W{i}, b{i}, 'approx');
    plot(RlxS, 'c');
    
    hold on;
    for j = 1:length(p)
        p{j} = point_layerReach(p{j}, W{i}, b{i});
        if length(p{j}) == 3
            plot3(p{j}(1),p{j}(2),p{j}(3), '*', 'color', 'black');
            xlabel('x1');
            ylabel('x2');
            zlabel('x3');
        else
            plot(p{j}(1),p{j}(2), '*', 'color', 'black');
            xlabel('x1');
            ylabel('x2');
        end
    end
    hold off;
    if i == 1
        title('Hidden Layer Reachable Set');
    elseif i == 2
        title('Output Layer Reachable Set');
    end
end


function r = point_layerReach(p, W, b)
    r = max(0,W*p + b);
end