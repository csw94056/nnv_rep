close all;
clear;
clc;

% steps = 3;
% for i = 1:steps
%     W{i} = rand(2,2);
%     b{i} = rand(2,1);

%     W{i} = randi([-5 5], 2,2);
%     b{i} = randi([-5 5], 2,1);
    
%     W{i} = -1 + 2*rand(2,2);
%     b{i} = -1 + 2*rand(2,1);
%     
%     W{i} = -10 + 20*rand(2,2);
%     b{i} = -10 + 20*rand(2,1);
%     
%     W{i} = randi([0 5], 2,2);
%     b{i} = randi([0 5], 2,1);
% end

% W{1} = [1 1; 1 -1];
% W{2} = [1 1; 1 -1];
% W{3} = [1 1; 1 0];
% 
% b{1} = [0; 0];
% b{2} = [0; 0];
% b{3} = [1; 0];

% W{1} = [0.299868459556700 0.382500389663171; 0.768333245781196 0.682814182234229];
% W{2} = [0.510477531753991 0.266520298232427; 0.726956495360460 0.374197862977184];
% W{3} = [0.722545463325952 0.639996828221631; 0.168575824568629 0.900263787130386];
% W{4} = [0.629738314417465 0.470847574408952; 0.358940343531266 0.162864077745110];
% 
% b{1} = [0.986340856903685; 0.765850765225542];
% b{2} = [0.332840367586042; 0.399614229990180];
% b{3} = [0.661881133589573; 0.656915139738211];
% b{4} = [0.533401905477953; 0.447770431174552];

%Ex 1 AbsS is conservative than Se
W{1} = [2 -2; 3 -3];
W{2} = [5 -2; 3 1];
% W{3} = [4 -5; 1 -3]; 
b{1} = [-2; -2]; 
b{2} = [3;  4];
% b{3} = [4; -1]
steps = length(W);




P = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);
S = Star(P);
R = RlxPoly(P, inf);        % RlxPoly with original constraints
RlxS = RlxStar(P, inf);
Ru = RlxPoly(P, inf);       % RlxPoly with only u >= -1 constraints
AbsS = AbsStar(P, inf);
Se = Star(P);

figure;
nexttile;
plot(R, 'r');
hold on;
plot(S, 'y');    
hold on;
plot(RlxS, 'g');
hold on;
plot(AbsS, 'g');
hold on;
plot(Se, 'c');
title('input');
legend('Sung Abs-Dom', 'NNV Abs-Dom','Sung Star with Abs-Dom bounds', 'NNV approx');

for i = 1:steps
    R = ReLU.layerReach(R, W{i}, b{i}, 'approx');
    S = ReLU.layerReach(S, W{i}, b{i}, 'abs_domain');
    RlxS = ReLU.layerReach(RlxS, W{i}, b{i}, 'approx');
    AbsS = ReLU.layerReach(AbsS, W{i}, b{i}, 'approx');
    Se = ReLU.layerReach(Se, W{i}, b{i}, 'approx');
    
    nexttile;
    plot(R, 'r');
    hold on;
    plot(S, 'y');
    hold on;
    plot(RlxS, 'g');
    plot(AbsS, 'g');
    hold on;
    plot(Se, 'c');
    if i < steps
        str = sprintf('Hidden Layer Rechable set (%d)', i);
    else 
        str = sprintf('Output Layer Rechable set');
    end
    title(str);
end

figure('Name', 'Output Layer Reachable Set');
% nexttile;
% plot(R, 'r');
% title('Sung Abs-Dom');
% nexttile;
% plot(S, 'y');
% title('NNV Abs-Dom');
% nexttile;
% plot(RlxS, 'g');
% title('Sung Star with Abs-Dom bounds');
nexttile
plot(AbsS, 'g');
title('AbsS');
nexttile;
plot(Se, 'c');
title('NNV approx');