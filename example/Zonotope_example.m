close all;
clear; 
clc; 

% P = ExamplePoly.randHrep('d',2);
P = Polyhedron('lb',[-1;-1],'ub',[1;1]);
I = Zonotope(P);

nN = [2 2 2 2];
% [W, b] =  ReLU.rand_Layers(nN);

W{1} = [1 1; 1 -1];
W{2} = [1 1; 1 -1];
W{3} = [1 1; 1 0];

b{1} = [0; 0];
b{2} = [0; 0];
b{3} = [1; 0];


R_z = ReLU.FNN_reach_BFS(I, W, b, 'approx');
% % method = 'approx';
% 
% 
% R = [];
% I1 = I;
% for i=1:length(W)
%     I1 = I1.affineMap(W{i}, b{i});
%     R = [R ReLU.reach_approx(I1, 'approx')];
% end
% R_z = R;                
                    
figure('Name','All methods');
% ReLU.plot_set_hold(I, 'Input set','r');
ReLU.plot_set_hold(R_z, 'Zonotope approx', 'c');
