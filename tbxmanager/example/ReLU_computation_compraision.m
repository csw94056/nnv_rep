close all;
clear
clc;

% input set
% random n-dimensional polyhedron
random = ExamplePoly.randHrep('d',2); 
rectangle = Polyhedron([1,1;0,0;0,2;2,0;2,2;]);

W2_1 = [2 0; 1 -1; 1 1];
b2 = [0.5;-1;-0.5];
W3_2 = [-1 0 1; 1 -1 0];
b3 = [-0.5;0.5];
%%
I = rectangle;
S = Star(I);
figure('Name','over-approximate analysis');
ReLU.plot_set(S, 'I star');
R = ReLU.layerReach(I, W2_1, b2, 'approx');
ReLU.plot_set(R,'R poly after W2_1');
R = ReLU.layerReach(S, W2_1, b2, 'approx');
ReLU.plot_set(R,'R star after W2_1');
%%
I = random;
S = Star(I);
W = rand(I.Dim, I.Dim);
b = rand(I.Dim,1);

figure('Name','over-approximate analysis 2');
ReLU.plot_set(I, 'I poly');
ReLU.plot_set(S, 'I star');

ReLU.plot_set(I.affineMap(W)+b, 'I affined poly');
ReLU.plot_set(S.affineMap(W,b), 'I affined star');

R = ReLU.layerReach(I, W, b, 'approx');
ReLU.plot_set(R,'R poly approx');
R = ReLU.layerReach(S, W, b, 'approx');
ReLU.plot_set(R,'R star approx');

%%
R_pe = ReLU.reach_exact(I);
R_pa = ReLU.reach_approx(I);
R_se = ReLU.reach_exact(S);
R_sa = ReLU.reach_approx(S);
figure('Name','ReLU computation');
ReLU.plot_set(I, 'I poly');
ReLU.plot_set(S, 'I star');

ReLU.plot_set(R_pe, 'R poly exact');
ReLU.plot_set(R_se, 'R star exact');
ReLU.plot_set(R_pa, 'R poly approx');
ReLU.plot_set(R_sa, 'R star approx');

figure('Name','ReLU computation');
ReLU.plot_set(I, 'I poly');
ReLU.plot_set(S, 'I star');

R = ReLU.layerReach(I, W2_1, b2, 'exact');
ReLU.plot_set(R,'R poly after W2_1 exact');
R = ReLU.layerReach(S, W2_1, b2, 'exact');
ReLU.plot_set(R,'R star after W2_1 exact');

R = ReLU.layerReach(I, W2_1, b2, 'approx');
ReLU.plot_set(R,'R poly after W2_1 approx');
R = ReLU.layerReach(S, W2_1, b2, 'approx');
ReLU.plot_set(R,'R star after W2_1 approx');

for i=1:4
    div = 2^i;
    tic
    R = ReLU.set_div(I, div, div);
    time = toc;
    ReLU.plot_set(R,'I poly div','rand');
    R = ReLU.layerReach(R, W2_1, b2, 'approx');
    str = sprintf('R poly div after W2_1 approx (%f sec)', time);
    ReLU.plot_set(R,str, 'rand');
    tic
    R = ReLU.set_div(S, div, div);
    ReLU.plot_set(R, 'I star div', 'rand');
    time = toc;
    R = ReLU.layerReach(R, W2_1, b2, 'approx');
    str = sprintf('R star div after W2_1 approx (%f sec)', time);
    ReLU.plot_set(R, str, 'rand');
end