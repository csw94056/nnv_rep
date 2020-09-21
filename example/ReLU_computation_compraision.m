clear all;
close all;
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
S = myStar(I);
figure('Name','over-approximate analysis');
myReLU.plot_set(S, 'I star');
R = myReLU.layerReach(I, W2_1, b2, 'approx');
myReLU.plot_set(R,'R poly after W2_1');
R = myReLU.layerReach(S, W2_1, b2, 'approx');
myReLU.plot_set(R,'R star after W2_1');
%%
I = random;
S = myStar(I);
W = rand(I.Dim, I.Dim);
b = rand(I.Dim,1);

figure('Name','over-approximate analysis 2');
myReLU.plot_set(I, 'I poly');
myReLU.plot_set(S, 'I star');

myReLU.plot_set(I.affineMap(W)+b, 'I affined poly');
myReLU.plot_set(S.affineMap(W,b), 'I affined star');

R = myReLU.layerReach(I, W, b, 'approx');
myReLU.plot_set(R,'R poly approx');
R = myReLU.layerReach(S, W, b, 'approx');
myReLU.plot_set(R,'R star approx');

%%
R_pe = myReLU.reach_exact(I);
R_pa = myReLU.reach_approx(I);
R_se = myReLU.reach_exact(S);
R_sa = myReLU.reach_approx(S);
figure('Name','ReLU computation');
myReLU.plot_set(I, 'I poly');
myReLU.plot_set(S, 'I star');

myReLU.plot_set(R_pe, 'R poly exact');
myReLU.plot_set(R_se, 'R star exact');
myReLU.plot_set(R_pa, 'R poly approx');
myReLU.plot_set(R_sa, 'R star approx');

figure('Name','ReLU computation');
myReLU.plot_set(I, 'I poly');
myReLU.plot_set(S, 'I star');

R = myReLU.layerReach(I, W2_1, b2, 'exact');
myReLU.plot_set(R,'R poly after W2_1 exact');
R = myReLU.layerReach(S, W2_1, b2, 'exact');
myReLU.plot_set(R,'R star after W2_1 exact');

R = myReLU.layerReach(I, W2_1, b2, 'approx');
myReLU.plot_set(R,'R poly after W2_1 approx');
R = myReLU.layerReach(S, W2_1, b2, 'approx');
myReLU.plot_set(R,'R star after W2_1 approx');

for i=1:4
    div = 2^i;
    tic
    R = set_div(I, div, div);
    time = toc;
    plot_set(R,'I poly div');
    R = layerReach(R, W2_1, b2, 'approx');
    str = sprintf('R poly div after W2_1 approx (%f sec)', time);
    plot_set(R,str);
    tic
    R = set_div(S, div, div);
    plot_set(R, 'I star div');
    time = toc;
    R = layerReach(R, W2_1, b2, 'approx');
    str = sprintf('R star div after W2_1 approx (%f sec)', time);
    plot_set(R, str);
end