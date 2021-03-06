close all;
clear;
clc;


P = ExamplePoly.randHrep('d',2);
Z = Zonotope(P);

W = rand(2, 2);
b = rand(3, 1);

Z1 = Z.affineMap(W);


%{
Rsig = Sigmoid.stepSigmoid_zono(Z, 1, 'logsig');
Rrelu = ReLU.approxStepReLU_zono(Z, 1);

figure;
nexttile
plot(Z, 'r');
hold on;
plot(Rsig, 'y');
hold on;
plot(Rrelu, 'b');
nexttile
plot(Rsig, 'y');
nexttile
plot(Rrelu, 'b');
%}




