close all;
clear;
clc;

P = ExamplePoly.randHrep('d',2);
I = Star(P);
W = [];
b = [];
T = Sigmoid.reach_approx_star(I, W, b, 'tansig');
L = Sigmoid.reach_approx_star(I, W, b, 'logsig');

figure;
ReLU.plot_set_hold(T,'','r');
ReLU.plot_set_hold(L,'','y');