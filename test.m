clear all;
close all;
clc;

I = ExamplePoly.randHrep('d',2);
S = myStar(I);

R = myReLU.reach_exact(I);
plot(I, 'r');
hold on;
polt(R, 'g');