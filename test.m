clear all;
close all;
clc;

I = ExamplePoly.randHrep('d',2);
S = myStar(I);
plot(S)
R = myReLU.reach_exact(S);
plot(I);
hold on;
plot(R);