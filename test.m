clear all;
close all;
clc;

I = ExamplePoly.randHrep('d',2);

S = myStar(I);
plot(S)
R1 = myReLU.reach_exact(S);
R2 = myReLU.reach_exact(I);

plot(I);
hold on;
plot(R1);