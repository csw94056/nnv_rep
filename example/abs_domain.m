 %abstract domain
close all;
clear;
clc;

%P = ExamplePoly.randHrep('d', 2);
%D = DeepPoly(P);
%R = ReLU.approxStepReLU_deeppoly(D);

%{
W{1} = [1 1; 1 -1];
W{2} = [1 1; 1 -1];
W{3} = [1 1; 0 1];
%}

%P = Polyhedron('ub', [1;1], 'lb', [-1;-1]);
%P1 = P.affineMap(W{1});

W1 = [1 1; 1 -1];
W2 = [1 1; 1 -1];
W3 = [1 1; 0 1];

b1 = [0; 0];
b2 = [0; 0];
b3 = [1; 0];


lower_a = [0 -1 0; 0  0 -1];
upper_a = [0 1 0; 0 0 1];
lb = [-1;-1];
ub = [1;1];
D = RelaxedPoly(lower_a, upper_a, lb, ub);

D_lower_a = D.lower_a
D_upper_a = D.upper_a
D_lb = D.lb
D_ub = D.ub


D1 = D.affineMap2(W1, b1);

D1_lower_a = D1.lower_a
D1_upper_a = D1.upper_a
D1_lb = D1.lb
D1_ub = D1.ub

D2 = D1;
for i = 1:D2.Dim
    D2 = ReLU.approxSingleStepReLU_rlxpoly(D2,i);
end

D2_lower_a = D2.lower_a 
D2_upper_a = D2.upper_a
D2_lb = D2.lb
D2_ub = D2.ub

D3 = D2.affineMap2(W2, b2);

D3_lower_a = D3.lower_a 
D3_upper_a = D3.upper_a
D3_lb = D3.lb

% D3.ub(end - 1) = 3;
D3_ub = D3.ub

D4 = D3;
for i = 1:D4.Dim
    D4 = ReLU.approxSingleStepReLU_rlxpoly(D4, i);
end

D4_lower_a = D4.lower_a 
D4_upper_a = D4.upper_a
D4_lb = D4.lb
D4_ub = D4.ub


D5 = D4.affineMap2(W3, b3);

D5_lower_a = D5.lower_a 
D5_upper_a = D5.upper_a
D5_lb = D5.lb
D5_ub = D5.ub

%{
x(1) = DeepPoly(-1, 1, -1, 1);
x(2) = DeepPoly(-1, 1, -1, 1);

x(3) = DeepPoly([x(1).lower_a x(2).lower_a], [x(1).upper_a x(2).upper_a], x(1).lb + x(2).lb, x(1).ub + x(2).ub);
x(4) = DeepPoly([x(1).lower_a -x(2).lower_a], [x(1).upper_a - x(2).upper_a], x(1).lb + x(2).lb, x(1).ub + x(2).ub);


x(3)
x(4)
%}
%{
P = ExamplePoly.randHrep('d',2);
S = Star(P);

R = ReLU.absDomainReachReLU_star(S);
Int = R.Intersect(S);

figure;
plot(S, 'r');
hold on;
plot(R, 'y');
hold on;
plot(Int, 'g');
%}