close all;
clear;
clc;
% 
% W{1} = [1 1; 1 -1];
% % W{2} = [1 1; 1 -1];
% % W{3} = [1 1; 1 0];
% % 
% b{1} = [0; 0];
% b{2} = [0; 1.5];
% b{3} = [1; 0];

for i = 1:5
%     W{i} = rand(2,2);
%     b{i} = rand(2,1);

%     W{i} = -1 + 2*rand(2,2);
%     b{i} = -1 + 2*rand(2,1);
    
%     W{i} = -10 + 20*rand(2,2);
%     b{i} = -10 + 20*rand(2,1);
    
%     W{i} = randi([0 5], 2,2);
%     b{i} = randi([0 5], 2,1);
end

% hidden layers have 3 nodes
% W{1} = rand(3,2);
% b{1} = rand(3,1);
% W{2} = rand(3,3);
% b{2} = rand(3,1);
% W{3} = rand(2,3);
% b{3} = rand(2,1);

% good example
% W{1} = [-5 -5; -2 -2];
% b{1} = [2; -3];
% W{2} = [2 5; 3 4];
% b{2} = [-2; 4];
% W{3} = [4 -2; 4 3];
% b{3} = [-2; 4];
% W{4} = [1 -5; -3 -1];
% b{4} = [0; -5];
% W{5} = [4 1; 3 5];
% b{5} = [2; 0];

% good example 2
% W{1} = [0.728214624538208 0.630757626740935; 0.542567725535353 0.269076315415997];
% W{2} = [0.925544898560559 0.066704479259003; 0.951119572011534 0.946137972467446];
% W{3} = [0.153641709762110 0.168748938145153; 0.296073114887304 0.735194614493405];
% b{1} = [0.591812006788566; 0.866305096716458];
% b{2} = [0.496633305608366; 0.068027386603233];
% b{3} = [0.879697065220770; 0.851423388162803];

% %good example 3
% W{1} = [0.465698898461883 0.343303486341360; 0.211809493190058 0.763338677110964];
% W{2} = [0.877803795161917 0.187516474910597; 0.151308867486537 0.949022459153999];
% W{3} = [0.424614495494384 0.852506608051607; 0.704501247596211 0.706492692448605];
% b{1} = [0.128698784708328; 0.015502586140148];
% b{2} = [0.948045881445746; 0.558770573875961];
% b{3} = [0.178788094416603; 0.305273656760658];

%example why AbsS is more conservative than exact-approx Star
% W{1} = [5 -1; 3 2];
% W{2} = [3 -2; 4 0];
% W{3} = [-1 1; 1 0];
% W{4} = [-3 -5; 4 -5];
% W{5} = [1 5; -3 5];
% b{1} = [3; -2];
% b{2} = [2; 2];
% b{3} = [0; 3];
% b{4} = [0; -5];
% b{5} = [1; 5];
%example2 why AbsS is more conservative than exact-approx Star
% W{1} = [-6.9348 9.5514; -2.4429 9.8134];
% W{2} = [2.5300 -2.6636; -6.8864 1.3701];
% W{3} = [-4.9032 7.8001; 8.3014 6.9620];
% W{4} = [-9.5327 -7.1667; 0.7005 8.5443];
% W{5} = [8.7960 -7.5174; 2.3997 -4.8236];
% b{1} = [-5.3579; -5.6636];
% b{2} = [-6.1321; 6.1095];
% b{3} = [-3.8838; -4.6882];
% b{4} = [-3.4347; -8.7425];
% b{5} = [1.5628; 9.7556];

%example3 why AbsS is more conservative than exact-approx Star
W{1} = [-6 9; -2 9];
W{2} = [2 -2; -6 1];
W{3} = [-4 7; 8 6];
W{4} = [-9 -7; 0 8];
W{5} = [8 -7; 2 -4];
b{1} = [-5; -5];
b{2} = [-6; 6];
b{3} = [-3; -4];
b{4} = [-3; -8];
b{5} = [1; 9];


% % W{1} = [1 -1; 1 -1];
% % W{2} = [1 1; 1 -1];
% W{3} = [1 1; 1 1];
% W{3} = [-1 1; 1 1];
%Ex 2 RlxS is more conservative than AbsS is conservative than Se
% W{1} = [1 -1; 1 -1];
% W{2} = [-1 1; -1 1];
% example AbsS more conservative
% W{1} = [1 -1; 1 -1];
% W{2} = [-1 1; -1 1];
% W{3} = [1 1; 1 -1];
% original example
% W{1} = [1 1; 1 -1];
% W{2} = [1 1; 1 -1];
% W{3} = [1 1; 0 1];
% b{1} = [0; 0]; 
% b{2} = [0; 0];
% b{3} = [0; 0];
% b{3} = [4; -1]
steps = length(W);


% P = Polyhedron('lb', [-1, 0], 'ub', [1, 0]);
P = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);
% P = Polyhedron('lb', [0; 0], 'ub', [1; 1]);
% P = Polyhedron('lb', [0; 0], 'ub', [1; 1]);
S = Star(P);
R = RlxPoly(P, inf);        % RlxPoly with original constraints
Ru = RlxPoly(P, inf);       % RlxPoly with only u >= -1 constraints
RlxS = RlxStar(P, inf);
Se = Star(P);
AbsS = AbsStar(P, inf);

figure;
nexttile;
% plot(R, 'r');
% hold on;
plot(Ru, 'b');
hold on;
% plot(S, 'y');    
% hold on;
plot(RlxS, 'g');
hold on;
plot(AbsS, 'm');
hold on;
plot(Se, 'c');
title('input');
% legend('RlxPoly','Abs-dom Star','upperRlxStar','Exact-approx Star');
% legend('RlxPoly', 'Abs-dom Star','RlxStar','Exact-approx Star');
legend('Sung Abs-Dom upper','Sung Abs-Dom Star', 'Sung Abs-Dom Star upper with 3 const','NNV Approx Star');

for i = 1:steps
    S = S.affineMap(W{i}, b{i});
    R = R.affineMap(W{i}, b{i});
    Ru = Ru.affineMap(W{i}, b{i});
    RlxS = RlxS.affineMap(W{i}, b{i});
    Se = Se.affineMap(W{i}, b{i});
    AbsS = AbsS.affineMap(W{i}, b{i});
    nexttile;
%     plot(R, 'r');
%     hold on;
    plot(Ru, 'b');
    hold on;
%     plot(S, 'y');    
%     hold on;
    plot(RlxS, 'g');
    hold on;
    plot(AbsS, 'm');
    hold on;
    plot(Se, 'c');
    
    
    title('affine map');
    R = ReLU.reach_approx(R, 'approx');
    Ru = ReLU.reach_approx(Ru, 'upper');
    RlxS = ReLU.reach_approx(RlxS);
    Se = ReLU.reach_approx(Se, 'approx');
    AbsS = ReLU.reach_approx(AbsS);
    
    nexttile;
%     plot(R, 'r');
%     hold on;

    S = ReLU.reach_approx(S, 'abs_domain');
    plot(Ru, 'b');
    hold on;
%     plot(S, 'y');
%     hold on;
    plot(RlxS, 'g');
    hold on;
    plot(AbsS, 'm');
    hold on;
    plot(Se, 'c');
    
    
    title('ReLU');
end

figure;
% nexttile;
% plot(R, 'r');
% title('RlxStar');
nexttile;
plot(Ru, 'b');
title('Sung Abs-Dom upper');
% nexttile;
% plot(S, 'y');
% title('abs-dom Star');
nexttile;
plot(RlxS, 'g');
title('Sung Abs-Dom Star');
nexttile;
plot(Se, 'c');
title('exact-approx Star');
nexttile;
plot(AbsS, 'm');
title('Sung Abs-Dom Star upper with 3 const');