%{
    Sung Woo Choi
    09142020
    Reachability Anylsis computation time of Star vs Polyhedron
%}
close all;
clear; 
clc;

%% computing reachable set of a relu network
dim = 2;
n_dim = 20;
out_dim = 2;

% random n-dimensional polyhedron
I = ExamplePoly.randHrep('d',dim); % input set
S = myStar(I);

R_pe = myReLU.reach_exact(I);
R_se = myReLU.reach_exact(S);
R_pa = myReLU.reach_approx(I);
R_sa = myReLU.reach_approx(S);


%% Reachability analysis of each method

figure('Name','ReLU computation');
myReLU.plot_set(I, 'I poly');
myReLU.plot_set(S, 'I star');
myReLU.plot_set(R_pe, 'R poly exact');
myReLU.plot_set(R_se, 'R star exact');
myReLU.plot_set(R_pa, 'R poly approx');
myReLU.plot_set(R_sa, 'R star approx');


%% FNN network of [2 n 2]
figure('Name','Computation of reachable set of a FNN relu network');
for i=dim:n_dim
    disp('------Computing reachable set of a FNN relu network------')   
    fprintf('network dimensions: [%d %d %d]\n',dim, i, out_dim');

    nN = [dim i out_dim]; % 2 inpus, n hidden layers, 2 outputs
    [W, b] = myReLU.rand_Layers(nN);

    disp('-Exact reachability analysis for Polyhedron')
    [R_pe, ct_pe] = myReLU.relu_FNN_reach_BFS(I, W, b, 'exact');
    fprintf('computation time: %f sec\n', ct_pe);

    disp('-Exact reachability analysis for Star Set')
    [R_se, ct_se] = myReLU.relu_FNN_reach_BFS(S, W, b, 'exact');
    fprintf('computation time: %f sec\n', ct_se);
    
    disp('-Approximate reachability analysis for Polyhedron')
    [R_pa, ct_pa] = myReLU.relu_FNN_reach_BFS(I, W, b, 'approx');
    fprintf('computation time: %f sec\n', ct_pa);

    disp('-Approximate reachability analysis for Star Set')
    [R_sa, ct_sa] = myReLU.relu_FNN_reach_BFS(S, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_sa);

    hold on
    plot(i,ct_pe,'x','Linewidth', 2, 'color', 'blue');
    plot(i,ct_se,'o','Linewidth', 2, 'color', 'red');
    plot(i,ct_pa,'+','Linewidth', 2, 'color', 'green');
    plot(i,ct_sa,'^','Linewidth', 2, 'color', 'magenta');
    xlim([dim n_dim]);
    hold off
    legend('exact polyhedron', 'exact star set', 'approx polyhedron', 'approx star set')
end
title('the verification time vs. the dimension of a hidden layer');
ylabel('verification time (sec)');
xlabel('dimension of a hidden layer');

%% FNN network of [2 n n n 2]
figure
for i=dim:n_dim
    disp('------Computing reachable set of a FNN relu network------')
    fprintf('network dimensions: [%d %d %d %d %d]\n',dim, i, i, i, out_dim');
    
    nN = [dim i i i out_dim]; % 2 inpus, n hidden layers, 2 outputs
    [W, b] = myReLU.rand_Layers(nN);

    disp('-Exact reachability analysis for Polyhedron')
    [R_pe, ct_pe] = myReLU.relu_FNN_reach_BFS(I, W, b, 'exact');
    fprintf('computation time: %f sec\n', ct_pe);

    disp('-Exact reachability analysis for Star Set')
    [R_se, ct_se] = myReLU.relu_FNN_reach_BFS(S, W, b, 'exact');
    fprintf('computation time: %f sec\n', ct_se);
    
    disp('-Approximate reachability analysis for Polyhedron')
    [R_pa, ct_pa] = myReLU.relu_FNN_reach_BFS(I, W, b, 'approx');
    fprintf('computation time: %f sec\n', ct_pa);

    disp('-Approximate reachability analysis for Star Set')
    [R_sa, ct_sa] = myReLU.relu_FNN_reach_BFS(S, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_sa);

    hold on
    plot(i,ct_pe,'x','Linewidth', 2, 'color', 'blue');
    plot(i,ct_se,'o','Linewidth', 2, 'color', 'red');
    plot(i,ct_pa,'+','Linewidth', 2, 'color', 'green');
    plot(i,ct_sa,'^','Linewidth', 2, 'color', 'magenta');
    xlim([dim n_dim]);
    hold off
    legend('exact polyhedron', 'exact star set', 'approx polyhedron', 'approx star set')
end
title('the verification time vs. the dimension of hidden layers in [2 n n n 2] network');
ylabel('verification time (sec)');
xlabel('dimension of hidden layers');