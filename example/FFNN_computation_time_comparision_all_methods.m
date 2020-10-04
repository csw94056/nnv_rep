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
n_dim = 15;
out_dim = 2;

% random n-dimensional polyhedron
I = ExamplePoly.randHrep('d',dim); % input set
S = Star(I);
Z = Zonotope(I);

%% FNN network of [2 n 2]
figure('Name','Computation of reachablility of a FNN network');
for i=dim:n_dim
    disp('------Computing reachable set of a FNN relu network------')   
    fprintf('network dimensions: [%d %d %d]\n',dim, i, out_dim');

    nN = [dim i out_dim]; % 2 inpus, n hidden layers, 2 outputs
    [W, b] = ReLU.rand_Layers(nN);

    disp('-Exact reachability analysis for Polyhedron')
    [R_pe, ct_pe] = ReLU.FNN_reach_BFS(I, W, b, 'exact');
    fprintf('computation time: %f sec\n', ct_pe);
    
    disp('-Approximate reachability analysis for Polyhedron')
    [R_pa, ct_pa] = ReLU.FNN_reach_BFS(I, W, b, 'approx');
    fprintf('computation time: %f sec\n', ct_pa);

    disp('-Exact reachability analysis for Star Set')
    [R_se, ct_se] = ReLU.FNN_reach_BFS(S, W, b, 'exact');
    fprintf('computation time: %f sec\n', ct_se);

    disp('-Approximate reachability analysis for Star Set')
    [R_sa, ct_sa] = ReLU.FNN_reach_BFS(S, W, b, 'approx');
    fprintf('computation time: %f sec\n', ct_sa);
    
    disp('-Abstract domain Approximate reachability analysis for Star Set')
    [R_sabs, ct_sabs] = ReLU.FNN_reach_BFS(S, W, b, 'abs_domain');
    fprintf('computation time: %f sec\n', ct_sabs);
    
    disp('-Approximate reachability analysis for Zonotope Set')
    [R_za, ct_za] = ReLU.FNN_reach_BFS(Z, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_za);

    hold on
    plot(i,ct_pe,'x','Linewidth', 1, 'color', 'blue');
    plot(i,ct_se,'x','Linewidth', 1, 'color', 'red');
    plot(i,ct_pa,'*','Linewidth', 1, 'color', 'green');
    plot(i,ct_sa,'*','Linewidth', 1, 'color', 'magenta');
    plot(i,ct_sabs,'*','Linewidth', 1, 'color', 'cyan');
    plot(i,ct_za,'square','Linewidth', 1, 'color', 'black');
    xlim([dim n_dim]);
    hold off
    legend('exact Polyhedron', 'approx Polyhedron', 'exact Star', 'approx Star', 'abs-domain Star','approx Zonotope');
end
title('the reachability time vs. the dimension of a hidden layer in [2 n 2] network');
ylabel('verification time (sec)');
xlabel('dimension of a hidden layer');

%% FNN network of [2 n n n 2]
figure('Name','Computation of reachablility of a FNN network');
for i=dim:n_dim
    disp('------Computing reachable set of a FNN relu network------')
    fprintf('network dimensions: [%d %d %d %d %d]\n',dim, i, i, i, out_dim');
    nN = [dim i out_dim]; % 2 inpus, n hidden layers, 2 outputs
    [W, b] = ReLU.rand_Layers(nN);

    disp('-Exact reachability analysis for Polyhedron')
    [R_pe, ct_pe] = ReLU.FNN_reach_BFS(I, W, b, 'exact');
    fprintf('computation time: %f sec\n', ct_pe);
    
    disp('-Approximate reachability analysis for Polyhedron')
    [R_pa, ct_pa] = ReLU.FNN_reach_BFS(I, W, b, 'approx');
    fprintf('computation time: %f sec\n', ct_pa);

    disp('-Exact reachability analysis for Star Set')
    [R_se, ct_se] = ReLU.FNN_reach_BFS(S, W, b, 'exact');
    fprintf('computation time: %f sec\n', ct_se);

    disp('-Approximate reachability analysis for Star Set')
    [R_sa, ct_sa] = ReLU.FNN_reach_BFS(S, W, b, 'approx');
    fprintf('computation time: %f sec\n', ct_sa);
    
    disp('-Abstract domain Approximate reachability analysis for Star Set')
    [R_sabs, ct_sabs] = ReLU.FNN_reach_BFS(S, W, b, 'abs_domain');
    fprintf('computation time: %f sec\n', ct_sabs);
    
    disp('-Approximate reachability analysis for Zonotope Set')
    [R_za, ct_za] = ReLU.FNN_reach_BFS(Z, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_za);

    hold on
    plot(i,ct_pe,'x','Linewidth', 1, 'color', 'blue');
    plot(i,ct_pa,'x','Linewidth', 1, 'color', 'red');
    plot(i,ct_se,'*','Linewidth', 1, 'color', 'green');
    plot(i,ct_sa,'*','Linewidth', 1, 'color', 'magenta');
    plot(i,ct_sabs,'*','Linewidth', 1, 'color', 'cyan');
    plot(i,ct_za,'square','Linewidth', 1, 'color', 'black');
    xlim([dim n_dim]);
    hold off
    legend('exact Polyhedron', 'approx Polyhedron', 'exact Star', 'approx Star', 'abs-domain Star','approx Zonotope');
end
title('the reachability time vs. the dimension of hidden layers in [2 n n n 2] network');
ylabel('verification time (sec)');
xlabel('dimension of hidden layers');
%}