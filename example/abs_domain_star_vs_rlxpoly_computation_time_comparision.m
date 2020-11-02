close all;
clear;
clc;

%% computing reachable set of a relu network
dim = 2;
n_dim = 10;
out_dim = 2;

iter = inf; %number of iteration for Relaxed Polyhedron

% random n-dimneesional polyhedron
I = ExamplePoly.randHrep('d', dim); % input set
S = Star(I);
R = RlxPoly(I,iter);

%% FNN network of [2 n 2]
figure('Name', 'Computation of reachability of a FNN ReLU network');
for i = dim : n_dim
    disp('------Computing reachable set of a FNN relu network------');
    fprintf('network dimensions: [%d %d %d]\n',dim, i, out_dim');

    nN = [dim i out_dim]; % 2 inpus, n hidden layers, 2 outputs
    [W, b] = ReLU.rand_Layers(nN);
    
    disp('-Abstract Domain reachability analysis for Star ');
    [S_ab, ct_sab] = ReLU.FNN_reach_BFS(S, W, b, 'abs_domain');
    fprintf('computation time: %f sec\n\n', ct_sab);
    
    disp('-Abstract Domain reachability analysis for Relaxed Polyhedron');
    [R_ab, ct_rab] = ReLU.FNN_reach_BFS(R, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rab);
    
    disp('-Abstract Domain reachability analysis for Relaxed Star');
    [RS_ab, ct_rsab] = ReLU.FNN_reach_BFS(R, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rsab);
    
    hold on
    plot(i,ct_sab,'x','Linewidth', 2, 'color', 'blue');
    plot(i,ct_rab,'o','Linewidth', 2, 'color', 'red');
    plot(i,ct_rab,'+','Linewidth', 2, 'color', 'yellow');
    xlim([dim n_dim]);
    hold off
    legend('abs-dom Star', 'RlxPoly', 'RlxStar');
end
title('the reachability time vs. the dimension of a hidden layer');
ylabel('verification time (sec)');
xlabel('dimension of a hidden layer');

%% FNN network of [2 n n n 2]
figure('Name','Computation of reachability of a FNN relu network');
for i=dim:n_dim
    disp('------Computing reachable set of a FNN relu network------')
    fprintf('network dimensions: [%d %d %d %d %d]\n',dim, i, i, i, out_dim');
    
    nN = [dim i i i out_dim]; % 2 inpus, n hidden layers, 2 outputs
    [W, b] = ReLU.rand_Layers(nN);
    
    disp('-Abstract Domain reachability analysis for Star ');
    [S_ab, ct_sab] = ReLU.FNN_reach_BFS(S, W, b, 'abs_domain');
    fprintf('computation time: %f sec\n\n', ct_sab);
    
    disp('-Abstract Domain reachability analysis for Relaxed Polyhedron');
    [R_ab, ct_rab] = ReLU.FNN_reach_BFS(R, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rab);
    
    disp('-Abstract Domain reachability analysis for Relaxed Star');
    [RS_ab, ct_rsab] = ReLU.FNN_reach_BFS(R, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rsab);
    
    hold on
    plot(i,ct_sab,'x','Linewidth', 2, 'color', 'blue');
    plot(i,ct_rab,'o','Linewidth', 2, 'color', 'red');
    plot(i,ct_rab,'+','Linewidth', 2, 'color', 'yellow');
    xlim([dim n_dim]);
    hold off
    legend('abs-dom Star', 'RlxPoly', 'RlxStar');
end
title('the reachability time vs. the dimension of a hidden layer');
ylabel('verification time (sec)');
xlabel('dimension of a hidden layer');