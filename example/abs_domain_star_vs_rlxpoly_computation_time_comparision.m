close all;
clear;
clc;

%% computing reachable set of a relu network
dim = 2;
n_dim = 30;
out_dim = 2;

iter = inf; %number of iteration for Relaxed Polyhedron

% random n-dimneesional polyhedron
I = ExamplePoly.randHrep('d', dim); % input set
S = Star(I);
Rp = RlxPoly(I,iter);
Rs = RlxStar(I,iter);
As = AbsStar(I,iter);

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
    
    disp('-Exact Approx reachability analysis for Star ');
    [S_e, ct_se] = ReLU.FNN_reach_BFS(S, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_se);
    
    disp('-Abstract Domain reachability analysis for Relaxed Polyhedron');
    [RP_ab, ct_rpab] = ReLU.FNN_reach_BFS(Rp, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rpab);
    
    disp('-Abstract Domain reachability analysis for Relaxed Star');
    [RS_ab, ct_rsab] = ReLU.FNN_reach_BFS(Rs, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rsab);
    
%     disp('-Abstract Domain reachability analysis for Exact Approx Star');
%     [S_eab, ct_seab] = ReLU.FNN_reach_BFS(As, W, b, 'approx');
%     fprintf('computation time: %f sec\n\n', ct_seab);
    
    hold on
    plot(i,ct_sab,'x','Linewidth', 2, 'color', 'blue');
    plot(i,ct_se,'+','Linewidth', 2, 'color', 'cyan');
    plot(i,ct_rpab,'o','Linewidth', 2, 'color', 'red');
    plot(i,ct_rsab,'^','Linewidth', 2, 'color', 'magenta');
%     plot(i,ct_seab,'s','Linewidth', 2, 'color', 'green');
    xlim([dim n_dim]);
    hold off
    legend('abs-dom Star', 'exact approx Star', 'RlxPoly', 'RlxStar');
%     legend('abs-dom Star', 'exact approx Star', 'RlxPoly', 'RlxStar', 'AbsStar');
%     legend('RlxPoly', 'RlxStar', 'AbsStar');
end
title('the reachability time vs. the dimension of a hidden layer in [2 n 2] network');
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
    
    disp('-Exact Approx reachability analysis for Star ');
    [S_e, ct_se] = ReLU.FNN_reach_BFS(S, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_se);
    
    disp('-Abstract Domain reachability analysis for Relaxed Polyhedron');
    [RP_ab, ct_rpab] = ReLU.FNN_reach_BFS(Rp, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rpab);
    
    disp('-Abstract Domain reachability analysis for Relaxed Star');
    [RS_ab, ct_rsab] = ReLU.FNN_reach_BFS(Rs, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rsab);
    
%     disp('-Abstract Domain reachability analysis for Exact Approx Star');
%     [S_eab, ct_seab] = ReLU.FNN_reach_BFS(As, W, b, 'approx');
%     fprintf('computation time: %f sec\n\n', ct_seab);
    
    hold on
    plot(i,ct_sab,'x','Linewidth', 2, 'color', 'blue');
    plot(i,ct_se,'+','Linewidth', 2, 'color', 'cyan');
    plot(i,ct_rpab,'o','Linewidth', 2, 'color', 'red');
    plot(i,ct_rsab,'^','Linewidth', 2, 'color', 'magenta');
%     plot(i,ct_seab,'s','Linewidth', 2, 'color', 'green');
    xlim([dim n_dim]);
    hold off
    legend('abs-dom Star', 'exact approx Star', 'RlxPoly', 'RlxStar');
%     legend('abs-dom Star', 'exact approx Star', 'RlxPoly', 'RlxStar', 'AbsStar');
%     legend('RlxPoly', 'RlxStar', 'AbsStar');
end
title('the reachability time vs. the dimension of hidden layers in [2 n n n 2] network');
ylabel('verification time (sec)');
xlabel('dimension of a hidden layer');