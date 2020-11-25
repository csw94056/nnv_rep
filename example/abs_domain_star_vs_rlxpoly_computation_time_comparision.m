close all;
clear;
clc;

%% computing reachable set of a relu network
dim = 2;
n_dim = 15;
out_dim = 2;

iter = inf; %number of iteration for Relaxed Polyhedron

% random n-dimneesional polyhedron
% P = ExamplePoly.randHrep('d', dim); % input set
P = Polyhedron('lb',[-1;-1],'ub',[1;1]);
I_z = Zonotope(P);
I_rlxp = RlxPoly(P, inf);
I_rlxs = RlxStar(P, inf);
I_s = Star(P);

%% FNN network of [2 n 2]
figure('Name', 'Computation of reachability of a FNN ReLU network');
for i = dim : n_dim
    disp('------Computing reachable set of a FNN relu network------');
    fprintf('network dimensions: [%d %d %d]\n',dim, i, out_dim');

    nN = [dim i out_dim]; % 2 inpus, n hidden layers, 2 outputs
    [W, b] = ReLU.rand_Layers(nN);
    
    disp('-Star exact');
    [S_exact, ct_sexact] = ReLU.FNN_reach_BFS(I_s, W, b, 'exact');
    fprintf('computation time: %f sec\n\n', ct_sexact);
    
    disp('-Star approx');
    [S_approx, ct_sapprox] = ReLU.FNN_reach_BFS(I_s, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_sapprox);
    
    disp('-Star with Abs-Dom bounds');
    [RlxS, ct_rsab] = ReLU.FNN_reach_BFS(I_rlxs, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rsab);
    
    disp('-Zonotope');
    [Z, ct_z] = ReLU.FNN_reach_BFS(I_z, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_z);
    
    hold on;
    plot(i,ct_sexact,'x','Linewidth', 2, 'color', 'blue');
    plot(i,ct_sapprox,'+','Linewidth', 2, 'color', 'magenta');
    plot(i,ct_rsab,'o','Linewidth', 2, 'color', 'red');
    plot(i,ct_z,'^','Linewidth', 2, 'color', 'green');
    xlim([dim n_dim]);
    hold off;
    legend('Star exact', 'Star approx','Star with Abs-Dom bounds', 'Zonotope');
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
    
    disp('-Star exact');
    [S_exact, ct_sexact] = ReLU.FNN_reach_BFS(I_s, W, b, 'exact');
    fprintf('computation time: %f sec\n\n', ct_sexact);
    
    disp('-Star approx');
    [S_approx, ct_sapprox] = ReLU.FNN_reach_BFS(I_s, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_sapprox);
    
    disp('-Star with Abs-Dom bounds');
    [RlxS, ct_rsab] = ReLU.FNN_reach_BFS(I_rlxs, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_rsab);
    
    disp('-Zonotope');
    [Z, ct_z] = ReLU.FNN_reach_BFS(I_z, W, b, 'approx');
    fprintf('computation time: %f sec\n\n', ct_z);
    
    hold on;
    plot(i,ct_sexact,'x','Linewidth', 2, 'color', 'blue');
    plot(i,ct_sapprox,'+','Linewidth', 2, 'color', 'magenta');
    plot(i,ct_rsab,'o','Linewidth', 2, 'color', 'red');
    plot(i,ct_z,'^','Linewidth', 2, 'color', 'green');
    xlim([dim n_dim]);
    hold off;
    legend('Star exact', 'Star approx','Star with Abs-Dom bounds', 'Zonotope');
end
title('the reachability time vs. the dimension of hidden layers in [2 n n n 2] network');
ylabel('verification time (sec)');
xlabel('dimension of a hidden layer');