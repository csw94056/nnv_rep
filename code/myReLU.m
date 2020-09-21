 classdef myReLU    
     
    properties
    end
    
    methods(Static)
        
        function [R, ct] = relu_FNN_reach_BFS(I, W, b, method)
            % @I: input set (Polyhedron or Star set)
            % @W: weight matrices
            % @b: bias vectors
            % @R: reachable set
            % @method: exact - exact reachability method
            %          approx - approximate reachability method 
            % @ct: computation time
            
            tic
            R = I;
            for i=1:length(W)
                R = myReLU.layerReach(R, W{i}, b{i}, method);
            end
            ct = toc;
        end

        function R = layerReach(I, W, b, method)
            % @I: input set
            % @W: weight matrices
            % @b: bias vectors
            % @R: reachable set of a layer
            % @method: exact - exact reachability method
            %          approx - approximate reachability method
            
            R = [];
            for i=1:length(I)
                if isa(I(i), 'Polyhedron')
                    I1 = I(i).affineMap(W) + b;                 
                elseif isa(I(i), 'myStar')
                    I1 = I(i).affineMap(W, b);
                else
                   error('Unknown type of input set');
                end
                
                if strcmp(method, 'approx')
                    R = [R myReLU.reach_approx(I1)];
                elseif strcmp(method, 'exact')
                    R = [R myReLU.reach_exact(I1)];
                else
                    error('foo:bar','Unknown method\nMethod options:\n\t''exact'' - exact reachbility analysis\n\t''approx'' - over-approximate rechability analysis\n'); 
                end
            end
        end        

        function R = reach_exact(I1)
            % @I1: intermediate input set
            if isa(I1, 'Polyhedron')
                I1.outerApprox;
                lb = I1.Internal.lb;
                ub = I1.Internal.ub;
            elseif isa(I1, 'myStar')
                [lb,ub] = I1.getRanges;  % get ranges of all input variables
            else
               error('Unknow input set'); 
            end
            
            map = find(lb<0);       % construct computation map
            In = I1;
            for i=1:length(map)
                if isa(I1, 'Polyhedron')
                    In = myReLU.stepReLU_poly(In, map(i), lb(map(i)), ub(map(i)));
                elseif isa(I1, 'myStar')
                    In = myReLU.stepReLU_star(In, map(i), lb(i), ub(i));
                end
            end
            R = In;
        end
        
        
        function R = reach_approx(I1)
            % @I1: intermediate input set
            if isa(I1, 'Polyhedron')
                R = myReLU.approxReLU_poly(I1);
            elseif isa(I1, 'myStar')
                R = myReLU.approxReachReLU_star(I1);
            end
        end
        
        %% Polyhedron 
        
        % Over-approximate method
        function B = approxReLU_poly(I1)
            % I: intermediate input set

            I1.outerApprox;
            lb = I1.Internal.lb;
            ub = I1.Internal.ub;
            for i=1:length(lb)
               if lb(i)<=0
                   lb(i)=0;
               end
               if ub(i)<0
                   ub(i)=0;
               end
            end
            B = Polyhedron('lb',lb,'ub',ub);
        end
        
        % Exact method
        function R = stepReLU_poly(I, i, lb, ub)
            % @I: input polyhedron set
            % @i: index of current neuron (x[i])
            % @lb: lower-bound of x[i]
            % @ub: upper-bound of x[i]
            % @R: array of output set

            dim = I.Dim;
            R = [];
            for j = 1:length(I)
                I1 = I(j);
                R1 = [];
                if lb >=0
                    R1 = I1;
                elseif ub < 0
                    % R1 = projection of I1 on x[i]=0
                    M = eye(dim);
                    M(i,i) = 0;
                    R1 = M*I1;
                elseif lb < 0 && ub >= 0
                    b  = vertcat(I1.b, [0]);
                    % Z1 = I1 && x_min>=0;
                    Z1 = zeros(1, dim);
                    Z1(1, i) = -1;
                    A1 = vertcat(I1.A, Z1);
                    Z1 = Polyhedron('A', A1, 'b', b, 'Ae', I1.Ae, 'be', I1.be);

                    % Z2 = I1 && x_max<0;
                    Z2 = zeros(1, dim);
                    Z2(1, i) = 1;
                    A2 = vertcat(I1.A, Z2);
                    Z2 = Polyhedron('A', A2, 'b', b, 'Ae', I1.Ae, 'be', I1.be);

                    M = eye(dim);
                    M(i, i) = 0;
                    Z2 = M*Z2;

                    e1 = Z1.isEmptySet;
                    e2 = Z2.isEmptySet;         
                    if e1 && e2
                        R1 = [];
                    elseif ~e1 && e2
                        R1 = Z1;
                    elseif e1 && ~e2
                        R1 = Z2;
                    else 
                        % if projection of Z2 is a subset of Z1
                        if Z1<=Z2
                            R1 = Z2;
                        else
                            R1 = [Z1, Z2];
                        end
                    end            
                end
                R = [R R1];
            end
        end
        %% Star Set
        
        % Over-approximate method
        function R = approxReachReLU_star(I)
            % @I: intermediate input set

            In = I;
            for i=1:I.Dim
                In = myReLU.approxStepReLU_star(In,i);        
            end
            R = In;
        end
        
        function R = approxStepReLU_star(I, i)
            % I: intermediate input set
            % i: index of current neuron
            % lb: lower-bound of x[i]
            % ub: upper-bound of x[i]
            % R: intermediate output set    
            [lb,ub] = I.getRanges;  % get ranges of all input variables

            lb = lb(i);
            ub = ub(i);
            if lb>=0
                R = I;
            elseif ub<=0
                % R1 = projection of I1 on x[i]=0
                star_V = I.V;
                star_V(i, :) = zeros(1, size(star_V,2));
                R = myStar(star_V, I.C, I.d);
            elseif lb<0 && ub>0
                % y[i] = ReLU(x[i]) = a_(m+1)
                V = I.get_V;
                c = I.get_c;
                % y[i] >= 0
                m = size(I.C,2);
                C1 = [zeros(1,m) -1];
                d1 = 0;
                % y[i] >= x[i]
                C2 = [V(i,:) -1];
                d2 = -c(i);
                % y[i] <= ub(x[i]-lb)/(ub-lb)
                C3 = [-ub*V(i,:)/(ub-lb) 1];
                d3 = ub*(c(i)-lb)/(ub-lb);

                C0 = [I.C zeros(size(I.C,1),1)];
                d0 = I.d;

                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];

                star_V = I.V;
                star_V(i, :) = zeros(1, size(star_V,2));
                e = zeros(I.Dim,1);
                e(i,1) = 1;
                star_V = [star_V e];

                R = myStar(star_V, new_C, new_d);
            end
        end

        % Exact method
        function R = stepReLU_star(I, i, lb, ub)
            % @I: input star set
            % @i: index of current neuron (x[i])
            % @lb: lower-bound of x[i]
            % @ub: upper-bound of x[i]
            % @R: array of output set
            
            dim = I.Dim;
            R = [];            
            for j = 1:length(I)
                I1 = I(j);
                R1 = [];
                if lb >= 0
                    R1 = I1;
                elseif ub < 0
                    % R1 = projection of I1 on x[i]=0
                    star_V = I1.V;
                    star_V(i, :) = zeros(1, size(star_V,2));
                    R1 = myStar(star_V, I.C, I.d);
                elseif lb <  0 && ub >= 0
                    c = I1.get_c;
                    V = I1.get_V;

                    % I1 && x_min>=0;
                    star_C = vertcat(I1.C, -V(i,:));
                    star_d = vertcat(I1.d, c(i,1));
                    Z1 = myStar(I1.V, star_C, star_d);

                    % I1 && x_max<0;
                    star_C = vertcat(I1.C, V(i,:));
                    star_d = vertcat(I1.d, -c(i,1));
                    star_V = I1.V;
                    star_V(i, :) = zeros(1, size(star_V,2));
                    Z2 = myStar(star_V, star_C, star_d);            

                    a = Z1.isEmptySet;
                    b = Z2.isEmptySet;
                    if a && ~b
                        R1 = Z2;
                    elseif ~a && b
                        R1 = Z1;
                    elseif ~a && ~b 
                        R1 = [Z1 Z2];
                    else
                        R1 = [];
                    end
                end
                R = [R R1];
            end
        end
        
        %% misc functions
        function [W, b] = rand_Layers(nN)
            % generate random weight matrices and bias vectors
            % example: nN = [2 3 2]; % 2 inpus, 1 hidden layers, 2 outputs

            W = cell(1, length(nN) - 1); % weight matrices
            b = cell(1, length(nN) - 1); % bias vectors

            for i=1:length(nN)-1
                W{i} = rand(nN(i+1), nN(i));
                b{i} = rand(nN(i+1), 1);
            end
        end

        function R = set_div(I, row, col)
            % @I: input set (Polyhedron or Star)
            % @row: number of rows divide to
            % @col: number of columns divide to
            % @R: array of sets
            if isa(I, 'Polyhedron') 
                I.outerApprox;
                lb = I.Internal.lb;
                ub = I.Internal.ub;
            elseif isa(I, 'myStar')
                [lb,ub] = I.getRanges;
            end
            col_div = (ub(1)-lb(1))/col;
            row_div = (ub(2)-lb(2))/row;
            k=1;
            for i=1:col
                for j=1:row
                    lb_1=lb(1)+(i-1)*col_div;
                    lb_2=lb(2)+(j-1)*row_div;
                    ub_1=lb(1)+i*col_div;
                    ub_2=lb(2)+j*row_div;
                    C = [eye(I.Dim);-eye(I.Dim)];
                    d = [ub_1;ub_2;-lb_1;-lb_2];
                    if isa(I, 'Polyhedron') 
                        R(k) = Polyhedron('A',C, 'b',d);
                    elseif isa(I, 'myStar')
                        V = [zeros(I.Dim,1) eye(I.Dim)];
                        S = myStar(V, C, d);
                        R(k) = Intersect(I,S);
                    end
                    k=k+1;
                end
            end
        end
        
        function plot_set(I, str)
            % I: input set
            % str: title
            % ploting polyhedron or star set with title
            nexttile
            plot(I);
            title(str);
            xlabel('x');
            ylabel('y');
            zlabel('z');
        end
    end
 end