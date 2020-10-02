 classdef ReLU    
     
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
                R = ReLU.layerReach(R, W{i}, b{i}, method);
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
                elseif isa(I(i), 'Star')
                    I1 = I(i).affineMap(W, b);
                elseif isa(I(i), 'Zonotope')
                    I1 = I(i).affineMap(W, b);
                else
                    error('Unrecognized input set');
                end
                
                if strcmp(method, 'exact') 
                    R = [R ReLU.reach_exact(I1)];
                elseif strcmp(method, 'approx') || strcmp(method, 'zono') || strcmp(method, 'abs_domain')
                    R = [R ReLU.reach_approx(I1, method)];
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
            elseif isa(I1, 'Star')
                [lb,ub] = I1.getRanges;  % get ranges of all input variables
            else
                error('Unrecognized input set for exact ReLU');
            end
            
            map = find(lb<0);       % construct computation map
            In = I1;
            for i=1:length(map)
                if isa(I1, 'Polyhedron')
                    In = ReLU.stepReLU_poly(In, map(i), lb(map(i)), ub(map(i)));
                elseif isa(I1, 'Star')
                    In = ReLU.stepReLU_star(In, map(i), lb(i), ub(i));
                else
                    error('Unrecognized input set for exact ReLU');
                end
            end
            R = In;
        end
        
        
        function R = reach_approx(I1, method)
            % @I1: intermediate input set
            if isa(I1, 'Polyhedron') && strcmp(method, 'approx')
                R = ReLU.approxReLU_poly(I1);
            elseif isa(I1, 'Star')
                if strcmp(method, 'approx')
                    R = ReLU.approxReachReLU_star(I1);
                elseif strcmp(method, 'zono')
                    R = ReLU.zonoApproxReachReLU_star(I1);
                elseif strcmp(method, 'abs_domain')
                    R = ReLU.absDomainReachReLU_star(I1);
                else
                    error('Unknown method for Star set for over-approximate ReLU');
                end
            elseif isa(I1, 'Zonotope')
                R = ReLU.approxReachReLU_zono(I1);
            else 
                error('Unrecognized input set/method for over-approximation ReLU');
            end
        end
        
        %% Polyhedron 
        
        % Over-approximate method
        function B = approxReLU_poly(I1)
            % I: intermediate input set
            
            if ~isa(I1, 'Polyhedron')
                error('Input set is not a Polyhedron object');
            end
            
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
            
            if ~isa(I, 'Polyhedron')
                error('Input set is not a Polyhedron object');
            end

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
        
        % Over-approximate method based on absstract domain
        
        function R = absDomainReachReLU_star(I)
            % @I: intermediate input set
            
            if ~isa(I, 'Star')
                error('Input set is not a Star object');
            end

            In = I;
            for i=1:I.Dim
                In = ReLU.absDomainStepReLU_star(In,i);        
            end
            R = In;
        end
        
        function R = absDomainStepReLU_star(I, i)
            % Zonotope approximation with decreasing function of lamda
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
                star_V(i,:) = 0;
                R = myStar(star_V, I.C, I.d);
            elseif lb<0 && ub>0
                % y[i] = ReLU(x[i]) = a[m+1]
                V = I.get_V;
                c = I.get_c;
                
                A1 = 0.5 * ub * (ub - lb);
                A2 = 0.5 * (-lb) * (ub - lb);
                
                % constraint 1: y[i] <= ub(x[i]-lb)/(ub-lb)
                C1 = [-ub*V(i,:)/(ub-lb) 1];
                d1 = ub*(c(i)-lb)/(ub-lb);
                     
                % constraint 2: y[i] >= 0
                m = size(I.C,2);
                C2 = [zeros(1,m) -1];
                d2 = 0;
                
                % constraint 3: y[i] >= x[i]
                C3 = [V(i,:) -1];
                d3 = -c(i);

                C0 = [I.C zeros(size(I.C,1),1)];
                d0 = I.d;

                star_V = I.V;
                star_V(i,:) = 0;
                e = zeros(I.Dim,1);
                e(i,1) = 1;
                star_V = [star_V e];
                
                if A1 <= A2
                    new_C = [C0; C1; C2];
                    new_d = [d0; d1; d2];
                    R = Star(star_V, new_C, new_d);
                else
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                    R = Star(star_V, new_C, new_d);
                end
            end
        end
        
        % Over-approximate method based on Zonotope
        function R = zonoApproxReachReLU_star(I)
            % @I: intermediate input set
            
            if ~isa(I, 'Star')
                error('Input set is not a Star object');
            end
            In = I;
            for i=1:I.Dim
                In = ReLU.zonoApproxStepReLU_star(In,i);        
            end
            R = In;
        end
        
        function R = zonoApproxStepReLU_star(I, i)
            % Zonotope approximation with decreasing function of lamda
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
                star_V(i,:) = 0;
                R = myStar(star_V, I.C, I.d);
            elseif lb<0 && ub>0
                % y[i] = ReLU(x[i]) = a[m+1]
                V = I.get_V;
                c = I.get_c;
                lamda = ub/(ub - lb);
                
                A1 = (1-lamda) * ub * (ub - lb);
                A2 = lamda * (-lb) * (ub - lb);
                
                % constraint 1: y[i] >= lamda*x[i]
                C1 = [lamda*V(i,:) -1];
                d1 = -lamda*c(i);

                % constraint 2: y[i] <= lamda*x[i] + ub*(1-lamda)
                C2 = [-lamda*V(i,:) 1];
                d2 = -lamda*c(i) + ub*(1-lamda);
               
                % constraint 3: y[i] <= lamda*(x[i] - lb)
                C3 = [-lamda*V(i,:) 1];
                d3 = lamda*(c(i)-lb);

                C0 = [I.C zeros(size(I.C,1),1)];
                d0 = I.d;
                
                star_V = I.V;
                star_V(i,:) = 0;
                e = zeros(I.Dim,1);
                e(i,1) = 1;
                star_V = [star_V e];
                
                if A1 <= A2
                    new_C = [C0; C1; C2];
                    new_d = [d0; d1; d2];
                    R = Star(star_V, new_C, new_d);

                else
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                    R = Star(star_V, new_C, new_d);
                end
            end
        end
        
        % Over-approximate method
        function R = approxReachReLU_star(I)
            % @I: intermediate input set
            
            if ~isa(I, 'Star')
                error('Input set is not a Star object');
            end

            In = I;
            for i=1:I.Dim
                In = ReLU.approxStepReLU_star(In,i);        
            end
            R = In;
        end
        
        function R = approxStepReLU_star(I, i)
            % I: intermediate input set
            % i: index of current neuron
            % lb: lower-bound of x[i]
            % ub: upper-bound of x[i]
            % R: intermediate output set    
            
            if ~isa(I, 'Star')
                error('Input set is not a Star object');
            end  
            
            [lb,ub] = I.getRanges;  % get ranges of all input variables
            
            lb = lb(i);
            ub = ub(i);
            
            if lb>=0
                R = I;
            elseif ub<=0
                % R1 = projection of I1 on x[i]=0
                star_V = I.V;
                star_V(i,:) = 0;
                R = Star(star_V, I.C, I.d);
            elseif lb<0 && ub>0
                % y[i] = ReLU(x[i]) = a[m+1]
                V = I.get_V;
                c = I.get_c;
                
                % constraint 1: y[i] >= 0
                m = size(I.C,2);
                C1 = [zeros(1,m) -1];
                d1 = 0;
                
                % constraint 2: y[i] >= x[i]
                C2 = [V(i,:) -1];
                d2 = -c(i);
                
                % constraint 3: y[i] <= ub(x[i]-lb)/(ub-lb)
                C3 = [-ub*V(i,:)/(ub-lb) 1];
                d3 = ub*(c(i)-lb)/(ub-lb);

                C0 = [I.C zeros(size(I.C,1),1)];
                d0 = I.d;

                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];

                star_V = I.V;
                star_V(i,:) = 0;
                e = zeros(I.Dim,1);
                e(i,1) = 1;
                star_V = [star_V e];

                R = Star(star_V, new_C, new_d);
            end
        end

        % Exact method
        function R = stepReLU_star(I, i, lb, ub)
            % @I: input star set
            % @i: index of current neuron (x[i])
            % @lb: lower-bound of x[i]
            % @ub: upper-bound of x[i]
            % @R: array of output set
            
            if ~isa(I, 'Star')
                error('Input set is not a Star object');
            end
            
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
                    R1 = Star(star_V, I.C, I.d);
                elseif lb <  0 && ub >= 0
                    c = I1.get_c;
                    V = I1.get_V;
                    
                    % I1 && x_min>=0;
                    star_C = vertcat(I1.C, -V(i,:));
                    star_d = vertcat(I1.d, c(i,1));
                    Z1 = Star(I1.V, star_C, star_d);

                    % I1 && x_max<0;
                    star_C = vertcat(I1.C, V(i,:));
                    star_d = vertcat(I1.d, -c(i,1));
                    star_V = I1.V;
                    star_V(i, :) = zeros(1, size(star_V,2));
                    Z2 = Star(star_V, star_C, star_d);            

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
        
        %% Zonotope
        
        function R = approxReachReLU_zono(I)
            % @I: intermediate input set
            
            if ~isa(I, 'Zonotope')
                error('Input set is not a Zonotope object');
            end

            In = I;
            for i=1:I.Dim
                In = ReLU.approxStepReLU_zono(In,i);        
            end
            R = In;
        end
        
        function R = approxStepReLU_zono(I, i)
            % @I: intermediate input set
            
            if ~isa(I, 'Zonotope')
                error('Input set is not a Zonotope object');
            end
            
            [lb, ub] = I.getRange(i);

            if lb > 0
                R = I;
            elseif ub <= 0
                c = I.c;
                X = I.X;

                c(i) = 0;
                X(i,:) = 0;
                R = Zonotope(c, X);
            elseif lb<0 && ub>0
                c = I.c;
                X = I.X;
                
                lamda = ub/(ub-lb);
                mu = -ub*lb/(2*(ub-lb));

                c(i) = lamda*c(i) + mu;
                X(i,:) = lamda*X(i,:);
                new_x = zeros(I.Dim,1);
                new_x(i) = mu;
                X = [X new_x];
                R = Zonotope(c,X);
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
            elseif isa(I, 'Star')
                [lb,ub] = I.getRanges;
            else 
                error('Unrecognized input set');
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
                    elseif isa(I, 'Star')
                        V = [zeros(I.Dim,1) eye(I.Dim)];
                        S = Star(V, C, d);
                        R(k) = Intersect(I,S);
                    end
                    k=k+1;
                end
            end
        end
        
        function plot_set(varargin)            
            str = [];
            color = 'red';
            
            switch nargin 
                case 1
                    I = varargin{1};
                case 2
                    I = varargin{1};
                    str = varargin{2};
                case 3
                    I = varargin{1};
                    str = varargin{2};
                    color = varargin{3};
            end
            
            ReLU.check_object(I);
            
            % I: input set
            % str: title
            % ploting polyhedron or star set with title
            nexttile;
            if isa(I, 'Polyhedron')
                if strcmp(color, 'rand')
                    plot(I);
                else
                    plot(I, 'color', color);
                end
            else
                plot(I, color);
            end
            title(str);
            xlabel('x');
            ylabel('y');
            zlabel('z');
        end
        
        function plot_set_hold(varargin)            
            str = [];
            color = 'red';
            
            switch nargin 
                case 1
                    I = varargin{1};
                case 2
                    I = varargin{1};
                    str = varargin{2};
                case 3
                    I = varargin{1};
                    str = varargin{2};
                    color = varargin{3};
            end
            
            ReLU.check_object(I);
            
            % I: input set
            % str: title
            % ploting polyhedron or star set with title
            hold on;
            if isa(I, 'Polyhedron') 
                if strcmp(color, 'rand')
                    plot(I);
                else
                    plot(I, 'color', color);
                end
            else
                plot(I, color);
            end
            hold off;
            title(str);
            xlabel('x');
            ylabel('y');
            zlabel('z');
        end
        
        function check_object(I)
            if ~isa(I, 'Polyhedron') && ~isa(I, 'Star') && ~isa(I, 'Zonotope')
                error('Unrecognized object');
            end
        end
    end
 end