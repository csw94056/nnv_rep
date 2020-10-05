classdef Sigmoid
   
    properties
    end
    
    methods(Static)
        function [R, ct] = FNN_reach_BFS(I, W, b, func)
            % @I: input set (Polyhedron or Star set)
            % @W: weight matrices
            % @b: bias vectors
            % @R: reachable set
            % @func: logsig - logsig reachability method
            %        tansig - tansig (tanh) reachability method 
            % @ct: computation time
            
            if ~isa(I, 'Star')
                error('Input set must be a Star set');
            end 

            if ~strcmp(func, 'tansig') && ~strcmp(func, 'logsig')
                error('foo:bar','Unknown function type:\nFunction options:\t ''logsig'' or ''tansig''');
            end
           
            tic
            R = I;
            for i=1:length(W)
                R = Sigmoid.reach_approx_star(R, W{i}, b{i}, func);
            end
            ct = toc;
        end
        
        function R = reach_approx_star(I, W, b, func)
            if ~isa(I, 'Star')
                error('Input set must be a Star set');
            end 

            if ~strcmp(func, 'tansig') && ~strcmp(func, 'logsig')
                error('foo:bar','Unknown function type:\nFunction options:\t ''logsig'' or ''tansig''');
            end
            
            if ~isempty(W)
                if ~isempty(b)
                    R = I.affineMap(W, b);
                else
                    R = I.affineMap(W);
                end
            else
                R = I;
            end
            
            n = R.Dim;
            for i = 1:n
               R = Sigmoid.stepSigmoid_star(R, i, func);
            end
        end

        function S = stepSigmoid_star(I, i, func)
            % @I: input star set
            % @i: index of the dimension we want to apply the ReLU activation
            % function

            if ~isa(I, 'Star')
                error('Input set is nto a star set');
            end

            if i <= 0 || i > I.Dim
                error('Invalid index');
            end
            [lb, ub] = I.getRanges;
            lb = lb(i);
            ub = ub(i);
            if strcmp(func, 'tansig')
                yl = tansig(lb);        % y(1) = sigmoid(l)
                yu = tansig(ub);        % y(u) = sigmoid(u)
                dyl = tansig('dn', lb); % y'(1): derivative at lower bound
                dyu = tansig('dn', ub); % y'(u): derivative at upper bound
            elseif strcmp(func, 'logsig')
                yl = logsig(lb);        % y(1) = sigmoid(l)
                yu = logsig(ub);        % y(u) = sigmoid(u)
                dyl = logsig('dn', lb); % y'(1): derivative at lower bound
                dyu = logsig('dn', ub); % y'(u): derivative at upper bound
            else
                error('foo:bar','Unknown function type:\nFunction options:\t ''logsig'' or ''tansig''');
            end
           
            if lb == ub % no approximation need

                V = I.get_V;
            	c = I.get_c;
                V(i,:) = 0;     % reset old center and basis vectors at neuro i to zero
                c(i,1) = yl;    % sigmoid(x(i)) = yl = yu
                S = Star(V, I.C, I.d);

            else % need approximation
                V = I.V;
                V(i, :) = 0; % reset old center and basis vectors at neuron i to zero
                v_new = zeros(I.Dim, 1);
                v_new(i) = 1;
                V_new = [V v_new]; % new center and basis vectors

                n = size(I.C, 1);
                C0 = [I.C zeros(n, 1)];
                d0 = I.d;

                % new constraints
                % case 1
                if lb >= 0
                    % y = a_new, a_new >= y1, a_new <= y2, a_new <= y3
                    % yl: line connect (1, y1) and (u, yu)
                    %     yl = ((yu - yl)/(u - 1)) * (x-1) + yl = a*(x-1) + yl

                    % y2: tangent line at (u, yu)
                    %     y2 = y'(u) * (x - u) + yu

                    % y3: tangent line at (1, yl)
                    %     y3 = y'(l) * (x - 1) + y1

                    % contraint (1): a_new >= yl <-> C1*a <= d1
                    % y1 >= ((yu - yl)/(u - 1)) * (x - 1) + y1 = a*(x-1) + y1
                    V = I.get_V;
                    c = I.get_c;
            
                    a = (yu-yl)/(ub-lb);
                    C1 = [a*V(i,:) -1];
                    d1 = -a*(c(i) - lb) - yl;

                    % constraint (2): a_new <= y2 <-> C2*a <= d2
                    % y2 <= y'(u) * (x - u) + yu
                    C2 = [-dyu*V(i,:) 1];
                    d2 = dyu*(c(i) - ub) + yu;

                    % constraint (3): a_new <= y3 <-> C3*a <= d3
                    % y3 <= y'(1) * (x - 1) + yl
                    C3 = [-dyl*V(i,:) 1];
                    d3 = dyl*(c(i) - lb) + yl;

                    C_new = [C0; C1; C2; C3];
                    d_new = [d0; d1; d2; d3];

                    S = Star(V_new, C_new, d_new);
                end

                % case 2
                if ub <= 0
                    % y = a_new, a_new <= y1, a_new >= y2, a_new >= y3
                    % y1: line connect (l, yl) and (u, yu)
                    %     y1 = ((yu - yl)/(u - l)) * (x - l) + yl = a*(x-1) + yl

                    % y2: tangent line at (u, yu)
                    %     y2 = y'(u) * (x - u) + yu

                    % y3: tangent line at (l, yl)
                    %     y3 = y'(l) * (x - l) + yl

                    % constraint (1): a_new <= y1 <-> C1*a <= d1
                    % y1 <= ((yu - yl)/(u - l)) * (x - l) + yl = a*(x-l) + yl
                    V = I.get_V;
                    c = I.get_c;
            
                    a = (yu - yl)/(ub - lb);
                    C1 = [-a*V(i,:) 1];
                    d1 = a*(c(i) - lb) + yl;

                    % constraint (2): a_new >= y2 <-> C2*a <= d2
                    % y2 >= y'(u) * (x - u) + yu
                    C2 = [dyu*V(i,:) -1];
                    d2 = -dyu*(c(i) - ub) - yu;

                    % constraint (3): a_new >= y3 <-> C3*a <= d3
                    % y3 >= y'(l) * (x - 1) + yl
                    C3 = [dyl*V(i,:) -1];
                    d3 = -dyl*(c(i) - lb) - yl;

                    C_new = [C0; C1; C2; C3];
                    d_new = [d0; d1; d2; d3];

                    S = Star(V_new, C_new, d_new);
                end

                % case 3
                if lb < 0 && ub > 0
                    % y = a_new, a_new <= y1, a_new <= y2, a_new >= y3, a_new >= y4
                    % y1: line connect (l, y1) and (0, a)
                    %     a = ? is the point of intersection of y2 and x=0

                    % y2: tangent line at (u, yu)
                    %     y2 = y'(u) * (x - u) + yu

                    % y3: tangent line at (l, yl)
                    %     y3 = y'(l) * (x - l) + yl

                    % y4: line connect (u, yu) and (0, b)
                    %     b = ? is the point of intersection of y3 and x=0


                    % constraint (1): a_new <= y1 <-> C1*a <= d1
                    % y1: line connect (l, yl) and (0, a)
                    %     a = ? is the point of interection of y2 and x=0
                    %     a = y2(x=0) = y'u*(-u) = yu;

                    % <y1 = ((a - yl)/(0 - 1))*x + a = gamma*x + a
                    V = I.get_V;
                    c = I.get_c;

                    a = -dyu*ub + yu;
                    gamma = (a - yl)/(-lb);
                    C1 = [-gamma*V(i,:) 1];
                    d1 = gamma*c(i) + a;
                    
                    % constraint (2): a_new <= y2 <-> C2*a <= d2
                    % y2 <= y'(u) * (x - u) + yu
                    C2 = [-dyu*V(i,:) 1];
                    d2 = dyu*(c(i) - ub) + yu;
                    
                    % constraint (3): a_new >= y3 <-> C3*a <= d3
                    % y3 >= y'(l) * (x - 1) + yl
                    C3 = [dyl*V(i,:) -1];
                    d3 = -dyl*(c(i) - lb) - yl;

                    % constraint (3): a_new >= y3 <-> C3*a <= d3
                    % y1: line connect (u, yu) and (0, b)
                    %     b = ? is the point of intersection of y3 and x=0
                    %     b = y3(x=0) = y'l*(-1) + yl;

                    % y1 >= ((b - yu)/(0 - u))*x +_ b = gamma*x + b

                    b = -dyl*lb + yl;
                    gamma = (b - yu)/(-ub);
                    C4 = [gamma*V(i,:) -1];
                    d4 = -gamma*c(i) - b;

                    C_new = [C0; C1; C2; C3; C4];
                    d_new = [d0; d1; d2; d3; d4];

                    S = Star(V_new, C_new, d_new);
                end
            end
        end
        
        
        function S = multiStepSigmoid(I, func)
            % @I: input star set
            % @i: index of the dimension we want to apply the ReLU activation
            % function

            if ~isa(I, 'Star')
                error('Input set is nto a star set');
            end

            [lb, ub] = I.getRanges;
            if strcmp(func, 'tansig')
                yl = tansig(lb);        % y(1) = sigmoid(l)
                yu = tansig(ub);        % y(u) = sigmoid(u)
                dyl = tansig('dn', lb); % y'(1): derivative at lower bound
                dyu = tansig('dn', ub); % y'(u): derivative at upper bound
            elseif strcmp(func, 'logsig')
                yl = logsig(lb);        % y(1) = sigmoid(l)
                yu = logsig(ub);        % y(u) = sigmoid(u)
                dyl = logsig('dn', lb); % y'(1): derivative at lower bound
                dyu = logsig('dn', ub); % y'(u): derivative at upper bound
            else
                error('foo:bar','Unknown sigmoid-type function:\nFunction options:\t ''logsig'' or ''tansig''');
            end
           
            if lb == ub % no approximation need

                V = I.get_V;
            	c = I.get_c;
                V(i,:) = 0;     % reset old center and basis vectors at neuro i to zero
                c(i,1) = yl;    % sigmoid(x(i)) = yl = yu
                S = Star(V, I.C, I.d);

            else % need approximation
                V = I.V;
                V(i, :) = 0; % reset old center and basis vectors at neuron i to zero
                v_new = zeros(I.Dim, 1);
                v_new(i) = 1;
                V_new = [V v_new]; % new center and basis vectors

                n = size(I.C, 1);
                C0 = [I.C zeros(n, 1)];
                d0 = I.d;

                % new constraints
                % case 1
                if lb >= 0
                    % y = a_new, a_new >= y1, a_new <= y2, a_new <= y3
                    % yl: line connect (1, y1) and (u, yu)
                    %     yl = ((yu - yl)/(u - 1)) * (x-1) + yl = a*(x-1) + yl

                    % y2: tangent line at (u, yu)
                    %     y2 = y'(u) * (x - u) + yu

                    % y3: tangent line at (1, yl)
                    %     y3 = y'(l) * (x - 1) + y1

                    % contraint (1): a_new >= yl <-> C1*a <= d1
                    % y1 >= ((yu - yl)/(u - 1)) * (x - 1) + y1 = a*(x-1) + y1
                    V = I.get_V;
                    c = I.get_c;
            
                    a = (yu-yl)/(ub-lb);
                    C1 = [a*V(i,1) -1];
                    d1 = -a*(c(i) - lb) - yl;

                    % constraint (2): a_new <= y2 <-> C2*a <= d2
                    % y2 <= y'(u) * (x - u) + yu
                    C2 = [-dyu*V(i,:) 1];
                    d2 = dyu*(c(i) - ub) + yu;

                    % constraint (3): a_new <= y3 <-> C3*a <= d3
                    % y3 <= y'(1) * (x - 1) + yl
                    C3 = [-dyl*V(i,:) 1];
                    d3 = dyl*(c(i) - lb) + yl;

                    C_new = [C0; C1; C2; C3];
                    d_new = [d0; d1; d2; d3];

                    S = Star(V_new, C_new, d_new);
                end

                % case 2
                if ub <= 0
                    % y = a_new, a_new <= y1, a_new >= y2, a_new >= y3
                    % y1: line connect (l, yl) and (u, yu)
                    %     y1 = ((yu - yl)/(u - l)) * (x - l) + yl = a*(x-1) + yl

                    % y2: tangent line at (u, yu)
                    %     y2 = y'(u) * (x - u) + yu

                    % y3: tangent line at (l, yl)
                    %     y3 = y'(l) * (x - l) + yl

                    % constraint (1): a_new <= y1 <-> C1*a <= d1
                    % y1 <= ((yu - yl)/(u - l)) * (x - l) + yl = a*(x-l) + yl
                    V = I.get_V;
                    c = I.get_c;
            
                    a = (yu - yl)/(ub - lb);
                    C1 = [-a*V(i,:) 1];
                    d1 = a*(c(i) - lb) + yl;

                    % constraint (2): a_new >= y2 <-> C2*a <= d2
                    % y2 >= y'(u) * (x - u) + yu
                    C2 = [dyu*V(i,:) -1];
                    d2 = -dyu*(c(i) - ub) - yu;

                    % constraint (3): a_new >= y3 <-> C3*a <= d3
                    % y3 >= y'(l) * (x - 1) + yl
                    C3 = [dyl*V(i,:) -1];
                    d3 = -dyl*(c(i) - lb) - yl;

                    C_new = [C0; C1; C2; C3];
                    d_new = [d0; d1; d2; d3];

                    S = Star(V_new, C_new, d_new);
                end

                % case 3
                if lb < 0 && ub > 0
                    % y = a_new, a_new <= y1, a_new <= y2, a_new >= y3, a_new >= y4
                    % y1: line connect (l, y1) and (0, a)
                    %     a = ? is the point of intersection of y2 and x=0

                    % y2: tangent line at (u, yu)
                    %     y2 = y'(u) * (x - u) + yu

                    % y3: tangent line at (l, yl)
                    %     y3 = y'(l) * (x - l) + yl

                    % y4: line connect (u, yu) and (0, b)
                    %     b = ? is the point of intersection of y3 and x=0


                    % constraint (1): a_new <= y1 <-> C1*a <= d1
                    % y1: line connect (l, yl) and (0, a)
                    %     a = ? is the point of interection of y2 and x=0
                    %     a = y2(x=0) = y'u*(-u) = yu;

                    % <y1 = ((a - yl)/(0 - 1))*x + a = gamma*x + a
                    V = I.get_V;
                    c = I.get_c;

                    a = -dyu*ub + yu;
                    gamma = (a - yl)/(-lb);
                    C1 = [-gamma*V(i,:) 1];
                    d1 = gamma*I.V(i,1) + a;
                    
                    % constraint (2): a_new <= y2 <-> C2*a <= d2
                    % y2 <= y'(u) * (x - u) + yu
                    C2 = [-dyu*V(i,:) 1];
                    d2 = dyu*(c(i) - ub) + yu;
                    
                    % constraint (3): a_new >= y3 <-> C3*a <= d3
                    % y3 >= y'(l) * (x - 1) + yl
                    C3 = [dyl*V(i,:) -1];
                    d3 = -dyl*(c(i) - lb) - yl;

                    % constraint (3): a_new >= y3 <-> C3*a <= d3
                    % y1: line connect (u, yu) and (0, b)
                    %     b = ? is the point of intersection of y3 and x=0
                    %     b = y3(x=0) = y'l*(-1) + yl;

                    % y1 >= ((b - yu)/(0 - u))*x +_ b = gamma*x + b

                    b = -dyl*lb + yl;
                    gamma = (b - yu)/(-ub);
                    C4 = [gamma*V(i,:) -1];
                    d4 = -gamma*c(i) - b;

                    C_new = [C0; C1; C2; C3; C4];
                    d_new = [d0; d1; d2; d3; d4];

                    S = Star(V_new, C_new, d_new);
                end
            end
        end
        
        function Z = stepSigmoid_zono(I, i, func)
            % @I: intermediate input set
            
            if ~isa(I, 'Zonotope')
                error('Input set is not a Zonotope object');
            end
            
            [lb, ub] = I.getRange(i);
            if strcmp(func, 'tansig')
                yl = tansig(lb);        % y(1) = sigmoid(l)
                yu = tansig(ub);        % y(u) = sigmoid(u)
                dyl = tansig('dn', lb); % y'(1): derivative at lower bound
                dyu = tansig('dn', ub); % y'(u): derivative at upper bound
            elseif strcmp(func, 'logsig')
                yl = logsig(lb);        % y(1) = sigmoid(l)
                yu = logsig(ub);        % y(u) = sigmoid(u)
                dyl = logsig('dn', lb); % y'(1): derivative at lower bound
                dyu = logsig('dn', ub); % y'(u): derivative at upper bound
            else
                error('foo:bar','Unknown function type:\nFunction options:\t ''logsig'' or ''tansig''');
            end

            if lb == ub
                X = I.X;
                c = I.c;
              
                X(i,:) = 0;
                c(i) = yl;
                Z = Zonotope(c, X);
            else
                X = I.X;
                c = I.c;
                
                lamda = min(dyl, dyu);
                mu1 = 0.5 * (yu + yl - lamda*(ub + lb));
                mu2 = 0.5 * (yu - yl - lamda*(ub - lb));
                
                c(i) = lamda*c(i) + mu1;
                X(i,:) = lamda*X(i,:);
                new_x = zeros(I.Dim,1);
                new_x(i) = mu2;
                X = [X new_x];
                Z = Zonotope(c,X);
            end
        end
    end
end