classdef Star
    properties
        V = []; % a set of basic vectors
        C = []; % predicate constraint matrix
        d = []; % predicate constraint vector
        Dim = 0;
        
        P_lb = []; % lower bound vector of predicate variable
        P_ub = []; % upper bound vector of predicate variable
    end
    methods
        function obj = Star(varargin)
            switch nargin
                case 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0
                    obj.V = [];
                    obj.C = [];
                    obj.d = [];
                    obj.Dim = 0;
                    
                    obj.P_lb = [];
                    obj.P_ub = [];
                case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1
                    P = varargin{1};
                    if ~isa(P, 'Polyhedron')
                        error('Input set is not a polyhedron object');
                    end
                    c = zeros(P.Dim,1);
                    Ve = eye(P.Dim);
                    V = [c Ve];
                    if ~isempty(P.Ae)
                        C = [P.A;P.Ae;-P.Ae];
                        d = [P.b;P.be;-P.be];
                    else
                        C = P.A;
                        d = P.b;
                    end
                    
                    P.outerApprox;
                    P_lb = P.Internal.lb;
                    P_ub = P.Internal.ub;
                    obj = Star(V, C, d, P_lb, P_ub);
                case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3
                    %V need to contain c and V
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    [nV, mV] = size(V);
                    [nC, mC] = size(C);
                    [nd, md] = size(d);
                    
                    if mV ~= mC + 1
                        error('Inconsistency between basic matrix and constraint matrix');
                    end

                    if nC ~= nd
                        error('Inconsistency between constraint matrix and constraint vector');
                    end

                    if md ~= 1
                        error('constraint vector should be one column');
                    end
                    
                    obj.V = V;
                    obj.C = C;
                    obj.d = d;
                    obj.Dim = nV;
                case 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4
                    c = varargin{1};
                    V = varargin{2};
                    C = varargin{3};
                    d = varargin{4};
                    
                    obj = Star([c V], C, d);
                case 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    P_lb = varargin{4};
                    P_ub = varargin{5};
                    
                    [nV, mV] = size(V);
                    [nC, mC] = size(C);
                    [nd, md] = size(d);
                    [n_lb, m_lb] = size(P_lb);
                    [n_ub, m_ub] = size(P_ub);
                    
                    if mV ~= mC + 1
                        error('Inconsistency between basic matrix and constraint matrix');
                    end

                    if nC ~= nd
                        error('Inconsistency between constraint matrix and constraint vector');
                    end

                    if md ~= 1
                        error('constraint vector should be one column');
                    end
                    
                    if n_lb ~= mC || n_ub ~= mC
                       error('dimensions of predicate bounds and constraints do not match');
                    end
                    
                    if m_lb ~= 1 || m_ub ~= 1
                       error('predicate lower and upper bound should be one column'); 
                    end
                    
                    obj.V = V;
                    obj.C = C;
                    obj.d = d;
                    obj.Dim = nV;
                    
                    obj.P_lb = P_lb;
                    obj.P_ub = P_ub;                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end switch    
            end
        end
            
        function S = affineMap(varargin)
            b = [];
            switch nargin
                case 3
                    obj = varargin{1};
                    W = varargin{2};
                    b = varargin{3};
                case 2
                    obj = varargin{1};
                    W = varargin{2};
            end
            
            [nW, mW] = size(W);
            [nb, mb] = size(b);
            
            if mW ~= obj.Dim
                error('Inconsistency between affine transformation mattrix and object dimension');
            end
            
            if mb > 1
                error('bias vector must have one column');
            end

            V = W * obj.V ;
            if mb ~= 0
                V(:, 1) = V(:, 1) + b;
            end
            if ~isempty(obj.P_lb)
                S = Star(V, obj.C, obj.d, obj.P_lb, obj.P_ub);
            else
                S = Star(V, obj.C, obj.d);
            end
        end
        
%         function [lb, ub] = getRanges(obj)
%             P = obj.toPolyhedron;
%             P.outerApprox;
%             lb = P.Internal.lb;
%             ub = P.Internal.ub;
%         end

%%%%%%%%%%%%%%%%%%%%%%%%% code from nnv ***********************************
        function [lb, ub] = getRanges(obj)
            
            % author: Dung Tran
            % date: 7/19/2019
            
            if ~obj.isEmptySet            
                n = obj.Dim;
                lb = zeros(n,1);
                ub = zeros(n,1);
                for i=1:n
                    %fprintf('\nGet range at index %d', i);
                    [lb(i), ub(i)] = obj.getRange(i);
                end
            else
                lb = [];
                ub = [];
            end

        end
        
        function [xmin, xmax] = getRange(obj, index)
            % @index: position of the state
            % range: min and max values of x[index]
            
            % author: Dung Tran
            % date: 11/16/2018
            
            if index < 1 || index > obj.Dim
                error('Invalid index');
            end
            
            f = obj.V(index, 2:end);
            if all(f(:)==0)
                xmin = obj.V(index,1);
                xmax = obj.V(index,1);
            else
                % **** linprog is much faster than glpk
                options = optimoptions(@linprog, 'Display','none'); 
                options.OptimalityTolerance = 1e-10; % set tolerance
                [~, fval, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], obj.P_lb, obj.P_ub, options);             
%                 [~, fval, exitflag, ~] = glpk(f, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
                if exitflag > 0
                    xmin = fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution, exitflag = %d', exitflag);
                end          
          
                [~, fval, exitflag, ~] = linprog(-f, obj.C, obj.d, [], [], obj.P_lb, obj.P_ub, options);   
                %[~, fval, exitflag, ~] = glpk(-f, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
                if exitflag > 0
                    xmax = -fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution');
                end

            end
                        
        end
%%%%%%%%%%%%%%%%%%%%%%%%% code from nnv ***********************************
        function V = get_V(obj)
            %V = obj.V(:,2:size(obj.V,2));
            V = obj.V(:, 2:end);
        end
        
        function c = get_c(obj)
            c = obj.V(:,1);
        end
        
        function P = toPolyhedron(obj)
%             if ~isempty(obj.P_lb) && ~isempty(obj.P_ub)
%                 m = size(obj.C, 2);
%                 C = [eye(m); -eye(m)];
%                 d = [obj.P_ub; -obj.P_lb];
%                 
%                 C = [obj.C; C];
%                 d = [obj.d; d];
%             else
                C = obj.C;
                d = obj.d;
%             end
            
            P1 = Polyhedron('A', [C], 'b', [d]);
            V = obj.get_V;
            c = obj.get_c;
            P = V*P1 + c;
        end
        
        function plot (varargin)
            color = 'red';
                   
            switch nargin
                case 1
                    obj = varargin{1};
                case 2
                    obj = varargin{1};
                    color = varargin{2};
            end
          
            n = length(obj);
            if ~strcmp(color, 'rand')
                c_rand = color;
            end
            hold on;
            for i=1:n
                I = obj(i);
                P = I.toPolyhedron;
                if strcmp(color, 'rand')
                    c_rand = rand(1,3);
                end
            
                plot(P, 'color', c_rand);
            end
            hold off
        end
        
        function S = Intersect(obj1, obj2)
          
            C1 = obj2.C * obj1.get_V;
            d1 = obj2.d - obj2.C * obj1.get_c;

            new_C = [obj1.C; C1];
            new_d = [obj1.d; d1];
            S = Star(obj1.V, new_C, new_d);     
            if isEmptySet(S)
                S = [];
            end
        end
            
            
        function a = isEmptySet(obj)
            a = 1;
            for i=1:length(obj)
                I = obj(i);
                P = Polyhedron('A', I.C, 'b', I.d);
                if P.isEmptySet==0
                    a = 0;
                end
            end
        end
%         function bool = isEmptySet(obj)
%             % author: Dung Tran
%             % date: 
%             % update: 6/16/2020
%             % update: 7/15/2020 The isEmptySet method in Polyhedron object
%             % has bug
%             
%             options = optimoptions(@linprog, 'Display','none'); 
%             options.OptimalityTolerance = 1e-10; % set tolerance
%             V = obj.get_V;
%             nVar = size(V,2);
%             f = zeros(1, nVar);
% %             [~, ~, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, options);
%             [~, ~, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], [], [], options);
%             if exitflag == 1
%                 bool = 0;
%             elseif exitflag == -2
%                 bool = 1;
%             else
%                 error('Error, exitflag = %d', exitflag);
%             end
%             
% %             P = Polyhedron('A', obj.C, 'b', obj.d, 'lb', obj.predicate_lb, 'ub', obj.predicate_ub);
% %             bool = P.isEmptySet;
%             
%         end
    end
end
