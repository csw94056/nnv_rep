classdef DeepPoly
    properties
        V = []; % c (center vector) and X (a set of basic vectors)
        lower_a = []; % matrix for lower polyhedral constraint (a[<=]) ((1 a[1] a[2] ... a[n])'
        upper_a = []; % matrix for upper polyhedral constraint (a[>=])
        lb = []; % matrix for lower bound
        ub = []; % matrix for upper bound
        Dim = 0;
        n = 0; % number of variables of each nodes
    end
    
    methods
        function obj = DeepPoly(varargin)
            switch nargin
                case 0
                case 1
                    I = varargin{1};
                    if isa(I, 'Polyhedron')
                      
                    elseif isa(I, 'Star')
                        
                    elseif isa(I, 'Zonotope')
                        
                        %{
                        A = [eye(I.Dim), -eye(I.Dim)];
                        b = ones(I.Dim,1);
                        [lb,ub] = I.getRanges;
                        obj = DeepPoly(I.V, A, b, lb, ub);
                        %}
                    else
                        error('Unkown imput set');
                    end
                case 4
                    lower_a = varargin{1};
                    upper_a = varargin{2};
                    lb = varargin{3};
                    ub = varargin{4};
                    
                    [nL, mL] = size(lower_a);
                    
                    obj.lower_a = lower_a;
                    obj.upper_a = upper_a;
                    obj.lb = lb;
                    obj.ub = ub;
                    obj.Dim = mL-1;
                case 5
                    V = varargin{1};
                    lower_a = varargin{2};
                    upper_a = varargin{3};
                    lb = varargin{4};
                    ub = varargin{5};
                    
                    [nL, mL] = size(lower_a);
                    
                    obj.V = V;
                    obj.lower_a = lower_a;
                    obj.upper_a = upper_a;
                    obj.lb = lb;
                    obj.ub = ub;
                    obj.Dim = mL-1;
            end
        end
        
        function c = get_c(obj)
           c = obj.V(:,1);
        end
        
        function X = get_X(obj)
            X = obj.V(:,2:size(obj.V,2));
        end
        
        function D = affineMap(varargin)
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
            
            V = W * obj.V;
            if mb ~= 0
                lower_a = W * obj.lower_a + b;
                upper_a = W * obj.upper_a + b;
                V(:, 1) = V(:, 1) + b;
            else
                lower_a = W * obj.lower_a;
                upper_a = W * obj.upper_a;
            end
            
            lb = abs(W) * obj.lb;
            ub = abs(W) * obj.ub;
            %}
            %{
            lb = zeros(obj.Dim, 1);
            ub = zeros(obj.Dim, 1);
            for i = 1:nW
                for j = 1:mW
                   if W(i,j) >= 0
                       lb(i) = lb(i) + W(i,j) * obj.lb(j);
                       ub(i) = ub(i) + W(i,j) * obj.ub(j);
                   else 
                       w = W(i,j)
                       i = i
                       j = j
                       ub = ub
                       lb = lb
                       obj_lb = obj.lb(j)
                       obj_ub = obj.ub(j)
                       ub(i) = ub(i) - W(i,j) * obj.lb(j);
                       lb(i) = lb(i) - W(i,j) * obj.ub(j);
                   end
                end
            end
            %}
            D = DeepPoly(V, lower_a, upper_a, lb, ub); 
        end
        
        function D = affineMap2(varargin)
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
            
%             lower_a = [zeros(obj.Dim,1) W*obj.lower_a(end - obj.Dim + 1: end, 2:end)];
%             upper_a = [zeros(obj.Dim,1) W*obj.upper_a(end - obj.Dim + 1: end, 2:end)];
            lower_a = [zeros(obj.Dim,1) W*eye(obj.Dim)];
            upper_a = [zeros(obj.Dim,1) W*eye(obj.Dim)];

            lower_a = [obj.lower_a; lower_a];
            upper_a = [obj.upper_a; upper_a];
            
            n = size(obj.lb,1) - obj.Dim + 1;
            
%             lb = -abs(W) * obj.lb(n:end);
%             ub = abs(W) * obj.ub(n:end);
            lb = zeros(obj.Dim, 1);
            ub = zeros(obj.Dim, 1);
            for i = 1:nW
                for j = 1:mW
                   if W(i,j) >= 0
                       lb(i) = lb(i) + W(i,j) * obj.lb(end - obj.Dim + j);
                       ub(i) = ub(i) + W(i,j) * obj.ub(end - obj.Dim + j);
                   else
                       ub(i) = ub(i) + W(i,j) * obj.lb(end - obj.Dim + j);
                       lb(i) = lb(i) + W(i,j) * obj.ub(end - obj.Dim + j);
                   end
                end
            end
            
            lb = [obj.lb; lb];
            ub = [obj.ub; ub];
            
            D = DeepPoly(lower_a, upper_a, lb, ub); 
        end
        
        function P = toPolyhedron(obj)
            P = Polyhedron('A', [obj.A], 'b', [obj.b]);
        end
        
        function [lb,ub] = getRange(obj, i)
            lb = obj.lb(i);
            ub = obj.ub(i);
        end
        
        function [lb,ub] = getRanges(obj)
            lb = I.lb;
            ub = I.ub;
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
    end
end