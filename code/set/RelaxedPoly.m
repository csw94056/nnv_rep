classdef RelaxedPoly
    properties
        lower_a = []; % matrix for lower polyhedral constraint (a[<=]) ((1 a[1] a[2] ... a[n])'
        upper_a = []; % matrix for upper polyhedral constraint (a[>=])
        lb = []; % matrix for lower bound
        ub = []; % matrix for upper bound

        nVar = []; % number of variables of each layer
        Dim = 0; % dimension of current relaxed polyhedron
    end
    
    methods
        function obj = RelaxedPoly(varargin)
            switch nargin
                case 0
                    obj.lower_a = [];
                    obj.upper_a = [];
                    obj.lb = [];
                    obj.ub = [];
                    obj.nVar = [];
                    obj.Dim = 0;
                case 1
                    I = varargin{1};
                    dim = I.Dim;
                    if isa(I, 'Polyhedron')
                        I.outerApprox;
                        lb = I.Internal.lb;
                        ub = I.Internal.ub;
                        
                        lower_a = zeros(dim, 1);
                        upper_a = zeros(dim, 1);
                        for i=1:dim
                            l = zeros(dim, 1); 
                            u = zeros(dim, 1);
                            l(i) = lb(i);
                            u(i) = ub(i);
                            lower_a = [lower_a l];
                            upper_a = [upper_a u];
                        end
                    elseif isa(I, 'Star')
                        [lb, ub] = I.getRanges();
                        
                        c = (ub + lb)*0.5;
                        x = (ub - lb)*0.5;

                        X = [];
                        for i=1:dim
                            X1 = zeros(dim, 1);
                            X1(i) = x(i);
                            X = [X X1];   
                        end
                        
                        lower_a = [c, -X];
                        upper_a = [c, X];
                    elseif isa(I, 'Zonotope')
                        [lb, ub] = I.getRanges();
                        lower_a = [I.c, I.X];
                        upper_a = [I.c, -I.X];
                    else
                        error('Unkown imput set');
                    end
                    
                    n_ub = size(ub,1);
                    nVar = dim * ones(n_ub,1);
                    
                    obj = RelaxedPoly(lower_a, upper_a, lb, ub, nVar);
                case 5
                    lower_a = varargin{1};
                    upper_a = varargin{2};
                    lb = varargin{3};
                    ub = varargin{4};
                    nVar = varargin{5};
                    
                    [n_nVar, m_nVar] = size(nVar);
                    [nL, mL] = size(lower_a);
                    [nU, mU] = size(upper_a);
                    
                    if mL ~= mU || nL ~= nU
                        error('sizes of matrix of upper and lower polyhedral constraints do not match');
                    elseif m_nVar ~= 1
                        error('nVar must be one column vector');   
                    end

                    obj.lower_a = lower_a;
                    obj.upper_a = upper_a;
                    obj.lb = lb;
                    obj.ub = ub;
                    obj.nVar = nVar;
                    obj.Dim = nVar(end);
            end
        end
        
        function R = affineMap(varargin)
            b = [];
            switch nargin
                case 3
                    obj = varargin{1};
                    W = varargin{2};
                    b = varargin{3};
                case 2
                    obj = varargin{1};
                    W = varargin{2};
                    b = zeros(size(W,1),1);
            end
            
            [nW, mW] = size(W);
            [nb, mb] = size(b);
            
            if mW ~= obj.Dim
                error('Inconsistency between affine transformation mattrix and object dimension');
            end
            
            if mb > 1
                error('bias vector must have one column');
            end
            
            if nW ~= nb
                error('inconsistency between affine transformation matrix and bias vector');
            end

            % update lower and upper polyhedral contraints
            lower_a = [b W];
            upper_a = [b W];
            lower_a = [obj.lower_a; lower_a];
            upper_a = [obj.upper_a; upper_a];
            
            % update number of variables of layer
            nVar = [obj.nVar; nW*ones(nW,1)];
            
            % get lower and upper bound using back substitution
            lb = obj.lb_backSub_while(lower_a, upper_a, nVar, inf);
            ub = obj.ub_backSub_while(lower_a, upper_a, nVar, inf);
            lb = [obj.lb; lb];
            ub = [obj.ub; ub];
            
            % create new relaxed polyhedron
            R = RelaxedPoly(lower_a, upper_a, lb, ub, nVar); 
        end

        function lb = lb_backSub_while(obj, lower_a, upper_a, nVar, maxIter)
            dim = nVar(end);
            [na, ma] = size(lower_a);
            alpha = lower_a(end-dim+1 : end, 2:end);
            lower_v = lower_a(na-dim+1 : na, 1); %zeros(dim,1);
            upper_v = zeros(dim,1); %upper_a(na-dim+1 : na, 1);
            
            na = na - dim;
            iter = 0;
            while(na > 2 && iter < maxIter)
                dim = nVar(na);
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);
                
                lower_v = max_a * lower_a(na-dim+1 : na, 1) + lower_v;
                upper_v = min_a * upper_a(na-dim+1 : na, 1) + upper_v;
                
                alpha = max_a * lower_a(na-dim+1 : na, 2:end) + ...
                        min_a * upper_a(na-dim+1 : na, 2:end);

                na = na-dim;
                iter = iter + 1;
            end
            
            dim = nVar(na);
            max_a = max(0, alpha);
            min_a = min(alpha, 0);
            
            if na == dim %iter <= maxIter && na == 0
                lb = max_a * lower_a(na-dim+1 : na, 2:end) * ones(2,1) + lower_v + ...
                     min_a * upper_a(na-dim+1 : na, 2:end) * ones(2,1) + upper_v;
            else
               lb1 = [obj.lb(na + 1 : na + dim)];
               ub1 = [obj.ub(na + 1 : na + dim)];
               
               lb = max_a * lb1 + lower_v + ...
                    min_a * ub1 + upper_v;
            end
        end
        
        function ub = ub_backSub_while(obj, lower_a, upper_a, nVar, maxIter)
            dim = nVar(end);
            [na, ma] = size(lower_a);
            alpha = upper_a(end-dim+1 : end, 2:end);
            lower_v = zeros(dim,1); %lower_a(na-dim+1 : na, 1);
            upper_v = upper_a(na-dim+1 : na, 1); %zeros(dim,1);
            
            na = na - dim;
            iter = 0;
            while(na > 2 && iter < maxIter)
                dim = nVar(na);
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);
                
                lower_v = min_a * lower_a(na-dim+1 : na, 1) + lower_v;
                upper_v = max_a * upper_a(na-dim+1 : na, 1) + upper_v; 
                
                alpha = min_a * lower_a(na-dim+1 : na, 2:end) + ...
                        max_a * upper_a(na-dim+1 : na, 2:end);

                na = na-dim;
                iter = iter + 1;
            end
            
            dim = nVar(na);
            max_a = max(0, alpha);
            min_a = min(alpha, 0);
            
            if na == dim              
                ub = min_a * lower_a(na-dim+1 : na, 2:end) * ones(2,1) + lower_v + ...
                     max_a * upper_a(na-dim+1 : na, 2:end) * ones(2,1) + upper_v;
            else
               ub1 = [obj.ub(na + 1 : na + dim)];
               lb1 = [obj.lb(na + 1 : na + dim)];
               
               ub = max_a * ub1 + lower_v + ...
                    min_a * lb1 + upper_v;
            end
        end
        
        
        function lb = lb_backSub_while_noIter(obj, lower_a, upper_a, nVar)
            dim = nVar(end);
            [na, ma] = size(lower_a);
            alpha = lower_a(end-dim+1 : end, 2:end); 
            lower_v = lower_a(na-dim+1 : na, 1); %zeros(dim,1);
            upper_v = zeros(dim,1); %upper_a(na-dim+1 : na, 1);
            
            na = na - dim;
            iter = 0;
            while(na > 0)
                dim = nVar(na);
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);
                
                lower_v = max_a * lower_a(na-dim+1 : na, 1) + lower_v;
                upper_v = min_a * upper_a(na-dim+1 : na, 1) + upper_v;
                
                alpha = max_a * lower_a(na-dim+1 : na, 2:end) + ...
                        min_a * upper_a(na-dim+1 : na, 2:end);

                na = na-dim;
                iter = iter + 1;
            end
            
            max_a = max(0, alpha);
            min_a = min(alpha, 0);
            
            lb = max_a * ones(2,1) + lower_v + ...
                 min_a * ones(2,1) + upper_v;
        end
        
        function ub = ub_backSub_while_noIter(obj, lower_a, upper_a, nVar)
            dim = nVar(end);
            [na, ma] = size(lower_a);
            alpha = upper_a(end-dim+1 : end, 2:end);
            lower_v = zeros(dim,1); %lower_a(na-dim+1 : na, 1);
            upper_v = upper_a(na-dim+1 : na, 1); %zeros(dim,1);
            
            na = na - dim;
            iter = 0;
            while(na > 0)
                dim = nVar(na);
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);
                
                lower_v = min_a * lower_a(na-dim+1 : na, 1) + lower_v;
                upper_v = max_a * upper_a(na-dim+1 : na, 1) + upper_v;
                
                
                alpha = min_a * lower_a(na-dim+1 : na, 2:end) + ...
                        max_a * upper_a(na-dim+1 : na, 2:end);

                na = na-dim;
                iter = iter + 1;
            end
            
            max_a = max(0, alpha);
            min_a = min(alpha, 0);
 
            ub = min_a * ones(2,1) + lower_v + ...
                 max_a * ones(2,1) + upper_v;
        end
        
        function [lb,ub] = getRange(obj, i)
            lb = obj.lb(end - obj.Dim + i);
            ub = obj.ub(end - obj.Dim + i);
        end
        
        function [lb,ub] = getRanges(obj)
           lb = obj.lb(end - obj.Dim + 1 : end);
           ub = obj.ub(end - obj.Dim + 1 : end);
        end
        
        function [lb,ub] = getAllRange(obj, i)
            lb = obj.lb(i);
            ub = obj.ub(i);
        end
        
        function [lb,ub] = getAllRanges(obj)
            lb = obj.lb;
            ub = obj.ub;
        end
        
        function P = toPolyhedron(obj)
            [lb, ub] = getRanges(obj);
            P = Polyhedron('lb', [lb], 'ub', [ub]);
        end
        
        function Z = toZonotope(obj)
            [lb,ub] = getRanges(obj);
            c = (ub + lb)*0.5;
            x = (ub - lb)*0.5;
            dim = obj.Dim;
            
            X = [];
            for i=1:dim
                X1 = zeros(dim, 1);
                X1(i) = x(i);
                X = [X X1];   
            end
            
            Z = Zonotope(c, X);
        end
        
        function S = toStar(obj)
            [lb,ub] = getRanges(obj);
            c = (ub + lb)*0.5;
            x = (ub - lb)*0.5;
            dim = obj.Dim;
            
            X = [];
            for i=1:dim
                X1 = zeros(dim, 1);
                X1(i) = x(i);
                X = [X X1];   
            end
            C = [eye(dim); -eye(dim)];
            d = ones(2*dim, 1); 
            
            S = Star([c X], C, d);
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
            
            if strcmp(color, 'rand')
                c_rand = rand(1,3);
            else
                c_rand = color;
            end
            
            P = toPolyhedron(obj);
            plot(P, 'color', c_rand);
        end
    end
end