classdef RlxPoly
    %RLXPOLY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lower_a = {[]}; % a set of matrix for lower polyhedral constraint (a[<=]) ((1 a[1] a[2] ... a[n])'
        upper_a = {[]}; % a set of matrix for upper polyhedral constraint (a[>=])
        lb = []; % a set of matrix for lower bound
        ub = []; % a set of matrix for upper bound

        iter = inf; % number of iterations for back substitution
        Dim = 0; % dimension of current relaxed polyhedron
    end
    
    methods
        function obj = RlxPoly(varargin)
            switch nargin
                case 0
                    obj.lower_a = {[]};
                    obj.upper_a = {[]};
                    obj.lb = [];
                    obj.ub = [];
                    obj.iter = inf;
                    obj.Dim = 0;
                case 2
                    I = varargin{1};
                    iter = varargin{2};
                    dim = I.Dim;
                    if isa(I, 'Polyhedron')
                        I.outerApprox;
                        l = I.Internal.lb;
                        u = I.Internal.ub;
                        
                        L = zeros(dim, 1);
                        U = zeros(dim, 1);
                        for i=1:dim
                            L1 = zeros(dim, 1); 
                            U1 = zeros(dim, 1);
                            L1(i) = l(i);
                            U1(i) = u(i);
                            L = [L L1];
                            U = [U U1];
                        end
                    elseif isa(I, 'Star')
                        [l, u] = I.getRanges();
                        
                        c = (u + l)*0.5;
                        x = (u - l)*0.5;

                        X = [];
                        for i=1:dim
                            X1 = zeros(dim, 1);
                            X1(i) = x(i);
                            X = [X X1];   
                        end
                        
                        L = [c, -X];
                        U = [c, X];
                    elseif isa(I, 'Zonotope')
                        [l, u] = I.getRanges();
                        L = [I.c, I.X];
                        U = [I.c, -I.X];
                    else
                        error('Unkown imput set');
                    end
                    
                    if iter <= 0
                        error('Iteration must be greater than zero');
                    end
                    
                    lower_a{1} = L;
                    upper_a{1} = U;
                    lb{1} = l;
                    ub{1} = u;
                    obj = RlxPoly(lower_a, upper_a, lb, ub, iter);                   
                case 5
                    lower_a = varargin{1};
                    upper_a = varargin{2};
                    lb = varargin{3};
                    ub = varargin{4};
                    iter = varargin{5};
                    
%                     [n_nVar, m_nVar] = size(nVar);
%                     [nL, mL] = size(lower_a);
%                     [nU, mU] = size(upper_a);
%                     
%                     if mL ~= mU || nL ~= nU
%                         error('sizes of matrix of upper and lower polyhedral constraints do not match');
%                     elseif m_nVar ~= 1
%                         error('nVar must be one column vector');   
%                     end
                    
                    obj.lower_a = lower_a;
                    obj.upper_a = upper_a;
                    obj.lb = lb;
                    obj.ub = ub;
                    
                    obj.iter = iter;
                    n = length(lower_a);
                    [nL, mL] = size(lower_a{n});
                    obj.Dim = nL;
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

            % new lower and upper polyhedral contraints
            lower_a = obj.lower_a;
            upper_a = obj.upper_a;
            len = length(lower_a);
            lower_a{len+1} = [b W];
            upper_a{len+1} = [b W];
            % new lower and uppper bounds
            lb = obj.lb;
            ub = obj.ub;
            lb{len+1} = obj.lb_backSub(lower_a, upper_a);
            ub{len+1} = obj.ub_backSub(lower_a, upper_a);
            
%             if len + 1 > obj.iter
%                 lower_a = {lower_a{2:end}};
%                 upper_a = {upper_a{2:end}};
%                 lb = {lb{2:end}};
%                 ub = {ub{2:end}};
%             end
                     
            % create new relaxed polyhedron
            R = RlxPoly(lower_a, upper_a, lb, ub, obj.iter); 
        end

        function lb = lb_backSub(obj, lower_a, upper_a)
            maxIter = obj.iter;
            len = length(lower_a);
            [nL, mL] = size(lower_a{len});
            alpha = lower_a{len}(:,2:end);
            lower_v = lower_a{len}(:,1);
            upper_v = zeros(nL, 1);
            % b[s+1] = v' + sum( max(0,w[j]')*lower_a[j] + min(w[j]',0)*upper_a[j}] ) for j is element of k and for k < i
            % iteration until lb' = b[s'] = v''
            len = len - 1;
            iter = 0;
            while (len > 1 && iter < maxIter)
                [nL, mL] = size(lower_a{len});
                dim = nL;
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);

                lower_v = max_a * lower_a{len}(:,1) + lower_v;
                upper_v = min_a * upper_a{len}(:,1) + upper_v;

                alpha = max_a * lower_a{len}(:,2:end) + ...
                        min_a * upper_a{len}(:,2:end);

                len = len - 1;
                iter = iter + 1;
            end
            
            max_a = max(0, alpha);
            min_a = min(alpha, 0);
            
            [lb1,ub1] = getRanges_L(obj,len);
            lb = max_a * lb1 + lower_v + ...
                 min_a * ub1 + upper_v;
        end
        
        function ub = ub_backSub(obj, lower_a, upper_a)
            maxIter = obj.iter;
            len = length(lower_a);
            [nL, mL] = size(lower_a{len});
            alpha = upper_a{len}(:,2:end);
            lower_v = zeros(nL, 1);
            upper_v = upper_a{len}(:,1);
            % c[t+1] = v' + sum( max(0,w[j]')*upper_a[j] + min(w[j]',0)*lower_a[j}] )  for j is element of k and for k < i
            % iteration until ub' = c[t'] = v''
            len = len - 1;
            iter = 1;
            while (len > 1 && iter < maxIter)
                [nL, mL] = size(lower_a{len});
                dim = nL;
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);
                
                lower_v = min_a * lower_a{len}(:,1) + lower_v;
                upper_v = max_a * upper_a{len}(:,1) + upper_v;
                
                alpha = min_a * lower_a{len}(:,2:end) + ...
                        max_a * upper_a{len}(:,2:end);

                len = len - 1;
                iter = iter + 1;
            end

            max_a = max(0, alpha);
            min_a = min(alpha, 0);
            
            [lb1,ub1] = getRanges_L(obj,len);
            ub = min_a * lb1 + lower_v + ...
                 max_a * ub1 + upper_v;
        end

%         function lb = lb_backSub(obj, lower_a, upper_a)
%             maxIter = obj.iter;
%             len = length(lower_a)
%             [nL, mL] = size(lower_a{len});
%             alpha = lower_a{len}(:,2:end)
%             lower_v = lower_a{len}(:,1)
%             upper_v = zeros(nL, 1)
% 
%             len = len - 1;
%             iter = 0;
%             while (len > 1 && iter < maxIter)
%                 [nL, mL] = size(lower_a{len});
%                 dim = nL;
%                 
%                 max_a = max(0, alpha)
%                 min_a = min(alpha, 0)
%                 
%                 prev_lower_v = lower_a{len}(:,1)
%                 prev_upper_v = upper_a{len}(:,1)
%                 lower_v = max_a * lower_a{len}(:,1) + lower_v
%                 upper_v = min_a * upper_a{len}(:,1) + upper_v
%                 
%                 prev_lower_a = lower_a{len}(:,2:end)
%                 prev_upper_a = upper_a{len}(:,2:end)
%                 alpha = max_a * lower_a{len}(:,2:end) + ...
%                         min_a * upper_a{len}(:,2:end)
% 
%                 len = len - 1;
%                 iter = iter + 1;
%             end
%             
%             max_a = max(0, alpha)
%             min_a = min(alpha, 0)
%             
%             if iter < maxIter
%                 [nL, mL] = size(lower_a{1});
%                 dim = nL;
%                 
%                 lb = max_a * lower_a{len}(:,2:end) * ones(dim,1) + lower_v + ...
%                      min_a * upper_a{len}(:,2:end) * ones(dim,1) + upper_v
%             else
%                 [lb1,ub1] = getRanges_L(obj,len);
%                 lb = max_a * lb1 + lower_v + ...
%                      min_a * ub1 + upper_v;
%             end
%         end
%         
%         function ub = ub_backSub(obj, lower_a, upper_a)
%             maxIter = obj.iter;
%             len = length(lower_a);
%             [nL, mL] = size(lower_a{len});
%             alpha = upper_a{len}(:,2:end);
%             lower_v = zeros(nL, 1);
%             upper_v = upper_a{len}(:,1);
% 
%             len = len - 1;
%             iter = 1;
%             while (len > 1 && iter < maxIter)
%                 [nL, mL] = size(lower_a{len});
%                 dim = nL;
%                 
%                 max_a = max(0, alpha);
%                 min_a = min(alpha, 0);
%                 
%                 lower_v = min_a * lower_a{len}(:,1) + lower_v;
%                 upper_v = max_a * upper_a{len}(:,1) + upper_v;
%                 
%                 alpha = min_a * lower_a{len}(:,2:end) + ...
%                         max_a * upper_a{len}(:,2:end);
% 
%                 len = len - 1;
%                 iter = iter + 1;
%             end
% 
%             max_a = max(0, alpha);
%             min_a = min(alpha, 0);
%             
%             if iter < maxIter
%                 [nL, mL] = size(lower_a{1});
%                 dim = nL;
%                 
%                 ub = min_a * lower_a{len}(:,2:end) * ones(dim,1) + lower_v + ...
%                      max_a * upper_a{len}(:,2:end) * ones(dim,1) + upper_v;
%             else
%                 [lb1,ub1] = getRanges_L(obj,len);
%                 ub = min_a * lb1 + lower_v + ...
%                      max_a * ub1 + upper_v;
%             end
%         end
        
        function [lb,ub] = getRange(obj, i)
            if i > obj.Dim
                error('i should not exceed dimnesion');
            end
            
            len = length(obj.lb);
            lb = obj.lb{len}(i);
            ub = obj.ub{len}(i);
        end
        
        function [lb,ub] = getRanges(obj)
            len = length(obj.lb);
            lb = obj.lb{len};
            ub = obj.ub{len};
        end
        
        function [lb,ub] = getRanges_L(obj, len)
            % get lower and upper bound of specific layer
            numL = length(obj.lower_a);
            if len > numL
                error('range request should be layers within iteration');
            end
            lb = obj.lb{len};
            ub = obj.ub{len};
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