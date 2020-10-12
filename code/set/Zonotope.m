classdef Zonotope
    properties
        c = []; % center coefficient
        X = []; % coefficient matrix; partial deviation around the center
        Dim = 0;
    end
    methods
        function obj = Zonotope(varargin)
            switch nargin
                case 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0
                    obj.c = [];
                    obj.X = [];
                    obj.Dim = 0;
                case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1
                    I = varargin{1};
                    if isa(I, 'Polyhedron')
                        I.outerApprox;
                        lb = I.Internal.lb;
                        ub = I.Internal.ub;
                    elseif isa(I, 'Star')
                        [lb, ub] = getRanges(I);
                    else
                        error('Unkown object');
                    end
                    
                    dim = I.Dim;
                    c = (ub + lb)*0.5;
                    x = (ub - lb)*0.5;

                    X = [];
                    for i=1:dim
                        X1 = zeros(dim, 1);
                        X1(i) = x(i);
                        X = [X X1];   
                    end

                    obj.c = c;
                    obj.X = X;
                    obj.Dim = dim;
                case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3
                    c = varargin{1};
                    X = varargin{2};
                    
                    if size(c,1) ~= size(X,1)
                        error('Inconsistent dimension between c and X')
                    end
                    
                    obj.c = c;
                    obj.X = X;
                    obj.Dim = size(X,1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end switch    
            end
        end
        
        function Z = affineMap(varargin)
            switch nargin
                case 3
                    obj = varargin{1};
                    W = varargin{2};
                    b = varargin{3};
                case 2
                    obj = varargin{1};
                    W = varargin{2};
                    b = [];
            end
            
            [nW, mW] = size(W);
            [nb, mb] = size(b);
            
            if mW ~= obj.Dim
                error('Inconsistency between affine transformation mattrix and object dimension');
            end
            
            if mb > 1
                error('bias vector must have one column');
            end
            
            if mb ~= 0
                c = W * obj.c + b;
            else
                c = W * obj.c;
            end
            
            X = W * obj.X;
            Z = Zonotope(c, X);
        end
        
        function [lb, ub] = getRanges(obj)
            ub = transpose(max(obj.X', [], 1));
            lb = -ub;
            
            ub = ub + obj.c;
            lb = lb + obj.c;
        end
        
        function [lb, ub] = getRange(obj, i)
            ub = transpose(max(obj.X(i,:)', [], 1));
            lb = -ub;
            
            ub = ub + obj.c(i);
            lb = lb + obj.c(i);
        end
        
        function P = toPolyhedron(obj)
            %{
            A = [eye(obj.Dim); -eye(obj.Dim)]
            b = ones(2*obj.Dim,1)
            P1 = Polyhedron('A', A, 'b', b);
            c = obj.c
            X = obj.X
            size(obj.X, 2)
            %P = obj.X' * P1 + obj.c;
            P = P1.affineMap(obj.X) + obj.c;
            %}

            n = size(obj.X, 2);
            A = [eye(n); -eye(n)];
            b = ones(2*n,1);
            P1 = Polyhedron('A', A, 'b', b);
            P = obj.X*P1 + obj.c;
        end
        
        function plot(varargin)
            color = 'red';
                   
            switch nargin
                case 1
                    obj = varargin{1};
                case 2
                    obj = varargin{1};
                    color = varargin{2};
            end
            
            if strcmp(color, 'rand')
                color = rand(1,3);
            end
            
            P = toPolyhedron(obj);
            plot(P, 'color', color);
        end
    end
end
