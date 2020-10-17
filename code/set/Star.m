classdef Star
    properties
        V = []; % a set of basic vectors
        C = []; % predicate constraint matrix
        d = []; % predicate constraint vector
        Dim = 0;
        
        predicate_lb = []; % lower bound vector of predicate variable
        predicate_ub = []; % upper bound vector of predicate variable
    end
    methods
        function obj = Star(varargin)
            switch nargin
                case 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0
                    obj.V = [];
                    obj.C = [];
                    obj.d = [];
                    obj.Dim = 0;
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
                    obj = Star(V, C, d);
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
                    %{
                    P = Polyhedron('A', obj.C, 'b', obj.d);
                    P.outerApprox;
                    obj.predicate_lb = P.Internal.lb;
                    obj.predicate_ub = P.Internal.ub;
                    %}
                case 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4
                    c = varargin{1};
                    V = varargin{2};
                    C = varargin{3};
                    d = varargin{4};
                    obj = Star([c V], C, d);
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
            S = Star(V, obj.C, obj.d);
        end
        
        function [lb, ub] = getRanges(obj)
            P = obj.toPolyhedron;
            P.outerApprox;
            lb = P.Internal.lb;
            ub = P.Internal.ub;
        end
        
        function V = get_V(obj)
            %V = obj.V(:,2:size(obj.V,2));
            V = obj.V(:, 2:end);
        end
        
        function c = get_c(obj)
            c = obj.V(:,1);
        end
        
        function P = toPolyhedron(obj)
            P1 = Polyhedron('A', [obj.C], 'b', [obj.d]);
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
    end
end
