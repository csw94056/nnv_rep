classdef myStar
    properties
        V = []; % a set of basic vectors
        C = []; % predicate constraint matrix
        d = []; % predicate constraint vector
        Dim = 0;
        
        predicate_lb = []; % lower bound vector of predicate variable
        predicate_ub = []; % upper bound vector of predicate variable
    end
    methods
        function obj = myStar(varargin)
            switch nargin
                case 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0
                    obj.V = 0;
                    obj.C = 0;
                    obj.d = 0;
                    obj.Dim = 0;
                case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1
                    P = varargin{1};
                    if ~isa(P, 'Polyhedron')
                        error('Input set is not a polyhedron object');
                    end
                    c = zeros(P.Dim,1);
                    V = eye(P.Dim);
                    obj.V = [c V];
                    if ~isempty(P.Ae)
                        obj.C = [P.A;P.Ae;-P.Ae];
                        obj.d = [P.b;P.be;-P.be];
                    else
                        obj.C = P.A;
                        obj.d = P.b;
                    end
                    obj.Dim = P.Dim;
                case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3
                    %V need to contain c and V
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    obj.V = V;
                    obj.C = C;
                    obj.d = d;
                    obj.Dim = size(V,1);
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

                    obj.V = [c V];
                    obj.C = C;
                    obj.d = d;
                    obj.Dim = size(V,2);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end switch    
            end
        end
            
        function S = affineMap(varargin)
            switch nargin
                case 3
                    obj = varargin{1};
                    W = varargin{2};
                    b = varargin{3};
                    V = W * obj.V ;
                    V(:, 1) = V(:, 1) + b;
                    S = myStar(V, obj.C, obj.d);
                case 2
                    obj = varargin{1};
                    W = varargin{2};
                    V = W * obj.V;
                    S = myStar(V, obj.C, obj.d);
            end
        end
        
        function [lb, ub] = getRanges(obj)
            P =  toPoly(obj);
            P.outerApprox;
            lb = P.Internal.lb;
            ub = P.Internal.ub;
        end
        
        function A = A(obj)
           A = obj.C; 
        end
        
        function b = b(obj)
            b = obj.d;
        end
        
        function V = get_V(obj)
            obj_V = obj.V;
            V = obj.V(:,2:size(obj.V,2));
        end
        
        function c = get_c(obj)
            c = obj.V(:,1);
        end
        
        function P = toPoly(obj)
            P1 = Polyhedron('A', obj.C, 'b', obj.d);
            V = obj.get_V;
            c = obj.get_c;
            P = V*P1 + c;
        end
        
        function plot(obj)
            n = length(obj);
            for i=1:n
                I = obj(i);
                hold on
                P = toPoly(obj);
                if n>1
                    plot(P, 'color', rand(1,3));
                else
                    plot(P);
                end
            end
            hold off
        end
        
        function S = Intersect(obj1, obj2)         
          C1 = obj2.C * obj1.get_V;
          d1 = obj2.d - obj2.C * obj1.get_c;
          
          new_C = [obj1.C; C1];
          new_d = [obj1.d; d1];
          S = myStar(obj1.V, new_C, new_d);     
          %if isEmptySet(S)
          %    S = [];
          %end
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
