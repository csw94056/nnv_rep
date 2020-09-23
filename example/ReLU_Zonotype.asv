function R = zonoDecApproxStepReLU_star(I, i)
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
                % y[i] = ReLU(x[i]) = a_(m+1)
                V = I.get_V;
                c = I.get_c;
                
                lamda = 
                
                % constraint 1: y
                
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

                R = myStar(star_V, new_C, new_d);
            end
        end