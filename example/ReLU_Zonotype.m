%close all;
clear;
clc;

P = ExamplePoly.randHrep('d',2);
%P = Polyhedron('lb',[-1;-1], 'ub', [1;1]);
S = Star(P);
Z = Zonotope(P);


Rs = ReLU.approxStepReLU_star(S,1);
Zd = zonoDecApproxStepReLU_star(S,1);
Zi = zonoIncApproxStepReLU_star(S,1);
Rz = ReLU.approxStepReLU_zono(Z,1);
%Rz_star = approxStepReLU_zono_star(S, 1);
Rabs = absDomainStepReLU_star(S, 1);
%{
plot(S);
nexttile
hold on;
plot(Rs,'g');
hold on;
plot(Zd, 'y');
hold on;
plot(Zi, 'b');
%}

figure;
ReLU.plot_set(S,'Input set');
ReLU.plot_set(Rs,'R star approx');
ReLU.plot_set(Zd,'R star zono dec');
ReLU.plot_set(Zi,'R star zono inc');
ReLU.plot_set(Rz,'R zono approx');

ReLU.plot_set(Rabs,'R abs domain');


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
        % y[i] = ReLU(x[i]) = a[m+1]
        V = I.get_V;
        c = I.get_c;
        lamda = ub/(ub - lb);

        % constraint 1: y[i] >= lamda*x[i]
        C1 = [lamda*V(i,:) -1];
        d1 = -lamda*c(i);

        % constraint 2: y[i] <= lamda*x[i] + ub*(1-lamda)
        C2 = [-lamda*V(i,:) 1];
        d2 = -lamda*c(i) + ub*(1-lamda);

        C0 = [I.C zeros(size(I.C,1),1)];
        d0 = I.d;

        new_C = [C0; C1; C2];
        new_d = [d0; d1; d2];

        star_V = I.V;
        star_V(i,:) = 0;
        e = zeros(I.Dim,1);
        e(i,1) = 1;
        star_V = [star_V e];

        R = Star(star_V, new_C, new_d);
    end
end


function R = zonoIncApproxStepReLU_star(I, i)
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
        R = Star(star_V, I.C, I.d);
    elseif lb<0 && ub>0
        % y[i] = ReLU(x[i]) = a[m+1]
        V = I.get_V;
        c = I.get_c;
        lamda = ub/(ub - lb);

        % constraint 1: y[i] >= lamda*x[i]
        C1 = [lamda*V(i,:) -1];
        d1 = -lamda*c(i);

        % constraint 2: y[i] <= lamda*(x[i] - lb)
        C2 = [-lamda*V(i,:) 1];
        d2 = lamda*(c(i)-lb);

        C0 = [I.C zeros(size(I.C,1),1)];
        d0 = I.d;

        new_C = [C0; C1; C2];
        new_d = [d0; d1; d2];

        star_V = I.V;
        star_V(i,:) = 0;
        e = zeros(I.Dim,1);
        e(i,1) = 1;
        star_V = [star_V e];

        R = Star(star_V, new_C, new_d);
    end
end

%{
function R = approxStepReLU_zono(I, i)

    if ~isa(I, 'Zonotope')
        error('Input set is not a Zonotope object');
    end
                    
    [lb, ub] = getRanges(I);
    
    lb = lb(i);
    ub = ub(i);
    
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
%}
function R = approxStepReLU_zono_star(I, i)
    [lb, ub] = getRanges(I);
    
    lb = lb(i);
    ub = ub(i);
    
    if lb > 0
        R = I;
    elseif ub <= 0
        c = I.get_c;
        X = I.get_V;
        
        c(i) = 0;
        X(i,:) = 0;
        R = Zonotope(c, X);
    else
        c = I.get_c;
        X = I.get_V;

        lamda = ub/(ub-lb);
        mu = -ub*lb/(2*(ub-lb));
        
        c(i) = lamda*c(i) + mu;
        X(i,:) = lamda*X(i,:);
        new_x = zeros(I.Dim,1);
        new_x(i) = mu;
        X = [X new_x];
        
        
        R = Star(c,X,I.C,I.d);
    end 
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

        % constraint 1: y[i] >= 0
        m = size(I.C,2);
        C1 = [zeros(1,m) -1];
        d1 = 0;
        
        % constraint 2: y[i] <= ub(x[i]-lb)/(ub-lb)
        C2 = [-ub*V(i,:)/(ub-lb) 1];
        d2 = ub*(c(i)-lb)/(ub-lb);

        C0 = [I.C zeros(size(I.C,1),1)];
        d0 = I.d;

        new_C = [C0; C1; C2];
        new_d = [d0; d1; d2];

        star_V = I.V;
        star_V(i,:) = 0;
        e = zeros(I.Dim,1);
        e(i,1) = 1;
        star_V = [star_V e];

        R = Star(star_V, new_C, new_d);
    end
end
