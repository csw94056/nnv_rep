close all;
clear;
clc;

W{1} = [1 1; 1 -1];
W{2} = [1 1; 1 -1];
W{3} = [1 1; 0 1];

b{1} = [0; 0];
b{2} = [0; 0];
b{3} = [1; 0];

lower_a{1} = [-1 0; 0 -1];
lower_a{2} = [0 0; 0 0];
lower_a{3} = [1 0; 0 0];

lower_b{1} = [0; 0];
lower_b{2} = [0; 0];
lower_b{3} = [1; 0];

upper_a{1} = [1 0; 0 1];
upper_a{2} = [0.5 0; 0 0.5];
upper_a{3} = [1 0; 0 0.5];

upper_b{1} = [0; 0];
upper_b{2} = [1; 1];
upper_b{3} = [1; 0];


%{
lower_beta = max(0,W{2}) * (max(0,lower_a{2}) * (max(0,W{1}) * max(0,lower_b{1}) + b{1}) + max(0,lower_b{2})) + b{2}
upper_beta = min(W{2},0) * (min(0,lower_a{2}) * (min(W{1},0) * min(0,upper_b{1}) + b{1}) + min(0,upper_b{2})) + b{2}
lower_m = max(0,W{2}) * max(lower_a{2},0) * max(0,W{1}) * min(lower_a{1},0)
upper_m = min(W{2},0) * min(0,upper_a{2}) * min(W{1},0) * min(upper_a{1},0)

lb_m = lower_m + upper_m
lb_b = lower_beta + upper_beta
%}

%{
lower_beta = max(0,W{2}) * (lower_a{2} * (W{1} * lower_b{1} + b{1}) + lower_b{2});
upper_beta = min(W{2},0) * (upper_a{2} * (W{1} * upper_b{1} + b{1}) + upper_b{2});
lower_m = max(0,W{2}) * lower_a{2} * W{1} * lower_a{1};
upper_m = min(W{2},0) * upper_a{2} * W{1} * upper_a{1};

lb_m = lower_m + upper_m;
lb_b = lower_beta + upper_beta + b{2};
lb = lb_m*[1; -1] + lb_b

lower_beta = min(W{2},0) * (lower_a{2} * (W{1} * lower_b{1} + b{1}) + lower_b{2});
upper_beta = max(0,W{2}) * (lower_a{2} * (W{1} * upper_b{1} + b{1}) + upper_b{2});
lower_m = min(W{2},0) * lower_a{2} * W{1} * lower_a{1};
upper_m = max(0,W{2}) * upper_a{2} * W{1} * upper_a{1};

ub_m = lower_m + upper_m ;
ub_b = lower_beta + upper_beta ;
ub = ub_m*[1; 1] + ub_b + b{2}
%}


%{
lower_beta = max(0,W{3}) * (lower_a{3} * (W{2} * (lower_a{2} * (W{1} * lower_b{1} + b{1}) + lower_b{2}) + b{2}) + lower_b{3});
upper_beta = min(W{3},0) * (upper_a{3} * (W{2} * (upper_a{2} * (W{1} * upper_b{1} + b{1}) + upper_b{2}) + b{2}) + upper_b{3});
lower_m = max(0,W{3}) * lower_a{3} * W{2} * lower_a{2} * W{1} * lower_a{1};
upper_m = min(W{3},0) * lower_a{3} * W{2} * upper_a{2} * W{1} * upper_a{1};

lb_m = lower_m + upper_m;
lb_b = lower_beta + upper_beta + b{3};
lb = -lb_m*ones(2,1) + lb_b;

lower_beta = min(W{3},0) * (lower_a{3} * (W{2} * (lower_a{2} * (W{1} * lower_b{1} + b{1}) + lower_b{2}) + b{2}) + lower_b{3})
upper_beta = max(0,W{3}) * (upper_a{3} * (W{2} * (upper_a{2} * (W{1} * upper_b{1} + b{1}) + upper_b{2}) + b{2}) + upper_b{3})
lower_m = min(W{3},0) * lower_a{3} * W{2} * lower_a{2} * W{1} * lower_a{1}
upper_m = max(0,W{3}) * upper_a{3} * W{2} * upper_a{2} * W{1} * upper_a{1}

ub_m = lower_m + upper_m
ub_b = lower_beta + upper_beta + b{3}
ub = ub_m*ones(2,1) + ub_b
%}

lower_beta = max(0,W{3}) * (lower_a{3} * (W{2} * (lower_a{2} * (W{1} * lower_b{1} + b{1}) + lower_b{2}) + b{2}) + lower_b{3});
upper_beta = min(W{3},0) * (upper_a{3} * (W{2} * (upper_a{2} * (W{1} * upper_b{1} + b{1}) + upper_b{2}) + b{2}) + upper_b{3});
lower_m = max(0,W{3}) * lower_a{3} * W{2} * lower_a{2} * W{1} * lower_a{1};
upper_m = min(W{3},0) * lower_a{3} * W{2} * upper_a{2} * W{1} * upper_a{1};

lb_m = lower_m + upper_m;
lb_b = lower_beta + upper_beta + b{3};
lb = -lb_m*ones(2,1) + lb_b;

disp('original')

lower_beta = min(W{3},0) * (lower_a{3} * (W{2} * (lower_a{2} * (W{1} * lower_b{1} + b{1}) + lower_b{2}) + b{2}) + lower_b{3})
upper_beta = max(0,W{3}) * (upper_a{3} * (W{2} * (upper_a{2} * (W{1} * upper_b{1} + b{1}) + upper_b{2}) + b{2}) + upper_b{3})
lower_m = min(W{3},0) * lower_a{3} * W{2} * lower_a{2} * W{1} * lower_a{1}
upper_m = max(0,W{3}) * upper_a{3} * W{2} * upper_a{2} * W{1} * upper_a{1}
ub_m = lower_m + upper_m
ub_b = lower_beta + upper_beta + b{3}
ub = ub_m*ones(2,1) + ub_b

disp('all max or min')

lower_beta = min(W{3},0) * (min(lower_a{3},0) * (min(W{2},0) * (min(lower_a{2},0) * (min(W{1},0) * lower_b{1} + b{1}) + lower_b{2}) + b{2}) + lower_b{3})
upper_beta = max(0,W{3}) * (max(0,upper_a{3}) * (max(0,W{2}) * (max(0,upper_a{2}) * (max(0,W{1}) * upper_b{1} + b{1}) + upper_b{2}) + b{2}) + upper_b{3})
lower_m = min(W{3},0) * min(lower_a{3},0) * min(W{2},0) * min(lower_a{2},0) * min(W{1},0) * min(lower_a{1},0)
upper_m = max(0,W{3}) * max(upper_a{3},0) * max(W{2},0) * max(0,upper_a{2}) * max(0,W{1}) * max(0,upper_a{1})

ub_m = lower_m + upper_m
ub_b = lower_beta + upper_beta + b{3}
ub = ub_m*[1;1] + ub_b

disp('partial')
lower_beta = min(W{3},0) * (lower_a{3} * (min(W{2},0) * (lower_a{2} * (min(W{1},0) * lower_b{1} + b{1}) + lower_b{2}) + b{2}) + lower_b{3})
upper_beta = max(0,W{3}) * (upper_a{3} * (max(0,W{2}) * (upper_a{2} * (max(0,W{1}) * upper_b{1} + b{1}) + upper_b{2}) + b{2}) + upper_b{3})
lower_m = min(W{3},0) * lower_a{3} * min(W{2},0) * lower_a{2} * min(W{1},0) * lower_a{1}
upper_m = max(0,W{3}) * upper_a{3} * max(W{2},0) * upper_a{2} * max(0,W{1}) * upper_a{1}

ub_m = lower_m + upper_m
ub_b = lower_beta + upper_beta + b{3}
ub = ub_m*[1;1] + ub_b

%{
lower_m = eye(2);
upper_m = eye(2)
beta = zeros(2,1);
for i = 1:3
    if i == 1
        lower_m = max(0,W{i})*lower_a{i}*lower_m;
        upper_m = max(0,W{i})*upper_a{i}*upper_m;
        beta = max(0,W{i})*(beta + lower_b{i}) + b{i};
    else
        lower_m = W{i}*lower_a{i}*lower_m;
        upper_m = W{i}*upper_a{i}*upper_m;
        beta = W{i}*(beta + lower_b{i}) + b{i};
    end
end

alpha = alpha
beta = beta
%}