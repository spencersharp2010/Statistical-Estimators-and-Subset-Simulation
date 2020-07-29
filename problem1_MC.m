clear;

%% setup
% define N such that coeff of variation is roughly 5%
N = 2e4;

% limit state function, derived by hand
g = @(X1,X2) -10*X1.^2 + 2*X1.*X2 + 2*X2.^2;

rng(1);

% inverted CDFs, used to generate samples from each distribution
cdf_x1_inv = @(u) 20./(2-u);
cdf_x2_inv = @(u) 20*sqrt(u);

max_iter = 1000;
%% run Monte Carlo simulation 1000 times
for i=1:max_iter
    u1 = rand(N,1);
    u2 = rand(N,1);
    
    x1 = cdf_x1_inv(u1);
    x2 = cdf_x2_inv(u2);
    
    g_eval = g(x1,x2);
    
    I = g_eval >= 0;
    Q(i) = sum(I)/N;
    Var_Pr(i)=var(I)/N;
    
    delta_Q(i)=sqrt(Var_Pr(i))/Q(i);
end

figure;
scatter(x1,x2);
title('x_1 and x_2 sampled using inverse transform')
xlabel('x_1')
ylabel('x_2')
% below are two methods of estimating the coefficient of variation
c_v1 = mean(delta_Q)
c_v2 = std(Q) / mean(Q)

% average probability of failure
Q_avg = mean(Q)

% since CI = 0.95, a = 1-0.95 = 0.05
a=0.05;

% calculate confidence intervals
ci_lower=Q_avg-norminv(1-a/2).*sqrt(mean(Var_Pr))
ci_upper=Q_avg+norminv(1-a/2).*sqrt(mean(Var_Pr))
