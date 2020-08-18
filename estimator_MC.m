clear;
%% setup
% define N such that coeff of variation is roughly 5%
N = 2e4;

% limit state function, derived by hand
g = @(X1,X2) -10*X1.^2 + 2*X1.*X2 + 2*X2.^2;

% set the seed
rng(1);

% inverted CDFs, used to generate samples from each distribution
cdf_x1_inv = @(u) 20./(2-u);
cdf_x2_inv = @(u) 20*sqrt(u);

max_iter = 1000;
%% run Monte Carlo simulation 1000 times
for i=1:max_iter
    % generate uniformly distributed samples in standard normal space
    u1 = rand(N,1);
    u2 = rand(N,1);
    
    % convert to physical space using inverse transform method by plugging
    % into the hand-derived inverted CDFs
    x1 = cdf_x1_inv(u1);
    x2 = cdf_x2_inv(u2);
    
    % evaluate g for each pair of x1, x2
    g_eval = g(x1,x2);
    
    % count how many times g exceeds 0
    I = g_eval >= 0;
    
    % Q is determined by dividing number of times g exceeds 0 by number of
    % sample evaluations
    Q(i) = sum(I)/N;
    
    % calculate the variance of the probability estimate
    Var_Pr(i)=var(I)/N;
    
    % calculate the coefficient of variation
    delta_Q(i)=sqrt(Var_Pr(i))/Q(i);
end

% plot the generated x1 and x2 values
figure;
scatter(x1,x2);
title('x_1 and x_2 sampled using inverse transform')
xlabel('x_1')
ylabel('x_2')

% below are two methods of estimating the coefficient of variation
c_v1 = mean(delta_Q);
c_v2 = std(Q) / mean(Q);

% average probability of failure
Q_avg = mean(Q);

% since CI = 0.95, a = 1-0.95 = 0.05
a=0.05;

% calculate confidence intervals
ci_lower=Q_avg-norminv(1-a/2).*sqrt(mean(Var_Pr));
ci_upper=Q_avg+norminv(1-a/2).*sqrt(mean(Var_Pr));

% print results to console
fprintf('Q_avg = %.4f \n', Q_avg);
fprintf('CoV = %.4f \n', c_v1);
fprintf('CI_lower = %.4f \n', ci_lower);
fprintf('CI_upper = %.4f \n', ci_upper);

