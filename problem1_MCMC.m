clear;
%% setup
% limit state function
g = @(X) -10*X(1).^2 + 2*X(1).*X(2) + 2*X(2).^2;

% joint pdf (derived by hand)
jointpdf = @(x) (10<=x(1) & x(1)<=20).*(0<=x(2) & x(2)<=20) ...
    .* (x(2)/(10*x(1)^2));

% Length of burn-in period
burn=1000;

% Sample size - same as MC portion of problem 1
N=2e4;

% % Seed value
% rng(1);

%% IMPORTANT: Select one of the following options by uncommenting it
%sigma_val=[2 5 10];
sigma_val = 5;

%% MCMC - random walk sampler
for iter=1:length(sigma_val)
    % Seed value
    rng(1);
    
    fprintf('sigma = %.1f\n',sigma_val(iter));
    x=zeros(burn+N+1,2);
    x(1,:)=[15,10];

    % Candidate sample
    sigma_mat = eye(2)*sigma_val(iter);
    %x_cand = mvnrnd([0,0],sigma_mat,burn+N);
    x_cand = normrnd(0,sigma_val(iter),burn+N,2);
    
    % Uniform sample
    u=rand(burn+N,1);

    % MCMC algorithm
    for j=1:burn+N

        % Generate proposal state of the Markov chain
        y=x(j,:)+x_cand(j,:);

        % evaluate joint pdf using proposal state
        eval = jointpdf(y)/jointpdf(x(j,:));
        alpha = min(1, eval );
        
        % Accept or reject sample by comparing random uniform sample to
        % alpha, store acceptance state
        if u(j)<=alpha
            x(j+1,:)=y;
            acceptance(j)=1;
        else
            x(j+1,:)=x(j,:);
            acceptance(j)=0;
        end

    end

    % Discard the samples and acceptance values from the burn-in period
    x_fin=x(burn+2:end,:);
    acceptance_fin=acceptance(burn+1:end);
    
    % evaluate g using final samples generated from random walk sampler
    for k=1:N
        g_eval(k) = g(x_fin(k,:));
    end
    
    % store x1 and x2 separately (used for plotting) to avoid 3D matrix
    x1_graph(:,iter) = x_fin(:,1);
    x2_graph(:,iter) = x_fin(:,2);
    
    % calculate acceptance rate for current iteration
    acceptance_mean(iter) = mean(acceptance_fin);
    
    % calculate probability that g >= 0 through simple summation divided by
    % number of samples
    I = g_eval >= 0;
    Q = sum(I)/N;

    % estimate the coefficient of variation of Q
    nlags = 50;
    [acf_Y, lags_Y] = autocorr(I',nlags);
    gamma = 2*sum((1-(1:nlags)/(length(I)))'.*acf_Y(2:end));
    Var_Q(iter) = var(I)/N*(1+gamma);
    delta_Q(iter)=sqrt(Var_Q(iter))/Q;
    
    % print results to console
    fprintf('\tP(g(X1,X2)>=0) = : %.4f\n',Q);
    fprintf('\tcoefficient of variation: %.4f\n',delta_Q(iter));
    fprintf('\tacceptance rate: %.4f\n',acceptance_mean(iter));
end

% plot evolution of sample generation 
figure;
scatter(x(:,1),x(:,2))
title("random walk sample evolution, start: [" + x(1,1) + ", " + x(1,2) + "]")
xlabel('x_1')
ylabel('x_2')

% plot autocorrelation of x1 samples
figure
hold on;
for i=1:length(sigma_val)
    [acf,Lags] = autocorr(x1_graph(:,i),50);
    plot(Lags,acf);
    legendinfo{i}=['\sigma = ', num2str(sigma_val(i))];
    %hold on;
end
title('x_1 sample correlation')
xlabel('Lag')
ylabel('Sample Autocorrelation')
legend(legendinfo)
hold off;

% plot autocorrelation of x2 samples using 
figure
hold on;
for i=1:length(sigma_val)
    [acf,Lags] = autocorr(x2_graph(:,i),50);
    plot(Lags,acf);
    legendinfo{i}=['\sigma = ', num2str(sigma_val(i))];
    %hold on;
end
title('x_2 sample correlation')
xlabel('Lag')
ylabel('Sample Autocorrelation')
legend(legendinfo)
hold off;
