clear;

% limit state function
g = @(X) -10*X(1).^2 + 2*X(1).*X(2) + 2*X(2).^2;

% joint pdf
jointpdf = @(x) (10<=x(1) & x(1)<=20).*(0<=x(2) & x(2)<=20) ...
    .* (x(2)/(10*x(1)^2));

% pdf_x1 = @(x1) 20/(x1^2);
% pdf_x2 = @(x2) x2/200;

% parameter sigma
sigma_val=[0.1 0.5 1 5 10 50];
sigma = 5;

% Length of burn-in period
burn=1000;

% Sample size
N=2e4;

% Seed value
rng(1,'twister');


for iter=1:length(sigma_val)
    
    x=zeros(burn+N+1,2);
    x(1,:)=[15,10];

    % Candidate sample
    sigma_mat = eye(2)*sigma_val(iter);
    x_cand = mvnrnd([0,0],sigma_mat,N+burn);
    % Uniform sample
    u=rand(burn+N,1);

    % MCMC algorithm
    for j=1:burn+N

        % Generate proposal state of the Markov chain
        y=x(j,:)+x_cand(j,:);

        % Acceptance probability
        %alpha=min(1,lognpdf(y,mu_log,si_log)/lognpdf(x(j),mu_log,si_log));
        eval = jointpdf(y)/jointpdf(x(j,:));
        alpha = min(1, eval );
        % Accept or reject sample
        if u(j)<=alpha
            x(j+1,:)=y;
            acceptance(j)=1;
        else
            x(j+1,:)=x(j,:);
            acceptance(j)=0;
        end

    end

    % Discard the samples from the burn-in period
    x_fin=x(burn:end,:);
    acceptance_fin=acceptance(burn+2:end);
    for k=1:N
        g_eval(k) = g(x_fin(k,:));
    end
    
    x_graph(:,iter) = x_fin(:,1);
    acceptance_mean(iter) = mean(acceptance_fin);
    I = g_eval >= 0;
    Q = sum(I)/N
    Var_Pr=var(I)/N;

    % Coefficient of variation of the estimator
    % estimate the coefficient of variation of Q
    nlags = 100;
    acf_Y = autocorr(I',nlags);
    gamma = 2*sum((1-(1:nlags)/(length(I)))'.*acf_Y(2:end));
    Var_Q(iter) = var(I)/N*(1+gamma);
    delta_Q(iter)=sqrt(Var_Q(iter))/Q
end

%figure;
%plot(x_fin(:,1))
%title('x1')
%figure;
%plot(x_fin(:,2))
%title('x2');
figure;
scatter(x(:,1),x(:,2))
title('random walk sampler, start: [15,10]')
xlabel('x_1')
xlim([-40,30])
ylim([-50,20])
ylabel('x_2')
% Mean of the estimator
%Q=mean(x_fin,);

% plot autocorrelation of samples
[acf,Lags] = autocorr(x_graph(:,1),50);
figure
plot(Lags,acf);
hold on;
[acf,Lags] = autocorr(x_graph(:,2),50);
plot(Lags,acf);
[acf,Lags] = autocorr(x_graph(:,3),50);
plot(Lags,acf);
xlabel('Lag')
ylabel('Sample Autocorrelation')



