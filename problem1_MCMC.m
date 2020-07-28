clear;

g = @(X) -10*X(1).^2 + 2*X(1).*X(2) + 2*X(2).^2;
jointpdf = @(x) (10<=x(1) & x(1)<=20).*(0<=x(2) & x(2)<=20) ...
    .* (x(2)/(10*x(1)^2));

pdf_x1 = @(x1) 20/(x1^2);
pdf_x2 = @(x2) x2/200;

% parameter sigma
%sigma=[0.1 0.5 1 5 10 50];
sigma = 5;

% Length of burn-in period
burn=1000;

% Sample size
N=2e4;
% Q=zeros(size(sigma,2),1);
% Var_Q=zeros(size(sigma,2),1);
% delta_Q=zeros(size(sigma,2),1);
% 
% f1=figure; hold on; xlabel('Lag','FontSize',14); ylabel('Sample autocorrelation','FontSize',14);
% set(f1, 'Units', 'normalized', 'Position', [0.2, 0.3, 0.5, 0.5])

% Target density parameters
% mu_log=0;
% si_log=1;

% Seed value
rng(1,'twister');


%for i=1:length(sigma)
    
x=zeros(burn+N+1,2);
x(1,:)=[5,10];

% Candidate sample
sigma_mat = eye(2)*sigma;
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
    else
        x(j+1,:)=x(j,:);
    end
   
end

% Discard the samples from the burn-in period
x_fin=x(burn+2:end,:);

for k=1:N
    g_eval(k) = g(x_fin(k,:));
end

I = g_eval >= 0;
Q = sum(I)/N
Var_Pr=var(I)/N;


%figure;
%plot(x_fin(:,1))
%title('x1')
%figure;
%plot(x_fin(:,2))
%title('x2');
figure;
scatter(x(:,1),x(:,2))
% Mean of the estimator
%Q=mean(x_fin,);


% Coefficient of variation of the estimator
% estimate the coefficient of variation of Q
nlags = 50;
acf_Y = autocorr(I',nlags);
gamma = 2*sum((1-(1:nlags)/(length(I)))'.*acf_Y(2:end));
Var_Q = var(I)/N*(1+gamma);
delta_Q=sqrt(Var_Q)/Q
% 
% plot(acf_x)
% legendInfo{i} = ['\sigma = ' num2str(sigma(i))];
% 
% %end
% 
% lgd=legend(legendInfo,'Location','eastoutside','Orientation','vertical','FontSize',14);
% xlim([0 nlags])
% 
% delta_Q
