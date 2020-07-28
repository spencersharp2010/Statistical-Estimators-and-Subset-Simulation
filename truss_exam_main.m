clear;
%% Setup
distributions(1) = ERADist('lognormal','MOM',[2e-3,2e-4]);

distributions(2) = ERADist('lognormal','MOM',[1e-3,1e-4]);

distributions(3:4) = ERADist('lognormal','MOM',[2.1e11,2.1e10]);

distributions(5:10) = ERADist('gumbel','MOM',[5e4,7.5e3]);

dim = 10;
corr = eye(dim);
Nataf = ERANataf(distributions,corr);

%% Limit state function
% limit state function
%gfun=@(r,s)s-r;
ulim = 0.12;
gfun = @(input) ulim - truss_exam(input);

%% Subset simulation
% threshold
gamma = 0;
rng(1)
%Q_exact = 0.0331

% Sample size
N_lev=2e3;

% quantile value
p = 0.1;

% correlation parameter for MCMC
rho = 0.8;

% various rho values for parameter study
rho_val = [0.4 0.5 0.6 0.7 0.8 0.9];

% various p values for parameter study
p_val = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
maxT = 100;

num_iter = 5;
for iter = 1:num_iter
    fprintf('iteration %d of %d \n',iter,num_iter);
    %for k = 1:size(rho_val,2)
    %for k = 1:size(p_val,2)
        % Generate samples at first level
        % samples at standard normal space
        u_sam = normrnd(0,1,dim,N_lev);
        
        % samples at original space
        x_sam = Nataf.U2X(u_sam)';
        
        % Evaluate the responses
        %g_sam = gfun(x_sam(:,1),x_sam(:,2));
        g_sam = gfun(x_sam);
        
        for t = 1:maxT
            
            %p = p_val(k);
            p=0.1; 
            % intermediate threshold
            gamma_t(t) = prctile(g_sam,p*100);
            gamma_t(t) = max(gamma_t(t),gamma);
            
            % Check if final level is reached
            if gamma_t(t) == gamma
                break
            end
            
            % Set seeds for MCMC
            I=g_sam<=gamma_t(t);
            
            u_seed = u_sam(:,I);
            g_seed = g_sam(I);
            
            % number of seeds
            N_seed = sum(I);
            
            % lenght of each chain
            N_chain = floor(N_lev/N_seed)*ones(N_seed,1);
            N_chain(1:mod(N_lev,N_seed)) = N_chain(1:mod(N_lev,N_seed))+1;
            
            % MCMC sampling
            count = 0;
            %rho = rho_val(k);
            for i = 1:N_seed
                
                count = count+1;
                u_sam(:,count) = u_seed(:,i);
                g_sam(count) = g_seed(i);
                
                for j = 1:N_chain(i)-1
                    
                    count = count+1;
                    
                    % generate the candidate state
                    z = normrnd(0,sqrt(1-rho^2),dim,1);
                    u_cand = rho*u_sam(:,count-1)+z;
                    
                    % sample at original space
                    x_cand = Nataf.U2X(u_cand)';
                    
                    % Evaluate the responses
                    %g_cand = gfun(x_cand(1),x_cand(2));
                    g_cand = gfun(x_cand);
                    
                    % accept or reject
                    if g_cand <= gamma_t(t)
                        u_sam(:,count) = u_cand;
                        g_sam(count) = g_cand;
                    else
                        u_sam(:,count) = u_sam(:,count-1);
                        g_sam(count) = g_sam(count-1);
                    end            
                                
                end %chains
                
            end %seeds
            
            
            
        end %subset levels
    
        % Final level
        T = t;
        
        % Evaluation of Pf
        I=g_sam<=gamma;
        %Q_SuS = (p^(T-1))*sum(I)/N_lev;
        Q_SuS(iter) = (p^(T-1))*sum(I)/N_lev;
    %end %rho/p values
        
end %iterations for averaging

%% Post processing

m = T-1;

Pf           = zeros(m,1);
Pf(1)        = p;
%Pf_line(1,:) = linspace(p0,1,Nc);
%b_line(1,:)  = prctile(gsort(1,:),Pf_line(1,:)*100);
for i = 2:m
   Pf(i)        = Pf(i-1)*p;
   %Pf_line(i,:) = Pf_line(i-1,:)*p0;
   %b_line(i,:)  = prctile(gsort(i,:),Pf_line(1,:)*100);
end
%Pf_line = sort(Pf_line(:));
%b_line  = sort(b_line(:));

% figure
% %semilogy(rho_val,mean(Q_SuS,1))
% semilogy(p_val,mean(Q_SuS,1))
% title('Q estimate vs p0')
% xlabel('p0')
% ylabel('Q')
% 
% figure
% plot(p_val,std(Q_SuS,1)./mean(Q_SuS,1))
% title('variance vs p0')
% xlabel('p0')
% ylabel('variance')

% figure
% %semilogy(mean(gamma_t,1), Q_SuS, 'r--');           % curve
% semilogy(mean(gamma_t(:,1:2),1), Pf,'ro', 'MarkerSize',6);   % points
% title('PoF vs \gamma')
% xlabel('\gamma')
% ylabel('PoF')
% grid
% hold on;
% semilogy(0, mean(Q_SuS), 'b*', 'MarkerSize',6);
% %semilogy(0, pf_ex, 'ro', 'MarkerSize',8);
% %set(gca,'yscale','log'); axis tight;
% hl = legend('Intermediate levels','PoF', 'Location', 'southeast');
% %set(hl,'Interpreter','latex'); set(gca,'FontSize',18);

%output = ['this is ', num2str(9), ' of ', num2str(10), ' runs'];
%disp(output)
%gamma_avg = mean(gamma_t,1);
%gamma_avg = reshape(gamma_avg,[8,8]);

% conv_steps = zeros(10,1);
% for i=1:size(p_val,2)
%     indices = gamma_t(i,:)==0;
%     
%     conv_steps(i) = sum(indices==0) + 1;
% end
%     
% figure
% plot(p_val,conv_steps,'-o')
% title('steps to convergence, \rho = 0.8')
% xlabel('p0')
% ylabel('steps to convergence')
