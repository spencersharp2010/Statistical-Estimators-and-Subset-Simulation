clear;
%% Setup
% define distributions as described in exam sheet
distributions(1) = ERADist('lognormal','MOM',[2e-3,2e-4]);
distributions(2) = ERADist('lognormal','MOM',[1e-3,1e-4]);
distributions(3:4) = ERADist('lognormal','MOM',[2.1e11,2.1e10]);
distributions(5:10) = ERADist('gumbel','MOM',[5e4,7.5e3]);

% ten random variables, so we have ten dimensional problem
dim = 10;

% define correlation. Identity matrix since distributions are independent
corr = eye(dim);

% define instance of ERANataf class to be able to transform samples in
% standard normal space back to physical space
Nataf = ERANataf(distributions,corr);

%% Limit state function
% note that with this definition, gfun<=0 means failure
ulim = 0.12;
gfun = @(input) ulim - truss_exam(input);

%% Subset simulation inputs
% threshold
rng(1)

% Sample size
N_lev=2e3;

% quantile value
p = 0.1;

% correlation parameter for MCMC
rho = 0.8;

% various rho values for parameter study
rho_val = linspace(0,1,11);

% various p values for parameter study
p_val = linspace(0.1,0.9,9);

%% IMPORTANT: Select one of the following options by uncommenting it
run_type = 'primary';   %estimates PoF and CoV of bridge
%run_type = 'p_study';  %studies convergence and variability as p varies
%run_type = 'rho_study'; %studies convergence and variability as rho varies

%% Run subset simulation according to run type selected
switch run_type
    case 'primary'
        disp('primary was selected')
        num_iter = 10;
        for iter = 1:num_iter
            fprintf('iteration %d of %d \n',iter,num_iter);
            [Q_SuS(iter),gamma_t(iter,:),T] = subsetSim(N_lev, p, rho, gfun, Nataf);
        end
        
        fprintf('P(u_max(X)>=u_lim) = %4.6f \n',mean(Q_SuS));
        fprintf('approximate coefficient of variation: %.6f \n', std(Q_SuS)./mean(Q_SuS));
        % creates vector of intermediate PoF values which correspond to
        % intermediate gamma values
        m = T-1;
        Pf           = zeros(m,1);
        Pf(1)        = p;
        for i = 2:m
           Pf(i)        = Pf(i-1)*p;
        end
        
        % plot PoF vs intermediate thresholds
        figure
        semilogy(mean(gamma_t(:,1:end-1),1), Pf,'ro', 'MarkerSize',8);   % points
        title("PoF vs \gamma, p = " + p + ", \rho = " + rho + ", " + num_iter + " iterations")
        xlabel('\gamma')
        ylabel('PoF')
        grid
        hold on;
        % add final PoF value at gamma=0
        semilogy(0, mean(Q_SuS), 'b*', 'MarkerSize',8);
        hl = legend('Intermediate levels','PoF', 'Location', 'southeast');
        set(hl,'Interpreter','latex'); 
        set(gca,'FontSize',14);

    case 'p_study'
        disp('p_study was selected')
        
        % calculate PoF for various p values
        for iter = 1:length(p_val)
            fprintf('iteration %d of %d, p = %.1f \n',iter,length(p_val), p_val(iter));
            [Q_SuS(iter),gamma_t{iter},T] = subsetSim(N_lev, p_val(iter), rho, gfun, Nataf);
            fprintf('\tP(u_max(X)>=u_lim) = %4.6f \n',Q_SuS(iter));
        end
        
        % calculate number of steps to convergence vs p value
        conv_steps = zeros(length(p_val),1);
        for i=1:length(p_val)
            conv_steps(i) = length(gamma_t{i});
        end
        
        % plot # of steps to convergence vs p value
        figure
        plot(p_val,conv_steps,'-o')
        title("steps to convergence vs p, \rho = " + rho)
        xlabel('p')
        ylabel('steps to convergence')
        grid
        
        % plot calculated PoF vs p value
        figure
        plot(p_val,Q_SuS,'-o')
        title("PoF vs p, \rho = " + rho)
        xlabel('p')
        ylabel('PoF')
        grid
             
    case 'rho_study'
        disp('rho_study was selected')
        
        % calculate PoF for various rho values
        for iter = 1:length(rho_val)
            fprintf('iteration %d of %d, rho = %.1f \n',iter,length(rho_val), rho_val(iter));
            [Q_SuS(iter),gamma_t{iter},T] = subsetSim(N_lev, p, rho_val(iter), gfun, Nataf);
            fprintf('\tP(u_max(X)>=u_lim) = %4.6f \n',Q_SuS(iter));
        end
        % calculate number of steps to convergence
        conv_steps = zeros(length(rho_val),1);
        for i=1:length(rho_val)
            conv_steps(i) = length(gamma_t{i});
        end
        
        % plot # of steps to convergence vs p value
        figure
        plot(rho_val,conv_steps,'-o')
        title("steps to convergence vs \rho, p = " + p)
        xlabel('\rho')
        ylabel('steps to convergence')
        grid
        
        % plot calculated PoF vs p value
        figure
        plot(rho_val,Q_SuS,'-o')
        title("PoF vs \rho, p = " + p)
        xlabel('\rho')
        ylabel('PoF')
        grid
        
    otherwise
        error('please enter a valid run type: either "primary", "p_study", or "rho_study" ');
end

