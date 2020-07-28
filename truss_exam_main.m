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
rho_val = linspace(0,1,3);

% various p values for parameter study
p_val = linspace(0.1,0.9,9);

%% IMPORTANT: Select one of the following options by uncommenting it
%run_type = 'primary';   %estimates PoF and CoV of bridge
%run_type = 'p_study';  %studies convergence and variability as p varies
run_type = 'rho_study'; %studies convergence and variability as rho varies

%% Run subset simulation according to run type selected
switch run_type
    case 'primary'
        disp('primary was selected')
        num_iter = 5;
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
        
        figure
        semilogy(mean(gamma_t(:,1:end-1),1), Pf,'ro', 'MarkerSize',6);   % points
        title("PoF vs \gamma, p = " + p + ", \rho = " + rho + ", " + num_iter + " iterations")
        xlabel('\gamma')
        ylabel('PoF')
        grid
        hold on;
        semilogy(0, mean(Q_SuS), 'b*', 'MarkerSize',6);
        %set(gca,'yscale','log'); axis tight;
        hl = legend('Intermediate levels','PoF', 'Location', 'southeast');
        %set(hl,'Interpreter','latex'); set(gca,'FontSize',14);

    case 'p_study'
        disp('p_study was selected')
        for iter = 1:length(p_val)
            fprintf('iteration %d of %d \n',iter,length(p_val));
            [Q_SuS(iter),gamma_t{iter},T] = subsetSim(N_lev, p_val(iter), rho, gfun, Nataf);
        end
        conv_steps = zeros(length(p_val),1);
        for i=1:length(p_val)
            conv_steps(i) = length(gamma_t{i});
        end
        
        figure
        plot(p_val,conv_steps,'-o')
        title("steps to convergence vs p, \rho = " + rho)
        xlabel('p')
        ylabel('steps to convergence')
        
        figure
        plot(p_val,Q_SuS,'-o')
        title("PoF vs p, \rho = " + rho)
        xlabel('p')
        ylabel('PoF')
        
    case 'rho_study'
        disp('rho_study was selected')
        for iter = 1:length(rho_val)
            fprintf('iteration %d of %d \n',iter,length(rho_val));
            [Q_SuS(iter),gamma_t{iter},T] = subsetSim(N_lev, p, rho_val(iter), gfun, Nataf);
        end
        conv_steps = zeros(length(rho_val),1);
        for i=1:length(rho_val)
            conv_steps(i) = length(gamma_t{i});
        end
        
        figure
        plot(rho_val,conv_steps,'-o')
        title("steps to convergence vs \rho, p = " + p)
        xlabel('\rho')
        ylabel('steps to convergence')
        grid
        
        figure
        plot(rho_val,Q_SuS,'-o')
        title("PoF vs \rho, p = " + p)
        xlabel('\rho')
        ylabel('PoF')
        grid
    otherwise
        error('please enter a valid run type: either "primary", "p_study", or "rho_study" ');
end

