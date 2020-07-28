function [Q_SuS,gamma_t,T] = subsetSim(N_lev, p, rho, g_fun, distribution)
%SUBSETSIM computes PoF of given random variables against given limit state
%function

n_dim   = length(distribution.Marginals);

% Generate samples at first level
% samples at standard normal space
u_sam = normrnd(0,1,n_dim,N_lev);

% samples at original space
x_sam = distribution.U2X(u_sam)';

% Evaluate the responses
g_sam = g_fun(x_sam);

gamma = 0;

% maximum number of subsets/levels
maxT = 100;
for t = 1:maxT

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

    % length of each chain
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
            z = normrnd(0,sqrt(1-rho^2),n_dim,1);
            u_cand = rho*u_sam(:,count-1)+z;

            % sample at original space
            x_cand = distribution.U2X(u_cand)';

            % Evaluate the responses
            g_cand = g_fun(x_cand);

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
Q_SuS = (p^(T-1))*sum(I)/N_lev;

end

