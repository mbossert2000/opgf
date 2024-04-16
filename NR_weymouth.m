function [p,q,C] = NR_weymouth(M,Qdot_star,p0,L,D,gas,f_curve,Re_vec,Parameters,idx)
%NR_WEYMOUTH Applies newton raphson method to the gas weymouth load flow equations
%   M: incidence matrix of the gas network
%   Load: load at each node (=-injection)
%   x0:initial vector of pressure
%   step_size:proportion of the standard newton step at each iteration
%   tol:tolerance of convergence (in MW)
%   max_iter:maximum number of iterations
%   L: vector of pipes length
%   D: vector of pipes diam§eter
%   gas: structure containing gas properties
%   f_curve: sample look up table of the AGA friction curve
%   Re_vec: Reynolds number for the sample points in f_curve
%   p: solution pressure at each node
%   q: solution flow in each pipe

[n_pipes,~] = size(M);
%idxs_slack = 1:n_nodes:Parameters.n_timesteps*n_nodes;
%idxs_demand = setdiff(1:Parameters.n_timesteps*n_nodes, idxs_slack);
%MCell = repmat({M}, 1, Parameters.n_timesteps);
%M_ext = blkdiag(MCell{:});
M2 = M;
M2(:,idx.slack_gas) = [];
p = p0;

%Friction factor:
Re = 1e3*ones(n_pipes,1);
fr = AGA(Re_vec,f_curve,Re);

%Weymouth coefficients:
%L_ext = repmat(L,Parameters.n_timesteps,1);
%D_ext = repmat(D,Parameters.n_timesteps,1);
C = gas.m3_to_MW*13.2989.*(gas.T_st/gas.p_st).*(1./(L.*gas.d_rel.*gas.T_avg.*gas.z_avg)).^0.5.*(1./fr).^0.5.*D.^(2.5)*1e5;

%% iteration loop gas
for k=1:Parameters.n_max
    n_iter = k;
    
    % Compute gas pipeline parameters
    sgn = sign(M*p);
    q = sgn.*C.*(sgn.*M*p.^2).^(1/2);
    %Friction factor:
    Re = max(abs(4*gas.rho_st/(pi*gas.mu)*(q./D)*gas.MW_to_m3), 1e3*ones(n_pipes,1));
    fr = AGA(Re_vec,f_curve,Re);
    %Weymouth coefficients:
    C = gas.m3_to_MW*13.2989.*(gas.T_st/gas.p_st).*(1./(L.*gas.d_rel.*gas.T_avg.*gas.z_avg)).^0.5.*(1./fr).^0.5.*D.^(2.5)*1e5;

    % Compute flux at every node
    Qdot = M2'*q;
    %% Mismatch calculation
    dQdot = Qdot - Qdot_star;    
    %% Convergence check
    if(max(abs(dQdot))<Parameters.tol)
        %fprintf('Power loadflow reached convergence at iteration: %i \n', n_iter);
        break;
    elseif(k==Parameters.n_max)
        disp('NR algorithm reached the maximum number of iterations!');
    end
    
    %% gas Jacobian
    dq_dp = M.*(C*p')./sqrt(sgn.*M*p.^2);
    J_gas = -M'*dq_dp;
    J_gas(idx.slack_gas,:) = [];
    J_gas(:,idx.slack_gas) = [];
    J_gas = [J_gas > 0].*min(J_gas,Parameters.max_jac*ones(size(J_gas))) + [J_gas < 0].*max(J_gas,-Parameters.max_jac*ones(size(J_gas)));
    
    %% Solution update

    % Solve
    
    dp = J_gas \ dQdot;
    %dp = dx(2*length(idx.pq)+1:end);
    
    % Update

    p(idx.demand) = p(idx.demand) + Parameters.step_size * dp;
end
