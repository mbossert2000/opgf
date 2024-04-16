function [J,E,S,p,q,n_iter] = NR_coupled(S_star,Y,E_0,idx,M,Qdot_star,p0,L,D,gas,f_curve,Re_vec,Parameters)
%
Yabs = abs(Y);
Yarg = angle(Y);
n_buses = length(E_0);

% Initialization
Eabs = abs(E_0);
Earg = angle(E_0);
J_el = [];

%Gas
[n_pipes,n_nodes] = size(M);
p = p0;

%Friction factor:
Re = 1e3*ones(n_pipes,1);
fr = AGA(Re_vec,f_curve,Re);
%Weymouth coefficients:
C = gas.m3_to_MW*13.2989.*(gas.T_st/gas.p_st).*(1./(L.*gas.d_rel.*gas.T_avg.*gas.z_avg)).^0.5.*(1./fr).^0.5.*D.^(2.5)*1e5;


%% iteration loop elec
for k=1:Parameters.n_max
    n_iter = k;
    
    % Compute nodal voltages/currents/power
    E = complex(Eabs.*cos(Earg),Eabs.*sin(Earg));
    I = Y*E;
    S = E.*conj(I);

    %% Mismatch calculation
    
    % Compute the mismatches for the entire network.
    dS = S_star-S;
    dP = real(dS);
    dQ = imag(dS);

    % Keep only the relevant mismatches.
    dP(idx.slack) = [];
    dQ(idx.slack) = [];
    
    % Compute gas pipeline parameters
    sgn = sign(M*p);
    q = sgn.*C.*(sgn.*M*p.^2).^(1/2);
    %Friction factor:
    Re = max(abs(4*gas.rho_st/(pi*gas.mu)*(q./D)*gas.MW_to_m3), 1e3*ones(n_pipes,1));
    fr = AGA(Re_vec,f_curve,Re);
    %Weymouth coefficients:
    C = gas.m3_to_MW*13.2989.*(gas.T_st/gas.p_st).*(1./(L.*gas.d_rel.*gas.T_avg.*gas.z_avg)).^0.5.*(1./fr).^0.5.*D.^(2.5)*1e5;

    % Compute flux at every node
    Qdot = M(:,2:end)'*q;
    dQdot = Qdot - Qdot_star;    
    %%
    dF = [dP;dQ;dQdot]; % mismatch of the power flow equations
    %% Convergence check
    if(max(abs(dF))<Parameters.tol)
        %fprintf('Power loadflow reached convergence at iteration: %i \n', n_iter);
        break;
    elseif(k==Parameters.n_max)
        disp('NR algorithm reached the maximum number of iterations!');
    end
    
    %% Jacobian construction
    
    % For the sake of simplicity, the blocks of J are constructed
    % for the whole network (i.e., with size n_buses x n_buses).
    % The unnecessary rows/columns are removed subsequently
    
    % Extract magnitude/angle
    Eabs = abs(E);
    Earg = angle(E);
    
    % Initialization
    J_PE = zeros(n_buses,n_buses); % derivative: P versus E_abs
    J_PT = zeros(n_buses,n_buses); % derivative: P versus E_arg (theta)
    J_QE = zeros(n_buses,n_buses); % derivative: Q versus E_abs
    J_QT = zeros(n_buses,n_buses); % derivative: Q versus E_arg (theta)
    
    % Construction
    for i=1:n_buses
        
        % Diagonal elements (terms outside the sum)
        J_PE(i,i) =  2*Yabs(i,i)*Eabs(i)*cos(Yarg(i,i));
        J_QE(i,i) = -2*Yabs(i,i)*Eabs(i)*sin(Yarg(i,i));
        
        for j=1:n_buses
            if i ~= j
                % Diagonal elements (terms inside the sum)
                J_PE(i,i) = J_PE(i,i) + Yabs(i,j)*Eabs(j)*cos(Earg(i)-Earg(j)-Yarg(i,j));
                J_QE(i,i) = J_QE(i,i) + Yabs(i,j)*Eabs(j)*sin(Earg(i)-Earg(j)-Yarg(i,j));
                J_PT(i,i) = J_PT(i,i) - Eabs(i)*Yabs(i,j)*Eabs(j)*sin(Earg(i)-Earg(j)-Yarg(i,j));
                J_QT(i,i) = J_QT(i,i) + Eabs(i)*Yabs(i,j)*Eabs(j)*cos(Earg(i)-Earg(j)-Yarg(i,j));

                % Offdiagonal elements
%                 J_PE(i,j) =  Y_abs(i,j)*E_abs(i)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j));
%                 J_QE(i,j) =  Y_abs(i,j)*E_abs(i)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j));
%                 J_PT(i,j) =  Y_abs(i,j)*E_abs(i)*E_abs(j)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j));
%                 J_QT(i,j) = -Y_abs(i,j)*E_abs(i)*E_abs(j)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j));
                
                J_PE(i,j) = J_PE(i,j)+Yabs(i,j)*Eabs(i)*cos(Earg(i)-Earg(j)-Yarg(i,j)); 
                J_QE(i,j) = J_QE(i,j)+Yabs(i,j)*Eabs(i)*sin(Earg(i)-Earg(j)-Yarg(i,j)); 
                J_PT(i,j) = J_PT(i,j)+ Yabs(i,j)*Eabs(i)*Eabs(j)*sin(Earg(i)-Earg(j)-Yarg(i,j));
                J_QT(i,j) = J_QT(i,j)-Yabs(i,j)*Eabs(i)*Eabs(j)*cos(Earg(i)-Earg(j)-Yarg(i,j)); 

            end
        end
    end
    
    % Remove extra rows (i.e., unnecessary equations)
    % slack bus: P & Q, PV buses: Q
    
    J_PE(idx.slack,:) = [];
    J_PT(idx.slack,:) = [];

    J_QE(idx.slack,:) = [];
    J_QT(idx.slack,:) = [];
    
    % Remove extra columns (i.e., variables)
    % slack bus: E_abs & E_arg, PV nodes: E_abs
    
    J_PE(:,idx.slack) = [];
    J_QE(:,idx.slack) = [];
    
    J_PT(:,idx.slack) = [];
    J_QT(:,idx.slack) = [];
    
    % Combination
    J_el = [J_PE,J_PT;J_QE,J_QT];
    

    %% gas Jacobian
    dq_dp = M.*(C*p')./sqrt(sgn.*M*p.^2);
    J_gas = -M'*dq_dp;
    J_gas = J_gas(2:end,2:end);
    J_gas = [J_gas > 0].*min(J_gas,Parameters.max_jac*ones(size(J_gas))) + [J_gas < 0].*max(J_gas,-Parameters.max_jac*ones(size(J_gas)));
    %% Combination
    J = [J_el zeros(2*(n_buses-1), n_nodes-1)
         zeros(n_nodes-1, 2*(n_buses-1)) J_gas];
    %% Solution update

    % Solve
    dx = J \ dF;

    % Reconstruct the solution
    
    dEabs = zeros(length(Eabs),1);
    dEabs(idx.pq,1) = dx(1:length(idx.pq));
    
    dEarg = zeros(length(Earg),1);
    dEarg(idx.pq,1) = dx((length(idx.pq)+1):2*length(idx.pq));
    
    % Update
    Eabs = Eabs + dEabs;
    Earg = Earg + dEarg;

    p(2:end) = p(2:end) + Parameters.step_size * dx(2*length(idx.pq)+1:end);
end

E = Eabs .* exp(1i*Earg);

end