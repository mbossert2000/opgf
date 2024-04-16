clear
clc
close all
%% run name:
run_name = "scheduling_64scenarios_decoupled";
is_scheduling = true; %activate the scheduling penalty
is_allscenarios = true; %choose betwne 1 or 64 scenarios for testing
rng(42) %for reproducibility
%% base parameters
%grid base parameter
Ab = 100e6; %100 MW
Vb = 20e3;
Zb = Vb^2/Ab;

%perfect gas properties
gas = load("gas_properties_struct.mat").gas_properties_struct;

%% Power grid data
%Lines and bus tables
lines = table2array(readtable("./data/lines_CIGRE_pandapower.csv"));
n_lines = height(lines);
table_buses = readtable("./data/nodes_CIGRE.csv");
n_buses = height(table_buses);
idx.slack = 1;
idx.pq = 2:n_buses;
idx.pv = [];
idx.src = [1];

%admittance matrix
is_pu = true; %are lines parameters in ohms/farad or pu.
[Y, YYL, YYT, Ib] = Ymatrix(lines, Ab, Vb, is_pu);

%% Gas network data
% Gas network tables
table_nodes = readtable("./data/node_list.csv","ReadRowNames",true);
table_pipes = readtable("./data/pipe_list.csv","ReadRowNames",true);

%Gas incidence matrix
[n_pipes,~] = size(table_pipes);
[n_nodes,~] = size(table_nodes);
idx.slack_gas = table_nodes{table_nodes.type == "slack","idx"};
idx.demand = table_nodes{table_nodes.type == "demand","idx"};

M = zeros(n_pipes,n_nodes);
for p=1:n_pipes
    from_idx = table_nodes{table_pipes{p,"From"}{1,1},"idx"};
    to_idx = table_nodes{table_pipes{p,"To"}{1,1},"idx"};
    M(p,from_idx) = 1;
    M(p,to_idx) = -1;
end

%pipe geometry
L = table_pipes{:,"Leq"};
D = table_pipes{:,"Diameter"}*1e-3;
roughness = table_pipes{:,"Roughness"}*1e-3;
A = pi.*D.^2./4;
%friction factor curve
[f_curve,Re_vec] = AGA_vec(D,roughness); %Friction factor
%semilogx(Re_vec,f_curve(:,:)') %to plot the curve

%pressure limits
pmin = table_nodes{:,"pmin"};
pmax = table_nodes{:,"pmax"};


%% daily profiles
%profiles = readtable("./data/load_profile.csv");
n_timesteps = 96;
resolution = 0.25; %resolution of timestep
profile.hours = [0:0.25:23.75];

%electrical load
profile.load = table2array(readtable("./data/scenarios/load_scenarios.csv")); %percentage of max value
profile.loadgas = table2array(readtable("./data/scenarios/gasload_scenarios.csv"));
profile.pv = table2array(readtable("./data/scenarios/pv_scenarios.csv")); %percentage of max value
profile.wind = table2array(readtable("./data/scenarios/wind_scenarios.csv")); %percentage of max value

if is_allscenarios
    profile.pv = profile.pv(:,[2,4,11,16]);
    profile.load = profile.load(:,[2,4,5,9]);
    profile.wind = profile.wind(:,[1,4]);
    profile.loadgas = profile.loadgas(:,[2,3]);
else
    profile.load = profile.load(:,1);
    profile.loadgas = profile.loadgas(:,1);
    profile.pv = profile.pv(:,1);
    profile.wind = profile.wind(:,1);
end

n_scenarios_load = width(profile.load);
n_scenarios_loadgas = width(profile.loadgas);
n_scenarios_pv = width(profile.pv);
n_scenarios_wind = width(profile.wind);

price_el_hour = table2array(readtable("./data/scenarios/elpricenextday.csv"));
%upscale from hourly prices to 1/4 hour prices
profile.price_el = [];
for i=1:24
    profile.price_el = [profile.price_el; price_el_hour(i); price_el_hour(i); price_el_hour(i); price_el_hour(i)];
end
profile.price_gas = table2array(readtable("./data/scenarios/gaspricenextday.csv"));
%build combination of scenarios
n_scenarios = n_scenarios_pv * n_scenarios_wind * n_scenarios_load * n_scenarios_loadgas;
s = 0;
scenarios.load = [];
scenarios.pv = [];
scenarios.wind = [];
scenarios.loadgas = [];
scenarios.price_el = [];
scenarios.price_gas = [];
for s_pv = 1:n_scenarios_pv
    for s_wind = 1:n_scenarios_wind
        for s_load = 1:n_scenarios_load
            for s_gas = 1:n_scenarios_loadgas
                s = s + 1;
                scenarios.load = [scenarios.load, profile.load(:,s_load)];
                scenarios.pv = [scenarios.pv profile.pv(:,s_pv)];
                scenarios.wind = [scenarios.wind profile.wind(:,s_wind)];
                scenarios.loadgas = [scenarios.loadgas profile.loadgas(:,s_gas)];
                scenarios.price_el = [scenarios.price_el profile.price_el(:,1)];
                scenarios.price_gas = [scenarios.price_gas profile.price_gas(:,1)];
            end
        end
    end
end
        

%% Distributed ressources
%Electrical distributed generation:
table_DG = readtable("./data/DG.csv");
%photovoltaics
idx.PV = table_DG{table_DG.type == "photovoltaic","Node_grid"}; %node with solar PV installed
Pmax.PV = table_DG{table_DG.type == "photovoltaic","P_max"}*1e3/Ab; %installed power of PV converted in pu.
n_PV = length(idx.PV);
%wind turbine
idx.WT = table_DG{table_DG.type == "wind_turbine","Node_grid"};
Pmax.WT = table_DG{table_DG.type == "wind_turbine","P_max"}*1e3/Ab;
n_WT = length(idx.WT);
%Fuel cell
idx.FC_elec = table_DG{table_DG.type == "fuel_cell","Node_grid"};
Pmax.FC = table_DG{table_DG.type == "fuel_cell","P_max"}*1e3/Ab*0;
%feedin
idx.feed_in_elec = table_DG{table_DG.type == "feed_in","Node_grid"};
Pmax.feedin = table_DG{table_DG.type == "feed_in","P_max"}*1e3/Ab*0;
%Battery
idx.BAT = table_DG{table_DG.type == "battery","Node_grid"};
Pmax.BAT = table_DG{table_DG.type == "battery","P_max"}*1e3/Ab;
BAT_max = table_DG{table_DG.type == "battery","capacity"}*1e3/Ab;
efficiency.BAT = table_DG{table_DG.type == "battery","efficiency"};%assuming same chargin and dischargin efficiency
n_BAT = length(idx.BAT);


%gas distributed ressources
idx.feed_in = table_DG{table_DG.type == "feed_in","Node_gas"};
n_feedin = length(idx.feed_in);
feed_in_max = table_DG{table_DG.type == "feed_in","P_max"}/1e3; %maximum feed-in power [MWth]
efficiency_feedin = table_DG{table_DG.type == "feed_in","efficiency"};
idx.FC_gas = table_DG{table_DG.type == "fuel_cell","Node_gas"}; %idx from gas network (different from power network)
n_FC = length(idx.FC_gas);
efficiency_FC = table_DG{table_DG.type == "fuel_cell","efficiency"};
idx.storage = table_DG{table_DG.type == "gas_storage","Node_gas"};
n_storage = length(idx.storage);
storage_max = table_DG{table_DG.type == "gas_storage","capacity"}/1e3;%in MWh_th
inj_storage_max = table_DG{table_DG.type == "gas_storage","P_max"}/1e3;%in MWh_th
efficiency.storage = table_DG{table_DG.type == "gas_storage","efficiency"};%assuming same chargin and dischargin efficiency


%% Power grid initial loadflow
P_load = zeros(n_buses,n_timesteps,n_scenarios);
Q_load = zeros(n_buses,n_timesteps,n_scenarios);

for s=1:n_scenarios
    P_load(:,:,s) = -table_buses{:,"P_max_house"}*scenarios.load(:,s)' - table_buses{:,"P_max_ind"}*scenarios.load(:,s)';
    Q_load(:,:,s) = -table_buses{:,"Q_max_house"}*scenarios.load(:,s)' - table_buses{:,"Q_max_ind"}*scenarios.load(:,s)';
end

S = P_load + j.*Q_load;
%adding PV and WT: (optionnal, might be unfeasible if too much PV)
for s=1:n_scenarios
    %S(idx.PV,:,s) = S(idx.PV,:,s) + Pmax.PV*scenarios.pv(:,s)';
    %S(idx.WT,:,s) = S(idx.WT,:,s) + Pmax.WT*scenarios.wind(:,s)';
end

E_star = []; %only used for pv nodes
E_start = ones(n_buses,1); %flat start
n_sources = length(idx.src);
Parameters.tol = 1e-6; %power mismatch tolerance, in pu.
Parameters.n_max = 1000; %maximum number of iteration
V0 = ones(n_buses,n_timesteps,n_scenarios);
S0 = ones(n_buses,n_timesteps,n_scenarios);
S0_losses = ones(1,n_timesteps,n_scenarios);
for s=1:n_scenarios
    for t=1:n_timesteps
        [J_t,V0_t,S0_t,n_iter_t] = NR_polar_sol(S(:,t,s),E_star,Y,E_start,idx,Parameters);
        V0(:,t,s) = V0_t;
        S0(:,t,s) = S0_t;
    
        %fprintf('Initial Power loadflow reached convergence at iteration: %i \n', n_iter);
        S0_losses(:,t,s) = sum(S0_t);
    end
end
Vmag_nosource = abs(V0); %for plotting latter
%% Gas network initial loadflow
%Gas loads:
load_gas = zeros(n_nodes,n_timesteps,n_scenarios);

for s=1:n_scenarios
    load_gas(:,:,s) = table_nodes{:,"demand"}*scenarios.loadgas(:,s)';
end

p_slack = 0.7; %pressure of the slack in bars
x0 = 0.6*ones(n_nodes-1,n_timesteps,n_scenarios)+rand(n_nodes-1,n_timesteps,n_scenarios)*0.1; %random initialisation necessary for good initial gradient
Parameters.tol_MW = 1e-5; %power mismatch tolerance, in MW
step_size = 0.1;
p_lf_init = ones(n_nodes-1,n_timesteps,n_scenarios);
q_lf_init = ones(n_pipes,n_timesteps,n_scenarios);
p0 = ones(n_nodes,n_timesteps,n_scenarios);
for s=1:n_scenarios
    for t=1:n_timesteps
        load_t = load_gas(2:end,t,s);
        [p_lf, q_lf] = NR_weymouth(M,load_t,p_slack,x0(:,t,s),step_size,Parameters.tol_MW,Parameters.n_max,L,D,gas,f_curve,Re_vec);
        p_lf_init(:,t,s) = p_lf;
        q_lf_init(:,t,s) = q_lf;
        p0(:,t,s) = [p_slack;p_lf];
    end
end
q0 = q_lf_init;
%% Linearization
Res_nodes = [2:n_buses];
nph = 1;
alphap = ones(n_buses,n_buses);
alphaq = ones(n_buses,n_buses);
Elin_0 = ones(n_buses,1);
K_p = zeros(n_buses,n_buses-1,n_timesteps,n_scenarios);
K_q = zeros(n_buses,n_buses-1,n_timesteps,n_scenarios);
C_rp = zeros(1,n_buses-1,n_timesteps,n_scenarios);
C_rq = zeros(1,n_buses-1,n_timesteps,n_scenarios);
C_xp = zeros(1,n_buses-1,n_timesteps,n_scenarios);
C_xq = zeros(1,n_buses-1,n_timesteps,n_scenarios);
KI_p = zeros(2*n_lines,n_buses-1,n_timesteps,n_scenarios);
KI_q = zeros(2*n_lines,n_buses-1,n_timesteps,n_scenarios);
I0 = zeros(2*n_lines,n_timesteps,n_scenarios);
%gas sensitivity coefficients
K_inj = zeros(n_nodes,n_nodes-1,n_timesteps,n_scenarios);

Re = max(abs(4*gas.rho_st/(pi*gas.mu)*(abs(q0)./D)*gas.MW_to_m3), 1e3*ones(n_pipes,1));
f = zeros(n_pipes,n_timesteps,n_scenarios);
for s=1:n_scenarios
    f(:,:,s) = AGA(Re_vec,f_curve,Re(:,:,s));
end

C = gas.m3_to_MW*13.2989.*(gas.T_st/gas.p_st).*(1./(L.*gas.d_rel.*gas.T_avg.*gas.z_avg)).^0.5.*(1./f).^0.5.*D.^(2.5)*1e5;
for s=1:n_scenarios
    for t=1:n_timesteps
        [K_p(:,:,t,s),K_com_p,K_q(:,:,t,s),K_comp_q] = Coeffs_Voltage_Alpha(Y,S0(:,t,s),V0(:,t,s),Res_nodes,idx.slack,nph,alphap,alphaq,Elin_0);
        [C_rp(:,:,t,s) , C_rq(:,:,t,s), C_xp(:,:,t,s), C_xq(:,:,t,s)] = Coeffs_Losses(Y, V0(:,t,s), K_com_p, K_comp_q, idx.pq);
        [KI_p(:,:,t,s), KI_q(:,:,t,s), I0(:,t,s)] = Coeffs_Currents(YYL,YYT,V0(:,t,s),K_com_p,K_comp_q,nph,lines);
        K_inj(:,:,t,s) = pressure_coeff(M,C(:,t,s),p0(:,t,s));
    end
end


%%
% costs grid:
cost_power_slack = reshape(scenarios.price_el*resolution,1,n_timesteps,n_scenarios);
% costs gas
cost_gas_slack = reshape(scenarios.price_gas*resolution,1,n_timesteps,n_scenarios);
%costs device utilisation
cost_FC = 1*resolution;
cost_feedin = 1*resolution;
cost_BAT = 100*resolution;
cost_STO = 1*resolution;

%initial point power grid
P0 = P_load;
Q0 = Q_load;
P0_losses = real(S0_losses);
Q0_losses = imag(S0_losses);
P0_FC = zeros(n_FC,n_timesteps,n_scenarios);
P0_feedin = zeros(n_feedin,n_timesteps,n_scenarios);
P0_PV = zeros(n_PV,n_timesteps,n_scenarios);
P0_WT = zeros(n_WT,n_timesteps,n_scenarios);
Parameters.max_iter_opti = 50;
Parameters.tol_v = 1e-6;%voltage tolerance in pu.
%BAT0 = BAT_max./2.*ones(n_BAT,1,n_scenarios); %initial state of the battery
BAT0 = zeros(n_BAT,1,n_scenarios); %initial state of the battery

%initial point gas network
feed_in0 = zeros(n_feedin,n_timesteps,n_scenarios);
load_FC0 = zeros(n_FC,n_timesteps,n_scenarios);
inj0 = -load_gas;
storage0 = zeros(n_storage,1,n_scenarios);%initiale state of the gas storage

%parameters for convergence
Parameters.tol_p = 1e-4; %pressure tolerance in bar
Parameters.tol_rel = 1e-3; % 0.1% relative eror
step_p = 0.05; %in bar
step_v = 0.05; %in pu


% iteration historic
cost_hist = [];
solution_hist = [];
lin_cost_hist =[];
%% optimisation loop
for iter=1:Parameters.max_iter_opti

    %Create optimization variables
    %power grid
    P = optimvar("P",[n_buses,n_timesteps,n_scenarios]); %Bus active power injection [pu.]
    Q = optimvar("Q",[n_buses,n_timesteps,n_scenarios]); %Bus reactive power injection [pu.]
    P_slack = optimvar("P_slack",[1,n_timesteps,n_scenarios]); %Active power drawn from the slack [pu.]
    Q_slack = optimvar("Q_slack",[1,n_timesteps,n_scenarios]); %Reactive power drawn from the slack [pu.]
    P_losses = optimvar("P_losses",[1,n_timesteps,n_scenarios]); %Active power losses in the grid [pu.]
    Q_losses = optimvar("Q_losses",[1,n_timesteps,n_scenarios]); %Reactive power losses in the grid [pu.]
    P_FC = optimvar("P_FC",[n_FC,n_timesteps,n_scenarios],"LowerBound",0,"UpperBound",Pmax.FC.*ones(1,n_timesteps,n_scenarios)); %Active power produced from Fuel Cells [pu.]
    P_PV = optimvar("P_PV",[n_PV,n_timesteps,n_scenarios],"LowerBound",0,"UpperBound",Pmax.PV.*ones(1,n_timesteps,n_scenarios)); %Active power produced from Photovoltaic pannels [pu.]
    P_WT = optimvar("P_WT",[n_WT,n_timesteps,n_scenarios],"LowerBound",0,"UpperBound",Pmax.WT.*ones(1,n_timesteps,n_scenarios)); %Active power produced from wind turbines pannels [pu.]
    P_feedin = optimvar("P_feedin",[n_feedin,n_timesteps,n_scenarios],"LowerBound",0,"UpperBound",Pmax.feedin.*ones(1,n_timesteps,n_scenarios)); %Active power drawned by the electrolyser [pu.]
    P_BAT = optimvar("P_BAT",[n_BAT,n_timesteps,n_scenarios],"LowerBound",-Pmax.BAT.*ones(1,n_timesteps,n_scenarios),"UpperBound",Pmax.BAT.*ones(1,n_timesteps,n_scenarios)); %Power injected/drawn by the battery [pu.]
    BAT = optimvar("BAT",[n_BAT,n_timesteps,n_scenarios],"LowerBound",0,"UpperBound",BAT_max.*ones(1,n_timesteps,n_scenarios)); %energy in the battery in pu.h
    I = optimvar("I",[2.*n_lines,n_timesteps,n_scenarios]); %Line currents [pu.]
    Vmag = optimvar("Vmag",[n_buses,n_timesteps,n_scenarios]); %Bus voltages magnitude [pu.]
    
    %gas network
    p = optimvar("p",[n_nodes,n_timesteps,n_scenarios],"LowerBound",pmin.*ones(1,n_timesteps,n_scenarios),"UpperBound",pmax.*ones(1,n_timesteps,n_scenarios)); %Node pressures [bar]
    %q = optimvar("q",[n_pipes,n_timesteps]); %Pipe flows [MWth]
    inj_slack = optimvar("inj_slack",[1,n_timesteps,n_scenarios],"LowerBound",0); %Gas injection from the external network [MWth]
    feed_in = optimvar("feed_in",[n_feedin,n_timesteps,n_scenarios],"LowerBound",0);%Local feed in of gas network [MWth]
    load_FC = optimvar("load_FC",[n_FC,n_timesteps,n_scenarios],"LowerBound",0);%Gas loads from Fuel Cells [MWth]
    inj = optimvar("inj",[n_nodes,n_timesteps,n_scenarios]); %Net nodal gas flow injection [MWth]
    storage = optimvar("storage",[n_storage,n_timesteps,n_scenarios],"LowerBound",0,"UpperBound",storage_max.*ones(1,n_timesteps,n_scenarios));
    Qcharge = optimvar("Qcharge",[n_storage,n_timesteps,n_scenarios],"LowerBound",0,"UpperBound",inj_storage_max*ones(1,n_timesteps,n_scenarios));
    Qdischarge = optimvar("Qdischarge",[n_storage,n_timesteps,n_scenarios],"LowerBound",0);

    %costs (objective function)
    c_elec = optimvar("c_elec",[1,n_timesteps,n_scenarios]);
    c_gas = optimvar("c_gas",[1,n_timesteps,n_scenarios]);
    c_BAT = optimvar("c_BAT",[n_BAT,n_timesteps,n_scenarios]);

    %Scheduling:
    if is_scheduling
        P_slack_scheduling = optimvar("P_slack_scheduling",[1,n_timesteps]); %Active power drawn from the slack [pu.]
        %penalty on daily total
        %daily_inj_scheduling = optimvar("daily_inj_scheduling",1,"LowerBound",0); %Gas injection from the external network [MWth]
        %penalty on each time step
        inj_scheduling = optimvar("inj_scheduling",[1,n_timesteps],"LowerBound",0);
        initialPoint.P_slack_scheduling = zeros(1,n_timesteps);
        %initialPoint.daily_inj_scheduling = 0;
        initialPoint.inj_scheduling = zeros(1,n_timesteps);
    end
    
    % Set initial starting point for the solver
    initialPoint.P = P0;
    initialPoint.Q = Q0;
    initialPoint.P_losses = P0_losses;
    initialPoint.Q_losses = Q0_losses;
    initialPoint.P_FC = P0_FC;
    initialPoint.P_feedin = P0_feedin;
    initialPoint.P_PV = P0_PV;
    initialPoint.P_WT = P0_WT;
    initialPoint.BAT = BAT0.*ones(n_BAT,n_timesteps,n_scenarios);
    initialPoint.P_BAT = zeros(n_BAT,n_timesteps,n_scenarios);
    initialPoint.P_slack = zeros(1,n_timesteps,n_scenarios);
    initialPoint.Q_slack = zeros(1,n_timesteps,n_scenarios);
    initialPoint.I = I0;
    initialPoint.Vmag = abs(V0);
    %gas network
    initialPoint.p = p0;
    %initialPoint.q = q0;
    initialPoint.feed_in = feed_in0;
    initialPoint.load_FC = load_FC0;
    initialPoint.inj_slack = q0(1,:,:);
    initialPoint.inj = inj0;
    initialPoint.Qcharge = zeros(n_storage,n_timesteps,n_scenarios);
    initialPoint.Qdischarge = zeros(n_storage,n_timesteps,n_scenarios);
    initialPoint.storage = zeros(n_storage,n_timesteps,n_scenarios);

    initialPoint.c_elec = zeros(1,n_timesteps,n_scenarios);
    initialPoint.c_gas = zeros(1,n_timesteps,n_scenarios);
    initialPoint.c_BAT = zeros(n_BAT,n_timesteps,n_scenarios);
    
    % Create problem
    problem = optimproblem('ObjectiveSense','min');
    
    % Define problem objective
    problem.Objective = sum(c_gas,'all') + sum(c_elec,'all') + cost_FC.*sum(P_FC,"all") + cost_feedin.*sum(P_feedin,"all") + cost_STO.*sum(Qcharge,"all")+ sum(c_BAT,"all");
    if is_scheduling
        for s=1:n_scenarios
            problem.Objective = problem.Objective + 1000*(Ab/1e6).^2*sum((P_slack_scheduling-P_slack(:,:,s)).^2,'all');
                %scheduling on daily total
            %problem.Objective = problem.Objective + 1000*(daily_inj_scheduling-sum(inj_slack(:,:,s),'all')).^2;
            %scheduling on timestep deviation
            problem.Objective = problem.Objective + 1000*sum((inj_scheduling-inj_slack(:,:,s)).^2,'all');
        end
    end
    % Define problem constraints
    %costs
    problem.Constraints.c_elec_constr = c_elec == P_slack.*cost_power_slack.*Ab/1e6;
    problem.Constraints.c_gas_constr = c_gas == inj_slack.*cost_gas_slack;
    problem.Constraints.c_BAT_constr_1 = cost_BAT.*P_BAT <= c_BAT;
    problem.Constraints.c_BAT_constr_2 = -cost_BAT.*P_BAT <= c_BAT;

    %linear approximation
    constr_Vmag = optimconstr(n_buses,n_timesteps,n_scenarios);
    constr_Ploss = optimconstr(1,n_timesteps,n_scenarios);
    constr_Qloss = optimconstr(1,n_timesteps,n_scenarios);
    constr_I = optimconstr(2*n_lines,n_timesteps,n_scenarios);
    constr_psens = optimconstr(n_nodes,n_timesteps,n_scenarios);
    for s=1:n_scenarios
        for t=1:n_timesteps
            constr_Vmag(:,t,s) = Vmag(:,t,s) - abs(V0(:,t,s)) == K_p(:,:,t,s) * (P(idx.pq,t,s) - P0(idx.pq,t,s)) + K_q(:,:,t,s) * (Q(idx.pq,t,s) - Q0(idx.pq,t,s));
            constr_Ploss(:,t,s) = P_losses(:,t,s) - P0_losses(:,t,s) == C_rp(:,:,t,s)*((P(idx.pq,t,s) - P0(idx.pq,t,s))) + C_rq(:,:,t,s) * (Q(idx.pq,t,s) - Q0(idx.pq,t,s));
            constr_Qloss(:,t,s) = Q_losses(:,t,s) - Q0_losses(:,t,s) == C_xp(:,:,t,s)*((P(idx.pq,t,s) - P0(idx.pq,t,s))) + C_xq(:,:,t,s) * (Q(idx.pq,t,s) - Q0(idx.pq,t,s));
            constr_I(:,t,s) = I(:,t,s) - I0(:,t,s) == KI_p(:,:,t,s)*(P(idx.pq,t,s) - P0(idx.pq,t,s)) + KI_q(:,:,t,s)*(Q(idx.pq,t,s) - Q0(idx.pq,t,s));
            constr_psens(:,t,s) = p(:,t,s)-p0(:,t,s) == K_inj(:,:,t,s)*(inj(2:end,t,s) - inj0(2:end,t,s));
        end
    end
    
    problem.Constraints.Vmagconstr = constr_Vmag;
    problem.Constraints.Ploss = constr_Ploss;
    problem.Constraints.Qloss = constr_Qloss;
    problem.Constraints.Iconstr = constr_I;
    %problem.Constraints.weymouth = C.*T.*(M*p) == q;
    problem.Constraints.psens = constr_psens;

    %local neighbourhood constraints:
    problem.Constraints.small_dev_leq = p - p0 <= step_p;
    problem.Constraints.small_dev_geq = p - p0 >= -step_p;
    problem.Constraints.small_dev_vmag_leq = Vmag - abs(V0) <= step_v;
    problem.Constraints.small_dev_vmag_geq = Vmag - abs(V0) >= -step_v;


    %nodal and branch limits
    problem.Constraints.Imax = I <= 0.05;
    problem.Constraints.Imin = I >= -0.05;
    problem.Constraints.Vmax = Vmag >= 0.95;
    problem.Constraints.Vmin = Vmag <= 1.05;
    
    %nodal power balance
    %Grid power
    P_bal = optimconstr(length(idx.pq),n_timesteps,n_scenarios);
    Q_bal = optimconstr(length(idx.pq),n_timesteps,n_scenarios);
    for n=1:length(idx.pq)
        bus_idx = idx.pq(n);
        [is_FC, loc_FC] = ismember(bus_idx,idx.FC_elec);
        [is_feedin, loc_feedin] = ismember(bus_idx,idx.feed_in_elec);
        [is_PV, loc_PV] = ismember(bus_idx,idx.PV);
        [is_WT, loc_WT] = ismember(bus_idx,idx.WT);
        [is_BAT, loc_BAT] = ismember(bus_idx,idx.BAT);
        P_bus = zeros(1,n_timesteps,n_scenarios);
        Q_bus = zeros(1,n_timesteps,n_scenarios);
        if is_FC
            P_bus = P_bus + P_FC(loc_FC,:,:);
        end
        if is_feedin
            P_bus = P_bus - P_feedin(loc_feedin,:,:);
        end
        if is_PV
            P_bus = P_bus + P_PV(loc_PV,:,:);
        end
        if is_WT
            P_bus = P_bus + P_WT(loc_WT,:,:);
        end
        if is_BAT
            P_bus = P_bus + P_BAT(loc_BAT,:,:);
        end
        P_bal(n,:,:) = P(bus_idx,:,:) == P_load(bus_idx,:,:) + P_bus;
        Q_bal(n,:,:) = Q(bus_idx,:,:) == Q_load(bus_idx,:,:) + Q_bus;
    end
    problem.Constraints.Pbal = P_bal;
    problem.Constraints.Qbal = Q_bal;

    flow_bal = optimconstr(n_nodes,n_timesteps,n_scenarios);
    for k=1:n_nodes
        [is_FC, loc_FC] = ismember(k,idx.FC_gas);
        [is_feedin, loc_feedin] = ismember(k,idx.feed_in);
        [is_storage, loc_storage] = ismember(k,idx.storage);
        [is_slack, loc_slack] = ismember(k,idx.slack_gas);
        balance_node = zeros(1,n_timesteps,n_scenarios);
        if is_FC
            balance_node = balance_node - load_FC(loc_FC,:,:);
        end
        if is_feedin
            balance_node = balance_node + feed_in(loc_feedin,:,:);
        end
        if is_storage
            balance_node = balance_node - Qcharge(loc_storage,:,:) + Qdischarge(loc_storage,:,:);
        end
        if is_slack
            balance_node = balance_node + inj_slack(loc_slack,:,:);
        end
        flow_bal(k,:,:) = inj(k,:,:) == balance_node - load_gas(k,:,:);
    end
    problem.Constraints.flow_bal = flow_bal;
    %problem.Constraints.kirchhoff_gas = -M'*q + inj == 0; %kirchhoff law for gas injections, positive flow leave the node hence minus sign

    %sum DG on slack
    P_DG_slack = zeros(1,n_timesteps,n_scenarios);
    Q_DG_slack = zeros(1,n_timesteps,n_scenarios);
    [is_FC, loc_FC] = ismember(idx.slack,idx.FC_elec);
    [is_feedin, loc_feedin] = ismember(idx.slack,idx.feed_in_elec);
    [is_PV, loc_PV] = ismember(idx.slack,idx.PV);
    [is_WT, loc_WT] = ismember(idx.slack,idx.WT);
    [is_BAT, loc_BAT] = ismember(idx.slack,idx.BAT);
    if is_FC
        P_DG_slack = P_DG_slack + P_FC(loc_FC,:,:);
    end
    if is_feedin
        P_DG_slack = P_DG_slack - P_feedin(loc_feedin,:,:);
    end
    if is_PV
        P_DG_slack = P_DG_slack + P_PV(loc_PV,:,:);
    end
    if is_WT
        P_DG_slack = P_DG_slack + P_WT(loc_WT,:,:);
    end
    if is_BAT
        P_DG_slack = P_DG_slack + P_BAT(loc_BAT,:,:);
    end
    %slack power balance
    problem.Constraints.Pslack = P_slack == -sum(P(idx.pq,:,:),1) - P_load(idx.slack,:,:) - P_DG_slack + P_losses;
    problem.Constraints.Qslack = Q_slack == -sum(Q(idx.pq,:,:),1) - Q_load(idx.slack,:,:) - Q_DG_slack + Q_losses;
    problem.Constraints.Pbus_slack = P(idx.slack,:,:) == P_slack + P_load(idx.slack,:,:) + P_DG_slack;
    problem.Constraints.Qbus_slack = Q(idx.slack,:,:) == Q_slack + Q_load(idx.slack,:,:) + Q_DG_slack;
    problem.Constraints.inj_slack_constr = inj_slack == sum(load_gas,1) + sum(load_FC,1) - sum(feed_in,1) + sum(Qcharge,1) - sum(Qdischarge,1);
    
    %DG ressources
    %problem.Constraints.P_FC = P_FC == 0;
    problem.Constraints.PV = P_PV <= (repmat(reshape(scenarios.pv,1,n_timesteps,n_scenarios),n_PV,1,1)).*(Pmax.PV.*ones(1,n_timesteps,n_scenarios));%PV generation with optionnal curtailment
    problem.Constraints.WT = P_WT <= repmat(reshape(scenarios.wind,1,n_timesteps,n_scenarios),n_WT,1,1).*(Pmax.WT.*ones(1,n_timesteps,n_scenarios));
    %storage constraints
    constr_soc = optimconstr(n_BAT,n_timesteps,n_scenarios);
    constr_storage = optimconstr(n_storage,n_timesteps,n_scenarios);
    constr_soc(:,1,:) = BAT(:,1,:) == BAT0 - resolution.*P_BAT(:,1,:);
    constr_storage(:,1,:) = storage(:,1,:) == storage0 + resolution*efficiency.storage*Qcharge(:,1,:) - resolution*Qdischarge(:,1,:);
    for t=2:n_timesteps
        constr_soc(:,t,:) = BAT(:,t,:) == BAT(:,t-1,:) - resolution.*P_BAT(:,t,:);
        constr_storage(:,t,:) = storage(:,t,:) == storage(:,t-1,:) + resolution.*efficiency.storage.*Qcharge(:,t,:) - resolution*Qdischarge(:,t,:);
    end
    problem.Constraints.soc = constr_soc;
    problem.Constraints.storage_constr = constr_storage;
    problem.Constraints.feedin = feed_in == efficiency_feedin.*ones(1,n_timesteps,n_scenarios).*P_feedin.*Ab./1e6;
    problem.Constraints.Load_FC = load_FC == Ab./1e6.*P_FC./(efficiency_FC.*ones(1,n_timesteps,n_scenarios));
    
    
    % Set nondefault solver options
    if is_scheduling
        options = optimoptions("quadprog",'Display','none');
    else
        options = optimoptions("linprog",'Display','none');
    end
    % Solve problem
    [solution, tot_cost] = solve(problem,initialPoint,"Options",options);
    %store history of solution
    solution_hist = [solution_hist solution];
    lin_cost_hist = [lin_cost_hist tot_cost];

    %electrical loadflow:
    S = solution.P + j*solution.Q;
    V_lf = ones(n_buses,n_timesteps,n_scenarios);
    S_lf = ones(n_buses,n_timesteps,n_scenarios);
    for s=1:n_scenarios
        for t=1:n_timesteps
            [J_t,V0_t,S0_t,n_iter_t] = NR_polar_sol(S(:,t,s),E_star,Y,E_start,idx,Parameters);
            V_lf(:,t,s) = V0_t;
            S_lf(:,t,s) = S0_t;
            
            %fprintf('Initial Power loadflow reached convergence at iteration: %i \n', n_iter);
            S0_losses(:,t,s) = sum(S0_t);
        end
    end

    %gas loadflow:
    p_lf = zeros(n_nodes-1,n_timesteps,n_scenarios);
    q_lf = zeros(n_pipes,n_timesteps,n_scenarios);
    x0 = solution.p(2:end,:,:);
    for s=1:n_scenarios
        for t=1:n_timesteps
            load_t = -solution.inj(2:end,t,s);
            [p_lf_t, q_lf_t] = NR_weymouth(M,load_t,p_slack,x0(:,t,s),step_size,Parameters.tol_MW,Parameters.n_max,L,D,gas,f_curve,Re_vec);
            p_lf(:,t,s) = p_lf_t;
            q_lf(:,t,s) = q_lf_t;
        end
    end
    %recompute ojective function
    P_slack_lf = -sum(real(S_lf(idx.pq,:,:)),1) - P_load(idx.slack,:,:) - P_DG_slack + real(S0_losses);
    exact_cost = sum(P_slack_lf.*cost_power_slack.*(Ab/1e6),"all") + sum(solution.inj_slack.*cost_gas_slack,"all") + cost_FC.*sum(solution.P_FC,"all") + cost_feedin.*sum(solution.P_feedin,"all") + cost_STO.*sum(solution.Qcharge,"all") + sum(solution.c_BAT,"all");
    if is_scheduling
        for s=1:n_scenarios
            exact_cost = exact_cost + 1000*(Ab/1e6).^2*sum((solution.P_slack_scheduling-P_slack_lf(:,:,s)).^2,'all');
            %exact_cost = exact_cost + 1000*(solution.daily_inj_scheduling-sum(solution.inj_slack(:,:,s),'all')).^2;
            exact_cost = exact_cost + 1000*sum((solution.inj_scheduling-solution.inj_slack(:,:,s)).^2,'all');
        end
    end
    cost_hist = [cost_hist exact_cost];
    

    %relative missmatch
    delta_p = max((abs(p_lf - x0)./p_slack),[],"all");
    delta_v = max(abs(abs(V_lf) - solution.Vmag)./1,[],"all"); %P slack is 1 pu.
    delta_slack = max(abs(abs(S_lf(idx.slack) - S(idx.slack))./Ab),[],"all");

    
    fprintf("Miss-match at iteration %i \n", iter)
    fprintf("pressure: %f \n", delta_p)
    fprintf("voltage: %f \n", delta_v)
    fprintf("Power at slack: %f \n", delta_slack)
    fprintf("\n")
    
    if (delta_v <= Parameters.tol_rel) && (delta_slack <= Parameters.tol_rel./10) && (delta_p <= Parameters.tol_rel)
        fprintf("optimisation convergence at iteration: %i\n", iter)
        break
    end
    
    %recover solution and update initial values
    %electrical grid
    V0 = V_lf;
    S0 = S_lf;
    P0 = real(S0);
    Q0 = imag(S0);
    S0_losses = sum(S0);
    P0_losses = real(S0_losses);
    Q0_losses = imag(S0_losses);
    P0_FC = solution.P_FC;
    P0_feedin = solution.P_feedin;
    P0_PV = solution.P_PV;
    P0_WT = solution.P_WT;
    %gas network
    p0 = [p_slack*ones(1,n_timesteps,n_scenarios);p_lf];
    q0 = q_lf;
    feed_in0 = solution.feed_in;
    load_FC0 = solution.load_FC;
    inj0 = solution.inj;
    %update C
    Re = max(abs(4*gas.rho_st/(pi*gas.mu)*(abs(q0)./D)*gas.MW_to_m3), 1e3*ones(n_pipes,1));
    f = zeros(n_pipes,n_timesteps,n_scenarios);
    for s=1:n_scenarios
        f(:,:,s) = AGA(Re_vec,f_curve,Re(:,:,s));
    end
    
    C = gas.m3_to_MW*13.2989.*(gas.T_st/gas.p_st).*(1./(L.*gas.d_rel.*gas.T_avg.*gas.z_avg)).^0.5.*(1./f).^0.5.*D.^(2.5)*1e5;
    for s=1:n_scenarios
        for t=1:n_timesteps
            [K_p(:,:,t,s),K_com_p,K_q(:,:,t,s),K_comp_q] = Coeffs_Voltage_Alpha(Y,S0(:,t,s),V0(:,t,s),Res_nodes,idx.slack,nph,alphap,alphaq,Elin_0);
            [C_rp(:,:,t,s) , C_rq(:,:,t,s), C_xp(:,:,t,s), C_xq(:,:,t,s)] = Coeffs_Losses(Y, V0(:,t,s), K_com_p, K_comp_q, idx.pq);
            [KI_p(:,:,t,s), KI_q(:,:,t,s), I0(:,t,s)] = Coeffs_Currents(YYL,YYT,V0(:,t,s),K_com_p,K_comp_q,nph,lines);
            K_inj(:,:,t,s) = pressure_coeff(M,C(:,t,s),p0(:,t,s));
        end
    end

    %step size decrease for convergence:
    step_p = 0.80 * step_p;
    step_v = 0.80 * step_v;

   
    % Clear variables
    clearvars P Q P_slack Q_slack P_losses Q_losses P_FC P_PV P_WT P_BAT BAT I Vmag p q inj_slack feed_in load_FC inj storage c_elec c_has initialPoint options
end

%% recover solution:
P = solution.P*Ab/1e6;
Q = solution.Q*Ab/1e6;
S = P + j*Q;
P_FC = solution.P_FC*Ab/1e6;
P_PV = solution.P_PV*Ab/1e6;
P_WT = solution.P_WT*Ab/1e6;
P_BAT = solution.P_BAT*Ab/1e6;
Q_sto = -solution.Qcharge + solution.Qdischarge;
BAT = solution.BAT*Ab/1e6;
STO = solution.storage;
P_feedin = solution.P_feedin*Ab/1e6;
P_slack = solution.P_slack*Ab/1e6;
Q_slack = solution.Q_slack*Ab/1e6;
if is_scheduling
    P_slack_scheduling = solution.P_slack_scheduling*Ab/1e6;
    %inj_daily = solution.daily_inj_scheduling;
    inj_slack_scheduling = solution.inj_scheduling;
end
line_losses = sum(S)*Ab/1e6;
Vmag = abs(V_lf);
p = [p_slack.*ones(1,n_timesteps,n_scenarios);p_lf];
q = q_lf;
feed_in = solution.feed_in;
load_FC = solution.load_FC;
inj = solution.inj;
inj_slack = solution.inj_slack;

%% Saving results:
[~,~] = mkdir("./results",run_name);
save('./results/' + run_name + "/" + run_name)

%% Loading results
%load('./results/' + run_name + "/" + run_name)
%% load previous run (for plotting)
%load('./results/' + "mincost_64scenarios/mincost_64scenarios")
%load('./results/' + "scheduling_64scenarios_timesteps/scheduling_64scenarios_timesteps")
%load('./results/' + "scheduling_64scenarios_decoupled/scheduling_64scenarios_decoupled")
%load('./results/' + "scheduling_64scenarios/scheduling_64scenarios")
%% recover costs
tbl_results = table();
tbl_results(1,"Power_slack") = {mean(sum(P_slack,2).*resolution)}; %avg. daily power consumption at the slack
tbl_results(1,"Injection_slack") = {mean(sum(inj_slack,2).*resolution)}; %avg. daily gas consumption at the slack
tbl_results(1,"c_elec") = {sum(solution.c_elec,'all')/n_scenarios}; %avg. daily power costs
tbl_results(1,"c_gas") = {sum(solution.c_gas,'all')/n_scenarios}; %avg. daily gas costs
operationnal_costs = sum(solution.c_gas,2) + sum(solution.c_elec,2) + sum(cost_FC.*sum(P_FC,1),2) + sum(cost_feedin.*sum(P_feedin,1),2) + sum(cost_STO.*sum(solution.Qcharge,1),2)+ sum(sum(solution.c_BAT,1),2);
tbl_results(1,"c_op") = {mean(operationnal_costs)}; %avg. daily operationnal costs


deviation_powerslack = 0;
deviation_gas_slack = 0;
deviation_gas_slack_daily = 0;
if is_scheduling == true
    for s=1:n_scenarios
        deviation_powerslack = deviation_powerslack + sum((P_slack_scheduling-P_slack(:,:,s)).^2);
        if exist('daily_inj_scheduling','var') == 1
            deviation_gas_slack = deviation_gas_slack_daily + abs(solution.daily_inj_scheduling-sum(inj_slack(:,:,s),'all'))*n_timesteps*resolution;
        else
            deviation_gas_slack = deviation_gas_slack + sum(abs(inj_slack_scheduling-inj_slack(:,:,s))*resolution);
        end
    end
    tbl_results(1,"dev_power") = {deviation_powerslack/n_scenarios}; %avg. deviation of the slack power
    tbl_results(1,"dev_gas") = {deviation_gas_slack/n_scenarios}; %avg. deviation of the slack gas injection at each timestep
end

%save above metrics in a csv
writetable(tbl_results,"./results/"+run_name+"/"+run_name+"_tbl_results.csv")
%% plotting results:
%% voltage magnitudes:
figure("Name","Mean Bus voltage magnitude")
boxplot(reshape(Vmag,n_buses,n_timesteps*n_scenarios)','Whisker',10)
hold on


yline(1.05,':k',"LineWidth",1)
yline(0.95,':k',"LineWidth",1)
hold off
ylim([0.92 1.07])
xlabel("Grid buse")
ylabel("|V| [pu.]")
title("Distribution of bus |V| across scenarios")
legend('boundaries', Location='northwest')

exportgraphics(gcf,'./results/'+ run_name + "/" +run_name+'_voltage.pdf','ContentType','vector')

%% pressure profiles:
figure("Name","Node pressure")

boxplot(reshape(p,n_nodes,n_timesteps*n_scenarios)')
hold on
yline(1,":k","LineWidth",1)
yline(min(pmin,[],1),":k","LineWidth",1)
xlim([0.5, n_nodes+0.5])
ylim([0.4 1.05])
legend('boundaries', Location='northwest')
xlabel('Gas network node')
ylabel('p [bar]')
title("Distribution of nodal pressure across scenarios")
exportgraphics(gcf,'./results/'+ run_name + "/" +run_name+'_pressure.pdf','ContentType','vector')

%% generation profiles
figure("Name","Power generation")
hours2 = [profile.hours, fliplr(profile.hours)];
plot(profile.hours,mean(P_WT,3),'Color','#A9A9A9')
hold on
plot(profile.hours,mean(sum(P_PV,1),3),"Color","#EDB120")
plot(profile.hours,mean(sum(P_FC,1),3),"b")
min_FC = min(sum(P_FC,1),[],3);
max_FC = max(sum(P_FC,1),[],3);
inBetween = [min_FC, fliplr(max_FC)];
%if we want to shade the region of variation across scenario (plot less
%readable)
%fill(hours2, inBetween,[0 1 1],'FaceAlpha',.1,'LineStyle','none');


plot(profile.hours,mean(sum(P_BAT,1),3),'Color',[0.4940 0.1840 0.5560])
min_BAT = min(sum(P_BAT,1),[],3);
max_BAT = max(sum(P_BAT,1),[],3);
inBetween = [min_BAT, fliplr(max_BAT)];
%fill(hours2, inBetween,[1 0 1],'FaceAlpha',.1,'LineStyle','none');

plot(profile.hours,-mean(sum(P_feedin,1),3),'Color','green')

if is_scheduling
    plot(profile.hours,P_slack_scheduling,'Color','red')
else
    plot(profile.hours,mean(P_slack,3),'Color','red')
end

min_slack = min(P_slack,[],3);
max_slack = max(P_slack,[],3);
inBetween = [min_slack, fliplr(max_slack)];
fill(hours2, inBetween,[0.6350 0.0780 0.1840],'FaceAlpha',.2,'LineStyle','none');

xlim([0,23.75])
ylim([-5,25])
ylabel("P [MW]")
xlabel("hour of day")
title("Distributed Generation (avg. of scenarios)")
legend("Wind","PV","FC","Bat.","P2G",'P slack','range');
%legend("Wind","PV","Bat.",'P slack','range');
exportgraphics(gcf,'./results/'+ run_name + "/" +run_name+'_power.pdf','ContentType','vector')
%% gas demand:
figure("Name","Nodal gas injections")
plot(profile.hours,-mean(sum(load_FC,1),3),'Color','b')
hold on
plot(profile.hours,mean(sum(feed_in,1),3),"Color","g")
plot(profile.hours,mean(sum(Q_sto,1),3),'Color',[0.4940 0.1840 0.5560])
plot(profile.hours,mean(inj_slack,3),"Color",'r')
plot(profile.hours,inj_slack_scheduling,"Color",'r')


min_slack = min(inj_slack,[],3);
max_slack = max(inj_slack,[],3);
inBetween = [min_slack, fliplr(max_slack)];
fill(hours2, inBetween,[0.6350 0.0780 0.1840],'FaceAlpha',.2,'LineStyle','none');

xlim([0,23.75])
ylim([-8 18])
ylabel("P [MWth]")
xlabel("hour of day")
title("Gas Network Injection (avg. of scenarios)")
legend("Fuel Cell","P2G","Storage","Slack scheduling","Slack range");
%legend("Fuel Cell","P2G","Storage","avg. Slack","Slack range");
%legend("Storage","avg. Slack","Slack range");

exportgraphics(gcf,'./results/'+ run_name + "/" +run_name+'_gas.pdf','ContentType','vector')
%% Battery charge
figure("Name","Battery charge")
plot(profile.hours,mean(sum(BAT,1),3))
hold on
min_bat = min(sum(BAT,1),[],3);
max_bat = max(sum(BAT,1),[],3);
inBetween = [min_bat, fliplr(max_bat)];
fill(hours2, inBetween,[0 1 1],'FaceAlpha',.15,'LineStyle','none');
plot(profile.hours,mean(STO,3))
min_STO = min(STO,[],3);
max_STO = max(STO,[],3);
inBetween = [min_STO, fliplr(max_STO)];
fill(hours2, inBetween,[0.9290 0.6940 0.1250],'FaceAlpha',.1,'LineStyle','none');
yline(sum(BAT_max*Ab/1e6),'--','Color',[0 0.4470 0.7410])
yline(10,'--r')
ylabel("E [MWh]")
xlabel("hour of day")
xlim([0,23.75])
ylim([0,11])
title("Charge Profile of the Battery (avg. of scenarios)")
legend('Battery charge','Battery range','Gas storage','Storage range','Max Battery', 'Max Tank',location = 'northwest')
exportgraphics(gcf,'./results/'+ run_name + "/"+run_name+'_bat.pdf','ContentType','vector')
%% profiles plotting
figure("Name", "load and generation profiles")
min_pv = min(profile.pv,[],2);
max_pv = max(profile.pv,[],2);
mean_pv = mean(profile.pv,2);
plot(profile.hours,mean_pv,"Color","#EDB120")
hold on
hours2 = [profile.hours, fliplr(profile.hours)];
inBetween = [max_pv', fliplr(min_pv')];
fill(hours2, inBetween,[1 1 0],'FaceAlpha',.2,'LineStyle','none');

min_load = min(profile.load,[],2);
max_load = max(profile.load,[],2);
mean_load = mean(profile.load,2);
plot(profile.hours,mean_load,"b")
inBetween = [min_load', fliplr(max_load')];
fill(hours2, inBetween,[0 1 1],'FaceAlpha',.2,'LineStyle','none');

min_wind = min(profile.wind,[],2);
max_wind = max(profile.wind,[],2);
mean_wind = mean(profile.wind,2);
plot(profile.hours,mean_wind,'Color','#A9A9A9')
inBetween = [min_wind', fliplr(max_wind')];
fill(hours2, inBetween,[0.5 0.5 0.5],'FaceAlpha',.2,'LineStyle','none');

min_loadgas = min(profile.loadgas,[],2);
max_loadgas = max(profile.loadgas,[],2);
mean_loadgas = mean(profile.loadgas,2);
plot(profile.hours,mean_loadgas,"Color",[0.4940, 0.1840, 0.5560])
inBetween = [min_loadgas', fliplr(max_loadgas')];
fill(hours2, inBetween,[0.75, 0, 0.75],'FaceAlpha',.2,'LineStyle','none');



xlim([0,23.75])
ylim([0,1])
xlabel("hour of day")
ylabel("% max")
title("Hourly Load and Generation Profiles (avg. & range of scenarios)")
legend("PV","","Electrical loads", "", "Wind", "","Gas load","",'Location','northeast')
exportgraphics(gcf,'./results/'+ run_name + "/"+run_name+'_profile.pdf','ContentType','vector')
%% price profiles
figure("Name", "costs")
plot(profile.hours,profile.price_el)
hold on
plot(profile.hours,profile.price_gas)
xlim([0,23])
ylim([0,150])
xlabel("hour of day")
ylabel("Eur/MWh")
title("Hourly Day-ahead Prices")
legend('Power spot price', 'Gas intraday price','Location','northeast')
exportgraphics(gcf,'./results/'+ run_name + "/"+run_name+'_costs.pdf','ContentType','vector')
%% evolution of objective function through iteration
figure("Name","cost iteration")
plot(cost_hist)
hold on
plot(lin_cost_hist)
xlim([1,iter])
xlabel("iteration")
ylabel("Objective function value")
title("Objective function evolution at each iteration")
legend("Exact loadflow", "Linear approximation")
exportgraphics(gcf,'./results/'+ run_name + "/"+run_name+'_costs_iter.pdf','ContentType','vector')
