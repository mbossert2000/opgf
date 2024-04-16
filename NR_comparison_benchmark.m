clear
clc
close all
%% run parameters:
run_name = "";
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
lines = table2array(readtable("./data/lines_benchmark.csv"));
n_lines = height(lines);
table_buses = readtable("./data/buses_benchmark.csv");
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
    from_idx = table_nodes{table_pipes{p,"From"},"idx"};
    to_idx = table_nodes{table_pipes{p,"To"},"idx"};
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
n_timesteps = 1;
resolution = 1; %resolution of timestep
profile.hours = [0:resolution:n_timesteps/resolution-resolution];

%electrical load
profile.load = table2array(readtable("./data/load_scenarios.csv")); %percentage of max value
profile.loadgas = table2array(readtable("./data/gasload_scenarios.csv"));
%profile.pv = table2array(readtable("./data/scenarios/pv_scenarios.csv")); %percentage of max value
%profile.wind = table2array(readtable("./data/scenarios/wind_scenarios.csv")); %percentage of max value

profile.price_elec = table2array(readtable("./data/elpricenextday.csv"));
profile.price_gas = table2array(readtable("./data/gaspricenextday.csv"));

P = zeros(n_buses*n_timesteps,1);
Q = zeros(n_buses*n_timesteps,1);
Qdot = zeros(n_nodes*n_timesteps,1);
for t=1:n_timesteps
    P((t-1)*n_buses+1:t*n_buses,:) = -table_buses{:,"P_load"}.*profile.load(t).*resolution*1000/Ab;
    Q((t-1)*n_buses+1:t*n_buses,:) = -table_buses{:,"Q_load"}.*profile.load(t).*resolution*1000/Ab;
    Qdot((t-1)*n_nodes+1:t*n_nodes,:) = -table_nodes{:,"demand"}.*profile.loadgas(t).*resolution;
end


%% Initial loadflow
S_star = P + j*Q;
E_0 = ones(n_buses,1);
Qdot_star = Qdot;
Qdot_star(1:n_nodes:n_timesteps) = [];
p_slack = 0.7;
p_rand = 0.6*ones(n_nodes-1,1)+rand(n_nodes-1,1)*0.1;
p0 = repmat([p_slack;p_rand],n_timesteps,1);
Parameters.tol = 1e-6;
Parameters.n_max = 50;
Parameters.step_size = 0.1;
Parameters.max_jac = 1e10;
Parameters.step_size = 0.4;
Parameters.n_timesteps = n_timesteps;

%% separate NR
tic
for i=1:100
    [J,E,S,n_iter] = NR_polar(S_star,Y,E_0,idx,Parameters);
    [p,q] = NR_weymouth(M,Qdot_star,p0,L,D,gas,f_curve,Re_vec,Parameters);
end
toc

%% joined NR
tic
for i=1:100
    [J,E,S,p,q,n_iter] = NR_coupled(S_star,Y,E_0,idx,M,Qdot_star,p0,L,D,gas,f_curve,Re_vec,Parameters);
end
toc
%% Conclusion:
% I will keep NR gas and power separates !

%% Comparison timestep together vs. separate

