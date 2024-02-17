% Main Script

clear 
clc

addpath('Channels/');
addpath('Differentials/');
addpath(genpath('SaveTools/')); 
addpath(genpath('Solvers/'));
addpath('PostProcess/');
addpath('data_management/');

rmpath('Channels/');
rmpath('Differentials/');
rmpath('SaveTools/'); addpath('SaveTools/EHL_Plots');
rmpath(genpath('Solvers/'));
rmpath('PostProcess/');
rmpath('data_management/');

%%

% EHL_Solver_Wi(ehl_solution, ehl_channel, velocity_list, wi_list, beta_list, load_list, "") 
EHL_Solver_De(velocity_list, de_list, beta_list, load_list, h_cent, "centh") 
% EHL_Solver_De_betadriven(velocity_list, de_list, beta_list, load_list, h_cent, "centh") 

%% Parabolic-Solver
clc
% clear 

solution_mesh = mesh_lubrication("N_x", 512, "N_y", 32);
parabolic_channel = parabolicslider("alfa", 0.1, "nodes", solution_mesh.N_x);
solution_mesh = solution_mesh.set_channel(parabolic_channel);
solution_mesh = solution_mesh.consturct_mesh_matrix(); solution_mesh = solution_mesh.set_Jacaboian();
 
fluid = fluid_domain("De", 0.01, "beta", 0.8);
sim_prob = simulation_proporties; %default simulation properties 
par_solution = HL_solution("simulation_proporties", sim_prob, "mesh", solution_mesh, "fluid", fluid);
par_solution = par_solution.LIN("cavitation", "on");
par_solution = par_solution.FieldProperties();
par_solution = par_solution.VR(6e-4, "CD2", "RK4", "cavitation", "on");
par_solution = par_solution.Calculate_Load();

%% Step Function
% clear

solution_mesh = mesh_lubrication("N_x", 512, "N_y", 32);
% step_channel = stepslider("a", 0.5, "e", 0.026, "mesh", solution_mesh);
step_channel = smoothstep("a", 0.5, "e", 0.026, "c", 400, "nodes", solution_mesh.N_x);

% step_channel  = pocket("type", "center", "c",  0.2, "curv",  300, "nodes", 512, "a", 0.9, "e", 0.02);
solution_mesh = solution_mesh.set_channel(step_channel);
solution_mesh = solution_mesh.consturct_mesh_matrix(); solution_mesh = solution_mesh.set_Jacaboian();
 
fluid = fluid_domain("De", 0.1, "beta", 0, "mesh", solution_mesh, "channel", step_channel);
sim_prob = simulation_proporties; %default simulation properties 
step_solution = solution("simulation_proporties", sim_prob, "mesh", solution_mesh, "fluid_domain", fluid);
step_solution = step_solution.LIN("cavitation", "off");
step_solution = step_solution.FieldProperties();
% step_solution = step_solution.VR(3e-3, "UW", "CN", "cavitation", "off");
% step_solution = step_solution.VR(1e-4, "UW", "RK4", "cavitation", "off");
% step_solution = step_solution.Calculate_Load();

%% Step Function 2
clear

a_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
e_list = [-0.025, -0.02, -0.015, -0.01, -0.005, 0.005, 0.010, 0.015, 0.02, 0.025];

solution_mesh = cell(length(a_list), length(e_list));
step_channel = cell(length(a_list), length(e_list));
fluid = cell(length(a_list), length(e_list));
step_solution = cell(length(a_list), length(e_list));

for i = 1:length(a_list)
    for j = 1:length(e_list)
solution_mesh{i, j} = mesh_lubrication("N_x", 512, "N_y", 64);
% step_channel{i, j} = stepslider("a", a_list(i), "e", e_list(j), "mesh", solution_mesh{i, j});
step_channel{i, j} = smoothstep("a", a_list(i), "e", e_list(j), "c", 400, "mesh", solution_mesh{i, j});
solution_mesh{i, j} = solution_mesh{i, j}.set_channel(step_channel{i, j});
solution_mesh{i, j} = solution_mesh{i, j}.consturct_mesh_matrix(); 
solution_mesh{i, j} = solution_mesh{i, j}.set_Jacaboian();
 
fluid{i, j} = fluid_domain("De", 0.1, "beta", 0, "mesh", solution_mesh{i, j}, "channel", step_channel{i, j});
sim_prob = simulation_proporties; %default simulation properties 
step_solution{i, j} = solution("simulation_proporties", sim_prob, "mesh", solution_mesh{i, j}, "fluid_domain", fluid{i, j});
step_solution{i, j} = step_solution{i, j}.LIN("cavitation", "off");
step_solution{i, j} = step_solution{i, j}.FieldProperties();
step_solution{i, j} = step_solution{i, j}.VR(3e-3, "UW", "CN", "cavitation", "off");
step_solution{i, j} = step_solution{i, j}.Calculate_Load();
    end
end
    

%%  Parabol 
clear
epsilon = [0.004];
Wi = [0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];

for i = 1:length(Wi)
    for j = 1:2
solution_mesh = mesh_lubrication("N_x", 512, "N_y", 32);
parabolic_channel = parabolicslider("alfa", 0.5, "nodes", solution_mesh.N_x);
solution_mesh = solution_mesh.set_channel(parabolic_channel);
solution_mesh = solution_mesh.consturct_mesh_matrix(); solution_mesh = solution_mesh.set_Jacaboian();
 
fluid{i, j} = fluid_domain("De", epsilon * Wi(i), "beta", 0.8, "mesh", solution_mesh, "channel", parabolic_channel);
sim_prob = simulation_proporties; %default simulation properties 
par_solution{i, j} = solution("simulation_proporties", sim_prob, "mesh", solution_mesh, "fluid_domain", fluid{i, j});
par_solution{i, j} = par_solution{i, j}.LIN("cavitation", "off");
par_solution{i, j} = par_solution{i, j}.FieldProperties();
if j == 1
    disp("----------------------RK4------------------------")
par_solution{i, j} = par_solution{i, j}.VR(1e-4, "CD2", "RK4", "cavitation", "off");
elseif j == 2
    disp("----------------------CN------------------------")
 par_solution{i, j} = par_solution{i, j}.VR(1e-3, "UW", "CN", "cavitation", "off");   
end
par_solution{i, j} = par_solution{i, j}.Calculate_Load();
    end
end


%% Rotine Load and Massflow Map

for i = 1:length(a_list)
    for j = 1:length(e_list)
LOAD_N(i, j) = step_solution{i, j}.load_dic("W_N");
LOAD_NN(i, j) = step_solution{i, j}.load_dic("W_NN");
LOAD_VR(i, j) = step_solution{i, j}.load_dic("W_VR");
p_ve = step_solution{i, j}.pressure_VR("p_ve");
p_n = step_solution{i, j}.pressure_VR("p_s") + step_solution{i, j}.pressure_VR("p_xy_N");
dpvedx = OneDcentraldiff(p_ve, step_solution{i, j}.mesh.dx, "CD6");
dpndx = OneDcentraldiff(p_n, step_solution{i, j}.mesh.dx, "CD6");
 q_n = -step_solution{i, j}.h.^3 .* dpndx; 
 q_ve =  -step_solution{i, j}.h.^3 .* dpvedx;
 Q_N(i, j) = q_n(16);
 Q_VE(i, j) = q_ve(16);
    end
end

LOAD_N_ext = interp2(LOAD_N, 3);
LOAD_NN_ext = interp2(LOAD_NN, 3);
LOAD_VR_ext = interp2(LOAD_VR, 3);
Q_N_ext = interp2(Q_N, 3);
Q_VE_ext = interp2(Q_VE, 3);
    

a_list_exten = linspace(0.1, 0.9, 65);
e_list_exten = linspace(-0.025, 0.025, 73);

figure(1);
surf(e_list_exten, a_list_exten, LOAD_N_ext, "LineStyle", "none")
xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
ylabel("\boldmath{$a$}","Interpreter", "latex")
set(figure(1),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
title("\boldmath{$L_{N}$}","Interpreter","latex")
colorbar; colormap("jet");

figure(2);
surf(e_list_exten, a_list_exten, LOAD_VR_ext - LOAD_N_ext, "LineStyle", "none")
xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
ylabel("\boldmath{$a$}","Interpreter", "latex")
set(figure(2),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
title("\boldmath{$L_{VR} - L_{N}$}","Interpreter","latex")
colorbar; colormap("jet");

figure(3);
surf(e_list_exten, a_list_exten, LOAD_VR_ext - LOAD_NN_ext, "LineStyle", "none")
xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
ylabel("\boldmath{$a$}","Interpreter", "latex")
set(figure(3),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
title("\boldmath{$L_{VR} - L_{LIN}$}","Interpreter","latex")
colorbar; colormap("jet");

figure(4);
surf(e_list_exten, a_list_exten, Q_VE_ext, "LineStyle", "none")
xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
ylabel("\boldmath{$a$}","Interpreter", "latex")
set(figure(4),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
title("\boldmath{$Q_{VE}$}","Interpreter","latex");
colorbar; colormap("jet");

figure(5);
surf(e_list_exten, a_list_exten, Q_N_ext, "LineStyle", "none")
xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
ylabel("\boldmath{$a$}","Interpreter", "latex")
set(figure(5),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
title("\boldmath{$Q_{N}$}","Interpreter","latex")
colorbar; colormap("jet");

%% Routine AllPressure for Diff. Loc.

for i = 0:2
    p_n = step_solution{1 + i*4, 10}.pressure_Newtonain;
    p_lin = step_solution{1 + i*4, 10}.pressure;
    p_vr = step_solution{1 + i*4, 10}.pressure_VR("p");
    figure(i+1); hold on; box on;
    set(figure(i+1),'Position',[10 90 850 720])
    set(gca,'FontSize',20,'TickLabelInterpreter','latex');
    xlabel("\boldmath{$x$}","Interpreter", "latex")
    plot(step_solution{1 + i*4, 10}.mesh.x, p_n, "LineWidth", 1.8)
    plot(step_solution{1 + i*4, 10}.mesh.x, p_lin, "LineWidth", 1.8)
    plot(step_solution{1 + i*4, 10}.mesh.x, p_vr, "LineWidth", 1.8)
end

%% Routine VRPressure for Diff. Loc.

for i = 0:2
    p_xx = step_solution{1 + i*4, 1}.pressure_VR("p_xx");
    p_xy = step_solution{1 + i*4, 1}.pressure_VR("p_xy_NN");
    p_ve = step_solution{1 + i*4, 1}.pressure_VR("p_ve");
    figure(i+1); hold on; box on;
    set(figure(1),'Position',[10 90 850 720])
    set(gca,'FontSize',20,'TickLabelInterpreter','latex');
    xlabel("\boldmath{$x$}","Interpreter", "latex")
    plot(step_solution{1 + i*4, 10}.mesh.x, p_xx, "LineWidth", 1.8)
    figure(2); hold on; box on;
    set(figure(2),'Position',[10 90 850 720])
    set(gca,'FontSize',20,'TickLabelInterpreter','latex');
    xlabel("\boldmath{$x$}","Interpreter", "latex")
    plot(step_solution{1 + i*4, 10}.mesh.x, p_xy, "LineWidth", 1.8)
    figure(3); hold on; box on;
    set(figure(3),'Position',[10 90 850 720])
    set(gca,'FontSize',20,'TickLabelInterpreter','latex');
    xlabel("\boldmath{$x$}","Interpreter", "latex")
    plot(step_solution{1 + i*4, 10}.mesh.x, p_ve, "LineWidth", 1.8)
end


%% Rotine LOAD w.r.t De

for i =1:13
   W_VR_RK(i) = par_solution{i, 1}.load_dic("W_VR");
   W_VR_CN(i) = par_solution{i, 2}.load_dic("W_VR");
   W_NN(i) =  par_solution{i, 2}.load_dic("W_NN");
   W_N(i) =  par_solution{i, 2}.load_dic("W_N");
end
Wi = [0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
figure(1); hold on; box on
plot(Wi, W_NN, "LineWidth", 1.8)
plot(Wi, W_VR_RK, "LineWidth", 1.8)
plot(Wi, W_VR_CN, "--","LineWidth", 1.8)
plot(Wi, W_N, "LineWidth", 1.8)
set(figure(1),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
legend("\bf{LIN}",  "\boldmath{$VR_{RK4}$}",  "\boldmath{$VR_{CN}$}", "\bf{RE}","box", "off", "Interpreter","latex")


%Humayuns Results

HUMAYUN = [0, 0; 
    1.9740315845053509, 0.010665048252712334;
12.925406806508901,0.02203906535231835;
15.50880814093025,0.024445843073787693;
17.839485431766903,0.026532248396571353;
22.262156194499106,0.029600491518312053;
25.28080449287188,0.03138792814210904;
28.369653914462624,0.03266531915605822;
31.45850333605337,0.034159024156326795;
34.54735275764411,0.03525881000256112;
37.63620217923486,0.036063548631432116;
40.7250516008256,0.03680488575741228;
43.81390102241635,0.03734680921787557;
46.9027504440071,0.03774528235056916;
49.99159986559784,0.03835096151226344;
53.080449287188586,0.038159694408570505;
56.16929870877934,0.038191572259186;
59.258148130370074,0.03806406085672405;
62.346997551960825,0.03785685482772337;
65.43584697355156,0.037617770948107226;
68.52469639514231,0.03729899244195235;
71.61354581673305,0.036948336085181985;
74.7023952383238,0.03651798510187289;
77.79124465991454,0.03616732874510252;
80.88009408150529,0.0357210388364857;
83.96894350309604,0.035306626778484354;
87.05779292468678,0.03490815364579075;
90.14664234627753,0.03470094761679009;
93.23549176786827,0.03438216911063521;
96.32434118945902,0.03433435233471198;
98.85158162530598,0.03433435233471198];

plot(HUMAYUN(:,1), HUMAYUN(:,2), "k:", "LineWidth", 1.7)

%% Routine VRStress 

T_xx =  step_solution{1,1}.stress_VR("tau_xx");
T_xy =  step_solution{1,1}.stress_VR("tau_xy");

figure(1); box on; hold on;
plot(step_solution{1, 10}.mesh.x, T_xx(1, :), "LineWidth", 1.8)
plot(step_solution{1, 10}.mesh.x, T_xx(8, :), "LineWidth", 1.8)
plot(step_solution{1, 10}.mesh.x, T_xx(24, :), "LineWidth", 1.8)
plot(step_solution{1, 10}.mesh.x, T_xx(32, :), "LineWidth", 1.8)
set(figure(1),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');

figure(2); box on; hold on;
plot(step_solution{1, 10}.mesh.x, T_xy(1, :), "LineWidth", 1.8)
plot(step_solution{1, 10}.mesh.x, T_xy(8, :), "LineWidth", 1.8)
plot(step_solution{1, 10}.mesh.x, T_xy(24, :), "LineWidth", 1.8)
plot(step_solution{1, 10}.mesh.x, T_xy(32, :), "LineWidth", 1.8)
set(figure(1),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');




