% Main Script

clear 
clc

addpath('Channels/');
addpath('Differentials/');
addpath(genpath('SaveTools/')); 
addpath(genpath('Solvers/'));
addpath('PostProcess/');
addpath('data_management/');
% 
% rmpath('Channels/');
% rmpath('Differentials/');
% rmpath('SaveTools/'); addpath('SaveTools/EHL_Plots');
% rmpath(genpath('Solvers/'));
% rmpath('PostProcess/');
% rmpath('data_management/');

%%

% EHL_Solver_Wi(ehl_solution, ehl_channel, velocity_list, wi_list, beta_list, load_list, "") 
EHL_Solver_De(velocity_list, de_list, beta_list, load_list, h_cent, "centh") 
% EHL_Solver_De_betadriven(velocity_list, de_list, beta_list, load_list, h_cent, "centh") 

%% Parabolic-Solver
clc
clear 

solution_mesh = mesh_lubrication("N_x", 512, "N_y", 32);
parabolic_channel = parabolicslider("alfa", 0.5, "nodes", solution_mesh.N_x);
% parabolic_channel = cosineslider("eccen", 0.2, "clearance", 0.8, "nodes", solution_mesh.N_x);
% parabolic_channel =    hinch("a", 0, "l", 1, "h_out", 0.9, "nodes", solution_mesh.N_x)                    
solution_mesh = solution_mesh.set_channel(parabolic_channel);
solution_mesh = solution_mesh.consturct_mesh_matrix(); solution_mesh = solution_mesh.set_Jacaboian();

phys_prob = physical_description("L", 0.015,"U", 1.5, "eta0", 0.00903 , "h0", 2.5e-5 , "lambda", 0.00134, "p0", 0 );
fluid = fluid_domain("De", 0.05, "beta", 0.4);
PIB_CI4 = fluid_domain("De",  phys_prob.lambda * phys_prob.U / phys_prob.L, "beta", 0.4, "name", "PIB/C14");
PEO_WG = fluid_domain("De",  phys_prob.lambda * phys_prob.U / phys_prob.L, "beta", 0.5, "name", "PEO_WG");


sim_prob = simulation_proporties; %default simulation properties 
par_solution = HL_solution("simulation_proporties", sim_prob, "mesh", solution_mesh, "fluid", PEO_WG, "PD",  phys_prob);
par_solution = par_solution.LIN("cavitation", "on");
par_solution = par_solution.FieldProperties();
% par_solution = par_solution.VR(6e-4, "UW", "RK4", "cavitation", "off", "semi_vel_correction", "on");
par_solution = par_solution.VR(3e-3, "UW", "CN", "cavitation", "on", "semi_vel_correction", "off");

par_solution = par_solution.Calculate_Load("on");
par_solution = par_solution.calculate_Friction();
par_solution = par_solution.convert_dimensional();

%% Step Function
% clear


cavitation_flag = "off";

solution_mesh = mesh_lubrication("N_x", 512, "N_y", 32);
% step_channel = stepslider("a", 0.5, "e", 0.026, "mesh", solution_mesh);
step_channel = smoothstep("a", 0.3, "e", 0.008, "c", 400, "nodes", solution_mesh.N_x);

% step_channel  = pocket("type", "center", "c",  0.2, "curv",  300, "nodes", 512, "a", 0.9, "e", 0.02);
solution_mesh = solution_mesh.set_channel(step_channel);
solution_mesh = solution_mesh.consturct_mesh_matrix(); solution_mesh = solution_mesh.set_Jacaboian();


phys_prob = physical_description("L", 0.015,"U", 1.5, "eta0", 0.00903 , "h0", 2.5e-5 , "lambda", 0.00134, "p0", 0 );
fluid = fluid_domain("De", 0.1 , "beta", 0.0);
PIB_CI4 = fluid_domain("De",  phys_prob.lambda * phys_prob.U / phys_prob.L, "beta", 0.5, "name", "PIB/C14");
PEO_WG = fluid_domain("De",  phys_prob.lambda * phys_prob.U / phys_prob.L, "beta", 0.5, "name", "PEO_WG");

sim_prob = simulation_proporties; %default simu +lation properties 

step_solution = HL_solution("simulation_proporties", sim_prob, "mesh", solution_mesh, "fluid", PEO_WG, "PD",  phys_prob);
step_solution.step_flag = "on"; 
step_solution = step_solution.LIN("cavitation", cavitation_flag);
step_solution = step_solution.FieldProperties();
step_solution = step_solution.VR(3e-3, "UW", "CN", "cavitation", cavitation_flag, "semi_vel_correction", "off");
% step_solution = step_solution.VR(1e-4, "UW", "RK4", "cavitation", "off", "semi_vel_correction", "off");
step_solution = step_solution.Calculate_Load(cavitation_flag);
step_solution = step_solution.calculate_Friction();
step_solution = step_solution.convert_dimensional();


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
solution_mesh{i, j} = mesh_lubrication("N_x", 512, "N_y", 32);
% step_channel{i, j} = stepslider("a", a_list(i), "e", e_list(j), "mesh", solution_mesh{i, j});
step_channel{i, j} = smoothstep("a", a_list(i), "e", e_list(j), "c", 400, "nodes", solution_mesh{i, j}.N_x);
solution_mesh{i, j} = solution_mesh{i, j}.set_channel(step_channel{i, j});
solution_mesh{i, j} = solution_mesh{i, j}.consturct_mesh_matrix(); 
solution_mesh{i, j} = solution_mesh{i, j}.set_Jacaboian();
 
fluid = fluid_domain("De", 0.4, "beta", 0);
sim_prob = simulation_proporties; %default simulation properties 
step_solution{i, j} = HL_solution("simulation_proporties", sim_prob, "mesh", solution_mesh{i, j}, "fluid", fluid);
step_solution{i, j}.step_flag = "off";
step_solution{i, j} = step_solution{i, j}.LIN("cavitation", "off");
step_solution{i, j} = step_solution{i, j}.FieldProperties();
disp("                   ")
disp(" i   =     " + string(i) + "           |            j   =     " + string(j)+ "   ")
disp("                   ")
step_solution{i, j} = step_solution{i, j}.VR(3e-3, "UW", "CN", "cavitation", "off", "semi_vel_correction", "off");
step_solution{i, j} = step_solution{i, j}.Calculate_Load();
    end
end

%% EHL-like

clc
clear

a_list = [0.2, 0.3, 0.4, 0.5];
l_list = [0.3, 0.4, 0.5];

channel_type = "sym";
midflat_solution = cell(length(a_list));

for i = 1:length(a_list)

a = a_list(i);
if channel_type == "sym"
    l  = (1 - a * 2);
elseif channel_type == "un-sym"
    l = l_list(i);
end

mesh = mesh_lubrication("N_x", 512, "N_y", 32);
midflat_channel = ehl_like("a", a, "l", (1 - a*2), "e", 0.02, "nodes", 512);                
mesh = mesh.set_channel(midflat_channel);
mesh = mesh.consturct_mesh_matrix(); mesh = mesh.set_Jacaboian();
 
fluid = fluid_domain("De", 0.1, "beta", 0.4);
sim_prob = simulation_proporties; %default simulation properties 
midflat_solution{i}= HL_solution("simulation_proporties", sim_prob, "mesh", mesh, "fluid", fluid);
midflat_solution{i} = midflat_solution{i}.LIN("cavitation", "off");
midflat_solution{i} = midflat_solution{i}.FieldProperties();
% par_solution = par_solution.VR(6e-4, "UW", "RK4", "cavitation", "off", "semi_vel_correction", "on");
midflat_solution{i} = midflat_solution{i}.VR(1e-3, "UW", "CN", "cavitation", "off", "semi_vel_correction", "off");
midflat_solution{i} = midflat_solution{i}.Calculate_Load(); 
    
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

a_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
e_list = [-0.025, -0.02, -0.015, -0.01, -0.005, 0.005, 0.010, 0.015, 0.02, 0.025];

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
colortype = "summer";

figure(1); hold on
contour_levels = linspace(min(LOAD_N_ext(:)), max(LOAD_N_ext(:)), 20);
colormap(colortype);
contourf(e_list_exten, a_list_exten, LOAD_N_ext, contour_levels,"Fill", "on" ,"LineColor", "k")
% surf(e_list_exten, a_list_exten, LOAD_N_ext, "LineStyle", "none")
xlim([-0.025 0.025])
xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
ylabel("\boldmath{$a$}","Interpreter", "latex")
set(figure(1),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
title("\boldmath{$L_{N}$}","Interpreter","latex")
colorbar; 

figure(2); hold on
% surf(e_list_exten, a_list_exten, LOAD_VR_ext - LOAD_N_ext, "LineStyle", "none")
colormap(colortype);
contourf(e_list_exten, a_list_exten, LOAD_VR_ext - LOAD_N_ext, 20,"Fill", "on" ,"LineColor", "k")
xlim([-0.025 0.025])
xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
ylabel("\boldmath{$a$}","Interpreter", "latex")
set(figure(2),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
title("\boldmath{$L_{VR} - L_{N}$}","Interpreter","latex")
colorbar; 


% figure(3); hold on
% surf(e_list_exten, a_list_exten, LOAD_VR_ext - LOAD_NN_ext, "LineStyle", "none")
% contour(e_list_exten, a_list_exten, LOAD_N_ext, 40, "LineColor", "k")
% xlim([-0.025 0.025])
% xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
% ylabel("\boldmath{$a$}","Interpreter", "latex")
% set(figure(3),'Position',[10 90 850 720])
% set(gca,'FontSize',20,'TickLabelInterpreter','latex');
% title("\boldmath{$L_{VR} - L_{LIN}$}","Interpreter","latex")
% colorbar; colormap(colortype);


figure(4); hold on
% surf(e_list_exten, a_list_exten, Q_VE_ext, "LineStyle", "none")
colormap(colortype);
contour_levels = linspace(min(Q_VE_ext(:)), max(Q_VE_ext(:)), 20);
contour(e_list_exten, a_list_exten, Q_VE_ext, contour_levels, "Fill", "on", "LineColor", "k", "LineStyle", "-")
xlim([-0.025 0.025])
xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
ylabel("\boldmath{$a$}","Interpreter", "latex")
set(figure(4),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
title("\boldmath{$Q_{VE}$}","Interpreter","latex");
colorbar; 


figure(5); hold on
% surf(e_list_exten, a_list_exten, Q_N_ext, "LineStyle", "none")
colormap(colortype);
contour_levels = linspace(min(Q_N_ext(:)), max(Q_N_ext(:)), 20);
contour(e_list_exten, a_list_exten, Q_N_ext, contour_levels, "Fill", "on", "LineColor", "k")
xlim([-0.025 0.025])
xlabel("\boldmath{$\Delta h$}","Interpreter", "latex")
ylabel("\boldmath{$a$}","Interpreter", "latex")
set(figure(5),'Position',[10 90 850 720])
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
title("\boldmath{$Q_{N}$}","Interpreter","latex")
colorbar; 

%% Routine TotalPressure for Diff. Loc.
n = 1;
ls = "-.";
colors = ["#0c2c84", "#990000", "#005824", "#6e016b"];

for i = 0:2
    p_n = step_solution{1 + i*4, n}.pressure_LIN("p0");
%     p_lin = step_solution{1 + i*4, n}.pressure;
    p_vr = step_solution{1 + i*4, n}.pressure_VR("p");
    p_ve = step_solution{1 + i*4, n}.pressure_VR("p_ve");
    figure(i+1); hold on; box on;
    set(figure(i+1),'Position',[0 0 850 720*1.4])
    set(gca,'FontSize',20,'TickLabelInterpreter','latex');
    xlabel("\boldmath{$x$}","Interpreter", "latex")
    plot(step_solution{1 + i*4, n}.mesh.x, p_n, ls, "LineWidth", 2.8, "color", colors(1))
%     plot(step_solution{1 + i*4, n}.mesh.x, p_ve, ls, "LineWidth", 1.8, "color", colors(2))
    plot(step_solution{1 + i*4, n}.mesh.x, p_vr, ls ,"LineWidth", 2.8, "color", colors(3))
    yline(0,':', "linewidth", 2)
end

%% Routine VE Pressure for Diff. Loc.
n = 10;
colors = ["#0c2c84", "#990000", "#005824", "#6e016b"];
ls = "-";

for i = 0:2
    p_xx = step_solution{1 + i*4, n}.pressure_VR("p_xx");
    p_xy = step_solution{1 + i*4, n}.pressure_VR("p_xy_NN");
    p_ve = step_solution{1 + i*4, n}.pressure_VR("p_ve");
    p_n = step_solution{1 + i*4, n}.pressure_LIN("p0");
    figure(1); hold on; box on;
    set(figure(1),'Position',[0 0 850 720])
    set(gca,'FontSize',23,'TickLabelInterpreter','latex');
    xlabel("\boldmath{$x$}","Interpreter", "latex")
    plot(step_solution{1 + i*4, n}.mesh.x, p_xx, "LineWidth", 2.8, "Color", colors(i+1))
    
    figure(2); hold on; box on;
    set(figure(2),'Position',[0 0 850 720])
    set(gca,'FontSize',23,'TickLabelInterpreter','latex');
    xlabel("\boldmath{$x$}","Interpreter", "latex")
    plot(step_solution{1 + i*4, n}.mesh.x, p_xy, ls,"LineWidth", 2.8, "Color", colors(i+1))
    
    figure(3); hold on; box on;
    set(figure(3),'Position',[0 0 850 720])
    set(gca,'FontSize',23,'TickLabelInterpreter','latex');
    xlabel("\boldmath{$x$}","Interpreter", "latex")
    plot(step_solution{1 + i*4, n}.mesh.x, p_ve, "LineWidth", 2.8,   "Color", colors(i+1))
end
legend("\bf{a = 0.1}",  "\bf{a = 0.5}",  "\bf{a = 0.9}", "box", "off", "Interpreter","latex")
%  ylabel("\boldmath{$p_{xx}$}","Interpreter", "latex")
% ylabel("\boldmath{$p_{xy}$}","Interpreter", "latex")
% ylabel("\boldmath{$p_{ve}$}","Interpreter", "latex")

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
n = 10; m = 5;
ls = "-";
colors = ["#0c2c84", "#990000", "#005824", "#6e016b"];
s = 3;

for i = 0:2
m = 1 + i*4;
s = i + 1;
T_xx =  step_solution{m, n}.stress_VR("tau_xx");
T_xy =  step_solution{m, n}.stress_VR("tau_xy");

figure(3); box on; hold on;
plot(step_solution{m, n}.mesh.x, trapz(T_xx) .* step_solution{m, n}.mesh.dy , ls, "LineWidth", 2.8, "Color", colors(s))
set(figure(3),'Position',[0 0 850 720])
set(gca,'FontSize',23,'TickLabelInterpreter','latex');
ylabel("\boldmath{$\tau_{xx}$}","Interpreter", "latex")
xlabel("\boldmath{$x$}","Interpreter", "latex")

figure(4); box on; hold on;
plot(step_solution{m, n}.mesh.x, trapz(T_xy) .* step_solution{m, n}.mesh.dy , ls, "LineWidth", 2.8, "Color", colors(s))
set(figure(4),'Position',[0 0 850 720])
set(gca,'FontSize',23,'TickLabelInterpreter','latex');
 ylabel("\boldmath{$\tau_{xy}$}","Interpreter", "latex")
 xlabel("\boldmath{$x$}","Interpreter", "latex")
end

%% Routine Mid-Flat-Channel

ls = "-";
colors = ["#0c2c84", "#990000", "#005824", "#6e016b"];
loads = [ ];
ms = 20;
for i = 0:3
s = i + 1;
i = i +1;
figure(1); hold on; box on;
plot(mesh.x, midflat_solution{i}.pressure_VR("p"), ls, "LineWidth", 2.8, "Color", colors(s))
set(figure(1),'Position',[0 0 850 720])
set(gca,'FontSize',23,'TickLabelInterpreter','latex');
ylabel("\boldmath{$p$}","Interpreter", "latex")
xlabel("\boldmath{$x$}","Interpreter", "latex")


figure(2); hold on; box on;
plot(mesh.x, midflat_solution{i}.mesh.channel.height, ls, "LineWidth", 2.8, "Color", colors(s))
set(figure(2),'Position',[0 0 850 720])
set(gca,'FontSize',23,'TickLabelInterpreter','latex');
ylabel("\boldmath{$h$}","Interpreter", "latex")
xlabel("\boldmath{$x$}","Interpreter", "latex")

T_xx =  midflat_solution{i}.stress_VR("tau_xx");
T_xy =  midflat_solution{i}.stress_VR("tau_xy");

figure(3); box on; hold on;
plot(midflat_solution{i}.mesh.x, trapz(T_xx) .* midflat_solution{i}.mesh.dy , ls, "LineWidth", 2.8, "Color", colors(s))
set(figure(3),'Position',[0 0 850 720])
set(gca,'FontSize',23,'TickLabelInterpreter','latex');
ylabel("\boldmath{$\tau_{xx}$}","Interpreter", "latex")
xlabel("\boldmath{$x$}","Interpreter", "latex")

figure(4); box on; hold on;
plot(midflat_solution{i}.mesh.x, trapz(T_xy) .* midflat_solution{i}.mesh.dy , ls, "LineWidth", 2.8, "Color", colors(s))
set(figure(4),'Position',[0 0 850 720])
set(gca,'FontSize',23,'TickLabelInterpreter','latex');
 ylabel("\boldmath{$\tau_{xy}$}","Interpreter", "latex")
 xlabel("\boldmath{$x$}","Interpreter", "latex")
 

loads = [loads midflat_solution{i}.load_dic("W_VR")];
       
end

 figure(5); box on; hold on;
plot(a_list, loads , "-square", "LineWidth", 2.8, "Color", colors(2), "Markersize", ms)
set(figure(5),'Position',[0 0 850 720])
set(gca,'FontSize',23,'TickLabelInterpreter','latex');
 ylabel("\boldmath{$W$}","Interpreter", "latex")
 xlabel("\boldmath{$Location of flatness$}","Interpreter", "latex")
 
 
%% Routine corrected U comparison

sol_int = step_solution;
mesh = solution_mesh;

uD = sol_int.velocity_field("uD");
% uC = sol_int.corrected_velocity_field("u");
uN = sol_int.velocity_field("u");

%Without correction, recalculation of u velocity
h = mesh.channel.height;
dpdx = OneDcentraldiff(sol_int.pressure_VR("p"), mesh.dx, "CD6");
beta = fluid.beta;
T_xx = sol_int.stress_VR("tau_xx");
gradTxx = TwoDcentraldiff(T_xx, mesh.dx, mesh.dy); DTxxDx = gradTxx{1} + gradTxx{2} .* mesh.Jacobian("dydx");
% T_xy =   +(sol_int.stress_VR("tau_xy") - (1 - beta) * sol_int.velocity_field("dudy"));
T_xy =   +(sol_int.stress_VR("tau_xy"));
[m, n] = size(T_xy);
dy = mesh.dy;
for i = 1:m
    dydy(i,:) = mesh.dy;
end
TD = mesh.Y ./ h; % Transformed Grid 
I11 = TD .* (dy .* trapz(dydy .*cumtrapz(flip(DTxxDx)))) - dydy .* cumtrapz(dydy .* cumtrapz(flip(DTxxDx)));
I12 = TD .* (dy .* trapz(flip(T_xy))) -  dydy .* cumtrapz((flip(T_xy)));

dpNdx = OneDcentraldiff(sol_int.pressure_LIN("p0"), mesh.dx, "CD2");

u_postprocess = 1 / 2 * dpdx .* (mesh.Y.^2  -  mesh.Y .* h)+ (1 - mesh.Y ./ h) ...
    + I11  + I12; 
u_postprocess_pois = 1 / 2 * dpdx .* (mesh.Y.^2  -  mesh.Y .* h) + (1 - mesh.Y ./ h);
%

figure(1); 
subplot(1,3, 1)
surf(u_postprocess - u_postprocess_pois, "linestyle", "none");
title("uC - uN")

subplot(1,3, 2)
surf(uN - u_postprocess, "linestyle", "none");
title("uN - upostprocess")

subplot(1, 3, 3)
surf((uN - u_postprocess_pois) ./ uN * 100, "linestyle", "none");
title("(uN - upostproces_p")


figure(2);
subplot(1, 1, 1); 
plot(mesh.x, sol_int.pressure_VR("p"), "--", "LineWidth" , 1.8)
hold on;
plot(mesh.x, sol_int.pressure_LIN("p"), "-.", "LineWidth" , 1.8)
plot(mesh.x, sol_int.pressure_LIN("p0"), "k:", "LineWidth" , 1.4)
legend("$p_{VR}$", "$p_{LIN}$", "$p_N$", "box", "off", "Interpreter", "latex")
% set(figure(2),'Position',[10 90 850 720])


figure(3);
subplot(1, 1, 1); 
plot(mesh.x, dpdx, "--", "LineWidth" , 1.8)
hold on;
plot(mesh.x, dpNdx, "k:", "LineWidth" , 1.4)
legend("$dpdx_{VR}$", "$dpdx_N$", "box", "off", "Interpreter", "latex")
title("Pressure Derivatives")
set(figure(2),'Position',[10 90 850 720])

Q_D = trapz(uD) .* dy;
q_d = trapz(uN) .* dy;

figure(4);
plot(mesh.x, Q_D, "--", "LineWidth" , 1.8); hold on
plot(mesh.x, q_d, "-.", "LineWidth" , 1.8)
title("Mass Flow Rates")
legend("$Q^D$", "$Q^N$", "box", "off", "Interpreter", "latex")

figure(5); hold on;
plot(mesh.x, u_postprocess(5,:), "--", "LineWidth" , 1.8); 

plot(mesh.x, u_postprocess_pois(5,:), ":", "LineWidth" , 1.8); 
plot(mesh.x, uN(5,:), "-.", "LineWidth" , 1.8)
plot(mesh.x, uN(5,:) + uD(5,:) * sol_int.fluid.De, "-.", "LineWidth" , 1.8)
legend("Post U with stress", "Post U without stress", "$U^N$", " $U^N + De \ U^D$", "box", "off", "Interpreter", "latex")


