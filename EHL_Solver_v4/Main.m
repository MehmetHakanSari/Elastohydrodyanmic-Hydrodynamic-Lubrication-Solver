% Main Script

clear 
clc

addpath('Differentials/');
addpath('SaveTools/'); addpath('SaveTools/EHL_Plots');
addpath('Solvers/');
addpath('PostProcess/');
addpath('data_management/');

% velocity_list = [0.00056234 ,0.0031623, 0.017783, 0.1, 1, 2.5];
% velocity_list = [0.0031623, 0.017783];
velocity_list = [0.00056234];
% velocity_list = [0.0001, 0.00056234, 0.0031623, 0.017783];
% velocity_list = [0.1];
% velocity_list = [1];
% velocity_list = [0.1, 1];
% wi_list = [5, 10];
de_list = [0, 0.01, 0.02, 0.05];
% beta_list = [0.8, 0.6, 0.4];
beta_list = [0.6];
% load_list = [10];
% load_list = [1];
load_list = [30, 50, 100, 200];
load_list = [1, 10, 30, 50];

% load_list = [100, 400];


%% Initialization of matrices

% h_cent = [1.8144e-07,5.4274e-07,1.8297e-06, 8.7596e-06,4.6202e-05]; %from U = 0.0001 to U = 0.1
h_cent = [4.6202e-05, 5.1512e-05]; %from U = 0.0001 to U = 0.1
% h_cent = [1.8144e-07];
% h_cent = [5.1512e-05];
% h_cent = [7.9147e-06, 7.9147e-06; 7.9147e-06, 7.9147e-06; 7.9147e-06, 7.9147e-06];
% h_cent = [7.9147e-06, 6.6786e-06, 5.4467e-06, 4.5244e-06; 5.0589e-05, 3.0545e-05, 2.3514e-05, 1.9840e-05];  %U = 1 and 0.1 
% h_cent = [4.3634e-07, 1.8436e-07, 1.3657e-07, 1.2129e-07; 2.3309e-06, 1.8309e-06, 1.2376e-06, 1.0778e-06];  %U = 0.0001 and 0.3162
h_cent = [4.0718e-07, 1.0778e-06, 8.7987e-07, 7.4273e-07; 3.8150e-06, 3.1834e-06, 2.6904e-06, 2.1860e-06];
h_cent = [1.6910e-06, 5.4922e-07, 4.0718e-07, 3.5394e-07];

% h_cent_F100_20bH_U01 = [7.9076e-06];
% h_cent_F200_20bH_U01 = [6.6786e-06];
% h_cent_F400_20bH_U01 = [5.4395e-06];
% h_cent_F800_20bH_U01 = [4.5244e-06];
% 
% h_cent_F100_25bH_U01 = [7.9147e-06];
% h_cent_F200_25bH_U01 = [6.6853e-06];
% h_cent_F400_25bH_U01 = [5.4467e-06];
% h_cent_F800_25bH_U01 = [4.5399e-06];
% 
% h_cent_F100_30bH_U1 = [5.0589e-05];
% h_cent_F200_30bH_U1 = [3.0545e-05];
% h_cent_F400_30bH_U1 = [2.3514e-05];
% h_cent_F800_30bH_U1 = [1.9840e-05];
% h_cent_F1600_30bH_U1 = [1.6155e-05];
% h_cent_F1600_30bH_U1 = [1.3665e-05];

% h_cent_F1_10bH_U1e-4 = [4.3634e-07];
% h_cent_F10_10bH_U1e-4 = [1.8436e-07];
% h_cent_F30_10bH_U1e-4 = [1.3657e-07];
% h_cent_F50_10bH_U1e-4 = [1.2129e-07];

% h_cent_F10_20bH_U3162e-3 = [1.8309e-06];
% h_cent_F30_20bH_U3162e-3 = [1.2376e-06];
% h_cent_F50_20bH_U3162e-3 = [1.0778e-06];
% h_cent_F100_20bH_U3162e-3 = [8.7987e-07];
% h_cent_F200_20bH_U3162e-3 = [7.4273e-07];

% h_cent_F1_12bH_U56e-4 = [1.6910e-06];
% h_cent_F10_12bH_U56e-4 = [5.4922e-07];
% h_cent_F30_12bH_U56e-4 = [4.0718e-07];
% h_cent_F50_12bH_U56e-4 = [3.5394e-07];
% h_cent_F10_12bH_U56e-4 = [2.9933e-07];
% h_cent_F200_12bH_U56e-4 = [2.6419e-07];

% h_cent_F30_22bH_U177e-2 = [3.8150e-06];
% h_cent_F50_22bH_U177e-2 = [3.1834e-06];
% h_cent_F10_22bH_U177e-2 = [2.6904e-06];
% h_cent_F200_22bH_U177e-2 = [2.1860e-06];


% h_cent = [1.9840e-05];

% h_cent = [];

%%

% EHL_Solver_Wi(ehl_solution, ehl_channel, velocity_list, wi_list, beta_list, load_list, "") 
EHL_Solver_De(velocity_list, de_list, beta_list, load_list, h_cent, "centh") 
% EHL_Solver_De_betadriven(velocity_list, de_list, beta_list, load_list, h_cent, "centh") 

%% Seperate-Solver
clc
clear 

parabolic_slider = domain;
parabolic_slider = parabolic_slider.set_length(1e-3);
parabolic_slider = parabolic_slider.set_hmax(4e-6);
parabolic_slider = parabolic_slider.set_viscosity(0.039);
parabolic_slider = parabolic_slider.set_initial_height("parabolicslider");
parabolic_slider = parabolic_slider.set_epsilon("max", "max");

parabol_solution = solution;
parabol_solution.domain = parabolic_slider;
parabol_solution.velocity = 4.57;
parabol_solution.wiessenberg_Number = 10;
parabol_solution.viscocity_ratio = 0.8;

parabol_solution = parabol_solution.initilizer();
parabol_solution = parabol_solution.LIN("cavitation", "on", 0);

parabol_solution = parabol_solution.FieldProperties(501, 64);
% parabol_solution = parabol_solution.VR(5e-4, 501, 64, "CD2", "RK4", 0.05, "on");


%% Step Function
clear
Nx = 256; Ny = 64;
De = 0.05;

step_slider = domain;
step_slider.N_x = Nx;
step_slider = step_slider.set_length(1e-3);
step_slider = step_slider.set_hmax(4e-6);
step_slider = step_slider.set_viscosity(0.039);
step_slider = step_slider.set_initial_height("step");
step_slider = step_slider.set_epsilon("max", "max");

step_solution = solution;
step_solution.domain = step_slider;
step_solution.velocity = 4.57;
step_solution.wiessenberg_Number = 1;
step_solution.viscocity_ratio = 0.0;
step_solution = step_solution.initilizer();
step_solution = step_solution.LIN("cavitation","off", De);
max(De * step_solution.pressure1_old) -  min(De * step_solution.pressure1_old)

step_solution = step_solution.FieldProperties(Nx,Ny);
step_solution = step_solution.VR(2e-4, Nx, Ny, "CD2", "RK4", De, "off");
max(step_solution.pressure_field{5}) -  min(step_solution.pressure_field{5})


%% EHL Solver-It is a function now!
ehl_channel(length(load_list), length(velocity_list)) = domain;


for i = 1:length(load_list)
    for j = 1:length(velocity_list)
        ehl_channel(i, j) = ehl_channel(i, j).set_domain_coeff(velocity_list(j));   
        ehl_channel(i, j) = ehl_channel(i, j).set_eq_radii_curvature();
        ehl_channel(i, j) = ehl_channel(i, j).set_eq_elastic_modules(); 
        ehl_channel(i, j) = ehl_channel(i, j).set_bH(load_list(i));
        ehl_channel(i, j) = ehl_channel(i, j).set_ehl_length(velocity_list(j));
        ehl_channel(i, j) = ehl_channel(i, j).elastic_map();
        ehl_channel(i, j) = ehl_channel(i, j).set_initial_height("punch");
        ehl_channel(i, j) = ehl_channel(i, j).set_epsilon(h_cent(j), "None");       
    end
end
         
ehl_solution(length(load_list), length(velocity_list), length(wi_list)) = solution;

for velocity = 1:length(velocity_list)   %Number of Velocities
    Ux = velocity_list(velocity);
for wi = 1:length(wi_list)
    Wi = wi_list(wi);
for bi = 1:length(beta_list)
    beta = beta_list(bi);
for force = 1:length(load_list) 
    F = load_list(force);
    
    loadflag = true;
    dispflag = true;
    
    ehl_solution(force, velocity, wi).domain = ehl_channel(force, velocity);
    ehl_solution(force, velocity, wi) = ehl_solution(force, velocity, wi).set_solution_parameters(Wi, beta, Ux, F);
    
    disp("Solving for U = " + string(Ux) + " Wi = " + string(Wi) + " beta = " + string(beta) + " load = " + string(F))
    
    if velocity == 0.0001 && F >= 100
        ehl_solution(force, velocity, wi).relaxD = 1e-5;
        ehl_solution(force, velocity, wi).relaxDmax = 1e-3;
        ehl_solution(force, velocity, wi).relaxDmin = 1e-4;
    end
    
    ehl_solution(force, velocity, wi) = ehl_solution(force, velocity, wi).initilizer();
    ehl_solution(force, velocity, wi) = ehl_solution(force, velocity, wi).set_hertizan_pressure(F);
    
    diplacement_counter = 0;
    load_counter = 0;
    
    tic
    while loadflag 
        load_counter = load_counter + 1;
        while dispflag
            diplacement_counter = diplacement_counter + 1;
            
%             [ehl_solution(force, velocity, wi)] = ehl_solution(force, velocity, wi).Reynolds();
            [ehl_solution(force, velocity, wi)] = ehl_solution(force, velocity, wi).LIN();
%             [ehl_solution(force, l, wi), pressure] = ehl_solution(force, l, wi).VR(domain, h, Wi, beta, Ux);
            
            [ehl_solution(force, velocity, wi), dispflag] = ehl_solution(force, velocity, wi).displacementsolver();
            
            if mod(diplacement_counter, 100) == 0
                disp("Error: " + string(ehl_solution(force, velocity, wi).errD) + " at itteration " + string(diplacement_counter) + " Load balance at: " + string(load_counter))
            end
            
%             if diplacement_counter == 1
%                 break
%             end
            
        end
        
        [ehl_solution(force, velocity, wi), loadflag] = ehl_solution(force, velocity, wi).loadsolver(F);
        dispflag = true;
        diplacement_counter = 0;
        
%         if diplacement_counter == 1
%             break
%         end
        disp("-------------     v   -----------")
        disp("Error: " + string(ehl_solution(force, velocity, wi).errL) + " Load balance at: " + string(load_counter) + " U = " + string(Ux) + " Wi = " + string(Wi) + " L = " + string(F))
        disp("--------------------------------------")

    end
    
    ehl_solution(force, velocity, wi) = ehl_solution(force, velocity, wi).FieldProperties(501,120);

    disp("the EHL is solved")
    toc
    disp("           ")
    disp("           ")
    
end
end
end
end


%Save the data
ehl_savedata(load_list, beta_list, wi_list, velocity_list, ehl_solution, "25bH")

%%

for i = 1
%    [ehl_solution(:, :, i)] = ehl_solution(:, :, i).LIN("cavitation", "on", 0.05);
   [ehl_solution(1, 1, 1)] = ehl_solution(1, 1, 1).VR(5e-4, 501, 150, "CD2", "RK4", 0.02, "on");
%    [ehl_solution(:, :, i)] = ehl._solution(:, :, i).VR(h , wi_list(i), 0.6, 0.1);
end

%%

for i = 1:7
   figure(i); hold on
   plot(linspace(0, 1, 501), ehl_solution(1, 1, i).pressure_VR, "--", "LineWidth", 1.7);
   plot(linspace(0, 1, 501), ehl_solution(1, 1, i).pressure, "-", "LineWidth", 1.7);
   legend("VR", "LIN", "Interpreter", "latex")
   title("Wi = " + string(ehl_solution(1, 1, i).wiessenberg_Number))
   ylabel("Pa", "Interpreter", "latex")
   xlabel("x","Interpreter", "latex")
   
   figure(8); hold on
   plot(linspace(0, 1, 501), ehl_solution(1, 1, i).pressure_VR, "-", "LineWidth", 1.5);
   legend("Wi = 0.001", "Wi = 0.01", "Wi = 0.1", "Wi = 0.5", "Wi = 1", "Wi = 4", "Wi = 5" )   
   ylabel("Pa", "Interpreter", "latex")
   xlabel("x", "Interpreter", "latex")
   title("U = 0.1, L = 100 N", "Interpreter", "latex")
   
   figure(9); hold on
   plot(linspace(0, 1, 501), ehl_solution(1, 1, i).pressure, "-", "LineWidth", 1.5);
   legend("Wi = 0.001", "Wi = 0.01", "Wi = 0.1", "Wi = 0.5", "Wi = 1", "Wi = 4", "Wi = 5", "Interpreter", "latex")
   title("U = 0.1, L = 100 N", "Interpreter", "latex")
   ylabel("Pa", "Interpreter", "latex")
   xlabel("x", "Interpreter", "latex")
end



% 
% external_datapath = 'E:\Bilkent Dökümanları\Masterımsı\Kodlar\EHL_Solver\Result';  
% savepath = string(pwd) + "\Figure_Results";
%     
    
    
    
    
    %%
   
% eccen = [0.1, 0.05, 0.01, 0.005, 0.001];
eccen = [0.1, 0.08, 0.05, 0.02];
list_of_values = [];

for i = 1:length(eccen)
    
%     load(string(eccen(i)) + "alfa.mat")
    load(string(eccen(i)) + "alfaDe0.05.mat")
    peak_pressure_diff =  max(parabol_solution.pressure) - max(parabol_solution.pressure_VR);
    peak_pressure_diff =  max(parabol_solution.pressure) / max(parabol_solution.pressure_VR);
%     list_of_values(i) = peak_pressure_diff / max(parabol_solution.pressure);
    list_of_values(i) = peak_pressure_diff;

end

figure(2); hold on
plot(eccen, list_of_values, "k.-", "LineWidth", 2.4, "MarkerSize", 30);
xlabel("\boldmath{$h_{min} / h_{max}$}","Interpreter","latex")
ylabel("\boldmath{$pLIN_{max} / pVR_{max}$}","Interpreter","latex")
set(gca,'fontsize',24)
set(gca,'TickLabelInterpreter','latex')
title("De = 0.02","Interpreter","latex")

% load("Results/EHL_Results/F_100_beta_0.6_Wi_5_10_U_0.1_.mat")

load("Results/EHL_Results/F_100_400_beta_0.6_De0.05_U_0.1_.mat")


% ehl_solution(1, 1, 1) = ehl_solution(1, 1, 1).VR(5e-4, 501, 150, "CD2", "RK4", 0);

ehl_diff = max(ehl_solution(1,1,1).pressure) / max(ehl_solution(1,1,1).pressure_VR);
eccen_ehl = min(ehl_solution(1,1,1).h) / max(ehl_solution(1,1,1).h);

figure(2);
scatter(eccen_ehl, ehl_diff, 75,"r", "filled")

legend("ParabolicSlider", "EHL Profile (L100,U01)","Interpreter","latex","edgecolor","none","Box","off")


    
(0.38 / step_solution.domain.href) * 4 + 57 * (0.02/step_solution.domain.href) / ((0.38 / step_solution.domain.href) * 12 + 99 * (0.02/step_solution.domain.href));
    
    
    
    
