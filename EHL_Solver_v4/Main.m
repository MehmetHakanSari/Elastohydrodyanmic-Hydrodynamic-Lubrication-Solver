% Main Script

clear 
clc

addpath('Differentials/');
addpath('SaveTools/'); addpath('SaveTools/EHL_Plots');
addpath('Solvers/');
addpath('PostProcess/');
addpath('data_management/');




%% EHL case
velocity_list = [0.00056234];
de_list = [0, 0.01, 0.02, 0.05];
beta_list = [0.6];
load_list = [30, 50, 100, 200];
% h_cent = [1.8144e-07,5.4274e-07,1.8297e-06, 8.7596e-06,4.6202e-05]; %from U = 0.0001 to U = 0.1

% EHL_Solver_Wi(ehl_solution, ehl_channel, velocity_list, wi_list, beta_list, load_list, "") 
EHL_Solver_De(velocity_list, de_list, beta_list, load_list, h_cent, "centh") 


%% Parabolic Slider
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


step_solution = step_solution.FieldProperties(Nx,Ny);
step_solution = step_solution.VR(2e-4, Nx, Ny, "CD2", "RK4", De, "off");



    
