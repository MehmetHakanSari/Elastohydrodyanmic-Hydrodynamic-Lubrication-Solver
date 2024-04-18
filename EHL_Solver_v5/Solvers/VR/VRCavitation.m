function [pressure_dic, theta, stress_dic, vel_dic] = VRCavitation(dt, Spatial_Scheme, Temporal_Scheme, solution, mesh, fluid, cavitation_flag, semi_vel_correction_flag, step_flag)

h = solution.h;
dhdx = OneDcentraldiff(h, mesh.dx);
d2hdx2 = OneDcentraldiff(dhdx, mesh.dx);

beta_t = (1 - fluid.beta);   % viscos = 1:full poylmer. viscos = 0: full solvent. viscos is reverse in main code. (1 - viscos) is needed.  

%Setting 
SSboundary = 10^(-7 - round(log10(dt^-1)));         %Steady-State Bound
Error_Boundary = 1000000;
% disp("-----------------------------------------")
velocity_field = solution.velocity_field;  vel_logic_map = solution.vel_logic_map;
stress_field = solution.stress_LIN;
u = velocity_field("u"); v = velocity_field("v");
DuDx = velocity_field("DuDx"); DvDx = velocity_field("DvDx");
dudy = velocity_field("dudy"); dvdy = velocity_field("dvdy");
T_xx = stress_field("tau_xx"); T_xy = stress_field("tau_xy"); T_yy = stress_field("tau_yy");
TD_xx = stress_field("tauD_xx"); TD_xy = stress_field("tauD_xy"); TD_yy = stress_field("tauD_yy");
vel_dic = "None";

if Temporal_Scheme == "RK4"
     if semi_vel_correction_flag == "on"
        disp("semi-step velocity correction approximation")
        [T_xx, T_xy, T_yy, DTxxDx, vel_dic] = ...
            StressSolver_RK4_vel_corr(mesh, fluid, Spatial_Scheme, Error_Boundary, SSboundary, dt, velocity_field, vel_logic_map, stress_field);
     else
        [T_xx, T_xy, T_yy, DTxxDx] = ...
            StressSolver_RK4(mesh, fluid, Spatial_Scheme, Error_Boundary, SSboundary, dt, velocity_field, vel_logic_map, stress_field); 
     end
elseif Temporal_Scheme == "IM_EU"
    [T_xx, T_xy, T_yy, DTxxDx] = ...
        StressSolver_IE_EU(mesh, fluid, Spatial_Scheme, Error_Boundary, SSboundary, dt, velocity_field, vel_logic_map, stress_field);
elseif Temporal_Scheme == "CN"
    if semi_vel_correction_flag == "on"
        disp("semi-step velocity correction approximation")
        [T_xx, T_xy, T_yy, DTxxDx, vel_dic] = ...
            StressSolver_CN_vel_corr(mesh, fluid, Spatial_Scheme, Error_Boundary, SSboundary, dt, velocity_field, vel_logic_map, step_flag);
    else
        [T_xx, T_xy, T_yy, DTxxDx] = ...
            StressSolver_CN(mesh, fluid, Spatial_Scheme, Error_Boundary, SSboundary, dt, velocity_field, vel_logic_map);
    end
elseif Temporal_Scheme == "EX_EU"
    [T_xx, T_xy, T_yy, DTxxDx] = ...
        StressSolver_EX_EU(mesh, fluid, Spatial_Scheme, Error_Boundary, SSboundary, dt, velocity_field, vel_logic_map, stress_field); 
end

if cavitation_flag == "on" 
    [pressure_dic, theta] = ...
        PressureSolver_Cavitated_v2(h, mesh, fluid, DTxxDx, T_xy);
elseif cavitation_flag == "off"
    pressure_dic = ...
        PressureSolver_v1(h, mesh, fluid, DTxxDx, T_xy);
    theta = ones(1, mesh.N_x);
end
stress_dic =  containers.Map({char("tau_xx"), char("tau_xy"), char("tau_yy")}, {T_xx, T_xy, T_yy});

end