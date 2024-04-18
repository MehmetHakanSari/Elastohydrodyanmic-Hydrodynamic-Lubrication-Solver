function [T_xx, T_xy, T_yy, DTxxDx, velocity_field] = ...
    StressSolver_RK4_vel_corr(mesh, fluid, Spatial_Scheme, Error_Boundary, SSboundary, dt, velocity_field, vel_logic_map, stress_field)

Scheme = Spatial_Scheme;
uD = velocity_field("uD");  vD = velocity_field("vD"); 

k = true;
if fluid.De == 0
    k = false;
end

message = 100;
counter = 0;

Stress = zeros(mesh.N_y, mesh.N_x, 3);
Stress(:,:,1) = stress_field("tau_xx");
Stress(:,:,2) = stress_field("tau_xy");
Stress(:,:,3) = stress_field("tau_yy");
% Stress(:,:,1) = zeros(mesh.N_y, mesh.N_x);%T_xx;
% Stress(:,:,2) = zeros(mesh.N_y, mesh.N_x);%T_xy;
% Stress(:,:,3) = zeros(mesh.N_y, mesh.N_x);%T_yy;

tic
while k
    counter = counter + 1;
    
    k1 = OlroydB_VR(Stress, mesh, fluid, Scheme, velocity_field, vel_logic_map);
    k2 = OlroydB_VR(Stress + dt * k1 / 2, mesh, fluid, Scheme, velocity_field, vel_logic_map);
    k3 = OlroydB_VR(Stress + dt * k2 / 2, mesh, fluid, Scheme, velocity_field, vel_logic_map);
    k4 = OlroydB_VR(Stress + dt * k3, mesh, fluid, Scheme, velocity_field, vel_logic_map);
   
    Stress_New = Stress + dt * 1/6 * (k1 + 2*k2 + 2*k3 + k4);

    Stress_New(:,1,:) = Stress_New(:,2,:);
    Stress_New(end,:,:) = Stress_New(end-1,:,:);
    
    gradTxx = TwoDcentraldiff(Stress_New(:,:,1), mesh.dx, mesh.dy);
    DTxxDx = gradTxx{1} + gradTxx{2} .* mesh.Jacobian("dydx");
    pressure_dic = PressureSolver_v1(mesh.channel.height, mesh, fluid, DTxxDx, Stress_New(:,:,2));
    velocity_field = velocity_updater(mesh.channel.height, mesh, fluid, pressure_dic("p"), DTxxDx, Stress_New(:,:,2), velocity_field("dudy"));
    vel_logic_map = containers.Map;
    vel_logic_map("u_uw") = velocity_field("u") > 0;
    vel_logic_map("v_uw") = velocity_field("v") > 0;
    vel_logic_map("u_dw") = velocity_field("u") <= 0;
    vel_logic_map("v_dw") = velocity_field("v") <= 0;
    
    residual_xx = max(abs(Stress_New(:,:,1) - Stress(:,:,1)),[],"all");
    residual_xy = max(abs(Stress_New(:,:,2) - Stress(:,:,2)),[],"all");
    residual_yy = max(abs(Stress_New(:,:,3) - Stress(:,:,3)),[],"all");
    
%     residual = max(abs(Stress_New - Stress),[],"all");
    residual_xx_nan = [residual_xx, 1];
    residual_xy_nan = [residual_xy, 1];
    residual_yy_nan = [residual_yy, 1];
%     disp(residual_nan)
    residual_xx_nan = isnan(residual_xx_nan);
    residual_xy_nan = isnan(residual_xy_nan);
    residual_yy_nan = isnan(residual_yy_nan);
    
    if residual_xx_nan(1) || residual_xy_nan(1) || residual_yy_nan(1)
       error("Gradient Explodes") 
    end
    
    if residual_xx < SSboundary && residual_xy < SSboundary && residual_yy < SSboundary
        k = false;
    end
    Stress = Stress_New;
    
    if mod(counter,message) == 0
        timer = toc;
        disp(string(counter) + " Itterations: "  + " " + string(residual_xx) + "   " + string(residual_xy) + "  " + string(residual_yy) + " ellapsed time: " + string(timer))     
%                 figure(1);
%         u = velocity_field("u");
%         T_yy = Stress_New(:,:,3);
% %         surf(u, "linestyle", "none")
%         figure(1); 
%         plot(T_yy(1, : ), "--" ,"LineWidth" , 1.8)
% %         plot(T_yy(17, :), "-." ,"LineWidth" , 1.8)
% %         plot(T_yy(end, :), ":" ,"LineWidth" , 1.8)
%         set(figure(1),'Position',[20 90 450 450])
%         figure(2);
%         plot(pressure_dic("p"), "LineWidth", 1.6)
%         pause(0.15)
%         set(figure(2),'Position',[500 90 450 450])
%         xlim([505 512])
    end
%     
%      if counter ==500
%         break
%     end
%     
end

T_xx = Stress(:,:,1);
T_xy = Stress(:,:,2);
T_yy = Stress(:,:,3);
gradTxx = TwoDcentraldiff(T_xx, mesh.dx, mesh.dy);
DTxxDx = gradTxx{1} + gradTxx{2} .* mesh.Jacobian("dydx");
velocity_field("uD") = uD;
velocity_field("vD") = vD;

end
