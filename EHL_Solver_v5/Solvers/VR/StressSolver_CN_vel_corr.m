function [T_xx, T_xy, T_yy, DTxxDx, velocity_field] = ...
    StressSolver_CN_vel_corr(mesh, fluid, Spatial_Scheme, Error_Boundary, SSboundary, dt, velocity_field, vel_logic_map, step_flag)

uD = velocity_field("uD");  vD = velocity_field("vD"); 
u_bench = velocity_field("u");

k = true;
if fluid.De == 0
    k = false;
end

message = 100;
counter = 0;

Stress = zeros(mesh.N_y, mesh.N_x, 3);
Stress(:,:,1) = zeros(mesh.N_y, mesh.N_x);%T_xx;
Stress(:,:,2) = zeros(mesh.N_y, mesh.N_x);%T_xy;
Stress(:,:,3) = zeros(mesh.N_y, mesh.N_x);%T_yy;

tic
while k
    counter = counter + 1;
    
    Stress_New = OlroydB_VR_CN_ADI(Stress, mesh, fluid, Spatial_Scheme, SSboundary, dt, velocity_field, vel_logic_map);
    
    Stress_New(:,1,:) = Stress_New(:,2,:);
    Stress_New(end,:,:) = Stress_New(end-1,:,:);
    
%     if step_flag == "on" 
%      uD_anal_Step = (1/2 * (1 - fluid.beta) * dhdx ./h - (1 - fluid.beta) / 8 * h.^3 .* dhdx .* (dpdx).^2 + 1/2 * h.^2 .* dpDdx ) .* (Y.^2 ./ h.^2 - Y ./ h);
%     else
    gradTxx = TwoDcentraldiff(Stress_New(:,:,1), mesh.dx, mesh.dy);
    DTxxDx = gradTxx{1} + gradTxx{2} .* mesh.Jacobian("dydx");
    pressure_dic = PressureSolver_v1(mesh.channel.height, mesh, fluid, DTxxDx, Stress_New(:,:,2));
    velocity_field = velocity_updater(mesh.channel.height, mesh, fluid, pressure_dic("p"), DTxxDx, Stress_New(:,:,2), velocity_field("dudy"));
%     end
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
    
    if residual_xx_nan(1) || residual_xx_nan(1) || residual_xx_nan(1)
       error("Gradient Explodes") 
    end
    
    if residual_xx < SSboundary && residual_xy < SSboundary && residual_yy < SSboundary
        k = false;
    end
    Stress = Stress_New;
    
    if mod(counter,message) == 0
        timer = toc;
        disp(string(counter) + " Itterations: "  + " " + string(residual_xx) + "   " + string(residual_xy) + "  " + string(residual_yy) + " ellapsed time: " + string(timer))     
%         figure(1);
%         u = velocity_field("u");
%         T_yy = Stress_New(:,:,3);
%         surf(u, "linestyle", "none")
%         figure(1); 
%         plot(T_yy(1, : ), "--" ,"LineWidth" , 1.8)
%         plot(T_yy(17, :), "-." ,"LineWidth" , 1.8)
%         plot(T_yy(end, :), ":" ,"LineWidth" , 1.8)
%         set(figure(1),'Position',[20 90 450 450])
%         figure(2);
%         plot(pressure_dic("p"), "LineWidth", 1.6)
%         pause(0.15)
%         set(figure(2),'Position',[500 90 450 450])
%         xlim([505 512])
    end
    
%      if counter == 1000
%         break
%     end
    

    
end


% u = velocity_field("u");
% 
% figure(1); hold on;
% plot(u(:,1), "--")
% plot(u_bench(:,1), "-.")
% figure(2); hold on;
% plot(u(:,17), "--")
% plot(u_bench(:,17), "-.")
% figure(3); hold on;
% plot(u(:,end), "--")
% plot(u_bench(:,end), "-.")





T_xx = Stress(:,:,1);
T_xy = Stress(:,:,2);
T_yy = Stress(:,:,3);
gradTxx = TwoDcentraldiff(T_xx, mesh.dx, mesh.dy);
% switch Spatial_Scheme
%     case "CD2"
%         gradTxx = TwoDcentraldiff(T_xx, mesh.dx, mesh.dy);
%     case "UW"
%         gradTxx = UpDownWind(T_xx, mesh.dx, mesh.dy, fluid_domain.vel_logic_map);
% end

DTxxDx = gradTxx{1} + gradTxx{2} .* mesh.Jacobian("dydx");
velocity_field("uD") = uD;
velocity_field("vD") = vD;

end