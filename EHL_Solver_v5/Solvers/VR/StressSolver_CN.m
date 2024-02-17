function [T_xx, T_xy, T_yy, DTxxDx] = ...
    StressSolver_CN(mesh, fluid, Spatial_Scheme, Error_Boundary, SSboundary, dt, velocity_field, vel_logic_map)

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
%     Stress_New = OlroydB_VR_Implicit(Stress, coeffs);
    
    Stress_New(:,1,:) = Stress_New(:,2,:);
    Stress_New(end,:,:) = Stress_New(end-1,:,:);
    
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
    end
    
end

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

end