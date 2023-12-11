function [T_xx, T_xy, T_yy, DTxxDx] = StressSolver_RK4(coeffs)

dx = coeffs{1};
De = coeffs{2};
T_xx = coeffs{3};
T_xy = coeffs{4};
T_yy = coeffs{5};
dy = coeffs{6};
u = coeffs{7};
% v = coeffs{8};
% DuDx = coeffs{9};
% DvDx = coeffs{10};
% dudy = coeffs{11};
% dvdy = coeffs{12};
dydx = coeffs{13};
% Scheme = coeffs{14};
% Error_Boundary = coeffs{15};
SSboundary = coeffs{16};
dt = coeffs{17};
% viscos = coeffs{18};
% x_MAT = coeffs{19};
% y_MAT = coeffs{20};

N_x = length(u(1,:));
N_y = length(u(:,1));

k = true;


if De == 0
    k = false;
end

message = 100;
counter = 0;

Stress = zeros(N_y, N_x, 3);
Stress(:,:,1) = T_xx;
Stress(:,:,2) = T_xy;
Stress(:,:,3) = T_yy;

tic
while k
    counter = counter + 1;
    
    k1 = OlroydB_VR(Stress, coeffs);
    k2 = OlroydB_VR(Stress + dt * k1 / 2, coeffs);
    k3 = OlroydB_VR(Stress + dt * k2 / 2, coeffs);
    k4 = OlroydB_VR(Stress + dt * k3, coeffs);
   
    Stress_New = Stress + dt * 1/6 * (k1 + 2*k2 + 2*k3 + k4);

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
gradTxx = TwoDcentraldiff(T_xx, dx, dy);
DTxxDx = gradTxx{1} + gradTxx{2} .* dydx;

end
