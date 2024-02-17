function [T_xx, T_xy, T_yy, DTxxDx] = StressSolver_v2(coeffs)

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

N = length(u(1,:));
J = length(u(:,1));



k = true;


if De == 0
    k = false;
end

message = 100;
counter = 0;

Stress = zeros(N,N,3);
Stress(:,:,1) = T_xx;
Stress(:,:,2) = T_xy;
Stress(:,:,3) = T_yy;

while k
    counter = counter + 1;
    
    if mod(counter,message) == 0
        timer = toc;
        disp(string(counter) + " Itterations: "  + " largest residual is: " + string(residual) + " ellapsed time: " + string(timer))     
    end
    
%     stress_coeffs = {N, dx, dy, dydx, DuDx, DvDx, dudy, dvdy, viscos, u, v, Scheme, De};
    
    k1 = OlroydB_VR(Stress, coeffs);
    k2 = OlroydB_VR(Stress + dt * k1 / 2, coeffs);
    k3 = OlroydB_VR(Stress + dt * k2 / 2, coeffs);
    k4 = OlroydB_VR(Stress + dt * k3, coeffs);
%   
    Stress_New = Stress + dt * 1/6 * (k1 + 2*k2 + 2*k3 + k4);
%     expl_euler = Stress + dt .* k1;
    
    
%     expl_euler(:,1,:) = expl_euler(:,2,:);
%     expl_euler(end,:,:) = expl_euler(end-1,:,:);
    
    Stress_New(:,1,:) = Stress_New(:,2,:);
    Stress_New(end,:,:) = Stress_New(end-1,:,:);
    
    residual = max(abs(Stress_New - Stress),[],"all");
%     residual_ = max(max(max(abs(expl_euler - Stress))));
    
    if residual < SSboundary
%         disp("Time step at: " + string(counter))
%         disp("Time: " + string(counter*dt) + " seconds")
        k = false;
    end
    Stress = Stress_New;
    

end


T_xx = Stress(:,:,1);
T_xy = Stress(:,:,2);
T_yy = Stress(:,:,3);
gradTxx = TwoDcentraldiff(T_xx, dx, dy);
DTxxDx = gradTxx{1} + gradTxx{2} .* dydx;

end
