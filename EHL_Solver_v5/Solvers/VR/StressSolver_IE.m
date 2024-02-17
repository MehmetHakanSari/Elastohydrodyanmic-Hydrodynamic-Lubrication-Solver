function [T_xx, T_xy, T_yy, DTxxDx] = StressSolver_IE(coeffs)

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

tic
while k
    counter = counter + 1;
    
    Stress_New = OlroydB_VR_Implicit(Stress, coeffs);
    
    Stress_New(:,1,:) = Stress_New(:,2,:);
    Stress_New(end,:,:) = Stress_New(end-1,:,:);
    
    residual = max(abs(Stress_New - Stress),[],"all");
    
    if residual < SSboundary
        k = false;
    end
    Stress = Stress_New;
    
    if mod(counter,message) == 0
        timer = toc;
        disp(string(counter) + " Itterations: "  + " largest residual is: " + string(residual) + " ellapsed time: " + string(timer))     
    end
    
end

T_xx = Stress(:,:,1);
T_xy = Stress(:,:,2);
T_yy = Stress(:,:,3);
gradTxx = TwoDcentraldiff(T_xx, dx, dy);
DTxxDx = gradTxx{1} + gradTxx{2} .* dydx;

end
