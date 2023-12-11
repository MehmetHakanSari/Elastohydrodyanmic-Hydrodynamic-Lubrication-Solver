function new_stress = OlroydB_VR_Implicit(Stress_Tensor, coeffs)

dx = coeffs{1};
De = coeffs{2};
dy = coeffs{6};
u = coeffs{7};
v = coeffs{8};
DuDx = coeffs{9};
DvDx = coeffs{10};
dudy = coeffs{11};
dvdy = coeffs{12};
dydx = coeffs{13};
Scheme = coeffs{14};
viscos = coeffs{18};
dt = coeffs{17};

N = length(u(1,:));

new_stress = zeros(N,N,3);

T_xx = Stress_Tensor(:,:,1);
T_xy = Stress_Tensor(:,:,2);
T_yy = Stress_Tensor(:,:,3);

if Scheme == "CD2"
    gradTxx = TwoDcentraldiff(T_xx, dx,dy);
    gradTxy = TwoDcentraldiff(T_xy, dx,dy);
    gradTyy = TwoDcentraldiff(T_yy, dx,dy);
elseif Scheme == "CD4"
    gradTxx = CD4(T_xx, dx,dy);
    gradTxy = CD4(T_xy, dx,dy);
    gradTyy = CD4(T_yy, dx,dy);
elseif Scheme == "UW"
    gradTxx = UpDownWind(T_xx, dx,dy);
    gradTxy = UpDownWind(T_xy, dx,dy);
    gradTyy = UpDownWind(T_yy, dx,dy);
end

DTxxDx = gradTxx{1} + gradTxx{2} .* dydx;
dTxxdy = gradTxx{2};

DTxyDx = gradTxy{1} + gradTxy{2} .* dydx;
dTxydy = gradTxy{2};

DTyyDx = gradTyy{1} + gradTyy{2} .* dydx;
dTyydy = gradTyy{2};

% 
% P_xx =  ((-T_xx/De) + 2*dudy.*T_xy + 2*T_xx.*DuDx - u.*DTxxDx - v.*dTxxdy);
% P_xy =  (((-T_xy + viscos * dudy) /De) + dudy.*T_yy + T_xx.*DvDx - u.*DTxyDx - v.*dTxydy);
% P_yy =  (((-T_yy + 2 * viscos * dvdy) /De) + 2*DvDx.*T_xy + 2*T_yy.*dvdy - u.*DTyyDx - v.*dTyydy);


for j = 1:N
    for i = 1:N
        
   
        


DM = [1 + De/dt - 2*De*DuDx(j,i), -2 * De * dudy(j,i), 0;
 -De * DvDx(j,i), 1 + De / dt, - De * dudy(j,i);
 0, -2 * De * DvDx(j,i), 1 + De / dt - 2 * De * dvdy(j,i)]; 

Sol_Vec = [De * T_xx(j,i) / dt - De * (u(j,i) * DTxxDx(j,i) + v(j,i) * dTxxdy(j,i));
            De * T_xy(j,i) / dt - De * (u(j,i) * DTxyDx(j,i) + v(j,i) * dTxydy(j,i)) + viscos * dudy(j,i);
            De * T_yy(j,i) / dt - De * (u(j,i) * DTyyDx(j,i) + v(j,i) * dTyydy(j,i)) + 2 * viscos * dvdy(j,i)];
        
stress_point = DM\Sol_Vec;

new_stress(j,i,:) = stress_point;

    end
end

% P_xx(:,1) = P_xx(:,2);
% P_xy(:,1) = P_xy(:,2);
% P_yy(:,1) = P_yy(:,2);
% 
% P_xx(end,:) = P_xx(end-1,:);
% P_xy(end,:) = P_xy(end-1,:);
% P_yy(end,:) = P_yy(end-1,:);

% new_stress(:,:,1) = P_xx;
% new_stress(:,:,2) = P_xy;
% new_stress(:,:,3) = P_yy;


end