function new_stress = OlroydB_VR(Stress_Tensor, coeffs)

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

N_x = length(u(1,:));
N_y = length(u(:,1));

new_stress = zeros(N_y, N_x, 3);

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
    gradTxx = UpDownWind(T_xx, dx,dy, u, v);
    gradTxy = UpDownWind(T_xy, dx,dy, u, v);
    gradTyy = UpDownWind(T_yy, dx,dy, u, v);
end

DTxxDx = gradTxx{1} + gradTxx{2} .* dydx;
dTxxdy = gradTxx{2};
 
DTxyDx = gradTxy{1} + gradTxy{2} .* dydx;
dTxydy = gradTxy{2};

DTyyDx = gradTyy{1} + gradTyy{2} .* dydx;
dTyydy = gradTyy{2};


P_xx =  ((-T_xx / De) + 2*dudy.*T_xy + 2*T_xx.*DuDx - u.*DTxxDx - v.*dTxxdy);
P_xy =  (((-T_xy + viscos * dudy) / De) + dudy.*T_yy + T_xx.*DvDx - u.*DTxyDx - v.*dTxydy);
P_yy =  (((-T_yy + 2 * viscos * dvdy) / De) + 2*DvDx.*T_xy + 2*T_yy.*dvdy - u.*DTyyDx - v.*dTyydy);

% P_xx(:,1) = P_xx(:,2);
% P_xy(:,1) = P_xy(:,2);
% P_yy(:,1) = P_yy(:,2);
% 
% P_xx(end,:) = P_xx(end-1,:);
% P_xy(end,:) = P_xy(end-1,:);
% P_yy(end,:) = P_yy(end-1,:);

new_stress(:,:,1) = P_xx;
new_stress(:,:,2) = P_xy;
new_stress(:,:,3) = P_yy;


end