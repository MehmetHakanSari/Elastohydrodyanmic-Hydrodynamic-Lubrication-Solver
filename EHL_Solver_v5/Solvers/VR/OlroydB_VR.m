function new_stress = OlroydB_VR(Stress_Tensor, mesh, fluid, Scheme, velocity_field, vel_logic_map)

dx = mesh.dx;
dy = mesh.dy;
De = fluid.De;

u = velocity_field("u"); v = velocity_field("v");
DuDx = velocity_field("DuDx"); DvDx = velocity_field("DvDx");
dudy = velocity_field("dudy"); dvdy = velocity_field("dvdy");
viscos = (1 - fluid.beta);

new_stress = zeros(mesh.N_y, mesh.N_x, 3);

T_xx = Stress_Tensor(:,:,1);
T_xy = Stress_Tensor(:,:,2);
T_yy = Stress_Tensor(:,:,3);

if Scheme == "CD2"
    gradTxx = TwoDcentraldiff(T_xx, mesh.dx, mesh.dy);
    gradTxy = TwoDcentraldiff(T_xy,mesh.dx, mesh.dy);
    gradTyy = TwoDcentraldiff(T_yy, mesh.dx, mesh.dy);
elseif Scheme == "CD4"
    gradTxx = CD4(T_xx, mesh.dx, mesh.dy);
    gradTxy = CD4(T_xy, mesh.dx, mesh.dy);
    gradTyy = CD4(T_yy, mesh.dx, mesh.dy);
elseif Scheme == "UW"
    gradTxx = UpDownWind(T_xx, mesh.dx, mesh.dy, vel_logic_map);
    gradTxy = UpDownWind(T_xy, mesh.dx, mesh.dy, vel_logic_map);
    gradTyy = UpDownWind(T_yy, mesh.dx, mesh.dy, vel_logic_map);
end

DTxxDx = gradTxx{1} + gradTxx{2} .* mesh.Jacobian("dydx");
dTxxdy = gradTxx{2};
 
DTxyDx = gradTxy{1} + gradTxy{2} .*mesh.Jacobian("dydx");
dTxydy = gradTxy{2};

DTyyDx = gradTyy{1} + gradTyy{2} .* mesh.Jacobian("dydx");
dTyydy = gradTyy{2};

P_xx =  ((-T_xx / De) + 2*dudy.*T_xy + 2*T_xx.*DuDx - u.*DTxxDx - v.*dTxxdy);
P_xy =  (((-T_xy + viscos * dudy) / De) + dudy.*T_yy + T_xx.*DvDx - u.*DTxyDx - v.*dTxydy);
P_yy =  (((-T_yy + 2 * viscos * dvdy) / De) + 2*DvDx.*T_xy + 2*T_yy.*dvdy - u.*DTyyDx - v.*dTyydy);

new_stress(:,:,1) = P_xx;
new_stress(:,:,2) = P_xy;
new_stress(:,:,3) = P_yy;

end