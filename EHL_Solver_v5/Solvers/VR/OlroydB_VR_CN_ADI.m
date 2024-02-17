function new_stress = OlroydB_VR_CN_ADI(Stress_Tensor, mesh, fluid, Spatial_Scheme, SSboundary, dt, velocity_field, vel_logic_map)

dx = mesh.dx;
dy = mesh.dy;
De = fluid.De;

u = velocity_field("u"); v = velocity_field("v");
DuDx = velocity_field("DuDx"); DvDx = velocity_field("DvDx");
dudy = velocity_field("dudy"); dvdy = velocity_field("dvdy");

dydx = mesh.Jacobian("dydx");
viscos = (1 - fluid.beta);

%
% T / De + d T dt + u d T dx + v d T dy - (.) (.) - (.) (.) = C
%                                                                               t + dt                                                                            t
% (T (t + dt) - T (t) / dt) + 1/2 [{T / De + u d T dx + v d T dy - (.) (.) - (.) (.)} + {u d T dx + v d T dy - (.) (.) - (.) (.)}]

new_stress = zeros(mesh.N_y, mesh.N_x, 3);

T_xx = Stress_Tensor(:,:,1);
T_xy = Stress_Tensor(:,:,2);
T_yy = Stress_Tensor(:,:,3);

inclusive_xx = -2 * DuDx / 2; inclusive_yy = -2 * dvdy / 2; inclusive_xy = 0 / 2;  % divide 2 due to CN implicit part

switch Spatial_Scheme
    case "UW" 
    gradTxx = UpDownWind(T_xx, mesh.dx, mesh.dy, vel_logic_map);
    gradTxy = UpDownWind(T_xy, mesh.dx, mesh.dy, vel_logic_map);
    gradTyy = UpDownWind(T_yy, mesh.dx, mesh.dy, vel_logic_map);
    case "CD2"
    gradTxx = TwoDcentraldiff(T_xx, dx,dy);
    gradTxy = TwoDcentraldiff(T_xy, dx,dy);
    gradTyy = TwoDcentraldiff(T_yy, dx,dy);
%     case "CD4"
%     gradTxx = TwoDcentraldiff(T_xx, dx,dy);
%     gradTxy = TwoDcentraldiff(T_xy, dx,dy);
%     gradTyy = TwoDcentraldiff(T_yy, dx,dy);
end



DTxxDx = gradTxx{1} + gradTxx{2} .* dydx;
DTxyDx = gradTxy{1} + gradTxy{2} .* dydx;
DTyyDx = gradTyy{1} + gradTyy{2} .* dydx;

dTxxdy = gradTxx{2};
dTxydy = gradTxy{2};
dTyydy = gradTyy{2};


C_xx = 0; C_xy = viscos * dudy / De; C_yy = 2 * viscos * dvdy / De;

C_xx = C_xx + 2 * T_xy .* dudy - 1/2 * (T_xx / De + u .* DTxxDx + v .* dTxxdy - 2 * T_xx .* DuDx);

T_xx = HalfStepCNADI(T_xx, inclusive_xx, C_xx, u, v, dx, dy, dydx, dt, De, Spatial_Scheme, vel_logic_map);
C_xy = C_xy + T_yy .* dudy +  T_xx .* DvDx - 1/2 * (T_xy / De + u .* DTxyDx + v .* dTxydy);
T_xy = HalfStepCNADI(T_xy, inclusive_xy, C_xy, u, v, dx, dy, dydx, dt, De, Spatial_Scheme, vel_logic_map);
C_yy = C_yy + 2 * T_xy .* DvDx - 1/2 * (T_yy / De + u .* DTyyDx + v .* dTyydy - 2 * T_yy .* dvdy); 
T_yy = HalfStepCNADI(T_yy, inclusive_yy, C_yy, u, v, dx, dy, dydx, dt, De, Spatial_Scheme, vel_logic_map);

% new_stress(:,:,1) = HalfStepCNADI(T_xx, inclusive_xx, C_xx, u, v, dx, dy, dydx, dt, De, Scheme);
% new_stress(:,:,2) = HalfStepCNADI(T_xy, inclusive_xy, C_xy, u, v, dx, dy, dydx, dt, De, Scheme);
% new_stress(:,:,3) = HalfStepCNADI(T_yy, inclusive_yy, C_yy, u, v, dx, dy, dydx, dt, De, Scheme);

new_stress(:,:,1) = T_xx;
new_stress(:,:,2) = T_xy;
new_stress(:,:,3) = T_yy;

end