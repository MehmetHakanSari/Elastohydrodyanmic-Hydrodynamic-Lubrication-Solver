function new_stress = OlroydB_VR_Implicit_ADI(Stress_Tensor, coeffs)

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

[Ny, Nx] = size(u);

new_stress = zeros(Ny, Nx, 3);

T_xx = Stress_Tensor(:,:,1);
T_xy = Stress_Tensor(:,:,2);
T_yy = Stress_Tensor(:,:,3);

inclusive_xx = -2 * DuDx; inclusive_yy = -2 * dvdy; inclusive_xy = 0;

C_xx = 0; C_xy = viscos * dudy / De; C_yy = 2 * viscos * dvdy / De;
C_xx = C_xx + 2 * T_xy .* dudy;
C_xy = C_xy + T_xx .* dvdy + T_yy .* dvdy;
C_yy = C_yy + 2 * T_xy .* DvDx; 

new_stress(:,:,1) = HalfStepADI(T_xx, inclusive_xx, C_xx, u, v, dx, dy, dydx, dt, De, Scheme);
new_stress(:,:,2) = HalfStepADI(T_xy, inclusive_xy, C_xy, u, v, dx, dy, dydx, dt, De, Scheme);
new_stress(:,:,3) = HalfStepADI(T_yy, inclusive_yy, C_yy, u, v, dx, dy, dydx, dt, De, Scheme);

end